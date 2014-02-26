use MooseX::Declare;
use Method::Signatures::Modifiers;
use feature qw(say);

class Biotools::CumulativeDensity {
    use MooseX::App::Command;            # important
    extends qw(Biotools);                
    use Bio::Moose::BedIO;
    use Bio::Moose::BedTools::Intersect;
    use Bio::Moose::BedTools::Complement;
    use Bio::Moose::BedTools::Slop;
    use Bio::Moose::BedTools::Flank;
    use Bio::Moose::BedTools::WindowMaker;
    use Moose::Util::TypeConstraints;
    use File::Temp;
    use File::Basename;
    use File::Copy;
    use File::Path qw(make_path rmtree);
    use AnyEvent;
    use AnyEvent::ForkManager;
    use Progress::Any::Output;
    use Sys::CPU;

    Progress::Any::Output->set('TermProgressBarColor');
 

    option 'input' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases   => 'i',
        required      => 1,
        documentation => q[Input genesBed File with 6 columns only!],
    );

    option 'output_file' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases   => 'o',
        required      => 1,
        documentation => q[Output filename],
    );

    option 'reads' => (
        is            => 'rw',
        isa           => 'Str',
        required      => 1,
        documentation => q[Input reads File],
    );

    option 'genome' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases   => 'g',
        required      => 1,
        documentation => q[genome],
    );

    option 'tss' => (
        is            => 'rw',
        isa           => 'Num',
        cmd_aliases   => 'l',
        default       => 2000,
        documentation => q[Amount to be subtracted from TSS in bp.],
    );

    option 'tts' => (
        is            => 'rw',
        isa           => 'Num',
        cmd_aliases   => 'r',
        default       => 2000,
        documentation => q[Amount to be subtracted from TTS in bp.],
    );

    option 'body_resolution' => (
        is            => 'rw',
        isa           => 'Int',
        cmd_aliases   => 'b',
        required      => 1,
        default       => 40,
        documentation => q[Number of body bins],
    );

    option 'window_size' => (
        is            => 'rw',
        isa           => 'Num',
        cmd_aliases   => 'w',
        default       => 100,
        documentation => q[Window size body],
    );

    option 'remove_overlapping_genes' => (
        is            => 'rw',
        isa           => 'Bool',
        default       => 0,
        documentation => q[Remove TSS+body+TTS overlapping genes form analysis],
    );

    option 'only_genes_with_reads' => (
        is            => 'rw',
        isa           => 'Bool',
        default       => 0,
        documentation => q[Only analyse genes with reads],
    );

    option 'n_forks' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases   => 'n',
        required      => 1,
        default       => 1,
        documentation => q[Number of forks to use],
    );

    has '_filter_input' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        lazy          => 1,
        builder       => 'build_filter_input',
        documentation => 'Hold filtered input',
    );

    has 'sense_genes' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        lazy => 1,
        default => sub{
            my $self = shift;
            my @aux;
            foreach my $f (@{$self->_filter_input}) {
                push @aux,$f if $f->strand eq '+';
            }
            return \@aux;
        },
        documentation => 'Positive filtered genes',
    );
    
    has 'antisense_genes' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        lazy => 1,
        default => sub{
            my $self = shift;
            my @aux;
            foreach my $f (@{$self->_filter_input}) {
                push @aux,$f if $f->strand eq '-';
            }
            return \@aux;
        },
        documentation => 'Negative filtered genes',
    );
   

    has 'relativeCoordSize' => (
        is            => 'rw',
        isa           => 'Int',
        required      => 1,
        default       => 10000,
        documentation => 'Body size to plot',
    );

    has 'bodyStep' => (
        is       => 'rw',
        isa      => 'Num',
        lazy => 1,
        default  => sub {
            my ($self) = @_;
            return $self->relativeCoordSize / $self->body_resolution;
        },
        documentation => 'Keep body step size',
    );
    
    has 'n_reads' => (
        is            => 'rw',
        isa           => 'Int',
        lazy      => 1,
        default => sub {
            my ($self) = @_;
            open( my $in, '<', $self->reads ) 
                || die "Cannot open/read file " . $self->reads . "!";
            my $n_reads=0;
            while ( my $row = <$in> ){
               $n_reads++  unless ($row =~ /chrM/i || $row =~ /random/i )
            }
            close( $in );

            return $n_reads;    
        },
        documentation => 'Read object',
    );
   
    has 'normFactor' => (
        is      => 'rw',
        isa     => 'Num',
        lazy    => 1,
        default => sub {
            my ($self) = @_;
            my $normfactor = (1e6 / $self->n_reads );
            return $normfactor;
        },
        documentation => 'Normalize by this factor',
    );


    method build_filter_input {
        $self->log->info( "Filtering " . $self->input );
        $self->log->info( "Increasing TSS and TTS  and removing overlap");
        
        # get slopped genes
        my $slopped_genes = $self->_get_slopped_input->stdout;
        my @slop = split "\n", $slopped_genes->[0];
        $self->log->info(
            "   - Number of genes before filter: " . scalar @slop );
        
        # intersect genes and remove genes with more than 1 intersection
        my $gene_i_gene = Bio::Moose::BedTools::Intersect->new(
            a => $slopped_genes,
            b => $slopped_genes,
            c => 1
        );
        $self->log->info( $gene_i_gene->show_cmd_line );
        $gene_i_gene->run;

        my $orig_input = Bio::Moose::BedIO->new(file=>$self->input)->features;

        my @genes_no_overlap;
        my $i =0;
        foreach my $f ( @{ $gene_i_gene->as_BedIO->features } ) {
            if ( $f->thickStart == 1){
                push @genes_no_overlap, $orig_input->[$i];
            }
            $i++;
        }

        $self->log->info( "   - Number of genes without overlap: "
                . scalar @genes_no_overlap );

        my $file;

        if ($self->remove_overlapping_genes){
            $self->log->info( "      Removing overlapping genes ( --remove_overlapping_genes = 1 )" );
            $file  = \@genes_no_overlap;
        }else{
            $self->log->info( "      Keeping overlapping genes ( --remove_overlapping_genes = 0 )" );
            $file = $orig_input;
        }

        my $in    = Bio::Moose::BedIO->new( file => $file );
        my $feats = $in->features;

        $self->log->info(
            "   - Number of genes before filter: " . scalar @{$feats} );

        my @aux;
        my $min_gene_size = ( $self->tss + $self->tts ) * 0.5 ;

        $self->log->info( "   - Removing genes smaller than " . $min_gene_size );
        $self->log->info("   - Removing chrM and chr*_random");

        foreach my $f ( @{$feats} ) {
            unless ( $f->chrom =~ /chrM/ || $f->chrom =~ /random/ ) {
                if ($f->size >= $min_gene_size){
                   #fix gene size
                   push @aux,$f;
                }
            }
        }

        $self->log->info( "   - Number of genes after filter: " . scalar @aux );
        return \@aux;
    }


    method _get_slopped_input {
        # Create slopBed object
        $self->log->info("   Slopping ");
        $self->log->info("   -TSS size: ".$self->tss);
        $self->log->info("   -TTS size: ".$self->tts);
        my $slop = Bio::Moose::BedTools::Slop->new(
            i => $self->input,
            g => $self->genome,
            l => $self->tss,
            r => $self->tts,
            s => 1,
        );
        $self->log->info($slop->show_cmd_line);
        # Run slopbed
        $slop->run();
        return $slop;
    }
 

=cut
    method _get_tss_genes (ArrayRef[Bio::Moose::Bed]|Bio::Moose::Bed $genes) {
        # Create slopBed object
        $self->log->trace("   Getting TSS ");
        $self->log->trace( "   -TSS size: " . $self->tss );
        my $flank = Bio::Moose::BedTools::Flank->new(
            i => $genes,
            g => $self->genome,
            l => $self->tss,
            r => 0,
            s => 1,
        );
        $self->log->trace( $flank->show_cmd_line );

        # Run slopbed
        $flank->run();

        # if there is a temporary file, remove it
        if ( -e $flank->i ){
            unlink $flank->i;
        }

        return $flank;
    }


    method _get_tts_genes (ArrayRef[Bio::Moose::Bed]|Bio::Moose::Bed $genes) {
        # Create slopBed object
        $self->log->trace("   Getting TTS ");
        $self->log->trace( "   -TTS size: " . $self->tts );
        my $flank = Bio::Moose::BedTools::Flank->new(
            i => $genes,
            g => $self->genome,
            l => 0,
            r => $self->tts,
            s => 1,
        );
        $self->log->trace( $flank->show_cmd_line );

        # Run slopbed
        $flank->run();

        # if there is a temporary file, remove it
        if ( -e $flank->i ){
            unlink $flank->i;
        }

        return $flank;
    }

        
    method build_body_bins($this_input) {
        $self->log->trace("   Divide genes body in bins");
        $self->log->trace("   - body resolution : ".$self->body_resolution.' bins');
        my $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $this_input,
            i => 'winnum',
            n => $self->body_resolution,
        );
        $self->log->trace($windows->show_cmd_line);
        $windows->run;

        # if there is a temporary file, remove it
        if ( -e $windows->b ){
            unlink $windows->b;
        }

        return $windows;
    }


    method intersect_bins ($this_input) {
        my $intersect = Bio::Moose::BedTools::Intersect->new(
            a => $this_input,
            b => $self->reads,
            c => 1
        );

        $intersect->run;

        # Removing temporary bed file
        unlink $intersect->a if -e $intersect->a;
        return $intersect;
    }


    method build_fixed_bins ($this_input) {
        $self->log->trace("   Divide in fixed bins");
        $self->log->trace("   -window size: ".$self->window_size);
        my $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $this_input,
            w => $self->window_size,
            i => 'winnum'
            #n => $self->body_resolution,
        );
        $self->log->trace($windows->show_cmd_line);
        $windows->run;
        # if there is a temporary file, remove it
        if ( -e $windows->b ){
            unlink $windows->b;
        }
        return $windows;
    }
=cut
=cut
    method process_gene (Object $pm, Bio::Moose::Bed $gene_interval, Int $i, Str $tmp_dir, Str $strand) {
        # FOR BODY
        # creating bins
        my $windows_body   = $self->build_body_bins($gene_interval)->stdout;
        my $intersect_body = $self->intersect_bins($windows_body);

        my $output = "$tmp_dir/$strand/body/" . $gene_interval->name . '_' . $i . '.bed';
        open( my $out, '>', $output )
            || die "Cannot open/read file " . $output . "!";
        print $out @{ $intersect_body->stdout };
        close($out);

        # FOR TSS
        # Create slopBed object
        my $tss           = $self->_get_tss_genes($gene_interval);
        my $windows_tss   = $self->build_fixed_bins( $tss->stdout )->stdout;
        my $intersect_tss = $self->intersect_bins($windows_tss);

        $output = "$tmp_dir/$strand/tss/" . $gene_interval->name . '_' . $i . '.bed';
        open( $out, '>', $output )
            || die "Cannot open/read file " . $output . "!";

        print $out @{ $intersect_tss->stdout };
        close($out);

        # FOR TSS
        # Create slopBed object
        my $tts           = $self->_get_tts_genes($gene_interval);
        my $windows_tts   = $self->build_fixed_bins( $tts->stdout )->stdout;
        my $intersect_tts = $self->intersect_bins($windows_tts);

        $output = "$tmp_dir/$strand/tts/" . $gene_interval->name . '_' . $i . '.bed';
        open( $out, '>', $output )
            || die "Cannot open/read file " . $output . "!";

        print $out @{ $intersect_tts->stdout };
        close($out);
    }
=cut    


    method process_gene_paralell (Object $pm, Bio::Moose::Bed $gene_interval, Int $i, Str $tmp_dir, Str $strand) {
        srand($i);
        # FOR BODY
        # creating bins
        my $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $gene_interval,
            i => 'winnum',
            n => $self->body_resolution,
        );
        $windows->run;

        # if there is a temporary file, remove it
        if ( -e $windows->b ) {
            unlink $windows->b;
        }

        my $windows_body = $windows->stdout;
        my $intersect    = Bio::Moose::BedTools::Intersect->new(
            a => $windows_body,
            b => $self->reads,
            c => 1,
            sorted => 1,
        );

        $intersect->run;

        # Removing temporary bed file
        unlink $intersect->a if -e $intersect->a;

        my $intersect_body = $intersect;

        my $output = "$tmp_dir/$strand/body/" . $gene_interval->name . '_' . $i . '.bed';
        open( my $out, '>', $output )
            || die "Cannot open/read file " . $output . "!";
        print $out @{ $intersect_body->stdout };
        close($out);

        # FOR TSS
        # Create slopBed object
        my $flank = Bio::Moose::BedTools::Flank->new(
            i => $gene_interval,
            g => $self->genome,
            l => $self->tss,
            r => 0,
            s => 1,
        );

        # Run slopbed
        $flank->run();

        # if there is a temporary file, remove it
        if ( -e $flank->i ) {
            unlink $flank->i;
        }
        my $tss = $flank;

        $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $tss->stdout,
            w => $self->window_size,
            i => 'winnum'

                #n => $self->body_resolution,
        );

        $windows->run;

        # if there is a temporary file, remove it
        if ( -e $windows->b ) {
            unlink $windows->b;
        }
        my $windows_tss = $windows->stdout;
        $intersect = Bio::Moose::BedTools::Intersect->new(
            a => $windows_tss,
            b => $self->reads,
            c => 1,
            sorted => 1,
        );

        $intersect->run;

        # Removing temporary bed file
        unlink $intersect->a if -e $intersect->a;

        my $intersect_tss = $intersect;

        $output = "$tmp_dir/$strand/tss/" . $gene_interval->name . '_' . $i . '.bed';
        open( $out, '>', $output )
            || die "Cannot open/read file " . $output . "!";

        print $out @{ $intersect_tss->stdout };
        close($out);

        # FOR TTS
        # Create slopBed object
        $flank = Bio::Moose::BedTools::Flank->new(
            i => $gene_interval,
            g => $self->genome,
            l => 0,
            r => $self->tts,
            s => 1,
        );
        $self->log->trace( $flank->show_cmd_line );

        # Run slopbed
        $flank->run();

        # if there is a temporary file, remove it
        if ( -e $flank->i ) {
            unlink $flank->i;
        }

        my $tts = $flank;

        $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $tts->stdout,
            w => $self->window_size,
            i => 'winnum'

                #n => $self->body_resolution,
        );

        $windows->run;

        # if there is a temporary file, remove it
        if ( -e $windows->b ) {
            unlink $windows->b;
        }
        my $windows_tts = $windows->stdout;
        $intersect = Bio::Moose::BedTools::Intersect->new(
            a => $windows_tts,
            b => $self->reads,
            c => 1,
            sorted => 1,
        );

        $intersect->run;

        # Removing temporary bed file
        unlink $intersect->a if -e $intersect->a;

        my $intersect_tts = $intersect;

        $output = "$tmp_dir/$strand/tts/" . $gene_interval->name . '_' . $i . '.bed';
        open( $out, '>', $output )
            || die "Cannot open/read file " . $output . "!";

        print $out @{ $intersect_tts->stdout };
        close($out);
    }


    method merge_files ($tmp_dir) {
        open( my $out, '>', $self->output_file . '_long' )
            || die "Cannot open/write file " . $self->output_file . "_long !";
        open( my $out_smooth, '>', $self->output_file )
            || die "Cannot open/write file " . $self->output_file . "!";

        my $n_genes        = 0;
        my $n_genes_signal = 0;
        my @sum_genes;

        foreach my $strand (qw/positive negative/) {

            my $tss_dir  = "$tmp_dir/$strand/tss/";
            my $body_dir = "$tmp_dir/$strand/body";
            my $tts_dir  = "$tmp_dir/$strand/tts";
            my @dirs     = ( $tss_dir, $body_dir, $tts_dir );

            opendir( my $TSS, $tss_dir ) or die $!;
            my %genes;
            while ( my $filename = readdir($TSS) ) {
                next if $filename !~ /\.bed/;
                if ( -e "$body_dir/$filename" && -e "$tts_dir/$filename" ) {
                    my @aux_gene;
                    my $signal = 0;
                    foreach my $dir (@dirs) {

                        open( my $in, '<', $dir . '/' . $filename )
                            || die "Cannot open/read file " . $dir . '/' . $filename . "!";

                        my @aux_region;
                        while ( my $row = <$in> ) {
                            chomp $row;

                            # cols: chr start end bin_num count
                            my ( $chr, $start, $end, $bin, $count ) = split /\t/, $row;

                            # Normalize count by bin size and library size
                            my $size = $end - $start;
                            $size = 1 if $size < 1;

                            push @aux_region, ( ( $count / $size ) * $self->normFactor );

                            $signal++ if $count > 0;
                        }
                        close($in);
                        @aux_region = reverse(@aux_region) if $strand =~ /negative/;
                        push @aux_gene, @aux_region;
                    }

                    # Add tss,body,tts array of bins to gene hash
                    say $out join "\t", ( $filename, @aux_gene );
                    my $i = 0;
                    foreach my $value (@aux_gene) {
                        $sum_genes[$i] += $value;
                        $i++;
                    }

                    $n_genes++;
                    $n_genes_signal++ if $signal;
                }
                else {
                    $self->log->error("Cannot find all files for $filename");
                }
            }
        }
        close($out);

        # normalize by library size perl million and number of genes
        foreach my $value (@sum_genes) {
            $value = ( $value / $n_genes );
        }

        say $out_smooth join "\t", ( 'bin', ( 1 .. ( scalar @sum_genes ) ) );
        say $out_smooth join "\t", ( 'smooth', @sum_genes );

    }


    method run {
        # Params
        my $MAX_WORKERS = 1;
        my $p_n         = Sys::CPU::cpu_count();
        if ( $self->n_forks > $p_n ) {
            $MAX_WORKERS = $p_n;
        }
        else {
            $MAX_WORKERS = $self->n_forks;
        }

        $self->log->info( "Using: " . $MAX_WORKERS . " workers!" );

        # Get filtered genes
        my $filtered_genes = $self->_filter_input;

        # Creating simulation count array
        my @genes = ( 1 .. scalar @{$filtered_genes} );

        my $tmp = File::Temp->new();
        my $tmp_dir = $tmp->newdir( DIR => '/dev/shm', CLEANUP => 1 )->dirname;

        $self->log->info("Using $tmp_dir as temporary directory for all files");
        foreach (qw/body tss tts/) {
            make_path( $tmp_dir . '/positive/' . $_ ) || die "Cannot create $_";
            make_path( $tmp_dir . '/negative/' . $_ ) || die "Cannot create $_";
        }

        my $progress = Progress::Any->get_indicator(
            task   => "coverage",
            target => ~~ @genes
        );

        my $pm = AnyEvent::ForkManager->new(
            max_workers => $MAX_WORKERS,
            on_start    => sub {
                my ( $pm, $pid, $gene_interval, $i, $tmp_dir, $strand ) = @_;
                $progress->update( message => "gene: $i" );
            },
            on_finish => sub {
                my ( $pm, $pid, $status, $gene_interval, $i, $tmp_dir, $strand ) = @_;
            }
        );

        my $i = 1;
        foreach my $gene_interval ( @{$filtered_genes} ) {
            my $strand = 'positive';
            $strand = 'negative' if $gene_interval->strand eq '-';

            $pm->start(
                cb => sub { $self->process_gene_paralell(@_) },
                args => [ $gene_interval, $i, $tmp_dir, $strand ]
            );

            $i++;
        }

        my $cv = AnyEvent->condvar;

        # wait with non-blocking
        $pm->wait_all_children(
            cb => sub {
                my ($pm) = @_;
                say $tmp_dir;
                $self->merge_files($tmp_dir);
                rmtree($tmp_dir);
                $cv->send;
            }
        );

        $cv->recv;
    }

}
