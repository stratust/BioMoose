use MooseX::Declare;
use Method::Signatures::Modifiers;
use feature qw(say);

class Biotools::CumulativeDensity {
    use MooseX::App::Command;            # important
    extends qw(Biotools);                   # purely optional
    use Bio::Moose::BedIO;
    use Bio::Moose::BedTools::Intersect;
    use Bio::Moose::BedTools::Complement;
    use Bio::Moose::BedTools::Slop;
    use Bio::Moose::BedTools::Flank;
    use Bio::Moose::BedTools::WindowMaker;
    use Moose::Util::TypeConstraints;
    use Data::Dumper;
    with 'Custom::Log';

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
    
    has '_reads' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        lazy      => 1,
        default => sub {
            my ($self) = @_;
            my $feats = Bio::Moose::BedIO->new(file=>$self->reads)->features;
            my @aux;
            foreach my $f (@{$feats}) {
                push @aux, $f unless ( $f->chrom =~ /chrM/ || $f->chrom =~ /random/ );
            }
            return \@aux;
        },
        documentation => 'Read object',
    );
    
    has 'normFactor' => (
        is      => 'rw',
        isa     => 'Num',
        lazy    => 1,
        default => sub {
            my ($self) = @_;
            my $normfactor = (1e6 / scalar @{ $self->_reads });
            return $normfactor;
        },
        documentation => 'Normalize by this factor',
    );

    has 'n_genes' => (
        is            => 'rw',
        isa           => 'Int',
        lazy => 1,
        default => sub{
            my ($self) = @_;
            my $n_gene=0;
            if ( $self->only_genes_with_reads ) {
                my $slop = Bio::Moose::BedTools::Slop->new(
                    i => $self->_filter_input,
                    g => $self->genome,
                    l => $self->tss,
                    r => $self->tts,
                    s => 1,
                );

                # Run slopbed
                $slop->run();

                my $intersect = Bio::Moose::BedTools::Intersect->new(
                    a => $slop->as_BedIO->features,
                    b => $self->_reads,
                    u => 1,
                );
                $intersect->run;

                $n_gene = scalar @{ $intersect->as_BedIO->features };
            }
            else {
                $n_gene = scalar @{ $self->_filter_input };
            }
            return $n_gene;
        }
    );

    has 'n_intergenic_regions' => (
        is            => 'rw',
        isa           => 'Int',
        lazy => 1,
        default => sub{
            my ($self) = @_;
            my $n_region = 0;
            if ( $self->only_genes_with_reads ) {
                my $intersect = Bio::Moose::BedTools::Intersect->new(
                    a => $self->intergenic_bed,
                    b => $self->_reads,
                    u => 1,
                );
                $intersect->run;
                $n_region = scalar @{ $intersect->as_BedIO->features };
            }
            else {
                $n_region=scalar @{ $self->intergenic_bed };                
            }
            return $n_region;
        }
    );

    
    has 'intergenic_bed' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        lazy      => 1,
        builder => '_get_complement_bed',
    );

    method build_filter_input {
        $self->log_info( "Filtering " . $self->input );
        $self->log_info( "Increasing TSS and TTS  and removing overlap");
        # get slopped genes
        #
        my $slopped_genes = $self->_get_slopped_input->as_BedIO->features;
        $self->log_info(
            "   - Number of genes before filter: " . scalar @{$slopped_genes} );
        
       
        # intersect genes and remove genes with more than 1 intersection
        my $gene_i_gene = Bio::Moose::BedTools::Intersect->new(
            a => $slopped_genes,
            b => $slopped_genes,
            c => 1
        );
        $self->log_info( $gene_i_gene->show_cmd_line );
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

        $self->log_info( "   - Number of genes without overlap: "
                . scalar @genes_no_overlap );

        my $file;

        if ($self->remove_overlapping_genes){
            $self->log_info( "      Removing overlapping genes ( --remove_overlapping_genes = 1 )" );
            $file  = \@genes_no_overlap;
        }else{
            $self->log_info( "      Keeping overlapping genes ( --remove_overlapping_genes = 0 )" );
            $file = $orig_input;
        }

        my $in    = Bio::Moose::BedIO->new( file => $file );
        my $feats = $in->features;

        $self->log_info(
            "   - Number of genes before filter: " . scalar @{$feats} );

        my @aux;
        my $min_gene_size = ( $self->tss + $self->tts ) * 0.5 ;

        $self->log_info( "   - Removing genes smaller than " . $min_gene_size );
        $self->log_info("   - Removing chrM and chr*_random");

        foreach my $f ( @{$feats} ) {
            unless ( $f->chrom =~ /chrM/ || $f->chrom =~ /random/ ) {
                if ($f->size >= $min_gene_size){
                   #fix gne size
                   push @aux,$f;
                }
            }
        }

        $self->log_info( "   - Number of genes after filter: " . scalar @aux );
        return \@aux;
    }

    method _get_slopped_input {
        # Create slopBed object
        $self->log_info("   Slopping ");
        $self->log_info("   -TSS size: ".$self->tss);
        $self->log_info("   -TTS size: ".$self->tts);
        my $slop = Bio::Moose::BedTools::Slop->new(
            i => $self->input,
            g => $self->genome,
            l => $self->tss,
            r => $self->tts,
            s => 1,
        );
        $self->log_info($slop->show_cmd_line);
        # Run slopbed
        $slop->run();
        return $slop;
    }
 
    method _get_tss_genes (ArrayRef[Bio::Moose::Bed] $genes) {
                            # Create slopBed object
        $self->log_info("   Getting TSS ");
        $self->log_info( "   -TSS size: " . $self->tss );
        my $flank = Bio::Moose::BedTools::Flank->new(
            i => $genes,
            g => $self->genome,
            l => $self->tss,
            r => 0,
            s => 1,
        );
        $self->log_info( $flank->show_cmd_line );

        # Run slopbed
        $flank->run();

        return $flank;
    }

    method _get_tts_genes (ArrayRef[Bio::Moose::Bed] $genes) {
        # Create slopBed object
        $self->log_info("   Getting TTS ");
        $self->log_info( "   -TTS size: " . $self->tts );
        my $flank = Bio::Moose::BedTools::Flank->new(
            i => $genes,
            g => $self->genome,
            l => 0,
            r => $self->tts,
            s => 1,
        );
        $self->log_info( $flank->show_cmd_line );

        # Run slopbed
        $flank->run();
        return $flank;
    }

    method _get_complement_bed {
        # Prepare
        $self->log_info("   Complement");
        $self->log_info("   -TSS: ".$self->tss);
        $self->log_info("   -TTS: ".$self->tts);
        my $complement = Bio::Moose::BedTools::Complement->new(
            i => $self->_get_slopped_input->as_BedIO->features,
            g => $self->genome,
        );
        $self->log_info($complement->show_cmd_line);
        # Run complement
        $complement->run();

        my $feats = $complement->as_BedIO->features;
        my @aux;

        $self->log_info("Filtering Intergenic regions (complement of slopped input");
        my $min_region_size = ( $self->tss + $self->tts ) * 0.5 ;
        $self->log_info("   Removing regions in chrM and chr*_random");
        $self->log_info("   Removing regions smaller than $min_region_size");
        foreach my $f ( @{$feats} ) {
            unless ( $f->chrom =~ /chrM/ || $f->chrom =~ /random/ ) {
                if ($f->size >= $min_region_size){
                   push @aux,$f;
                }
            }
        }
        return \@aux; 
    }
        
    method build_body_bins($this_input) {
        $self->log_info("   Divide genes body in bins");
        $self->log_info("   - body resolution : ".$self->body_resolution.' bins');
        my $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $this_input,
            i => 'winnum',
            n => $self->body_resolution,
        );
        $self->log_info($windows->show_cmd_line);
        $windows->run;
        return $windows;
    }

    method build_fixed ($this_input) {
        $self->log_info("   Divide in fixed bins");
        $self->log_info("   -window size: ".$self->window_size);
        my $windows = Bio::Moose::BedTools::WindowMaker->new(
            b => $this_input,
            w => $self->window_size,
            i => 'winnum'
            #n => $self->body_resolution,
        );
        $self->log_info($windows->show_cmd_line);
        $windows->run;
        return $windows;
    }

    method get_intergenicD (
        ArrayRef[Bio::Moose::Bed] $bed_windows,
        ) {
        
        $self->log_info("Calculating Density...");
        my %bodyD;
        my $bodyStart = $self->relativeCoordSize +
        ($self->tts * 1.5) - $self->bodyStep;

        foreach my $f ( @{$bed_windows} ) {

            # Index by relative position
            my $key = ( $f->name * $self->bodyStep ) + $bodyStart;

            # Normalize reads by bin size and add to relatie bin
            $bodyD{$key}{rpb} += ($f->score / $f->size);
            $bodyD{$key}{index} = $f->name;
        }

        $self->log_info("Smoothing...");
        
        my $n_intergenic = $self->n_intergenic_regions;
        # Normalizing by gene number * factor
        foreach my $k ( keys %bodyD ) {
            $bodyD{$k}{smooth} = ($bodyD{$k}{rpb}
                / $n_intergenic)  * $self->normFactor;
        }
        #say "$_ => $bodyD{$_}{smooth} ($bodyD{$_}{index})"
        #    for ( sort { $a <=> $b } keys %bodyD );

        return \%bodyD;
    }

    method get_bodyD (
        ArrayRef[Bio::Moose::Bed] $bed_windows_sense,
        ArrayRef[Bio::Moose::Bed] $bed_windows_antisense
        ) {

        # Reverse bins in antisense array
        foreach my $f ( @{ $bed_windows_antisense } ) {
            $f->name($self->body_resolution - $f->name + 1);
        }
        # Combine sense and antisense
        my @bed_windows = (@{$bed_windows_sense},@{$bed_windows_antisense
            });
        
        $self->log_info("Calculating Density...");
        my %bodyD;
        my $bodyStart = $self->bodyStep / 2;

        foreach my $f ( @bed_windows ) {

            # Index by relative position
            my $key = ( $f->name * $self->bodyStep ) - $bodyStart;

            # Normalize reads by bin size and add to relatie bin
            $bodyD{$key}{rpb} += ($f->score / $f->size);
            $bodyD{$key}{index} = $f->name;

            # kee Bed object for earch bin
            #push @{$bodyD{$key}{bed}},$f;
        }

        $self->log_info("Smoothing...");

        # Normalizing by gene number * factor
        
        foreach my $k ( keys %bodyD ) {
            $bodyD{$k}{smooth} = ($bodyD{$k}{rpb}
                / $self->n_genes) * $self->normFactor ;
        }
        #say "$_ => $bodyD{$_}{smooth} ($bodyD{$_}{index})"
        #    for ( sort { $a <=> $b } keys %bodyD );

        return \%bodyD;
    }

    method get_fixedD (
     ArrayRef[Bio::Moose::Bed] $bed_windows_sense, 
     ArrayRef[Bio::Moose::Bed] $bed_windows_antisense, 
     Str $region) {

        my $nbin = 0;
        if ( $region =~ /tss/i ) {
            $nbin = int( $self->tss / $self->window_size );
        }
        elsif ( $region =~ /tts/i ) {

            $nbin = int( $self->tts / $self->window_size );
        }

        # Reverse bins in antisense array
        foreach my $f ( @{ $bed_windows_antisense } ) {
            $f->name($nbin - $f->name + 1);
        }
        # Combine sense and antisense
        my @bed_windows = (@{$bed_windows_sense},@{$bed_windows_antisense
            });

        my %fixedD;
        $self->log_info("Calculating Density...");
        foreach my $f ( @bed_windows ) {
            my $key;

            # Index by relative position (name holds bin number)
            if ( $region =~ /tts/i ) {
                $key
                    = ( $f->name * $self->window_size )
                    - ( $self->window_size / 2 )
                    + $self->relativeCoordSize;
            }
            elsif ( $region =~ /tss/i ) {
                $key = ( -1 * abs( $self->tss ) )
                    + ( ( $f->name * $self->window_size ) );
            }

            # Normalize reads by bin size and add to relatie bin
            $fixedD{$key}{rpb} += ($f->score / $f->size);
            $fixedD{$key}{index} = $f->name;

            # kee Bed object for earch bin
            #push @{$fixedD{$key}{bed}},$f;
        }

        $self->log_info("Smoothing...");

        # Normalizing by gene number * factor
        foreach my $k ( keys %fixedD ) {
            $fixedD{$k}{smooth} = ($fixedD{$k}{rpb}
                / $self->n_genes)  * $self->normFactor ;
        }
        #say "$_ => $fixedD{$_}{smooth} ($fixedD{$_}{index})"
        #    for ( sort { $a <=> $b } keys %fixedD );

        return \%fixedD;
    }

    method intersect_genes {
        # Body Density
        # =====================================================================
        $self->log_info("Calculating gene body bins");

        # FOr positivve genes
        $self->log_info(" Positive Genes");      
        my $gene_body_sense_bins = $self->build_body_bins($self->sense_genes);

        $self->log_info( "Intersecting Positive gene body bins with " . $self->reads );
        my $body_sense_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $gene_body_sense_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );
        $self->log_info($body_sense_intersected->show_cmd_line);        
        $body_sense_intersected->run;

        # For negative genes
        $self->log_info(" Positive Genes");      
        my $gene_body_antisense_bins = $self->build_body_bins($self->antisense_genes);

        $self->log_info( "Intersecting Negative gene body bins with " . $self->reads );
        my $body_antisense_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $gene_body_antisense_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );
        $self->log_info($body_antisense_intersected->show_cmd_line);        
        $body_antisense_intersected->run;

        # calcualte body density
        my $bodyD = $self->get_bodyD( 
            $body_sense_intersected->as_BedIO->features,
            $body_antisense_intersected->as_BedIO->features 
        );
 
        # TSS Density
        # =====================================================================
        $self->log_info("Calculating TSS bins");
        
        # For positive
        $self->log_info(
            "Intersecting Positive TSS bins with " . $self->reads );
        my $tss_sense_genes = $self->_get_tss_genes( $self->sense_genes );
        my $gene_tss_sense_bins
            = $self->build_fixed( $tss_sense_genes->as_BedIO->features );

        # Get TTS intersecton
        my $tss_sense_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $gene_tss_sense_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );

        $self->log_info( $tss_sense_intersected->show_cmd_line );
        $tss_sense_intersected->run;
        
        # For negative
        $self->log_info(
            "Intersecting Negative TSS bins with " . $self->reads );
        my $tss_antisense_genes = $self->_get_tss_genes( $self->antisense_genes );
        my $gene_tss_antisense_bins
            = $self->build_fixed( $tss_antisense_genes->as_BedIO->features );


        # Get TSS intersecton
        my $tss_antisense_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $gene_tss_antisense_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );

        $self->log_info( $tss_antisense_intersected->show_cmd_line );
        $tss_antisense_intersected->run;

        my $tssD
            = $self->get_fixedD( $tss_sense_intersected->as_BedIO->features,
            $tss_antisense_intersected->as_BedIO->features, 'TSS' );

        # TTS Density
        # =====================================================================
         $self->log_info("Calculating TTS bins");
        
        # For positive
        $self->log_info(
            "Intersecting Positive TTS bins with " . $self->reads );
        my $tts_sense_genes = $self->_get_tts_genes( $self->sense_genes );
        my $gene_tts_sense_bins
            = $self->build_fixed( $tts_sense_genes->as_BedIO->features );

        # Get TTS intersecton
        my $tts_sense_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $gene_tts_sense_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );

        $self->log_info( $tts_sense_intersected->show_cmd_line );
        $tts_sense_intersected->run;
        
        # For negative
        $self->log_info(
            "Intersecting Negative TTS bins with " . $self->reads );
        my $tts_antisense_genes = $self->_get_tts_genes( $self->antisense_genes );
        my $gene_tts_antisense_bins
            = $self->build_fixed( $tts_antisense_genes->as_BedIO->features );

        # Get TTS intersecton
        my $tts_antisense_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $gene_tts_antisense_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );

        $self->log_info( $tts_antisense_intersected->show_cmd_line );
        $tts_antisense_intersected->run;

        my $ttsD
            = $self->get_fixedD( $tts_sense_intersected->as_BedIO->features,
            $tts_antisense_intersected->as_BedIO->features, 'TTS' );
     
        # create hash with density
        my %geneD = (%{$tssD}, %{$bodyD}, %{$ttsD});
        return \%geneD;
    }
    
    method intersect_intergenic {
        # Body Density
        # =====================================================================
        $self->log_info("Calculating Intergenic bins");

        # FOr positivve genes
        my $intergenic_bins =
        $self->build_body_bins($self->intergenic_bed);

        $self->log_info( "Intersecting intergenic bins with " . $self->reads );
        my $intergenic_intersected = Bio::Moose::BedTools::Intersect->new(
            a => $intergenic_bins->as_BedIO->features,
            b => $self->reads,
            c => 1
        );
        $self->log_info($intergenic_intersected->show_cmd_line);        
        $intergenic_intersected->run;

        # calcualte body density
        my $intergenicD = $self->get_intergenicD( 
            $intergenic_intersected->as_BedIO->features,
        );

        return $intergenicD;
    }
 
    method run {
        my $geneD_file = $self->output_file . '.gene';
        open( my $out, '>', $geneD_file )
            || die "Cannot open/write file " . $geneD_file . "!";
        
        my $geneD = $self->intersect_genes;
        foreach my $pos (sort {$a <=> $b} keys %{ $geneD }) {
            say $out join "\t", ($pos, $geneD->{$pos}->{smooth});
        }
        close( $out );

        my $intergenicD_file = $self->output_file . '.intergenic';
        open( $out, '>', $intergenicD_file )
            || die "Cannot open/write file " . $intergenicD_file . "!";
        
        my $intergenicD = $self->intersect_intergenic;
        foreach my $pos (sort {$a <=> $b} keys %{ $intergenicD }) {
            say $out join "\t", ($pos, $intergenicD->{$pos}->{smooth});
        }
        close( $out );

    }
}
