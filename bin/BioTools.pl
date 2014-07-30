#!/usr/bin/env perl
use Moose;
use feature qw(say);
use MooseX::Declare;
use Method::Signatures::Modifiers;


# Store the log_path
our $log_path;
our $logfile_path;

sub logfile {
    return $logfile_path;
}

# Log Role
role Custom::Log {
    use Log::Log4perl qw(:easy);
    with 'MooseX::Log::Log4perl::Easy';
 
    use Cwd 'abs_path';
    use File::Basename;
    use File::Path;

    # Configuring log 
    BEGIN {
        my $logconf_file    = 'log4perl.conf';
        my $log_conf_path   = '';
        my $full_path       = abs_path($0);
        my $script_path     = dirname($full_path);
        my $current_path    = &Cwd::cwd();
        my $script_filename = basename($full_path);
        my $script_name     = $script_filename;

        # Removing extension
        $script_name =~ s/\.\S+$//;
        
        $log_path = &Cwd::cwd().'/logs/';
        unless (-e $log_path){
            mkpath($log_path);
        }
        $logfile_path = $log_path . $script_name . '.log';
        my $logtracefile_path = $log_path . $script_name . '_trace.log';


        # Verifify conf path
        if ( -d $current_path . '/conf' ) {
            $log_conf_path = $current_path.'/conf/';
        }
        elsif ( -d $current_path . '/../conf' ) {
            $log_conf_path = $current_path. '/../conf/';
        }

        # Name of the custom file: "script_name"_log4perl.conf
        my $personal_logconf_file = $script_name . '_log4perl.conf';

        if ( -e $current_path . $personal_logconf_file ) {
            $logconf_file = $personal_logconf_file;
        }

        $log_conf_path .= $logconf_file
          if ( -e $log_conf_path . $logconf_file );
       
        if ($log_conf_path){
            Log::Log4perl->init($log_conf_path);
        }
        else {
            Log::Log4perl->init(
                \qq{

                log4perl.rootLogger = TRACE, LOGFILE, Screen, AppTrace

                # Filter to match level ERROR
                log4perl.filter.MatchError = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchError.LevelToMatch  = ERROR
                log4perl.filter.MatchError.AcceptOnMatch = true
 
                # Filter to match level DEBUG
                log4perl.filter.MatchDebug = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchDebug.LevelToMatch  = DEBUG
                log4perl.filter.MatchDebug.AcceptOnMatch = true
 
                # Filter to match level WARN
                log4perl.filter.MatchWarn  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchWarn.LevelToMatch  = WARN
                log4perl.filter.MatchWarn.AcceptOnMatch = true
 
                # Filter to match level INFO
                log4perl.filter.MatchInfo  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchInfo.LevelToMatch  = INFO
                log4perl.filter.MatchInfo.AcceptOnMatch = true
 
                # Filter to match level TRACE
                log4perl.filter.MatchTrace  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.MatchTrace.LevelToMatch  = TRACE
                log4perl.filter.MatchTrace.AcceptOnMatch = true
 
                # Filter to match level TRACE
                log4perl.filter.NoTrace  = Log::Log4perl::Filter::LevelMatch
                log4perl.filter.NoTrace.LevelToMatch  = TRACE
                log4perl.filter.NoTrace.AcceptOnMatch = false


                log4perl.appender.LOGFILE=Log::Log4perl::Appender::File
                log4perl.appender.LOGFILE.filename= $logfile_path
                log4perl.appender.LOGFILE.mode=append
                log4perl.appender.LOGFILE.layout=Log::Log4perl::Layout::PatternLayout
                log4perl.appender.LOGFILE.layout.ConversionPattern=%d %p> %F{1}:%L %M%n%m%n%n
                log4perl.appender.LOGFILE.Filter = NoTrace

                # Error appender
                log4perl.appender.AppError = Log::Log4perl::Appender::File
                log4perl.appender.AppError.filename = $logfile_path
                log4perl.appender.AppError.layout   = SimpleLayout
                log4perl.appender.AppError.Filter   = MatchError
 
                # Warning appender
                log4perl.appender.AppWarn = Log::Log4perl::Appender::File
                log4perl.appender.AppWarn.filename = $logfile_path
                log4perl.appender.AppWarn.layout   = SimpleLayout
                log4perl.appender.AppWarn.Filter   = MatchWarn

                # Debug  appender
                log4perl.appender.AppDebug = Log::Log4perl::Appender::File
                log4perl.appender.AppDebug.filename = $logfile_path
                log4perl.appender.AppDebug.layout   = SimpleLayout
                log4perl.appender.AppDebug.Filter   = MatchDebug

                # Trace  appender
                log4perl.appender.AppTrace = Log::Log4perl::Appender::File
                log4perl.appender.AppTrace.filename = $logtracefile_path
                log4perl.appender.AppTrace.layout   = SimpleLayout
                log4perl.appender.AppTrace.Filter   = MatchTrace

                # Screen Appender (Info only)
                log4perl.appender.Screen = Log::Log4perl::Appender::ScreenColoredLevels
                log4perl.appender.Screen.stderr = 0
                log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
                log4perl.appender.Screen.layout.ConversionPattern = %d %m %n
                log4perl.appender.Screen.Filter = MatchInfo


            });
        }
    }
}

class MyApp is dirty {
    use MooseX::App qw(Color);
}

class MyApp::Bedpe2Bed12 {
    use MooseX::App::Command;            # important
    extends qw(MyApp);                   # purely optional
    use Bio::Moose::HydraBreaksIO;

    option 'input_file' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(i)],
        required      => 1,
        documentation => q[Bedpe File],
    );  

    option 'genome' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(g)],
        required      => 1,
        documentation => q[Chromosome size file],
    );  

    option 'dist' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(d)],
        required      => 1,
        default      => 1000000,
        documentation => q[Maxium distance between pairs before split in two  entries],
    );  

    option 'name' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(n)],
        required      => 1,
        documentation => q[Track Name],
    );  

    option 'description' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(e)],
        required      => 1,
        documentation => q[Track Description],
    );  

    method parse_genome {
        open( my $in, '<', $self->genome ) 
            || die "Cannot open/read file " . $self->genome . "!";
        
        my %hash;
        while ( my $row = <$in> ){
            chomp $row;
            my ($chr,$size);
            ($chr,$size) = split /\s+/,$row;
            $hash{$chr} = $size;
        }
        
        close( $in );
        return \%hash;        
    }

    method run {
        my $in = Bio::Moose::HydraBreaksIO->new( file => $self->input_file );
        #say $in->count_features;
        #say $in->all_features;
        say "track name='".$self->name."' description='".$self->description."' itemRgb='On'";
        while (my $feat = $in->next_feature){
            say $feat->write_bed12($self->dist,$self->parse_genome);
        }
    }   
}

class MyApp::Cluster {
    use MooseX::App::Command;            # important
    extends qw(MyApp);                   # purely optional
    use Bio::Moose::BedIO;
    use Math::CDF;
    use File::Basename;
    use Digest::SHA qw(sha256_hex);
    with 'Custom::Log';

    option 'input_file' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases     => [qw(i)],
        required      => 1,
        documentation => q[Bedpe File],
    );  

    option 'genome' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases     => [qw(g)],
        required      => 0,
        documentation => q[Chromosome size file],
    );  

    option 'human' => (
        is            => 'ro',
        isa           => 'Bool',
        required      => 0,
        documentation => q[Use hardcoded human genome size],
    );  

    option 'mouse' => (
        is            => 'ro',
        isa           => 'Bool',
        required      => 0,
        documentation => q[Use hardcoded mouse genome size],
    );  

    option 'cutoff' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases     => [qw(t)],
        required      => 1,
        default      => 0.00000001,
        documentation => q[Cutoff],
    );  

    option 'cutoff_pair' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases     => [qw(p)],
        required      => 1,
        default      => 0.01,
        documentation => q[Cutoff for pairs],
    );  

    option 'minority' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases     => [qw(m)],
        required      => 1,
        default      => 0.1,
        documentation => q[Minimum of both primers to be accepted as a hotspot],
    );  

    option 'min_cluster_number' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases     => [qw(n)],
        required      => 1,
        default      => 3,
        documentation => q[Minimum of translocations to form a "hotspot" cluster],
    );

    has 'transloc_bed' => (
        is            => 'ro',
        isa           => 'Bio::Moose::BedIO',
        lazy          => 1,
        builder       => '_builder_transloc_bed',
        documentation => 'Hold parsed BED file',
    );

    has 'genome_size' => (
        is            => 'ro',
        isa           => 'Int',
        lazy          => 1,
        builder       => '_builder_genome_size',
        documentation => 'Keep the genome size',
    );
    
    option 'index_file' => (
        is            => 'ro',
        isa           => 'Str',
        cmd_aliases   => [qw(f)],
        lazy          => 1,
        builder        => '_builder_index_file',
        documentation => q[Index filename],
    );
  

   method _builder_transloc_bed {
        my $in = Bio::Moose::BedIO->new( file => $self->input_file );
        return $in;
   } 

   method _builder_genome_size {
        my $genome_size;
        # Sum chromosome size
        if ($self->genome){
        	$genome_size += $self->parse_genome->{$_} foreach keys %{$self->parse_genome};
	}
	else{

        	$genome_size=2861343702 if $self->human;
        	$genome_size=2123929214 if $self->mouse;
	}
        die "Select genome size" unless $genome_size;

        #$self->log_info( "genome size " . $genome_size );
        return $genome_size;
   } 

   method _builder_index_file {
       my $filename = $self->input_file;

       if ($filename =~ /\.bed$/){
            $filename =~ s/\.bed/\.hotspots_index/g;
       }
       else{
            $filename.=".hotspots_index";    
       }
       return $filename;
   } 

    method parse_genome {
        open( my $in, '<', $self->genome )
            || die "Cannot open/read file " . $self->genome . "!";
        my %hash;
        while ( my $row = <$in> ) {
            chomp $row;
            my ( $chr, $size );
            ( $chr, $size ) = split /\s+/, $row;
            next if $size !~ /\d+/;
            $hash{$chr} = $size;
        }
        close($in);
        return \%hash;
    }

    method probability($dist,$n) {
        # Probability of success is number of translocations/size of genome
        my $p = $self->transloc_bed->count_features/$self->genome_size;
        
        #  computes the negative binomial cdf at each of the values in $dist using
        #  the corresponding number of successes, $n and probability of success
        #  in a single trial, $p
        my $prob = &Math::CDF::pnbinom($dist,$n,$p);

        return $prob;
    }

    method get_clusters {
        # Hold current and last bed feature
        my ($this, $last);
        # Keep the number of features in each chrom
        my $i = 1;
        # Keep number of clusters created
        my $j = 1;
        # Hash of Hash of Array to keep all cluster of bed features
        my %cluster;

        foreach my $feat ( @{ $self->transloc_bed->features_sorted } ) {
            # receive current object
            $this = $feat;

            if ($last) {
                if ( $this->chrom eq $last->chrom ) {
                    # Calculate distance between current and last feature
                    my $dist = $this->chromStart - $last->chromStart;

                    # Calculate the probability of the distance be within the
                    # expected by a random uniform distribution
                    my $p = $self->probability( $dist, 1 );
 
                    if ( $p < $self->cutoff_pair && $cluster{$j} ) {
                        push @{ $cluster{$j}->{features} }, $this;
                        push @{ $cluster{$j}->{pvalue} }, $p;
                    }
                    elsif ( $p < $self->cutoff_pair && !$cluster{$j} ) {
                        push @{ $cluster{$j}->{features} }, ( $last, $this );
                        push @{ $cluster{$j}->{pvalue} }, $p;
                    }
                    elsif ( $p >= $self->cutoff_pair && $cluster{$j} ) {
                        #say $_->chromStart for @{$cluster{$j}->{features}};
                        #say join " ", @{$cluster{$j}->{pvalue}};
                        $j++;
                    }
                }
                else {
                    # Reset chromosome count
                    $i = 1;
                    if ( $cluster{$j} ) {
                        $j++;
                    }
                }
            }
            # current object becomes the last
            $last = $this;
        }
        return \%cluster;
    }

    method get_filtered_clusters {
        # Keep filtered clusters
        my %filtered_cluster;

        # Get all clusters
        my $c = $self->get_clusters;

        foreach my $key ( sort { $a <=> $b } keys %{$c} ) {

            # Filter by minimun of features in a cluster (default: 3)
            if ( scalar @{ $c->{$key}->{features} } >= $self->min_cluster_number ) {

               # Calculate P-value for each hotspot
               # ($hotspot_len - $n_total) is the number of failures in the negative
               # binomial; n_total is the number of success
               my @cluster = @{ $c->{$key}->{features} };
               my $hotspot_len
                    = $cluster[$#cluster]->chromEnd - $cluster[0]->chromStart;
               my $n_total = scalar @cluster;
               my $diff = ( $hotspot_len - $n_total ) + 1;
               my $p = $self->probability( $diff, $n_total );
                
               #$p = 1 unless defined($p);
                
                # If $hotspot_len - $n_total is negative $p is null
                if ($diff < 0){
                    $p = 0;
                }
                #$p = 1 if $hotspot_len < 0;
                die "Hotspot_lengh is negative" if $hotspot_len < 0;
 
                # Counting left and right
                my ( $n_left, $n_right ) = ( 0, 0 );
                
                foreach my $feat (@cluster) {
                    $n_left++  if $feat->name =~ /left/;
                    $n_right++ if $feat->name =~ /right/;
                }
                
                my $has_minority = 0;
                if ( $self->minority ) {
                    if ( $n_left == 0 && $n_right == 0 ) {
                        die "Houston, we have a problem! No left and no right primer found!";
                    }

                    if (   ( $n_left / $n_total >= $self->minority )
                        && ( $n_right / $n_total >= $self->minority ) )
                    {
                        $has_minority = 1;
                    }

                }
                else {
                    $has_minority = 1;
                }

                if ( $p < $self->cutoff  && $has_minority ) {
                    $filtered_cluster{$key}->{features} = \@cluster;
                    $filtered_cluster{$key}->{pvalue} = $p;
                    $filtered_cluster{$key}->{hotspot_len} = $hotspot_len;
                    $filtered_cluster{$key}->{chr} = $cluster[0]->chrom;
                    $filtered_cluster{$key}->{start} = $cluster[0]->chromStart;
                    $filtered_cluster{$key}->{end} = $cluster[$#cluster]->chromEnd;
                    $filtered_cluster{$key}->{n_left} = $n_left;
                    $filtered_cluster{$key}->{n_right} = $n_right;
                }
            }
        }
        return \%filtered_cluster;
    }

    method run {
       my %cluster = %{$self->get_filtered_clusters};
       #say "Total clusters: ",scalar keys %cluster;
       my $i=1;
       open( my $out_index, '>', $self->index_file ) 
           || die "Cannot open/read file " . $self->index_file . "!";
      
       foreach (sort {$a <=> $b } keys %cluster){
            my @feat = @{$cluster{$_}->{features}};
            my $n_feat = scalar @feat;
            my $p = $cluster{$_}->{pvalue};
            my $hotspot_len = $cluster{$_}->{hotspot_len};
            my $chr = $cluster{$_}->{chr};
            my $start = $cluster{$_}->{start};
            my $end = $cluster{$_}->{end};

            my $ht_id = "hotspot".$i++;

            say join "\t",
            ($chr,$start,$end,$ht_id,$hotspot_len,$n_feat,$cluster{$_}->{n_left},$cluster{$_}->{n_right},$p);

            my @shear_names;
            foreach my $f (@feat) {
                push @shear_names, $f->name;
            }
            
            my $sorted_shear_names_string = join( ",", sort { $a cmp $b } @shear_names );
            my $hotspot_sha256_id = sha256_hex($sorted_shear_names_string);

            say $out_index join( "\t", ( $ht_id, $hotspot_sha256_id ,$n_feat, $sorted_shear_names_string));

       }
       close( $out_index );
 
    }
}

class MyApp::Cumulative_Density {
    use MooseX::App::Command;            # important
    extends qw(MyApp);                   # purely optional
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

class Main {
    import MyApp;
    MyApp->new_with_command->run();
}


=head1 NAME 

    MyApp

=head1 SYNOPSIS
  This application requires Perl 5.10.0 or higher   
  This application requires, at least, the following modules to work:
    - Moose
    - MooseX::App::Command

  Here, you want to concisely show a couple of SIMPLE use cases.  You should describe what you are doing and then write code that will run if pasted into a script.  

  For example:

  USE CASE: PRINT A LIST OF PRIMARY IDS OF RELATED FEATURES

    my $gene = new Modware::Gene( -feature_no => 4161 );

    foreach $feature ( @{ $gene->features() } ) {
       print $feature->primery_id()."\n";
    }

=head1 DESCRIPTION

   Here, AT A MINIMUM, you explain why the object exists and where it might be used.  Ideally you would be very detailed here. There is no limit on what you can write here.  Obviously, lesser used 'utility' objects will not be heavily documented.

   For example: 

   This object attempts to group together all log_information about a gene
   Most of this log_information is returned as references to arrays of other objects.  For example
   the features array is one such association.  You would use this whenever you want to read or write any 
   properties of a gene.


=head1 AUTHOR

Thiago Yukio Kikuchi Oliveira E<lt>stratust@gmail.comE<gt>

Copyright (c) 2012 Rockefeller University - Nussenzweig's Lab

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html

=head1 METHODS

=cut

