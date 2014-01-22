use MooseX::Declare;
use Method::Signatures::Modifiers;
use feature qw(say);

class Biotools::Cluster {
    use MooseX::App::Command;            # important
    extends qw(Biotools);                   # purely optional
    use Bio::Moose::BedIO;
    use Math::CDF;
    use File::Basename;

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

        #$self->log->info( "genome size " . $genome_size );
        return $genome_size;
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
       foreach (keys %cluster){
            my @feat = @{$cluster{$_}->{features}};
            my $p = $cluster{$_}->{pvalue};
            my $hotspot_len = $cluster{$_}->{hotspot_len};
            my $chr = $cluster{$_}->{chr};
            my $start = $cluster{$_}->{start};
            my $end = $cluster{$_}->{end};

            say join "\t",
            ($chr,$start,$end,"hotspot".$i++,$hotspot_len,scalar
                @feat,$cluster{$_}->{n_left},$cluster{$_}->{n_right},$p);     
       }
    }
}
