#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class Bio::Moose::BedIO {
    with 'Bio::Moose::Role::BEDFile';
    with 'Bio::Moose::Role::TrackI';
    use Bio::Moose::Bed;
    use IO::ScalarArray;
    use Class::MOP;
    use POSIX;
    use Storable 'dclone';

    has 'file' => (
        is            => 'ro',
        isa           => 'BEDFileI',
        required      => 1,
        coerce        => 1,
        documentation => 'Bed file to be open',
    );

    has 'outfile' => (
        is            => 'rw',
        isa           => 'BEDFileO',
        coerce        => 1,
        documentation => 'Bed file to be written',
    );

    has 'track' => (
        is            => 'ro',
        isa           => 'TrackI',
        required      => 0,
        coerce        => 1,
        documentation => 'Track Object',
    );

    has 'features' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        traits        => ['Array'],
        lazy          => 1,
        builder       => '_build_features',
        documentation => 'ArrayRef of features',
        handles       => {
            all_features   => 'elements',
            add_feature    => 'push',
            next_feature   => 'shift',
            map_features   => 'map',
            count_features => 'count',
        },
    );

    has 'features_sorted' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::Bed]',
        traits        => ['Array'],
        lazy          => 1,
        builder       => '_build_features_sorted',
        documentation => 'ArrayRef of sorted features by chrom, chromStar and chromEnd',
#        handles       => {
            #all_features   => 'elements',
            #add_feature    => 'push',
            #next_feature   => 'shift',
            #map_features   => 'map',
            #count_features => 'count',
        #},
    );

    method _create_bed_object (Str $row, Str $track_row, Int $init_pos) {
        my @column = split /\s+/, $row;
        my $column_number = scalar @column;

        #Check minimum number of columns
        die "File " . $self->file . " has < than 3 columns"
            if ( $column_number < 3 );

        my $feat = Bio::Moose::Bed->new(
            chrom      => $column[0],
            chromStart => $column[1],
            chromEnd   => $column[2],
            init_pos   => $init_pos,
        );
        $feat->track_line($track_row) if $track_row;

        my @attr = (
            qw(
                name score strand thickStart thickEnd
                itemRgb blockCount blockSizes blockStarts
                )
        );

        my $i = 0;

        foreach my $value ( @column[ 3 .. $#column ] ) {
            my $attribute = $attr[$i];
            $feat->$attribute($value);
            $i++;
        }

        return $feat;
    }

    method _build_features {
        my @objects;
        my $track_row = 0;
        my $init_pos  = 1;

        my $in = $self->file;
        while ( my $row = <$in> ) {
            chomp $row;

            if ( $row =~ /^track/ ) {
                $track_row = $row;
                next;
            }

            push( @objects,
                $self->_create_bed_object( $row, $track_row, $init_pos ) );

            $init_pos++;
        }
        close($in);
        return \@objects;
    }

    method _build_features_sorted {
        my @sorted_features = sort {
                   $a->chrom cmp $b->chrom
                || $a->chromStart <=> $b->chromStart
                || $a->chromEnd <=> $b->chromEnd

        } @{ $self->features };
        return \@sorted_features;
    }
    
    method write (Bio::Moose::Bed $bed_object) {
        if ( $self->outfile ) {
            my $out = $self->outfile;
            print $out $bed_object->row;
        }
        else {
            print $bed_object->row;
        }
    }

    method write_all_features {
        
        foreach my $f ( @{ $self->features } ) {
            $self->write($f);
        }
    }

    # return BEDIO object with center of intervals
    method get_middle_intervals {
        my $middleIO = dclone $self->features;
        foreach my $f (@{$middleIO}) {
            my $middle_start =  floor( ( ($f->chromEnd + $f->chromStart ) ) / 2 );
            $f->chromStart($middle_start );
            $f->chromEnd( $middle_start + 1);
        }
        return $middleIO;
    }
}
