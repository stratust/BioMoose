#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class Bio::Moose::BedIO {
    use Bio::Moose::Bed;

    has 'file' => (
        is            => 'ro',
        isa           => 'Str',
        required      => 1,
        documentation => 'Bed file to be open',
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
        chomp $row;
        my @column = split /\s+/, $row;
        my $column_number = scalar @column;

        #Check minimum number of columns
        die "File " . $self->file . " has < then 3 columns"
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

        open( my $in, '<', $self->file )
            || die "Cannot open/read file " . $self->file . "!";

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


}
