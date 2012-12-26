#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class Bio::Moose::BedPEIO {
    use Bio::Moose::BedPE;

    has 'file' => (
        is            => 'ro',
        isa           => 'Str',
        required      => 1,
        documentation => 'Bed file to be open',
    );

    has 'features' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::BedPE]',
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

    method _create_bedpe_object (Str $row, Str $track_row, Int $init_pos) {
        chomp $row;
        my @column = split /\s+/, $row;
        my $column_number = scalar @column;

        #Check minimum number of columns
        die "File " . $self->file . " has < then 6 columns"
            if ( $column_number < 6 );

        my $feat = Bio::Moose::BedPE->new(
            chrom1      => $column[0],
            chromStart1 => $column[1],
            chromEnd1   => $column[2],
            chrom2      => $column[3],
            chromStart2 => $column[4],
            chromEnd2   => $column[5],
            init_pos    => $init_pos,
        );

        $feat->track_line($track_row) if $track_row;

        my @attr = (qw(name mate1End strand1 strand2 editDist1 editDist2));

        my $i = 0;

        foreach my $value ( @column[ 6 .. $#column ] ) {
            if ( $attr[$i] ) {
                my $attribute = $attr[$i];
                $feat->$attribute($value);
            }
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

            if ( $row =~ /^[\#]track/ ) {
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
}
