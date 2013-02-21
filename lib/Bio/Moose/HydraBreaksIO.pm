#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class Bio::Moose::HydraBreaksIO {
    use Bio::Moose::HydraBreaks;

    has 'file' => (
        is            => 'ro',
        isa           => 'Str',
        required      => 1,
        documentation => 'File to be open',
    );
    
    has 'features' => (
        is            => 'rw',
        isa           => 'ArrayRef[Bio::Moose::HydraBreaks]',
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

    method _create_break_object (Str $row, Str $track_row, Int $init_pos) {
        chomp $row;
        my @column = split /\s+/, $row;
        my $column_number = scalar @column;

        #Check minimum number of columns
        die "File " . $self->file . " has < then 10 columns"
            if ( $column_number < 10 );

        my $feat = Bio::Moose::HydraBreaks->new(
            chrom1           => $column[0],
            chromStart1      => $column[1],
            chromEnd1        => $column[2],
            chrom2           => $column[3],
            chromStart2      => $column[4],
            chromEnd2        => $column[5],
            breakpointId     => $column[6],
            numDistinctPairs => $column[7],
            strand1          => $column[8],
            strand2          => $column[9],
            init_pos         => $init_pos,
        );

        $feat->track_line($track_row) if $track_row;

        my @attr = (
            qw(
                meanEditDist1
                meanEditDist2
                meanMappings1
                meanMappings2
                breakpointSize
                numMappings
                allWeightedSupport
                finalSupport
                finalWeightedSupport
                numUniquePairs
                numAnchoredPairs
                numMultiplyMappedPairs
                )
        );

        my $i = 0;

        foreach my $value ( @column[ 10 .. $#column ] ) {
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
                $self->_create_break_object( $row, $track_row, $init_pos ) );

            $init_pos++;
        }
        close($in);
        return \@objects;
    }

    method summary_breaks {
        while ( my $feat = $self->next_feature ){
            
        }
    }
}
