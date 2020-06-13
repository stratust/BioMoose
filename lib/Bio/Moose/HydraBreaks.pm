#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

=cut
Hydra breakpoint output file format (.final and .all)
Field Nbr.  Field Name  Description
1.   chrom1  Chromosome for end 1 of the breakpoint.
2.   start1  Start position for end 1 of the breakpoint.
3.   end1    End position for end 1 of the breakpoint.
4.   chrom2  Chromosome for end 2 of the breakpoint.
5.   start2  Start position for end 2 of the breakpoint.
6.   end2    End position for end 2 of the breakpoint.
7.   breakpointId    Unique Hydra breakpoint identifier.
8.   numDistinctPairs    Number of distinct pairs in breakpoint.
9.   strand1     Orientation for the first end of the breakpoint.
10.  strand2     Orientation for the second end of the breakpoint.
11.  meanEditDist1   Mean edit distance observed in end1 of the breakpoint pairs.
12.  meanEditDist2   Mean edit distance observed in end2 of the breakpoint pairs.
13.  meanMappings1   Mean number of mappings for end1 of all pairs in the breakpoint.
14.  meanMappings2   Mean number of mappings for end2 of all pairs in the breakpoint.
15.  breakpointSize  Size of the breakpoint.
16.  numMappings     Total number of mappings included in the breakpoint.
17.  allWeightedSupport  Amount of weighted support from the mappings in the breakpoint.
18.  finalSupport    Amount of final support from the mappings in the breakpoint.
19.  finalWeightedSupport    Amount of final weighted support from the mappings in the breakpoint.
20.  numUniquePairs  Number of pairs in the breakpoint that were uniquely mapped to the genome.
21.  numAnchoredPairs    Number of pairs in the breakpoint that were mapped to the genome in an "anchored" fashion (i.e. 1xN).
22.  numMultiplyMappedPairs  Number of pairs in the breakpoint that were multiply mapped to the genome in fashion (i.e. NxN).

Hydra breakpoint output detail file format
Field Nbr.  Field Name  Description
1.   chrom1  Chromosome for end 1
2.   start1  Start position for end 1
3.   end1    End position for end 1
4.   chrom2  Chromosome for end 1
5.   start2  Start position for end 2
6.   end1    End position for end 2
7.   name    Name of the pair
8.   mate1End    To which mate of the pair do fields 1,2,3,9,11 correspond? (values: 1 or 2)
9.   strand1     Orientation for end 1 (+ or -)
10.  strand2     Orientation for end 2 (+ or -)
11.  editDist1   Alignment edit distance for end 1 (can be extracted from NM tag in SAM)
12.  editDist2   Alignment edit distance for end 2 (can be extracted from NM tag in SAM)
13.  numMappings1    Number of mappings for end 1 of this pair
14.  numMappings2    Number of mappings for end 2 of this pair
15.  mappingType     What type of mapping is this? (1=unique, 2=anchored, 3-multiply)
16.  includedInBreakpoint    Was this pair ultimately included in this breakpoint?
17.  breakpointId    Unique Hydra breakpoint identifier.
=cut

class Bio::Moose::HydraBreaks {
    use Bio::Moose::DependencyTypes;
    use MooseX::Attribute::Dependent;

    # BEDPE Fields
    has 'chrom1' => ( is => 'rw', isa => 'Str', required => 1 );
    has 'chromStart1' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        dependency => SmallerThan ['chromEnd1']
    );
    has 'chromEnd1' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        dependency => BiggerThan ['chromStart1']
    );
    has 'chrom2' => ( is => 'rw', isa => 'Str', required => 1 );
    has 'chromStart2' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        dependency => SmallerThan ['chromEnd2']
    );
    has 'chromEnd2' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        dependency => BiggerThan ['chromStart2']
    );
    has 'breakpointId'     => ( is => 'rw', isa => 'Str', required => 1 );
    has 'numDistinctPairs' => ( is => 'rw', isa => 'Int', required => 1 );
    has 'strand1'          => ( is => 'rw', isa => 'Str', required => 1 );
    has 'strand2'          => ( is => 'rw', isa => 'Str', required => 1 );
    has 'meanEditDist1'    => ( is => 'rw', isa => 'Num' );
    has 'meanEditDist2'    => ( is => 'rw', isa => 'Num' );
    has 'meanMappings1'    => ( is => 'rw', isa => 'Num' );
    has 'meanMappings2'    => ( is => 'rw', isa => 'Num' );
    has 'breakpointSize'   => ( is => 'rw', isa => 'Num' );
    has 'numMappings'      => ( is => 'rw', isa => 'Int' );
    has 'allWeightedSupport'     => ( is => 'rw', isa => 'Num' );
    has 'finalSupport'           => ( is => 'rw', isa => 'Num' );
    has 'finalWeightedSupport'   => ( is => 'rw', isa => 'Num' );
    has 'numUniquePairs'         => ( is => 'rw', isa => 'Int' );
    has 'numAnchoredPairs'       => ( is => 'rw', isa => 'Int' );
    has 'numMultiplyMappedPairs' => ( is => 'rw', isa => 'Int' );
    has 'track_line'             => ( is => 'rw', isa => 'Str' );

    # Attributes used to get gene names
    has 'genome' => ( is => 'ro', isa => 'Str', required => 0 );
    has 'table_name' => ( is => 'ro', isa => 'Str', );
    has 'init_pos'   => ( is => 'ro', isa => 'Int', );

    method get_space_end ($chr_size,$chr_end) {
        my $space = 500;
        if ( $chr_end + $space > $chr_size ) {
            $space = $chr_size - $chr_end;
        }
        return $space;
    }

    method get_space_start ($chr_start) {
        my $space = 500;
        if ( $chr_start - $space < 0 ) {
            $space = 0;
        }
        return $space;
    }

    method write_bed12 ($dist=10000,$hash='') {
        my $color;
        my $feature_name;
        my @rows;

        if ( $self->strand1 eq '+' && $self->strand2 eq '-' ) {
            $color = "153,0,0";    # deletion breakpoints are red
        }

        #elsif ( $self->strand1 eq '-' && $self->strand2 eq '+' ) {
        #    $color = "0,102,0";    # duplication breakpoints are green
        #}
        elsif ( $self->strand1 eq '-' && $self->strand2 eq '+' ) {
            $color = "153,0,0";    # duplication breakpoints are red
        }
        elsif ( $self->strand1 eq '+' && $self->strand2 eq '+' ) {
            $color = "0,51,204";    # inversion breakpoints are blue
        }
        elsif ( $self->strand1 eq '-' && $self->strand2 eq '-' ) {
            $color = "0,51,204";    # inversion breakpoints are blue
        }

        $feature_name
            = $self->breakpointId . ","
            . $self->numDistinctPairs . ":"
            . $self->strand1 . "/"
            . $self->strand2
            . ":intra:"
            . $self->chrom1 . ":"
            . $self->chromStart1 . "-"
            . $self->chromEnd1 . ","
            . $self->chrom2 . ":"
            . $self->chromStart2 . "-"
            . $self->chromEnd2;

        # space to avoid go over chromosome end
        my $space1_pos = $self->get_space_end( $hash->{ $self->chrom1 },
            $self->chromEnd1 );
        my $space2_pos = $self->get_space_end( $hash->{ $self->chrom2 },
            $self->chromEnd2 );

        my $space1_neg = $self->get_space_start( $self->chromStart1 );
        my $space2_neg = $self->get_space_start( $self->chromStart2 );


        # intrachromosomals
        if (   ( $self->chrom1 eq $self->chrom2 )
            && ( abs( $self->chromEnd2 - $self->chromStart1 ) <= $dist ) )
        {
            push @rows, join "\t",
                (
                $self->chrom1,
                $self->chromStart1,
                $self->chromEnd2,
                $feature_name,
                abs( $self->chromEnd2 - $self->chromStart1 ),
                "+",
                $self->chromStart1,
                $self->chromEnd2,
                $color,
                2,
                ( $self->chromEnd1 - $self->chromStart1 ) . ","
                    . ( $self->chromEnd2 - $self->chromStart2 ),
                "0," . ( $self->chromStart2 - $self->chromStart1 )
                );
        }

        # intrachromosomals that exceed dist
        elsif (( $self->chrom1 eq $self->chrom2 )
            && ( abs( $self->chromEnd2 - $self->chromStart1 ) > $dist ) )
        {
            if ( $self->strand1 eq "+" ) {
                push @rows, join "\t",
                    (
                    $self->chrom1,
                    $self->chromStart1,
                    ( $self->chromEnd1 + $space1_pos ),
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "+",
                    $self->chromStart1,
                    ( $self->chromEnd1 + $space1_pos ),
                    $color, 2,
                    ( $self->chromEnd1 - $self->chromStart1 ) . ",1",
                    "0,"
                        . (
                        $self->chromEnd1 - $self->chromStart1 + $space1_pos - 1
                        )
                    );
            }
            if ( $self->strand1 eq "-" ) {
                push @rows, join "\t",
                    (
                    $self->chrom1,
                    ( $self->chromStart1 - $space1_neg ),
                    $self->chromEnd1,
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "-",
                    ( $self->chromStart1 - $space1_neg ),
                    $self->chromEnd1,
                    $color,
                    2,
                    "1," . ( $self->chromEnd1 - $self->chromStart1 ),
                    "0,$space1_neg"
                    );
            }
            if ( $self->strand2 eq "+" ) {
                push @rows, join "\t",
                    (
                    $self->chrom2,
                    $self->chromStart2,
                    ( $self->chromEnd2 + $space2_pos ),
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "+",
                    $self->chromStart2,
                    ( $self->chromEnd2 + $space2_pos ),
                    $color, 2,
                    ( $self->chromEnd2 - $self->chromStart2 ) . ",1",
                    "0,"
                        . (
                        $self->chromEnd2 - $self->chromStart2 + $space2_pos - 1
                        )
                    );
            }
            if ( $self->strand2 eq "-" ) {
                push @rows, join "\t",
                    (
                    $self->chrom2,
                    ( $self->chromStart2 - $space2_neg ),
                    $self->chromEnd2,
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "-",
                    ( $self->chromStart2 - $space2_neg ),
                    $self->chromEnd2,
                    $color,
                    2,
                    "1," . ( $self->chromEnd2 - $self->chromStart2 ),
                    "0,$space2_neg"
                    );
            }
        }

        # interchromosomals:
        elsif ( $self->chrom1 ne $self->chrom2 ) {

            $feature_name =~ s/intra/inter/g;

            if ( $self->strand1 eq "+" ) {
                push @rows, join "\t",
                    (
                    $self->chrom1,
                    $self->chromStart1,
                    ( $self->chromEnd1 + $space1_pos ),
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "+",
                    $self->chromStart1,
                    ( $self->chromEnd1 + $space1_pos ),
                    $color, 2,
                    ( $self->chromEnd1 - $self->chromStart1 ) . ",1",
                    "0,"
                        . (
                        $self->chromEnd1 - $self->chromStart1 + $space1_pos - 1
                        )
                    );
            }
            if ( $self->strand1 eq "-" ) {
                push @rows, join "\t",
                    (
                    $self->chrom1,
                    ( $self->chromStart1 - $space1_neg ),
                    $self->chromEnd1,
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "-",
                    ( $self->chromStart1 - $space1_neg ),
                    $self->chromEnd1,
                    $color,
                    2,
                    "1," . ( $self->chromEnd1 - $self->chromStart1 ),
                    "0,$space1_neg"
                    );
            }
            if ( $self->strand2 eq "+" ) {
                push @rows, join "\t",
                    (
                    $self->chrom2,
                    $self->chromStart2,
                    ( $self->chromEnd2 + $space2_pos ),
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "+",
                    $self->chromStart2,
                    ( $self->chromEnd2 + $space2_pos ),
                    $color, 2,
                    ( $self->chromEnd2 - $self->chromStart2 ) . ",1",
                    "0,"
                        . (
                        $self->chromEnd2 - $self->chromStart2 + $space2_pos - 1
                        )
                    );
            }
            if ( $self->strand2 eq "-" ) {
                push @rows, join "\t",
                    (
                    $self->chrom2,
                    ( $self->chromStart2 - $space2_neg ),
                    $self->chromEnd2,
                    $feature_name,
                    abs( $self->chromEnd2 - $self->chromStart1 ),
                    "-",
                    ( $self->chromStart2 - $space2_neg ),
                    $self->chromEnd2,
                    $color,
                    2,
                    "1," . ( $self->chromEnd2 - $self->chromStart2 ),
                    "0,$space2_neg"
                    );
            }

        }

        return join "\n", @rows;
    }


}

