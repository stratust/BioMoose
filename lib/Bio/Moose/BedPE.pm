#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

=cut
Field Nbr.  Field Name  Description
1.   chrom1  Chromosome for end 1
2.   start1  Start position for end 1
3.   end1    End position for end 1
4.   chrom2  Chromosome for end 2
5.   start2  Start position for end 2
6.   end1    End position for end 2
7.   name    Name of the pair
8.   mate1End    To which mate of the pair do fields 1,2,3,9,11 correspond? (values: 1 or 2)
9.   strand1     Orientation for end 1 (+ or -)
10.  strand2     Orientation for end 2 (+ or -)
11.  editDist1   Alignment edit distance for end 1 (can be extracted from NM tag in SAM)
12.  editDist2   Alignment edit distance for end 2 (can be extracted from NM tag in SAM)
=cut

class Bio::Moose::BedPE {
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
    has 'name'       => ( is => 'rw', isa => 'Str' );
    has 'mate1End'   => ( is => 'rw', isa => 'Int' );
    has 'strand1'    => ( is => 'rw', isa => 'Str' );
    has 'strand2'    => ( is => 'rw', isa => 'Str' );
    has 'editDist1'  => ( is => 'rw', isa => 'Str' );
    has 'editDist2'  => ( is => 'rw', isa => 'Str' );
    has 'track_line' => ( is => 'rw', isa => 'Str' );

    # Attributes used to get gene names 
    has 'genome' =>  ( is => 'ro', isa => 'Str', required => 0 );
    has 'table_name' => ( is => 'ro', isa => 'Str', );
    has 'init_pos' => ( is => 'ro', isa => 'Int', );

}

