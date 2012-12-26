#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class Bio::Moose::Bed {
    use Bio::Moose::DependencyTypes;
    use MooseX::Attribute::Dependent;

    # BED Fields
    has 'chrom' => ( is => 'rw', isa => 'Str', required => 1 );
    has 'chromStart' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        dependency => SmallerThan['chromEnd']
    );
    has 'chromEnd' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        dependency => BiggerThan ['chromStart']
    );
    has 'name'        => ( is => 'rw', isa => 'Str' );
    has 'score'       => ( is => 'rw', isa => 'Int' );
    has 'strand'      => ( is => 'rw', isa => 'Str' );
    has 'thickStart'  => ( is => 'rw', isa => 'Int' );
    has 'thickEnd'    => ( is => 'rw', isa => 'Int' );
    has 'itemRgb'     => ( is => 'rw', isa => 'Str' );
    has 'blockCount'  => ( is => 'rw', isa => 'Str' );
    has 'blockSizes'  => ( is => 'rw', isa => 'Str' );
    has 'blockStarts' => ( is => 'rw', isa => 'Str' );
    has 'track_line'  => ( is => 'rw', isa => 'Str' );

    # Attributes used to get gene names 
    has 'genome' =>  ( is => 'ro', isa => 'Str', required => 0 );
    has 'table_name' => ( is => 'ro', isa => 'Str', );
    has 'init_pos' => ( is => 'ro', isa => 'Int', );
   
    method make_windows (Int :$number_of_windows ) {
        
    }
    

}

