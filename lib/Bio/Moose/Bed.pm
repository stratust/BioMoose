#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::Bed {
    use Bio::Moose::DependencyTypes;
    use MooseX::Attribute::Dependent;

    # BED Fields
    has 'chrom' =>
        ( is => 'rw', isa => 'Str', required => 1, predicate => 'has_chrom' );

    has 'chromStart' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        predicate => 'has_chromStart',
        dependency => SmallerThan ['chromEnd']
    );
    has 'chromEnd' => (
        is         => 'rw',
        isa        => 'Int',
        required   => 1,
        predicate => 'has_chromEnd',
        dependency => BiggerThan ['chromStart']
    );
    has 'name'        => ( is => 'rw', isa => 'Str' );
    has 'score'       => ( is => 'rw', isa => 'Num' );
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
  
    has size => ( is => 'rw', isa => 'Int', lazy => 1, builder => 'build_size' );
    
    # store any information you want in this part
    has misc => ( is => 'rw', isa => 'Any' );


    method make_windows (Int :$number_of_windows ) {
        
    }
    
    method build_size {
        return $self->chromEnd - $self->chromStart;
    }

    method row {
        my @keys = qw/chrom chromStart chromEnd name score strand thickStart
            thickEnd itemRgb blockCount blockSizes blockStarts/;
        my @k_print;
        foreach my $k (@keys) {
            push @k_print, $k if ( defined $self->$k );
        }
        my $str = join "\t", @{$self}{@k_print};
        return $str."\n";

    }
   
    method remove_upstream (Int $size) {
        if ($self->strand eq '+') {
            $self->chromStart( $self->chromStart + $size);
        }
        elsif ($self->strand eq '-') {
            $self->chromEnd($self->chromEnd - $size);
        }
        else{
            die "Cannot remove upstream without strand information";
        }
    }
    
    method remove_downstream (Int $size) {
        if ($self->strand eq '+') {
            $self->chromEnd($self->chromEnd - $size);
        }
        elsif ($self->strand eq '-') {
            $self->chromStart($self->chromStart + $size);
        }
        else{
            die "Cannot remove downstream without strand information";
        }
    }

}

