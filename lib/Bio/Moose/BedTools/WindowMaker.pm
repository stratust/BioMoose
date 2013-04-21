#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::BedTools::WindowMaker {
    extends 'Bio::Moose::BedTools';
    use MooseX::Attribute::Dependent;

    sub build_bin_name {'windowMaker'}

    has 's' => (
        is  => 'rw',
        isa => 'Num',
        documentation =>
            q[Step size: i.e., how many base pairs to step before creating a new window. Used to create "sliding" windows. 
            - Defaults to window size (non-sliding windows)]
    );

    has 'w' => (
        is            => 'rw',
        isa           => 'Num',
        documentation => q[Window size],
        dependency    => None [qw/n/]
    );

    has 'n' => (
        is            => 'rw',
        isa           => 'Num',
        documentation => q[number_of_windows>],
        dependency    => None [qw/w/]
    );

    has 'i' => ( is => 'rw', isa => 'Str', default => 'srcwinnum', );

    # Inherit without modifying
    has [qw/+g/];

    # Inherit modifying
    has '+b' => (
        required => 1,
        documentation =>
            q[BED file (with chrom,start,end fields). Windows will be created for each interval in the file.]
    );

}
