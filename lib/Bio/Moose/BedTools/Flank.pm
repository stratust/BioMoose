#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::BedTools::Flank {
    extends 'Bio::Moose::BedTools';

    sub build_bin_name { 'flankBed'}
    
    has 's' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => q[Define -l and -r based on strand.E.g. if used, -l
        500 for a negative-stranded feature, it will add 500 bp downstream.
        Default = false.]
    );

    has 'l' => (
        is            => 'rw',
        isa           => 'Num',
        documentation => q[The number of base pairs to subtract from the start
        coordinate.
        - (Integer) or (Float, e.g. 0.1) if used with -pct.]
    );

    has 'r' => (
        is  => 'rw',
        isa => 'Num',
        documentation =>
            q[The number of base pairs to add to the end coordinate.
            - (Integer) or (Float, e.g. 0.1) if used with -pct.]
    );

    has 'pct' => (
        is            => 'rw',
        isa           => 'Num',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => q[Define -l and -r as a fraction of the feature's
        length. E.g. if used on a 1000bp feature, -l 0.50, will add 500 bp "upstream".Default =false.]
    );

    has 'header' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => q[Print the header from the input file prior to
        results.]
    );
   
    # Inherit without modifying
    has [qw/+i +g/];

    # Inherit modifying
    has '+b' => ( isa => 'Num');

}
