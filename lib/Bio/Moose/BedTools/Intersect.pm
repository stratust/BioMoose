#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class Bio::Moose::BedTools::Intersect {
    extends 'Bio::Moose::BedTools';

    sub build_bin_name { 'intersectBed'}

    has [qw/+abam +ubam +bed +loj +a +wa +wb +wo +wao/];
    has '+b' => ( required => 1); 
    has 'c' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => 'Add a count column',
    );
    
}
