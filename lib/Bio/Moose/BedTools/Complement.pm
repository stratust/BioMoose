#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::BedTools::Complement {
    extends 'Bio::Moose::BedTools';

    sub build_bin_name { 'complementBed'}
    
    # Inherit without modifying
    has [qw/+i +g/];

}
