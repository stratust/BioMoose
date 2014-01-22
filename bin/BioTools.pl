#!/usr/bin/env perl
use MooseX::Declare;

class Main {
    use Biotools;
    Biotools->new_with_command->run();
}
