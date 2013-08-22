use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;

class Bio::Moose::Track {
    has 'name' => (
        is            => 'rw',
        isa           => 'Str',
        default => 'User Track',
    );

    has 'description' => (
        is            => 'rw',
        isa           => 'Str',
        default => '"User Supplied Track"' 
    );
    
    has 'type' => (
        is            => 'rw',
        isa           => 'Str',
        required      => 1,
        default      => 'bedDetail',
        documentation => 'BAM, bedDetail, bedGraph, bigBed, bigWig, broadPeak, narrowPeak, Microarray, VCF, wig',
    );

    has 'visibility' => (
        is            => 'rw',
        isa           => 'Int',
        required      => 1,
        default       => 1,
        documentation => '0:hide, 1:dense, 2:full , 3:pack , 4:squish',
    );

    has 'color' => (
        is            => 'rw',
        isa           => 'Str',
        traits        => ['Getopt'],
        cmd_aliases   => '',
        required      => 1,
        documentation => 'describe',
    );
   
    has 'itemRgb' => (
        is            => 'rw',
        isa           => 'Str',
        default       => 'On',
        documentation => 'describe',
    );
}
