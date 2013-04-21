#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
class Bio::Moose::BedTools {
    with 'MooseX::Role::Cmd';
    with 'Bio::Moose::Role::BEDFile';
    use MooseX::Attribute::Dependent;
    use Bio::Moose::BedIO;

    sub build_bin_name { 'bedtools'}

    has 'i' => (
        is            => 'rw',
        isa           => 'BEDFile',
        documentation => 'Input BED file',
        coerce        => 1,
        dependency    => None [qw/abam ubam a/],
    );

    #Implementing all options from bedtools
    has 'a' => (
        is            => 'rw',
        isa           => 'BEDFile',
        documentation => 'Input BED file A',
        coerce        => 1,
        dependency    => None[qw/abam ubam/],
    );

    has 'b' => (
        is  => 'rw',
        isa => 'BEDFile',
        coerce        => 1,
        documentation => 'Input BED file B',
    );

    has 'abam' => (
        is  => 'rw',
        isa => 'Str',
        documentation =>
            "The A input file is in BAM format.  Output will be BAM as well.",
        dependency => None[qw/a/],
    );

    has 'ubam' => (
        is            => 'rw',
        isa           => 'Str',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation =>
            'Write uncompressed BAM output. Default writes compressed BAM.',
        dependency => None[qw/a/],
    );

    has 'bed' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => "When using BAM input (-abam), write output as BED. The
        default is to write output in BAM when using -abam.",
    );

    has 'wa' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => "Write the original entry in A for each overlap.",
        dependency => None[qw/wb wao wo/],
    );

    has 'wb' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => "Write the original entry in B for each overlap.",
        dependency => None[qw/wa wao wo/],
    );

    has 'loj' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation =>
            "Perform a \"left outer join\". That is, for each feature
        in A report each overlap with B.  If no overlaps are found, report a
        NULL feature for B.",
    );

    has 'wo' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => "Write the original A and B entries plus the number of
        base pairs of overlap between the two features.",
        dependency => None[qw/wa wao wb/],
    );

    has 'wao' => (
        is            => 'rw',
        isa           => 'Bool',
        traits        => ['CmdOpt'],
        cmdopt_prefix => '-',
        documentation => "Write the original A and B entries plus the number of
        base pairs of overlap between the two features.",
        dependency => None[qw/wa wo wb/],
    );

    has 'g' =>
        ( is => 'rw', isa => 'Str', documentation => 'Genome file', );


    has 'show_cmd_line' =>
        ( is => 'rw', isa => 'Str', lazy => 1, builder => '_build_cmd_line', );
    
    method _build_cmd_line {
        my @args = $self->cmd_args;
        my $arg_line = join " ", @args;
        return 'command: "'.$self->bin_name .' '.$arg_line.'"';
    }
    
    
    method as_BedIO {
        return ( Bio::Moose::BedIO->new( file => $self->stdout ) );
    }


}
