#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl;

class MyApp {
    use MooseX::App qw(Color);
}

class MyApp::Bedpe2Bed12 {
    use MooseX::App::Command;            # important
    extends qw(MyApp);                   # purely optional
    use Bio::Moose::HydraBreaksIO;
    use Data::Dumper;

    option 'input_file' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(i)],
        required      => 1,
        documentation => q[Bedpe File],
    );  

    option 'genome' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(g)],
        required      => 1,
        documentation => q[Chromosome size file],
    );  

    option 'dist' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(d)],
        required      => 1,
        default      => 1000000,
        documentation => q[Maxium distance between pairs before split in two  entries],
    );  

    option 'name' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(n)],
        required      => 1,
        documentation => q[Track Name],
    );  

    option 'description' => (
        is            => 'rw',
        isa           => 'Str',
        cmd_aliases     => [qw(e)],
        required      => 1,
        documentation => q[Track Description],
    );  

    method parse_genome {
        open( my $in, '<', $self->genome ) 
            || die "Cannot open/read file " . $self->genome . "!";
        
        my %hash;
        while ( my $row = <$in> ){
            chomp $row;
            my ($chr,$size);
            ($chr,$size) = split /\s+/,$row;
            $hash{$chr} = $size;
        }
        
        close( $in );
        return \%hash;        
    }

    method run {
        my $in = Bio::Moose::HydraBreaksIO->new( file => $self->input_file );
        #say $in->count_features;
        #say $in->all_features;
        say "track name='".$self->name."' description='".$self->description."' itemRgb='On'";
        while (my $feat = $in->next_feature){
            say $feat->write_bed12($self->dist,$self->parse_genome);
        }
    }   
}

class Main {
    import MyApp;
    MyApp->new_with_command->run();
}


=head1 NAME 

    MyApp

=head1 SYNOPSIS
  This application requires Perl 5.10.0 or higher   
  This application requires, at least, the following modules to work:
    - Moose
    - MooseX::App::Command

  Here, you want to concisely show a couple of SIMPLE use cases.  You should describe what you are doing and then write code that will run if pasted into a script.  

  For example:

  USE CASE: PRINT A LIST OF PRIMARY IDS OF RELATED FEATURES

    my $gene = new Modware::Gene( -feature_no => 4161 );

    foreach $feature ( @{ $gene->features() } ) {
       print $feature->primery_id()."\n";
    }

=head1 DESCRIPTION

   Here, AT A MINIMUM, you explain why the object exists and where it might be used.  Ideally you would be very detailed here. There is no limit on what you can write here.  Obviously, lesser used 'utility' objects will not be heavily documented.

   For example: 

   This object attempts to group together all information about a gene
   Most of this information is returned as references to arrays of other objects.  For example
   the features array is one such association.  You would use this whenever you want to read or write any 
   properties of a gene.


=head1 AUTHOR

Thiago Yukio Kikuchi Oliveira E<lt>stratust@gmail.comE<gt>

Copyright (c) 2012 Rockefeller University - Nussenzweig's Lab

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html

=head1 METHODS

=cut

