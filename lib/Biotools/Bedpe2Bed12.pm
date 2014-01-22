use MooseX::Declare;
use Method::Signatures::Modifiers;
use feature qw(say);

class Biotools::Bedpe2Bed12 {
    use MooseX::App::Command;            # important
    extends qw(Biotools);                   # purely optional
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
