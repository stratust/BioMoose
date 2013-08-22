use MooseX::Declare;
use Method::Signatures::Modifiers;
 
role Bio::Moose::Role::TrackI {
   use Moose::Util::TypeConstraints;

   subtype 'TrackI'
        => as 'Bio::Moose::Track'
        => where { -e $_ };

    coerce 'TrackI' => from 'Str' => via { _set_input_string($_) };

    sub _set_input_string {
        my $value = shift;
        my @F = split /\s+/, $value;
        my $track = Bio::Moose::Track->new();

        foreach my $f (@F) {
            my ( $key, $value );
            if ( $f =~ /(.*)=["'](.*)["']/ ) {
                ( $key, $value ) = ( $1, $2 );
            }
            else{
                die "Cannot parse values in $f";
            }
            
            $track->$key($value);
        }
        return $track;
    }
}

