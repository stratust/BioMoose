use MooseX::Declare;
use Method::Signatures::Modifiers;
use feature qw(say);
BEGIN { $ENV{LOG_LEVEL} = 'info' unless $ENV{LOG_LEVEL}  }

class Biotools is dirty {
    use MooseX::App qw(Color);
    use Log::Any::App '$log',
        -screen => { pattern_style => 'script_long' },
        -file   => { path          => 'logs/', level => 'debug' };

    has 'log' => (
        is            => 'ro',
        isa           => 'Object',
        required      => 1,
        default       => sub { return $log },
        documentation => 'Keep Log::Any::App object',
    );
}
