#!/usr/bin/env perl
use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;
use Modern::Perl '2011';

class Bio::Moose::DependecyTypes {
    use MooseX::Attribute::Dependency;
    use List::MoreUtils ();

    MooseX::Attribute::Dependency::register(
        {
            name       => 'SmallerThan',
            message    => 'The value must be smaller than %s',
            constraint => sub {
                my ( $attr_name, $params, @related ) = @_;
                return List::MoreUtils::all {
                    $params->{$attr_name} < $params->{$_};
                }
                @related;
            },
        }
    );

    MooseX::Attribute::Dependency::register(
        {
            name       => 'BiggerThan',
            message    => 'The value must be Bigger than %s',
            constraint => sub {
                my ( $attr_name, $params, @related ) = @_;
                return List::MoreUtils::all {
                    $params->{$attr_name} > $params->{$_};
                }
                @related;
            },
        }
    );
}


