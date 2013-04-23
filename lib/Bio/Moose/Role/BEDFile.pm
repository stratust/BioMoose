#!/usr/bin/env perl

=head1 NAME 

    Modware::Some_module

=head1 SYNOPSIS

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

Thiago Yukio Kikuchi Oliveira E<lt>stratus@gmail.comE<gt>

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html

=head1 METHODS

=cut

use Moose;
use MooseX::Declare;
use Method::Signatures::Modifiers;

role Bio::Moose::Role::BEDFile {
    use Moose::Util::TypeConstraints;
    use File::Temp;
    use IO::ScalarArray;
    use IO::Scalar;
    
    subtype 'BEDFile'
        => as 'Str'
        => where { -e $_ };

    coerce 'BEDFile' 
        => from 'Str' 
        => via {_set_input_bed($_)};

    coerce 'BEDFile' 
        => from 'ArrayRef[Bio::Moose::Bed]' 
        => via {_set_input_bed_arrayobj($_)};

    coerce 'BEDFile' 
        => from 'ArrayRef' 
        => via {_set_input_bed($_)};

    # Check if is a filename or a arrey ref
    our $tmp;    # let temp file exist after loosing sub subscope
    our $tmpi=0;#indexer of object (if overwritten, object loose scope too)

    sub _set_input_bed {
        my $value = shift;
        $tmpi++;
        $tmp->{$tmpi} = File::Temp->new( UNLINK => 1, SUFFIX => '.bed' );
        my $out = $tmp->{$tmpi};
        print $out @{$value};
        return $out->filename;
    }

    sub _set_input_bed_arrayobj {
        my $value = shift;
        $tmpi++;
        $tmp->{$tmpi} = File::Temp->new( UNLINK => 1, SUFFIX => '.bed' );
        my $out = $tmp->{$tmpi};
        foreach my $r (@{$value}) {
            print $out $r->row;
        }
        return $out->filename;
    }

    subtype 'BEDFileI' => as 'FileHandle';

    coerce 'BEDFileI' 
        => from 'ArrayRef[Bio::Moose::Bed]' 
        => via {_return_scalar_io($_) };

    coerce 'BEDFileI' 
        => from 'ArrayRef' 
        => via {_set_input_bedI_arrayref($_)};

    sub _set_input_bedI_arrayref {
        my $arrayref = shift;
        my $in       = IO::ScalarArray->new($arrayref);
        return ($in);
    }

    coerce 'BEDFileI' 
        => from 'Str' 
        => via { _set_input_bedI_str($_) };

    sub _set_input_bedI_str {
        my $string = shift;
        if ( -e $string ) {
            open( my $in, '<', $string )
                || die "Cannot open/read file " . $string . "!";
            return ($in);            
        }
        else {
            my $in = IO::Scalar->new($string);
            return ($in);
        }
    }
 
    sub _return_scalar_io {
        my $value = shift;
        my $string;
        foreach my $f (@{$value}) {
            $string .= $f->row;
        }
        my $in = IO::Scalar->new(\$string);
        return ($in);
    }

    subtype 'BEDFileO' => as 'FileHandle';
    coerce 'BEDFileO' => from 'Str' => via { _set_input_bedO_str($_) };

    sub _set_input_bedO_str {
        my $string = shift;
        open( my $out, '>', $string )
            || die "Cannot open/read file " . $string . "!";
        return ($out);
    }

}

