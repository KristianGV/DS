#!/usr/bin/perl -w
#
# Script to add static libraries (or .o files) together to one new library.
# This is needed to make sure we have all the objects referring to each other
# in one file.
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: 2015-12-07

if (scalar @ARGV <=1) {
print <<END;
add_libs adds libraries (or .o files) together to one library.
Note 1: If the same .o file exists in several libraries the one that exists
in the first one in the list will be the one that goes into the new library
(the old libraries are extracted from the last one to the first one).
This is used for the concept of replaceable functions where we want
to replace .o files with user replaced versions of them.

Note 2: If the new library ends in .a a static library will be created.
If it ends in .so a shared library will be created (requires that the
object files in the libraries are compiled with the -fPIC option).    

Usage: There are two ways to call this script:
    
1) add_libs <new_library> <file-with-libraries, one per line>
2) add_libs <new_library> <old_library_1> <old_library_2> etc    
END
die "\n";  
}

$newlib=shift;
unlink("$newlib");
if ($newlib =~ /\.so/) {
    $static=0;
    print "Will create shared library $newlib from given libraries.\n";
} else {
    $static=1;
    print "Will create static library $newlib from given libraries.\n";
}

$tmpdir="tmp_lib";

if (-d "$tmpdir") {
    die "ERROR: Directory tmp_lib already existing. Exiting.\n";
} 

# Figure out if we have list of libraries and .o files in a file or given
# on the command line. We do this by checking if the next argument ends
# in .a or .o in which case we assume that we get libraries and .o filed
# on the command line

if (($ARGV[0] =~ /\.a/) || ($ARGV[0] =~ /\.o/)) { # on command line
    @libs=@ARGV;
} else { # read from file
    $liblistfile=shift;
    open(LIBLIST,$liblistfile) or die "ERROR: $liblistfile no such file. Exiting.\n";
    chomp (@libs = <LIBLIST>);
    close(LIBLIST);
}

mkdir("$tmpdir");
chdir("$tmpdir");

# Note: 'reverse' is crucial for replaceable function concept to work correctly
foreach $lib (reverse @libs) { 
    print "   $lib \n";
    if ($lib =~ /^\//) { # absolute path
        if ($lib =~ /\.o/) { # .o file
	    system("cp -p $lib ./");
	} else {
	    system("ar -x $lib");
	}
    } else { # relative path
        if ($lib =~ /\.o/) { # object file
	    system("cp -p ../$lib ./");
	} else {
	    system("ar -x ../$lib");
	}
    }
}

if ($static == 1) {
    if ($newlib =~ /^\//) { # absolute path
	system("ar r $newlib *.o");
    } else {
	system("ar r ../$newlib *.o");
    }
} else {
    if ($newlib =~ /^\//) { # absolute path
	system("gfortran -shared -fopenmp *.o -o $newlib");
    } else {
	system("gfortran -shared -fopenmp *.o -o ../$newlib");
    }
}
chdir("..");
system("rm -fr $tmpdir");
print "Done.\n";
