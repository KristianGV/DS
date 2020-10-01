#!/usr/bin/perl -w
#
# Configure script to configure DarkSUSY.
# NOTE. Run this script before compiling DarkSUSY. This is done automatically
# either at the configure stage or at the make stage.
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: April 12, 2001
# Modified: 
#   20080107 Paolo Gondolo split lines in dsdir.h
#   20080226 Joakim Edsjo, rewrite to use parameters instead
#   20130814 Paolo Gondolo, split dsdir.h and dsver.h


$dsroot = shift;
$dsroot =~ s#/$##;   # take away final / if any
$dsversion = $dsroot;
$dsversion =~ s#^.*/##;

$dsinstall=shift;
$dsinstall =~ s#/$##;   # take away final / if any
$dsinstall .= "/"; # put it back so that we know we have it.


# Add revision number to dsversion
# The logic is to extract revision info from svn or git
# It is then added to the folder name as the version tag.
# However, if a file 'README.RELEASE' exists, the info in the first line of
# that file will be used instead. 
# (the rest of that file contains a list of main updates wrt the previous version)
$rev="";
if (open(IN,"svn info|")) { # if we are on subversion
    while(defined($line=<IN>)) {
	if ($line =~ /^Revision:\s+(\d+)$/) {
	    $rev=" (rev $1)";
	}
    }
    close(IN);
} else {
    #print "Couldn't invoke svn. No revision number added to version tag.\n";
}
if (open(IN,"git status|")) { # check if we are on git and can get info
    while(defined($line=<IN>)) {
	if ($line =~ /^On branch (.+)$/) {
	    $branch=$1;
	    $rev=" (branch $branch)";
	}
    }
    close(IN);
    $found=0;
    if (open(IN,"git log|")) { # check if we are on git and can get info
	while(defined($line=<IN>)) {
	    if ($line =~ /^Date:\s+(\S.+)$/) {
		$gitdate=$1;
		$rev=" (branch $branch from $gitdate)";
		$found=1;
	    }
	    last if ($found==1);
	}
	close(IN);
    }
} else {
    print "Couldn't invoke git. No revision info added to version tag.\n";
}
$dsversion .= $rev;

$dsvfile="README.RELEASE";
if (-f $dsvfile) {
    open(IN,"$dsvfile") || die "Can't open $dsvfile for reading.\n";
    $line=<IN>;
    chomp($line);
    if ($line =~ /^darksusy/i) { # check if raw number or with Darksusy text
	$dsversion=$line;
    } else {
	$dsversion="darksusy " . $line;
    }
    print "Using DarkSUSY version tag from file $dsvfile.\n";
    close(IN);
}


print "DarkSUSY version: $dsversion\n";

# Define subdirectory where data files are stored
$dsdatapath=$dsinstall . "data/";


# OK, Now dsversion, dsdatapath contain the information needed.
# Now make the include files.

### dsver.h

$defverfile="src/include/dsver.h";

open(OUT,">$defverfile") || die "Can't open $defverfile for writing.\n";
print OUT <<END
*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsver.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsver.h              ***
c----------------------------------------------------------------------c
END
;

$date=`date`;
print OUT "*** This file is created by config2.pl on $date\n";

$len = length ($dsversion);
print OUT "      character*${len} dsver\n";
$line=" " x 6 . "parameter(dsver='${dsversion}')\n";
$line=contline_string($line);
print OUT $line;

$len2=int(length($dsversion)/10)*10+50;
print OUT " " x 6 . "character*${len2} dsversion\n";
print OUT " " x 6 . "common /dsv/dsversion\n";
print OUT " " x 6 . "save /dsv/\n";

print OUT "***" . " " x 66 . "***\n";
print OUT "*" x 26 . " end of dsver.h " . "*" x 27 . "\n";

### dsdir.h

$defdirfile="src/include/dsdir.h";

open(OUT,">$defdirfile") || die "Can't open $defdirfile for writing.\n";
print OUT <<END
*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsdir.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsdir.h              ***
c----------------------------------------------------------------------c
END
;

$date=`date`;
print OUT "*** This file is created by config2.pl on $date\n";

$len = length ($dsdatapath);
print OUT "      character*${len} dsdatapath\n";
#$len = length ($dsinstall);
#print OUT "      character*${len} dsinstall\n";

#$line=" " x 6 . "parameter(dsdatadir='${dsdatadir}')\n";
#$line=contline_string($line);
#print OUT $line;

#$line=" " x 6 . "parameter(dsinstall='${dsinstall}')\n";
$line=" " x 6 . "dsdatapath='${dsdatapath}'\n";
$line=contline_string($line);
print OUT $line;

print OUT "***" . " " x 66 . "***\n";
print OUT "*" x 26 . " end of dsdir.h " . "*" x 27 . "\n";

### Split long lines

sub contline {
    my $line = $_[0];
    my $out;
    my $i;
 
    $out="";
    if (length($line) >= 71) {
        print "*** A line is too long. I will split it.\n";
        print "Line before: \n$line";
        while (length($line) >= 71) {
            for ($i=70; $i==20; $i--) {
                if (substr($line,$i,1) =~ m@(\*|\+|-|/)@ ) {
                    last;
                }
            }
            $out=$out . substr($line,0,$i-1) . "\n";
            $line = " " x 5 . "&  " . substr($line,$i-1)
        }
        $out = $out . $line;
        print "Line after: \n$out";
    } else {
        $out=$line;
    }
 
    return $out;
}


# The following only works for strings

sub contline_string {
    my $line = $_[0];
    my $out;
    my $max;
 
    $max=70;
    $out="";
    while (length($line) > $max) {
#        print "*** A line is too long. I will split it.\n";
#        print "Line before: \n$line";
	$max -=4 if (length($line) <$max+4); # to avoid cutting at wrong place
	$out .= substr($line,0,$max) . "'\n";
	$line= "     & //'" . substr($line,$max);
#        print "Line after: \n$out";
    }
    $out .= $line;
 
    return $out;
}
