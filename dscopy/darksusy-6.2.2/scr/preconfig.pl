#!/usr/bin/perl -w
#
# Script to set up makefile.in:s and configure.ac and run auto-conf in one go.
# This script is intended only for DarkSUSY developers. It will
#   - run makemf.pl to generate all makefile.in:s
#   - update configure.ac with all the makefile.in:s that makemf.pl generated
#   - run autoconf to generate a new configure script
#     [Note: autoconf is usually not intalled by default, get with e.g.
#     apt-get, homebrew, MacPorts or similar package managers before running.]
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: 2015-10-29

if (not(-d "src" and -d "src_models")) {
    print "preconfig.pl needs to be run from the DarkSUSY root directory.\n";
    die "Please run it as scr/preconfig.pl from the DarkSUSY root.\n";
}

# Check if autoconf is available
if (not(`which autoconf`)) {
    print "ERROR: autoconf is not found on your system. You need to install it\n";
    die "prior to running preconfig.pl.\n";
}

# Set up makefiles
print "Running makemf.pl...\n";
system("scr/makemf.pl");

# Update main makefile.in
print "Running update_mainmf.pl...\n";
system("scr/update_mainmf.pl");

# Now we have a list of makefiles in tmp/makefiles-to-make.txt, add them
# configure.ac

open(IN,"<configure.ac") || die "Can't open configure.ac for reading.\n";
open(OUT,">configure.ac.tmp") || die "Can't open configure.ac.tmp for writing.\n";
while (defined($line=<IN>)) {
    if ($line =~ /^AC_CONFIG_FILES/) {
	print OUT $line;
	while(defined($line=<IN>)) {
	    if (not(($line =~ /^\s*src\//) || ($line =~ /^\s*src_models\//) || ($line=~ /^\s*src_halos\//) || ($line =~ /\]\)/))) {
		print OUT $line;
	    }
	    last if $line =~ /\]\)/;
	}

	# check if last line should be printed or not
	if (not(($line =~ /^\s*src\//) || ($line =~ /^\s*src\_models\//) || ($line =~ /^\s*src\_halos\//) )) {
	    $line =~ s/\]\)//;
	    print OUT $line;
	}

	open(MTM,"<tmp/makefiles-to-make.txt") || die "Can't open tmp/makefiles-to-make.txt for reading.\n";
	$mno=0;
	while(defined($mline=<MTM>)) {
	    if ($mno>0) {
		print OUT "\n";
	    }
	    chomp($mline);
	    print OUT "                 $mline";
	    $mno++;
	}
	close(MTM);
	print OUT "])\n";
	$line="";
    }
    print OUT $line;
}
close(IN);
close(OUT);

rename("configure.ac.tmp","configure.ac");

print "configure.ac is updated.\n";
print "Running autoconf...\n";
system("autoconf");
print "Done. You can now configure and compile DarkSUSY.\n";

