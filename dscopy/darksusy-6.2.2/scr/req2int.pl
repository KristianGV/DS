#!/usr/bin/perl -w
#
# Script to change text in source files from required to interface functions.
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: 2015-10-29

while (defined($file=shift)){
    $outfile="$file-tmp";
    open(IN,$file) || die "Can't open $file for reading.\n";
    open(OUT,">$outfile") || die "Can't open $outfile for writing.\n";
    $ch=0;
    while(defined($line=<IN>)) {
	if ($line =~ /^\S+\s*type\s*:\s*required/i) {
	    $line =~ s/REQUIRED/INTERFACE/;
	    $line =~ s/Required/Interface/;
	    $line =~ s/required/interface/i;
	    $ch++;
	}
	$ch += ($line =~ s/Required function/Interface function/);
	$ch += ($line =~ s/required function/interface function/i);
	print OUT $line;
    }
    close(OUT);
    close(IN);
    if ($ch>0) {
	print "File $file changed in $ch places.\n";
	rename($outfile,$file);
    } else {
	unlink($outfile);
    }	
}
print "Done.\n";

