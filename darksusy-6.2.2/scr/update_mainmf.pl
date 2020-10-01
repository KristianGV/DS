#!/usr/bin/perl -w
#
# Script to update the main makefile.in in the DarkSUSY root.
# It will check which modules exist in src_models and add them (with build
# instructions in makefile.in).
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: 2016-02-07

if (not(-d "src" and -d "src_models")) {
    print "update_mainmf.pl needs to be run from the DarkSUSY root directory.\n";
    die "Please run it as scr/update_mainmf.pl from the DarkSUSY root.\n";
}

$date=`date +"%b %d, %Y"`;
chomp($date);

# First figure out which modules we have
$modules=();
foreach $dir (<src_models/*>) {
    if (-d $dir) {
	$add=1;
	$add=0 if ($dir =~ /docs/);
	$add=0 if ($dir =~ /include/);
	@srcfiles = <$dir/*/*.{f,f90,F,F90,c,C,cc,CC}>;
        $add=0 if ($#srcfiles < 0);     
	if ($add == 1) {
	    push(@modules,$dir);
	}
    }
}
print "Found the following modules: @modules\n";

# Update main makefile.in
print "Updating main makefile.in... \n";

open(IN,"<makefile.in") || die "Can't open makefile.in for reading.\n";
open(OUT,">makefile.in-tmp") || die "Can't open makefile.in-tmp for writing.\n";
while(defined($line=<IN>)) {

    if ($line =~ /^darksusy\s*:/) { # main target line
	$line="darksusy : ds_core";
	foreach $mod (@modules) {
	    $shortmod=$mod;
	    $shortmod =~ s/src_models\///;
	    $line .= " ds_$shortmod";
	}
	$line .= " inst_tab_if_loc\n";
    }

    if ($line =~ /\[MODSTART\]/) {
	$newline=$line;
	$newline .= "\n";
	# Read through old module build instructions and trash them
	while (defined($line=<IN>)) {
	    last if ($line =~ /\[MODEND\]/);
	}

	# Add new build instructions
	foreach $mod (@modules) {
	    next if $mod =~ /common/;
	    $shortmod=$mod;
	    $shortmod =~ s/src_models\///;
	    $newline .= "# Module ds_$shortmod. Code below auto-generated by update_mainmf.pl\n";
	    $newline .= "ds_$shortmod :\n";
	    $newline .= "\tmkdir -p lib\n";
	    $newline .= "\tmkdir -p tmp/build-$shortmod\n";
	    $newline .= "\tmkdir -p tmp/build-$shortmod-user\n";
	    $newline .= "\tcd \$(DS_ROOT)/$mod; make all\n";
	    $newline .= "\n";
	}

#	# Add instructions to build executables in examples (pg20160412)
#	$newline .= "# Code below auto-generated by update_mainmf.pl\n";
#	$newline .= "rmmodx :\n";
#	foreach $mod (@modules) {
#         $mod =~ s/src_models\///;
#	  $newline .= "\trm -f \$(DS_ROOT)/examples/ds_$mod\n";
#	}
#	$newline .= "\n";

	# Add instructions to remove executables in examples (pg20160412)
	$newline .= "rmmodx :\n";
	foreach $mod (@modules) {
	  $shortmod=$mod;
          $shortmod =~ s/src_models\///;
	  $newline .= "\trm -f \$(DS_ROOT)/examples/ds_$shortmod\n";
	}
	$newline .= "\n";

	$newline .= $line; # add last comment as well
	$line = $newline;
    }
    print OUT $line;
}
close(IN);
close(OUT);

rename("makefile.in-tmp","makefile.in");
print "Done.\n\n";


# Update module_defs.F in src_models
print "Updating module_defs.F... \n";

$infile="src_models/include/module_defs.F";
$outfile="src_models/include/module_defs.F-tmp";
open(IN,"<$infile") || die "Can't open $infile for reading.\n";
open(OUT,">$outfile") || die "Can't open $outfile for writing.\n";
while(defined($line=<IN>)) {

    if ($line =~ /\[MODSTART\]/) {
	$newline=$line;
	# Read through old module build instructions and trash them
	while (defined($line=<IN>)) {
	    last if ($line =~ /\[MODEND\]/);
	}

	# Add new build instructions
	$i=0;
	foreach $mod (@modules) {
	    if ($mod =~ /src_models/) {
		$shortmod=$mod;
		$shortmod =~ s/src_models\///;
		$i++;
		$newline .= "#define MODULE_$shortmod $i\n";
	    }
	  }
	$newline .= $line; # add last comment as well
	$line = $newline;
    }
    print OUT $line;
}
close(IN);
close(OUT);

rename("$outfile","$infile");
print "Done.\n\n";

