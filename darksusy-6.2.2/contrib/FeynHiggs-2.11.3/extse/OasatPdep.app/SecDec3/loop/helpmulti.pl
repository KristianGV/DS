#! /usr/bin/perl -w
#
$dirbase=`pwd`;
chomp $dirbase;
$graph="P126";
$paramfile="${graph}.input";
$pointname="p";
#
$s1="0.";
$s2="0.";
$ms1="1.0";
#
$destination=$dirbase;
$filename="multi${graph}.input";
$s12min = 0.4;
$s12max = 20;
$Nsteps=6;
$lines=eval($Nsteps+1);
#
$lmin=0;
$lmax=$Nsteps;
#
open(EWRITE, ">","$destination/$filename") || die "cannot open $destination/$filename\n";
print EWRITE "paramfile=$paramfile\n";
print EWRITE "pointname=$pointname\n";
print EWRITE "# only lines 1 to [lines] will be calculated\n";
print EWRITE "lines=$lines\n";
print EWRITE "#---------------------------------------------\n";
print EWRITE "# each line below contains numsij sij's, then numpi2 pi^2, then numms2 mi^2\n";
print EWRITE "numsij=1\n";
print EWRITE "numpi2=3\n";
print EWRITE "numms2=1\n";
print EWRITE "xplot=1\n";
print EWRITE "#---------------------------------------------\n";
for ($l=$lmin;$l<=$lmax;$l++) { 
    $s12 = $s12min + $l*($s12max - $s12min)/$Nsteps;
    print "s12= $s12\n";
    print EWRITE "$s12,$s1,$s2,$s12,$ms1\n";
}
close(EWRITE)








 

