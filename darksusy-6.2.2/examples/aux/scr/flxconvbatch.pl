#!/usr/bin/perl
#
# Script to run flxconv for a set of input parameters and produce one set of
# conversions. This script does not take any input arguments, instead you have
# to modify the setup below
#
# Author: Joakim Edsjö, edsjo@fysik.su.se
# Date: December 1, 2018

$prog='./flxconv';
$wh='su';
$dens=0.92;
$kind=1;
$emin=1.0;
$thmax=29;
$start_type=0;
$end_type=27; # 27=neutrino conversion rate, 28=mu flux

# These are the input parameters looped over, the entries need to match (i.e.
# entry i in each array is the one use for the i:th iteration in the loop)
@mass=qw(20
35
35
50
50
100
100
100
250
250
250
500
500
500
1000
1000
1000
3000
3000
3000
5000
5000
5000
10000
10000
10000);

@ch=qw(11
5
11
5
11
5
8
11
5
8
11
5
8
11
5
8
11
5
8
11
5
8
11
5
8
11);

$len=@mass;
print "Number of cases to go through: ",$len,"\n";
print "Conversion factors:\n";
$i=0;
while ($i<$len) {
# write input script, this might need to be modified depending on
# input options    
    open(TMP,">tmp.flxconv.in") || die "Can't open tmp.flxconv.in for writing.\n";
    print TMP $mass[$i], "\n";
    print TMP "$wh\n";
    print TMP $ch[$i],"\n";
    print TMP "$dens\n";
    print TMP "$start_type\n";
    print TMP "$end_type\n";
    print TMP "$kind\n";
    print TMP "$emin\n";
    print TMP "$thmax\n";
    print TMP "0\n";   
    close(TMP);
    $cf2f="ERROR";
    if (open(IN,"$prog < tmp.flxconv.in|")) {
	while (defined($line=<IN>)) {
	    if ($line =~ /cf2f = (.*)$/) {
		$cf2f=$1;
	    }
	}
	close(IN);
    } else {
	die "Could not run flxconv.\n";
    }
#    print $mass[$i]," ",$ch[$i],"\n";
    printf("%.6e\n",$cf2f);
    $i++;
    unlink("tmp.flxconv.in");
}

    
