#!/usr/bin/perl -w
#
# This script goes through the files given as arguments and changes names
# according to the table in the hash %names below. The names on the left
# are changed to the names on the right. If the newname is set to - this
# has a special meaning; the new name is 'ds' plus the old name appended to
# it.

# First define list of names to change
# If the new name is '-', the new name is set to 'ds' plus the old name
%names=qw(dsntsundens dssem_sundens
dsntsuncdens dssem_suncdens
dsntsuncdensint dssem_suncdensint
dsntsuncdfunc dssem_suncdfunc
dsntsunmass dssem_sunmass
dsntsundenscomp dssem_sundenscomp
dsntsunmfrac dssem_sunmfrac
dsntsunne dssem_sunne
dsntsunne2x dssem_sunne2x
dsntsunpot dssem_sunpot
dsntsunpotint dssem_sunpotint
dsntsunread dssem_sunread
dsntsunset dssem_sunset
dsntsunvesc dssem_sunvesc
dsntsunx2z dssem_sunx2z
dsntsunz2x dssem_sunz2x
dsntearthdens dssem_earthdens
dsntearthdenscomp dssem_earthdenscomp
dsntearthdensint dssem_earthdensint
dsntearthmass dssem_earthmass
dsntearthmassint dssem_earthmassint
dsntearthne dssem_earthne
dsntearthpot dssem_earthpot
dsntearthpotint dssem_earthpotint
dsntearthvesc dssem_earthvesc
dsntsun.h dssem_sun.h
dsntearth.h dssem_earth.h
dsntepfunc dssem_epfunc
dsntedfunc dssem_edfunc
dsntspfunc dssem_spfunc
dsntcapsun dssenu_capsun
dsntcapsunnum dssenu_capsunnum
dsntcapsunnumff dssenu_capsunnumff
dsntcapsunnumffi dssenu_capsunnumffi
dsntcapsunnumi dssenu_capsunnumi
dsntcapsuntab dssenu_capsuntab
dsntcapearth dssenu_capearth
dsntcapearth2 dssenu_capearth2
dsntcapearthfull dssenu_capearthfull
dsntcapearthnum dssenu_capearthnum
dsntcapearthnumi dssenu_capearthnumi
dsntcapearthtab dssenu_capearthtab
dsntcapcom dssenu_capcom
dsntcsint dssenu_csint
dsntcsint2 dssenu_csint2
dsntcsintff dssenu_csintff
dsntcsintff2 dssenu_csintff2
dsntcsintff3 dssenu_csintff3
dsntceint dssenu_ceint
dsntceint2 dssenu_ceint2
dsntctabcreate dssenu_ctabcreate
dsntctabget dssenu_ctabget
dsntctabread dssenu_ctabread
dsntctabwrite dssenu_ctabwrite
dsntse dssenu_se
dsntsefull dssenu_sefull
dsntss dssenu_ss
dsntlitlf_e dssenu_litlf_e
dsntlitlf_s dssenu_litlf_s
dsntveoutjupiter dssenu_veoutjupiter
dsntmoderf dsmoderf
dsntdqagse dssenu_dqagse
dsntdqagseb dssenu_dqagseb
dsntdqk21 dssenu_dqk21
dsntdqk21b dssenu_dqk21b
dsntfoveru dssenu_foveru
dsntfoveruearth dssenu_foveruearth
dsntset  dssenu_set
dsff dssea_ff
dsff2 dssea_ff2
dsff3 dssea_ff3
dsfff2 dssea_fff2
dsfff3 dssea_fff3
dsgauss1 dssea_gauss1
dsntismbkg dssea_ismbkg
dsntsunbkg dssea_sunbkg
dsntismrd dssea_ismrd
dshonda dssea_honda
dslnff dssea_lnff
dsatm_mu dssea_atmmu
dsntnuism dssea_nuism
dsntnusun dssea_nusun
dshiprecint dssenu_hiprecint
dshiprecint2 dssenu_hiprecint2
dsntICanglike dsseIC_anglike
dsntICangres dsseIC_angres
dsntICbgangpdf dsseIC_bgangpdf
dsntICbginit dsseIC_bginit
dsntICbglikeprecomp dsseIC_bglikeprecomp
dsntICbgpredinit dsseIC_bgpredinit
dsntICbgspec dsseIC_bgspec
dsntICbounds dsseIC_bounds
dsntICdP_SdE dsseIC_dP_SdE
dsntICdPhi_SdE dsseIC_dPhi_SdE
dsntICeainit dsseIC_eainit
dsntICedisp dsseIC_edisp
dsntICedispcheckout dsseIC_edispcheckout
dsntICedispinit dsseIC_edispinit
dsntICeffarea dsseIC_effarea
dsntICeventinit dsseIC_eventinit
dsntICinit dsseIC_init
dsntICnlike dsseIC_nlike
dsntICpsf dsseIC_psf
dsntICpval dsseIC_pval
dsntICsigintegrand dsseIC_sigintegrand
dsntICsignal dsseIC_signal
dsntICspecintegrand dsseIC_specintegrand
dsntICspeclike dsseIC_speclike
dsntannrate dssenu_annrate
dsntannrateff dssenu_annrateff
);  
# Now go through name changes and add ds where appropriate
foreach $name (keys %names) {
    $newname=$names{$name};
    if ($newname eq "-") {
        $newname = "ds$name";
        $names{$name}=$newname;
    }
}


while(defined($file=shift)) {
    print "Now fixing file $file\n";
    $outfile="tmp.f";

    open(FILE,$file) || die "Can't open $file\n";
    open(OUT,">$outfile") || die "Can't open $outfile for writing.\n";
    $ch=0;
    $i=0;
    while(defined($line=<FILE>)) {
        $i++;
        chomp($line);
#	$ch += ($line =~ s/c     include 'dssm.h'/      include 'dsmssm.h'/);
#        $ch += ($line =~ s/edsjo\@cfpa.berkeley.edu/edsjo\@fysik.su.se/);
#        $ch += ($line =~ s/edsjo\@teorfys.uu.se/edsjo\@fysik.su.se/);
#        $ch += ($line =~ s/edsjo\@physto.se/edsjo\@fysik.su.se/);
#	$ch += ($line =~ s/include 'dssun.h'/include 'dsntsun.h'/);
#	$ch += ($line =~ s/include 'dsearth.h'/include 'dsntearth.h'/);
#	$ch += ($line =~ s/include 'dsmssm.h'/include 'dssm.h'\n      include 'dsmssm.h'/);
#	 $ch += ($line =~ s/'halo.h'/'hmcom.h'/);
#	 $ch += ($line =~ s/'ddcom.h'/'ddset.h'/);
#        $ch += ($line =~ s/double precision/real*8/i);
#        $ch += ($line =~ s/(\W)gef_int(\W|$)/$1dsf_int$2/ig);
#        $ch += ($line =~ s/(\W)gef_int2(\W|$)/$1dsf_int2$2/ig);
#        $ch += ($line =~ s/(\W)absq/$1dsabsq/gi);
#        $ch += ($line =~ s/\s+$//);
#Now change names to ds...
        foreach $name (keys %names) {
	    $newname=$names{$name};
            $ii=1;
            $jj=0;
            while ($ii > 0) {
		$ii = ($line =~ s/(\W|^)$name(\W|$)/$1$newname$2/gi);
		$jj += $ii;
	    }
	    $ch += $jj;
        }
#        $line=lc($line);   # uncomment if you want lowercase
        if (length($line)>72 and substr($line,0,1)=~ /\s/
           and not(substr($line,0,72) =~ /!/)) {
	    print "***Line $i is too long\n";
	}
   	print OUT "$line\n";
    }

    close(FILE);
    close(OUT);
    if ($ch>0) {
	print "  $file changed ($ch place", 
        ($ch==1) ? '' : "s",
        ") - replacing with new version\n";
        rename($outfile,$file);
        $newfile=$file;
        $ch=0;
        foreach $name (keys %names) {
	    $newname=$names{$name};
            $ch += ($newfile =~ s/(\W|^)$name(\W|$)/$1$newname$2/gi);
        }
        if ($ch>0) {
            print "  filename changed to $newfile\n";
#	    system("svn mv $file ../../../src/se_nu/$newfile");
#            rename($file,$newfile);
	}

    } else {
	unlink($outfile);
    }
}

