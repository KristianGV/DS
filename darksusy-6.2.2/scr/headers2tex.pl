#!/usr/bin/perl -w
#
# Script to go through the files in DarkSUSY and produce a tex file
# with the headers of all the routines and additional documentation
# distributed with DarkSUSY
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: March 8, 2006
# Modified: April 4, 2014
#

# This is the list of subdirectories of src/ from which headers and
# documentation will be read.

$texfile="docs/Manual.tex";
$srcfile="docs/src-dirs-to-include.txt";

$include_headers=1; # include routine headers by default
foreach $arg (@ARGV) {
    if ($arg =~ /\-\-no\-headers/) {
	$include_headers=0;
        $texfile="docs/Manual-short.tex";
    }
}

# Now read in which directories to include and their names
open(IN,$srcfile) || die "Can't open $srcfile for reading.\n";
while(defined($line=<IN>)) {
    next if $line =~ /^#/;
    ($dir)=split(/\s/,$line);
    $name = $line;
    $name =~ s/^\S+\s//;
    chomp($name);
    push(@dirlist,$dir);
    $ch_names{$dir}=$name;
}


# Determine name of current DarkSUSY directory
$dsver=`pwd`;
chomp($dsver);
$dsver =~ s#^.*/##;
$dsver =~ s#^ds-##;

# Add revision number to dsver
$rev="";
if (open(IN,"svn info|")) {
    while(defined($line=<IN>)) {
	if ($line =~ /^Revision:\s+(\d+)$/) {
	    $rev=" (rev $1)";
	}
    }
    close(IN);
} else {
    print "Couldn't invoke svn. No revision number added to version tag.\n";
}

$dsver .= $rev;

$date=localtime;

open(TEX,">$texfile") || die "Can't open $texfile for writing.\n";

print_texheader($dsver,$date,$include_headers);

print_part("Prelude");

foreach $file (<docs/I*tex>) {
    include_texfile($file);
}


# Start going through the directories and files

$root=`pwd`;
chomp($root);

# Include chapter definition for src and general text for src structure

print_part("General DarkSUSY routines in src/");

#include_texfile("docs/Src-Intro.tex"); Now above in general I-file instead.

$model="";
foreach $dir (@dirlist) {    
    chdir($root) || die "Can't cd to root directory $root\n";
    chdir($dir) || die "Can't cd to $dir\n";
    print "Directory: $dir\n";
    $newtex=1;
    $newfortran=1;
    if ($dir =~ /src_models\/docs/) { # Now we reach the modules in src_models
	    print_part("Particle physics modules in src\\_models");
    } 
    if ($dir =~ /docs$/) {
	$model = $dir;
	$model =~ s#src_models/##;
	$model =~ s#/docs##;
	print "Documentation directory, using tex files directly.\n";
    } else {
	if ($dir =~ /src_models/) { # create sections, otherwise chapters
           print_texdir_sec($dir,$ch_names{$dir});
	} else {
           print_texdir_chap($dir,$ch_names{$dir});
	}
    }
    if ($include_headers){
	@files=<*.tex *.f>;
    } else {
	@files=<*.tex>;
    }
    foreach $file (@files) {
        print "   File: $file\n";
        if ($newtex and $file =~ /\.tex/) {
#            print TEX "\\section{Overview (from tex files)}\n";
            $newtex=0;
	}
        if ($newfortran and $file =~ /\.f/) {
	    if ($dir =~ /src_models/) {
		print TEX "\\subsection{Routine headers -- fortran files}\n";
	    } else {
		print TEX "\\section{Routine headers -- fortran files}\n";
	    }
            $newfortran=0;
	}

	print_texfile($file);
    }
#    chdir("..") || die "Can't cd to ..\n";
}
#chdir("..") || die "Can't cd to ..\n";
chdir($root) || die "Can't cd to root directory $root\n";



# Now insert some files at the end, start with appencies
$first=1;
foreach $file (<docs/A*tex>) {
    if ($first) {
	print TEX "\\appendix\n";
	$first=0;
    }
    include_texfile($file);
}

# Now insert the rest of the files at the end
foreach $file (<docs/E*tex>) {
    include_texfile($file);
}


print_texend();

# Copy figures to the docs/fig/ directory (not needed anymore)
# system ('mkdir docs/fig') unless -d 'docs/fig';
# system ('cp -p docs/fig/* docs/fig/');

print "Done.\n";
exit;

##########
sub print_texheader {
    my $ver=$_[0];   # DarkSUSY version
    my $date=$_[1];
    my $inc_head=$_[2];

    $template="./docs/Template.tex";
    open(IN,$template) or die "Can't open $template for reading.\n";
    if ($inc_head) {
	$title="Manual and long description of routines";
    } else {
	$title="Manual and short description of routines";
    }
    while(defined($line=<IN>)) {
        $line=~ s/\[Date\]/$date/;
        $line=~ s/\[DarkSUSYVersion\]/$ver/;
        $line=~ s/\[Title\]/$title/;
        if ($line=~ /\[Macros\]/) {
            foreach $file (<docs/headers/*.tex>) {
                include_texfile($file);
            }
	    next;
	}
        if ($line=~ /\[HeaderFiles\]/) {
            foreach $file (<docs/H*.tex>) {
                include_texfile($file);
            }
	    next;
	}
	print TEX $line;
    }
    close(IN);
}


##########
sub print_texdir_sec{
    my $dir=$_[0];
    my $name=$_[1];

    $printdir = $dir;
    $printdir =~ s/_/\\_/g;

print TEX <<END;

\\newpage
\\section\[$printdir: $name\]{\\codeb{src/$printdir}:\\\\ $name}
\\label{sec:$dir}

END
}

##########
sub print_texdir_chap{
    my $dir=$_[0];
    my $name=$_[1];

    $printdir = $dir;
    $printdir =~ s/_/\\_/g;

print TEX <<END;

\\chapter\[$printdir: $name\]{\\codeb{src/$printdir}:\\\\ $name}
\\label{ch:$dir}

END
}


##########
sub print_texfile{
    my $file=$_[0];
    my $tfile=$file;
    my $line;
    my $tline;
    my $i;
    $tfile =~ s#\_#\\\_#g;

    $verb=1;  # 0=print as tex, 1=print as verbatim (to look like code)
    $type="fortran";
    if ($file =~ /\.txt/) {
	$type="text";
    }
    if ($file =~ /\.tex/) { # print file directly
        $type="latex";
        $verb=0;
	open (IN,"$file") || die "Can't open $file for reading.\n";
        while(defined($line=<IN>)) {
	    print TEX $line;
	}
        return;
    }


print TEX <<END;

%%%%% routine $file %%%%%
\\begin{routine}{$tfile}
END
    if ($verb) {
	print TEX "\\begin{verbatim}\n";
    }

# Now go through the file

    open(IN,"$file") || die "Can't open $file for reading.\n";
    $i=0;
    $tex=0;
    $tex=1 if $type eq "latex";
    $clean=0; # 0=print line as is, 1=take away leading comment signs and spaces
    $defcont=0; # continuation line for subroutine / function found or not
    while(defined($line=<IN>)) {
        if (substr($line,0,1) ne " " or $type eq "text") {
            if ($line =~ /BeginTex/i) {
		print TEX "\\end{verbatim}\n" if ($verb);
                $tex=1;
                next;
	    }
            if ($line =~ /EndTex/i) {
		print TEX "\\begin{verbatim}\n" if ($verb);
                $tex=0;
                next;
	    }
            if ($clean==1 or $tex==1) {
		$line =~ s/^\S+\s*// if $type eq "fortran";
	    }
	    print TEX $line;
            $i++;
            next;
	}
        if ($line =~ /subroutine/i or $line =~ /function/i) {
            print TEX $line;
            if ($line =~ /\(/) {
		$defcont=1 unless $line =~ /\)/;
	    }
	    next;
	}
        if ($defcont) {
            print TEX $line;
            $defcont=0;
            $defcont=1 unless $line =~ /\)/;
	    next;
	}
        $tline=chomp($line);
	$tline =~ s/\s*//g;
        next if length($tline)==0;
        next if substr($line,5,1) ne " ";
        last;  # OK, now we have come to the code
    }
    close(IN);
    print TEX "No header found.\n" if $i==0;

    if ($verb==0 or $tex==1) {  # still in tex mode
	print TEX "\\end{routine}\n";
    } else {
	print TEX "\\end{verbatim}\n \\end{routine}\n";
    }

}


##########
sub print_texend{

print TEX <<'END';
\end{document}
END

}



##########
sub include_texfile{
    my $file =$_[0];
    print "Including tex file: $file\n";
    print TEX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    print TEX "%%% Here comes $file %%%\n";
    print TEX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    open (IN2,"$file") || die "Can't open $file for reading.\n";
    while (defined($line=<IN2>)) {
	print TEX $line;
    }
    close(IN2);
}

##########
sub print_part{
    my $ptitle=$_[0];
    print "Part division: $ptitle\n";
    print TEX "\n";
    print TEX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    print TEX "%%% Here comes part $ptitle\n";
    print TEX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    print TEX "\\part{$ptitle}\n\n";
}
