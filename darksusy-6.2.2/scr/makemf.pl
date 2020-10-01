#!/usr/bin/perl -w
#
# Script to create makefile.in's in source directories. It includes all
# *.f files and properly defines all included files as dependencies.
# Use this script with the directory path to the source directory as 
# argument. It will include slightly different include directories
# depending on if the source directory is in src/ or in src_models/xxx/
#
# As of October, 2015 the script will no longer have a hard-coded
# list of directories, instead it will (if no arguments are given)
# go through ALL directories in src/ and src_models/ and add makefile.in's
#
# It will also figure out if there are main programs in a directory, like
# in src_models/mssm/examples and then create the makefile.in for a main
# program instead of a library
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: August 31, 2000.
# Changed by Paolo Gondolo 2011 to separately compile every object file
# into a build directory.
# Updated to new diretory structure by Joakim Edsjö, 2013.
# Updated to be more general and handle main programs by Joakim Edsjö, 2015.

use Cwd;
#use Text::Diff;

$pwd=cwd();
@errorlog=();

if (not(-d "src" and -d "src_models")) {
    print "makemf.pl needs to be run from the DarkSUSY root directory.\n";
    print "Please run it as scr/makemf.pl [directories] from the DarkSUSY root.\n";
}

if (@ARGV >0) {
    @dirlist=@ARGV;
    $auto=0;
} else {    
    # Build dirlist automatically
    $auto=1;
    # create tmp directory if it does not exist and open file for list of makefiles
    mkdir("tmp") if (not(-d "tmp"));
    unlink("tmp/makefiles-to-make.txt");
    open(MK,">tmp/makefiles-to-make.txt");

    @dirlist=();
    print "Automatically determining directories to make makefile.in's in.\n";
    print "Directories in which to generate a regular or mains makefile.in in:\n";
    foreach $file (<src/* src_models/*{}/*>) {  # {} just to avoid being treated as comment by editor
       if (-d $file) {
	   $file="" if $file =~ /include/;
	   $file="" if $file =~ /docs/;
	   if ($file =~ /src_models\/(.*)\//) { # check if whole module is empty
	       $dir=$1;
	       @srcfiles = <src_models/$dir/*/*.{f,f90,F,F90,c,C,cc,CC}>;
               $file="" if ($#srcfiles < 0);
           }
	   if ($file ne "") { # Check if directory is empty
	       @srcfiles=<$file/*.f $file/*.F $file/*.f90 $file/*.F90 $file/*.c>;
	       if ((scalar @srcfiles) == 0) {
		   push (@errorlog,"\nDirectory $file does not contain source files.\n");
		   push (@errorlog,"Please remove directory, or add source files.\n");
		   $file="";
	       }
	   }
	   if ($file ne "") {
 	      print "   Directory: $file\n";
              push(@dirlist,$file);
	      print MK "$file/makefile\n";
	   }
       } 
    }

    # Add to dirlist directories where to make upper-level makefiles
    print "Automatically determine where to add steering makefile.in's: \n";
    foreach $file (<src src_models/*>) {  
       if (-d $file) {
	   $file="" if $file =~ /include/;
	   $file="" if $file =~ /docs/;
	   $file="" if $file =~ /common/;
	   if ($file =~ /src_models\/(.*)/) { # check if whole module is empty
	       $dir=$1;
	       @srcfiles = <src_models/$dir/*/*.{f,f90,F,F90,c,C,cc,CC}>;
               $file="" if ($#srcfiles < 0);
           }
#	   $file="" if dir_with_main($file);
	   if ($file ne "") {
	       print "   Directory: $file\n";
	       push(@dirlist,$file);
 	       print MK "$file/makefile\n";
	   }
       } 
    }
    close(MK);
    print "Let's get started creating makefile.in's.\n\n";
}

$date=`date +"%b %d, %Y"`;
chomp($date);

print "Going through subdirectories and creating regular makefile.in's...\n";

if ($auto) {
print "Also creating list of directories from which to build documentation...\n";
# Documentation header
open(DOC,">docs/src-dirs-to-include.txt.auto") || die "Can't open docs/src-dirs-to-include.txt.new for writing.\n";
print DOC<<END;
# This file includes the directories in src/ and src_models/,
# which should be searched by the script headers2tex.pl.
# The script will look in these
# directories for documentation tex-files and for headers in the
# Fortran routines.
#
# The format below is that each directory resides on one row with the first
# element being the directory name and the rest being the title
# of that chapter as it will appear in the manual. The order in which 
# the directories appear below is the order in which the chapters/sections will
# appear in the manual. Note that for particle physics modules, the first entry
# should be to the subdirectory docs where general text and documentation 
# about the model should be given.
# This version of the file is automatically created by makemf.pl
#
END
}

foreach $dir (@dirlist) {
    $errno=scalar @errorlog;
    print "Now taking care of $dir...";
    chdir($dir) || die "Can't cd to $dir\n";

### Steering makefiles
    if (make_steering("makefile.in")==1) { # steering makefile
	print " ...steering type...";
	open(IN,"<makefile.in") || die "Can't open makefile.in for reading in $dir\n";
	open(FILE,">makefile.in.new") || die "Can't open makefile.in.new for writing in $dir\n";
	while(defined($line=<IN>)) {
	    if ($line =~ /^DIRS/) { # now update dirs
		while($line=~ /\\/) { # take care of continuation lines
		    $line=<IN>;
	    	}
		@dirs_to_build=dirs_without_main(".");
		if ($dir =~ /src\_models/) {
		    @cdirs=dirs_without_main("../common/");
		    foreach $cdir (@cdirs) {
			$cdir="../common/$cdir";
			push(@dirs_to_build,$cdir);
		    }
		}
		$line = "DIRS = " . join(" ",@dirs_to_build) . "\n";
	    }
	    print_line_no_nl($line);
	}
	close(IN);
	close(FILE);
	
	rename("makefile.in.new","makefile.in");
    } elsif (dir_with_main(".")==1) {
### Directory with main programs, make makefile.in for building executables
	print " ...main program type...";
	check_auto("makefile.in"); # check if OK to generate
	$modname=$dir;
	$modname =~ s/^.*src_models\///;
	$modname =~ s/\/.*$//;
	open(FILE,">makefile.in.new") || die "Can't open makefile.in.new in $dir.\n";
	print_header_mains($modname);
	@mains=list_of_mains(".");

	print FILE "all:";
	foreach $file (@mains) {
	    $base=$file;
	    $base =~ s/\.f//;
	    $base =~ s/\.F//;
	    $base =~ s/\.f90//;
	    $base =~ s/\.F90//;
	    print FILE " $base";
	}
	print FILE "\n\n";
	    
	foreach $file (@mains) {
	    $base=$file;
	    $base =~ s/\.f//;
	    $base =~ s/\.F//;
	    $base =~ s/\.f90//;
	    $base =~ s/\.F90//;
	    print FILE "$base: $file \$(LIBDEPS)\n";
	    print FILE "\t \$(FF) \$(FOPT) \$(INC) \$(INC_MODULE) -L\$(LIB) -o $base $file \\\n";
	    print FILE "\t \$(LIBS)\n\n";
	}
	close(FILE);
	rename("makefile.in.new","makefile.in");
	    	
    } else {
### Regular automatically generated makefiles	    
	print " ...regular library type...";
	@tmpfiles=<*.f *.F *.f90 *.F90 *.c>;
	@files=();
	$suf=".f";
	foreach $file (@tmpfiles) {
	    push(@files,$file) unless ($file eq "ds$dir.f" or $file eq "ds$dir.F" or $file eq "ds$dir.f90" or $file eq "ds$dir.F90" or $file eq "ds$dir.c");
	    $suf=".F" if $file =~ /\.F$/;
	    $suf=".f90" if $file =~ /\.f90$/;
	    $suf=".F90" if $file =~ /\.F90$/;
	    $suf=".c" if $file =~ /\.c$/;
	}
	@files=sort(@files);
	@checkfiles=@files;
	if ($dir=~ "src_models") { # Add src_models/<module>/include/*.h to check for dependencies (if they include other files)
	    foreach $file (<../include/*.h>) {
		push(@checkfiles,$file);
	    }
	}
	@deps=getdeps(@checkfiles);
	check_auto("makefile.in");
	open(FILE,">makefile.in.new") || die "Can't open makefile.in.new in src/$dir.\n";
	if ($dir=~ "src_models") {
	    @incdirs=qw(../include ../../include ../../../src/include ../../../contrib/include);
	    check_deps(\@deps,\@incdirs);
	    $modname=$dir;
	    $modname =~ s/^.*src_models\///;
	    $modname =~ s/\/.*$//;
	    if ($dir =~ 'user_replaceable'){
		print_header_src_models($modname,'user');
		$libname="libds_${modname}_user.a";
		$libdir="../../../lib";
		$fdep="$libname files-to-include.txt";
	    } else {
		if ($auto) {
		    $modno{$modname}++;
		    if ($modno{$modname}==1) {
			if ($dir =~ "src_models") {
			    print DOC "src_models/$modname/docs Module $modname\n";
			}
		    }
		    $modtitle=$dir;
		    $modtitle =~ s/_/\\_/g;
		    print DOC "$dir $modtitle\n";
		}
		print_header_src_models($modname,'');
		$fdep='';
	    }
	    $back="../../..";
	} else { # /src
	    @incdirs=qw(../include ../../contrib/include);
	    check_deps(\@deps,\@incdirs);
	    if ($dir =~ 'user_replaceable'){
		print_header_src('user');
		$libname="libds_core_user.a";
		$libdir="../../lib";
		$fdep="$libname files-to-include.txt";
	    } else {
		if ($auto) {
		    $modtitle=$dir;
		    $modtitle =~ s/_/\\_/;
		    print DOC "$dir $modtitle\n";
		}
		print_header_src('');
		$fdep='';
	    }
	    $back="../..";
	}
	$line="INC_DEP = " . join(" ",@deps);
	print_line($line);
	print FILE "vpath %.h \$(INC)\n\n";
	if ($dir =~ "user_replaceables") {
	    $line = "include files-to-include.txt";
	    $line2= "# To avoid creating an empty archive\n";
	    $line2 .= "SRC += empty_dummy.f\n\n";
	} else {
	    $line="SRC = " . join(" ",@files);
	    $line2="";
	}
	print_line($line);
	print FILE $line2;
	print FILE "OBJ1 = \$(patsubst %.f,\$(DOBJ)/%.o,\$(SRC))\n\n";
	print FILE "OBJ2 = \$(patsubst %.F,\$(DOBJ)/%.o,\$(OBJ1))\n\n";
	print FILE "OBJ3 = \$(patsubst %.f90,\$(DOBJ)/%.o,\$(OBJ2))\n\n";
	print FILE "OBJ4 = \$(patsubst %.F90,\$(DOBJ)/%.o,\$(OBJ3))\n\n";
	print FILE "OBJ = \$(patsubst %.c,\$(DOBJ)/%.o,\$(OBJ4))\n\n";
	print FILE "all : \$(OBJ) $fdep mods\n\n";
	if ($dir =~ 'user_replaceable') {
	    $line = "LIB=$libdir\n\n";
	    $line .= "$libname: always\n";
	    $line .= "\trm -f \$(LIB)/$libname\n";
	    $line .= "\tar rS \$(LIB)/$libname \$(OBJ)\n";
	    $line .= "\tranlib \$(LIB)/$libname\n\n";
	    $line .= ".PHONY: always\n\n";
	    print FILE $line;
	}
	print_compile('');
	if ($dir =~ 'user_replaceable') {
	    print FILE ".NOTPARALLEL: \n";
	}
	print FILE "\nmods : \n";
	print FILE "\tif \[ -f *.mod \]; then\\\n";
	print FILE "\t\tcp -p *.mod ../include/ ;\\\n";
	print FILE "\tfi\n";
	close(FILE);
	rename("makefile.in.new","makefile.in");
    }
    #    chdir($back) || die "Can't cd to $back\n";
    chdir($pwd) || die "Can't cd to root directory $pwd\n";
    if ($errno == scalar @errorlog) {
	print " done.\n";
    } else {
	print " done, but errors found (see below)\n";
    }
}
if (scalar @errorlog >0) { # errors found
    print "\n**************************************************\n";
    print "Some errors were found when making makefile.in:s:\n";
    print "**************************************************\n";
    foreach $line(@errorlog) {
	print $line;
    }
    print "**************************************************\n";	
} else {
    print "makemf.pl finished successfully.\n";
}

close(DOC);

### getdeps ###
sub getdeps{
    %depstmp=();
    my @files_tmp = @_;
    foreach $file (@files_tmp) {
        open(DFILE,"$file");
        while(defined($line=<DFILE>)) {
	    if ($line =~ /^ .*include\s+'(.+)'/i) {
	        $depstmp{$1}++;
    	    }
        }
        close(DFILE);
    }
    return (sort(keys %depstmp));
}


### print_header ###
sub print_header_src{
    my ($user) = shift @_;
    if ($user =~ 'user') {
	$builddir='../../tmp/build-src-user';
    } else {
	$builddir='../../tmp/build-src';
    }
print FILE <<END;
# Makefile for $dir directory
# Author: Joakim Edsjo, edsjo\@fysik.su.se
# Changed by Paolo Gondolo (2011), Joakim Edsjo (2013, 2014, 2015, 2016, 2017)
# This file is automatically created by makemf.pl.

# Define fortran compiler and options (set when ./configure is run
# in the DarkSUSY root directory).
FF=\@F77\@
FOPT=\@FOPT\@

FC=\$(FF)
FFLAGS=\$(FOPT) -c \$(DINC)

CC=\@CC\@
COPT=\@CFLAGS\@
CCFLAGS=\$(COPT) -c \$(DINC)

# Dependencies and libraries
INC=../include ../../contrib/include
DINC=-I../include -I../../contrib/include
DOBJ=$builddir

END
}

### print_header_src_models ###
sub print_header_src_models{
    my ($modname) = shift @_;
    my ($user) = shift @_;
    my $locinc="";
    my $locdinc="";
    if ($user =~ 'user') {
	$builddir="../../../tmp/build-$modname-user";
    } else {
	$builddir="../../../tmp/build-$modname";
    }
    # Check if module has its own include directory
    if (-d '../include') {
	$locinc="../include ";
	$locdinc="-I../include ";
    } 
print FILE <<END;
# Makefile for $dir directory
# Author: Joakim Edsjo, edsjo\@fysik.su.se
# Changed by Paolo Gondolo (2011), Joakim Edsjo (2013)
# This file is automatically created by makemf.pl.

# Define fortran compiler and options (set when ./configure is run
# in the DarkSUSY root directory).
FF=\@F77\@
FOPT=\@FOPT\@

FC=\$(FF)
FFLAGS=\$(FOPT) -c \$(DINC)

CC=\@CC\@
COPT=\@CFLAGS\@
CCFLAGS=\$(COPT) -c \$(DINC)

# Dependencies and libraries
INC=${locinc}../../include ../../../src/include ../../../contrib/include
DINC=${locdinc}-I../../include -I../../../src/include -I../../../contrib/include
DOBJ=$builddir

END
}

### print_header_mains ###
# Prints a header for a folder where we have main programs
sub print_header_mains{
    my $module=$_[0];
print FILE <<END;
# Makefile for $dir directory
# Author: Joakim Edsjo, edsjo\@fysik.su.se
# Changed by Paolo Gondolo (2011), Joakim Edsjo (2013, 2014, 2015)
# This file is automatically created by makemf.pl.
# To aviod missing libraries this automatically created makefile.in     
# links to all libraries that could typically be needed. 
# If you use this makefile.in as a template for your own makefile
# you can of course delete the linking to not needed libraries.
#
# Define fortran compiler and options (set when ./configure is run
# in the DarkSUSY root directory).
FF=\@F77\@
FOPT=\@FOPT\@

CC=\@CC\@
COPT=\@CFLAGS\@

# Determine where DarkSUSY is installed
prefix=\@prefix\@
DS_INSTALL=\${prefix}

LIB=\$(DS_INSTALL)/lib
INC=-I./ -I\$(DS_INSTALL)/src/include -I\$(DS_INSTALL)/src_models/include -I\$(DS_INSTALL)/contrib/include -I\$(DS_INSTALL)/src/templates
INC_MODULE=-I\$(DS_INSTALL)/src_models/$module/include    

LIBDEPS=\$(LIB)/libds_core.a \$(LIB)/libds_${module}.a \$(LIB)/libisajet.a
LIBS=-lds_core -lds_${module} -lds_${module}_user -lFH -lHB -lisospin -lisajet

END
}

### print_line ###
sub print_line{
    my $line=$_[0];
    $cols=72;
    while(length($line) != 0) {
	if (length($line)>$cols) {
            $i=rindex($line," ",$cols);
	    print FILE substr($line,0,$i);
	    print FILE " \\\n";
            substr($line,0,$i+1)="";
	} else {
	    print FILE "$line\n";
	    $line="";
	}
    }
    print FILE "\n";
}

### print_line ###
sub print_line_no_nl{ # as print_line but no \n are added at the end
    my $line=$_[0];
    $cols=72;
    while(length($line) != 0) {
	if (length($line)>$cols and not($line =~ /^#/)) {
            $i=rindex($line," ",$cols);
	    if (substr($line,0,$i) =~/#/) { # comment found
		print FILE "$line";
		$line="";
	    } else {
		print FILE substr($line,0,$i);
		print FILE " \\\n";
		substr($line,0,$i+1)="";
	    }
	} else {
	    print FILE "$line";
	    $line="";
	}
    }
}

### print_compile ###
sub print_compile{
    my ($fdep) = shift @_;
    print FILE "\$(DOBJ)/%.o : %.F \$(INC_DEP) $fdep\n";
    print FILE "\t\$(FC) \$(FFLAGS) \$< -o \$@\n\n";
    print FILE "\$(DOBJ)/%.o : %.f \$(INC_DEP) $fdep\n";
    print FILE "\t\$(FC) \$(FFLAGS) \$< -o \$@\n\n";
    print FILE "\$(DOBJ)/%.o : %.F90 \$(INC_DEP) $fdep\n";
    print FILE "\t\$(FC) \$(FFLAGS) \$< -o \$@\n\n";
    print FILE "\$(DOBJ)/%.o : %.f90 \$(INC_DEP) $fdep\n";
    print FILE "\t\$(FC) \$(FFLAGS) \$< -o \$@\n\n";
    print FILE "\$(DOBJ)/%.o : %.c \$(INC_DEP) $fdep\n";
    print FILE "\t\$(CC) \$(CCFLAGS) \$< -o \$@\n\n";
}

### dir_with_main ###
# Figure out if a directory contains main programs or not
# Returns 1 if it does, 0 if it doesn't
sub dir_with_main{
    my $dir=$_[0];
    $main=0;
    foreach $file (<$dir/*.f $dir/*.F $dir/*.f90 $dir/*.F90>) {
	open(TMP,"<$file") || die "Can't open file $file for reading.\n";
	while(defined($line=<TMP>)) {
	    $main=1 if $line =~ /^\s*program/i;
	}
	close(TMP);
    }
    return $main;
}

### list_of_mains ###
# Get a list of main programs back for a given directory
sub list_of_mains{
    my $dir=$_[0];
    @list_dirs=();
    foreach $file (<$dir/*.f $dir/*.F $dir/*.f90 $dir/*.F90>) {
	open(TMP,"<$file") || die "Can't open file $file for reading.\n";
	while(defined($line=<TMP>)) {
	    if($line =~ /^\s*program/i) {
		$file =~ s/^.*\///;
		push(@list_dirs,$file);
	    }
	}
	close(TMP);
    }
    @list_dirs=sort(@list_dirs);
    return @list_dirs;
}


### dirs_without_main ###
# Returns the a list of all the subdirectories of the given directory
# which does not contain main programs
sub dirs_without_main{
    my $dir=$_[0];
    @dirs=();
#    print "Directory: $dir\n";
    foreach $file (<$dir/*>) {
	if (-d $file && dir_with_main($file)==0) {
	    $file="" if $file =~ /include/;
	    $file="" if $file =~ /docs/;
	    $file="" if $file =~ /user_replaceables/;
	    $file="" if $file =~ /templates/;
	    $file =~ s/.*\///;
            # Check if fortran 90 modules exist, the put it first on list
	    if (contain_modules($file)) {
		@dirs=($file,@dirs);
	    } else {
		push(@dirs,$file) if $file ne "";
	    }
	}
    }
#    print "Subdirectories: @dirs\n";
    return @dirs;
}

	    
### make_auto ###
# make_auto checks if a given makefile is automatically generated or not
# it returns 1 if it is either automatically generated or does not exist
# otherwise it returns 0
sub make_auto{
    my $file=$_[0];
    $ret=0;
    if (-f $file) {
	open(TMP,"<$file") || die "Can't open $file for reading.";
	$af=0;
	while(defined($line=<TMP>)) {
	    $af++ if $line =~ /This file is automatically created by makemf\.pl/;
	}
	close(TMP);
	$ret=1 if $af>0;
    } else {
	$ret=1;
    }
    return $ret
}

### make_steering ###
# make_steering checks if a given makefile.in is a steering makefile where
# only the variable DIRS should be updated.
# It returns 1 if it is a steering makefile, 0 if not and -1 if
# the makefile.in does not exist
sub make_steering{
    my $file=$_[0];
    $ret=0;
    if (-f $file) {
	open(TMP,"<$file") || die "Can't open $file for reading.";
	$af=0;
	while(defined($line=<TMP>)) {
	    $af++ if $line =~ /DIRS/;
	}
	close(TMP);
	$ret=1 if $af>0;
    } else {
	$ret=-1;
    }
    return $ret
}


### check_auto ###
# check_auto checks if the makefile.in is automatically generated as expected
sub check_auto{
    my $file=$_[0];
    if (make_auto($file)==0) { # File not automatically created
	print "\n\n************************************************************\n";
	print "ERROR in makemf.pl:\n";
	print "You have tried to run makemf.pl on the directory $dir,\n";
	print "but there already exists a makefile.in in this directory\n";
	print "that is NOT automatically created.\n";
	print "makemf.pl cowardly refuses to overwrite and destroy the manual\n";
	print "work that has gone into this makefile.in.\n";
	print "NB: you can delete this makefile.in by hand, and run makemf.pl again,\n";
	print "but only if you know what you are doing!\n";
	print "Script stopping.\n";
	print "************************************************************\n";
	die "\n";
    }
}

### check_deps ###
# Checks if all dependecies (.h files) exist
sub check_deps{
    my $depsref=shift;
    my $incref=shift;
    my $dep;
    my $idir;
    return if ($dir =~ /templates/);
    foreach $dep(@$depsref) {
        $found=0;
	foreach $idir (@$incref) {
	    if (-f "$idir/$dep") {
#		print "Found file $dep in $idir\n";
		$found++;
	    } elsif (($dep =~ /dsdir\.h/) || ($dep =~ /dsver\.h/)) {
		print "File $dep will be automatically created later.\n";
		$found++;
	    }
	}
	if ($found==0) {
	    push @errorlog,"\nERROR for directory $dir\n";
	    push @errorlog,"Did not find include file $dep in @$incref\n";
	    push @errorlog,"You need to make sure this include file exists.\n";

	}
    }
}

### contain_modules ###
# check if a directory contain fortran code with modules
sub contain_modules{
    my $dir=$_[0];
    my $ret=0;
    my $file;
    foreach $file (<$dir/*.f $dir/*.F $dir/*.f90 $dir/*.F90>) {
	open(TMPIN,"<$file") || die "Can't open $file for reading.\n";
	while(defined($line=<TMPIN>)) {
	    if ($line =~ /^\s*module\s+/) {
		$ret=1;
	    }
	}
	close(TMPIN);
    }
    return $ret; 
}
