#!/usr/bin/perl -w
#
# Script to create a module by copying another module, either by copying the
# full module (default) or by only keeping interface functions. 
# The script makes sure that the module name is used consistently.
# This script is useful when setting up a new module, both when starting
# more or less from scratch (staring e.g. from the empty or generic wimp
# module) or when using a more complex module as staring point, like the
# mssm module.
#
# This script works for both particle physics modules and halo modules
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: 2015-10-29
# Updated 2016-02-07 to use update_mainmf instead of manually upating 
# makefile.in

if (not(-d "src" and -d "src_models")) {
    print "make_module.pl needs to be run from the DarkSUSY root directory.\n";
    die "Please run it as scr/make_module.pl <old> <new> from the DarkSUSY root.\n";
}

# Check if autoconf is available
if (not(`which autoconf`)) {
    print "ERROR: autoconf is not found on your system. You need to install it\n";
    die "prior to running make_module.pl.\n";
}

if (scalar @ARGV <=1) {
    print <<END;
make_module takes an existing module in src_models or src_halos and
creates a copy out of it.
Usage:
make_module [OPT] <name-of-old-module> <name-of-new-module>

The script will take the old module and create the new one by copying it.
The old module is expected to be in src_models or src_halos (will determine
this automatically) and the new one will be put in the same module directory.
If the option '-i' is given, only interface functions are copied, otherwise
the full module is copied. The -i option is useful to have a cleaner starting
point even if it most likely will not compile as is.
END
die "\n";  
} else {
    if ($ARGV[0] eq "-i") {
	$int=1;
	shift;
	print "Will only copy interface functions...\n";
    } else {
	$int=0;
    }
    $modfrom=shift;
    $modto=shift;
}

$date=`date +"%b %d, %Y"`;
chomp($date);

if (-d "src_models/$modfrom") {
    $moddir="src_models";
} else {
    print "Starting module $modfrom does not exist in src_models.\n"; 
    die "Stopping.\n";
}

if (-d "$moddir/$modto") {
    die "New module $moddir/$modto already exists. Stopping.\n";
}

# Start copying

print "Start copying files from module $moddir/$modfrom to module $moddir/$modto\n\n";

# First set up which include files to copy for interface functions
@files=<$moddir/$modfrom/*{}/*.{f,F,f90,F90}>; # {} to avoid being comment
@includes = interface_includes(@files); # list of include files in int. func.
# Create hash for quick access
%include_hash = map {$_, 1} @includes;

@files=<$moddir/$modfrom/*{}/*.{f,F,f90,F90,c,cc,h,tex,txt}>; # {} to avoid being comment
push(@files,"$moddir/$modfrom/makefile.in");
for $file (@files) {
    print "File: $file";
    $newfile=$file;
    $newfile =~ s/$modfrom/$modto/g;
    $newfile=~ /(.*\/)/;
    $dir=$1;
    $basename=$file;
    $basename =~ s/.*\///;

# Regular source files
    if ($file =~ /\.(f|F|f90|F90|c|cc)/) { # regular file
	if ((($int) && interface($file)) || $int==0) { # copy it
	    system("mkdir -p $dir") if (not(-d "$dir"));
	    system("cp -p $file $newfile");
	    print " ...copied to $newfile\n";
	} else {
	    print " ...not copied\n";
	}
    }

# Header files
    if ($file =~ /\.h/) { # header file
        if ((($int) && $include_hash{$basename}) || $int==0) {
	    system("mkdir -p $dir") if (not(-d "$dir"));
	    system("cp -p $file $newfile");
	    print " ...copied\n";
        } else {
            print " ...not copied\n";
        }	    
    }

# other files that should be copied
    if ($file =~ /\.(in|tex|txt)/) { # regular file
	system("mkdir -p $dir") if (not(-d "$dir"));
	system("cp -p $file $newfile");
        print " ...copied to $newfile\n";
    }
}
print "Done.\n\n";

### Now we have copied all files
### Modify files to change module name
# The module name will be changed both in filenames and in the files

print "Going through files in new module and changing module name...\n";
@files = <$moddir/$modto/*{}/*.{f,F,f90,F90,h,txt} $moddir/$modto/makefile.in>;
for $file (@files) {
    change_name($file,$modfrom,$modto);
}
print "Done.\n\n";

# Update main makefile.in
# This is not needed as this is done in preconfig.pl
#print "Updating main makefile.in... ";
#system("scr/update_mainmf.pl");
#print "done.\n\n";

# Run preconfigure scripts etc
print "Running preconfigure scripts...";
system("scr/preconfig.pl");
print " done\n\n";

print "************************************************************\n";
print "Your new module $modto is now set up.\n";
print "The module is in $moddir/$modto.\n";
print "The main make system is set up to automatically build the new module.\n";
if ($int) {
    print "As you have only included interface functions, it is unlikely that\n";
    print "it will compile without you modifying the module to your liking.\n";
}
print "Now run configure again and make DarkSUSY to make DarkSUSY and your new module.\n";
print "The main programs have not been updated, if you want to run e.g.\n";
print "examples/dsmain_wimp.F on your new module, modify examples/makefile to\n";
print "choose your new module, make dsmain_wimp and run it.\n";
print "Documentation (except txt-files) have NOT been updated.\n";
print "If you are working on this in a subversion repository,\n";
print "no files have been added, you have to do that manually\n";
print "(on e.g. the whole new module directory -- do it before you configure\n";
print "to avoid adding derived files to the repository).\n";
print "Please note that if you create a new module and want to share it\n";
print "send an e-mail to edsjo\@fysik.su.se so that we can link to it\n";
print "from the DarkSUSY web page.\n";
print "************************************************************\n";
if ($moddir =~ /src_models/) {
    print "New particle-physics module created: $moddir/$modto\n";
} else {
    print "New halo module created: $moddir/$modto\n";
}

##################################################
########## Subroutines
##################################################

# Returns 1 if interface function, otherwise 0
sub interface{
    my $file = shift @_;
    open(IN,"<$file") || die "Can't open $file for reading.\n";
    $found=0;
    while (defined($line=<IN>)) {
	$found=1 if ($line =~ /^\S+\s*type\s*:\s*interface/i);
    }
    close(IN);
    return $found;
}


######
# returns a list of all include files in all interface fuctions
# from the list of files supplied
sub interface_includes{
    my @files = @_;
    my @includes=();
    my $file;
    
    for $file (@files) {
	if (interface($file)) {
	    open(IN,"<$file") || die "Can't open $file for reading.\n";
	    while (defined($line=<IN>)) {
		if ($line =~ /^\s+include\s+\'(.+)\'/) {
		    push(@includes,$1);
		}
		if ($line =~ /^#include\s+\"(.+)\"/) {
		    push(@includes,$1);
		}
	    }
	    close(IN);
	}
    }
    return @includes;
}
	

######
# This routine changes module name in the supplied file
# Both instances in the file and in the filename are changed.
sub change_name{
    my $file=shift @_;
    my $oldname=shift @_;
    my $newname=shift @_;
   
    print "Now going through file $file\n";
    $outfile="$file-tmp";
    open(FILE,"<$file") || die "Can't open $file for reading.\n";
    open(OUT,">$outfile") || die "Can't open $outfile for writing.\n";
    $ch=0;
    $i=0;
    while(defined($line=<FILE>)) {
        $i++;
        chomp($line);
	$ch+= ($line =~ s/(\W|^)$oldname(\W|$)/$1$newname$2/gi);
	$ch+= ($line =~ s/ds$oldname/ds$newname/gi);
	$ch+= ($line =~ s/ds_$oldname/ds_$newname/gi);
	$ch+= ($line =~ s/dsgivemodel_$oldname/dsgivemodel_$newname/gi);
#        if (length($line)>72 and substr($line,0,1)=~ /\s/
#           and not(substr($line,0,72) =~ /!/)) {
#	    print "***Line $i is too long\n";
#	}
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
        $ch += ($newfile =~ s/ds$oldname/ds$newname/gi);                       
        if ($ch>0) {
            print "  filename changed to $newfile\n";
	    system("mv $file $newfile");
	}

    } else {
	unlink($outfile);
    }
}
