#!/usr/bin/perl -w
#
# Script to take any source file and make a user_replaceable file out of it.
# It will
#    - copy the file to the correct location
#    - add some statements and comments that it needs to be modified, 
#    - update makefile.in to include the new file
# There is also an option (-d) to delete a user_replaceable file, which will
#    - delete the file
#    - update makefile.in to not include the deleted file
#
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: October 29, 2015

#use Cwd;
#$pwd=cwd();

if (not(-d "src" and -d "src_models")) {
    print "make_replaceable.pl needs to be run from the DarkSUSY root directory.\n";
    die "Please run it as scr/make_replaceable.pl [files] from the DarkSUSY root.\n";
}

if (scalar @ARGV == 0) {
    print <<END;
make_replaceable takes any source file and makes a user replaceable file
out of it.
Usage:
make_replaceable [OPT] file(s)
If no option is given, make_replaceable will 
  - copy the file(s) to the correct location (the folder user_replaceables
    in the given module)					 
  - add some statements and comments that it needs to be modified, 
  - update makefile.in to include the new file
  Note: you need to add the file to files-to-include.txt in the same directory
  as the file to make it compile and be included in the library. This you do
  by hand.
If the option '-d' is given, make_replaceable will
  - delete the given user replaceable file(s)
  - update makefile.in to not include the deleted file(s)
END
die "\n";  
} else {
    if ($ARGV[0] eq "-d") {
	$delete=1;
	shift;
	print "Will delete user replaceable files...\n";
    } else {
	$delete=0;
    }
    @files=@ARGV;
}

$date=`date +"%b %d, %Y"`;
chomp($date);

foreach $file (@files) {
    # Figure out directory and filenames
    $fname=$file;
    $fname =~ s/^.*\///;
    $base=substr($file,0,length($file)-length($fname)-1);
    $dir=$base;
    $dir =~ s/^.*\///;
    $base=substr($base,0,length($base)-length($dir)-1);
    $userfile="$base/user_replaceables/$fname";
    print "File: $file\n";
#    print "File name: $fname\n";
#    print "Base directory: $base\n";
#    print "User replaceable file: $userfile\n";
    if ($delete==0) { # create user replaceable file
	if (-f $userfile) {
	    print "************************************************************\n";
	    print "ERROR: User replaceable file already exists: \n  $userfile\n";
	    print "Please remove it if you want to create it from source again.\n";
	    print "Stopping.\n";
	    print "************************************************************\n";
	    die "\n";
	}
	print "Creating file $userfile\n";
	create_replaceable($file,$userfile);
	print "Running makemf.pl to update makefile.in\n";
	system("scr/makemf.pl $base/user_replaceables");

	# Check if the file should be added to files-to-include
	# First delete it if it is there
	open(IN,"<$base/user_replaceables/files-to-include.txt") || die "Can't opne files-to-include.txt for reading.\n";
	open(OUT,">$base/user_replaceables/files-to-include.txt.new") || die "Can't opne files-to-include.txt.new for writing.\n";
	$ch=0;
	while(defined($line=<IN>)) {
	    $ch +=($line =~s/ $fname//);
	    $ch +=($line =~s/$fname //);
	    $ch +=($line =~s/$fname//);
	    print OUT "$line";
	}
	close(IN);
	close(OUT);
	if ($ch>0) {
	    rename("$base/user_replaceables/files-to-include.txt.new","$base/user_replaceables/files-to-include.txt");
	} else {
	    unlink("$base/user_replaceables/files-to-include.txt.new");
	}

        # Now add it
	open(IN,"<$base/user_replaceables/files-to-include.txt") || die "Can't opne files-to-include.txt for reading.\n";
	open(OUT,">$base/user_replaceables/files-to-include.txt.new") || die "Can't opne files-to-include.txt.new for writing.\n";
	while(defined($line=<IN>)) {
	    if ($line=~ /^SRC/) {
		while ($line =~ /\\/) {
		    print OUT $line;
		    $line=<IN>;
		}
		chomp($line);
		$line = "$line $fname";
		print OUT "$line\n";
	    }
	}
	close(IN);
	close(OUT);
	rename("$base/user_replaceables/files-to-include.txt.new","$base/user_replaceables/files-to-include.txt");
	
	
	print "**********************************************************************\n";
	print "File created: $userfile\n";
	print "The file is also added to the list in\n";
	print "$base/user_replaceables/files-to-include.txt\n";
	print "You now need to run configure, modify your file and compile.\n";
	print "**********************************************************************\n";
    } elsif ($delete==1) { # delete file
	if (not($file =~ /user_replaceable/)) {
	    print "************************************************************\n";
	    print "ERROR: You have run make_replaceable.pl and asked it to delete the file\n";
	    print "$file\n";
	    print "but this file is NOT in a user_replaceable directory. I will not delete it.\n";
	    print "Script stopping.\n";
	    print "************************************************************\n";
	    die "\n";
	}
	if (not(-f $file)) {
	    print "File $file does not exist.\n";
	    die "Stopping.\n";
	}
	unlink($file);
	print "File deleted.\n";
	print "Running makemf.pl to update makefile.in\n";
	system("scr/makemf.pl $base/user_replaceables");

	# Check if the file should be deleted from files-to-include
	open(IN,"<$base/user_replaceables/files-to-include.txt") || die "Can't opne files-to-include.txt for reading.\n";
	open(OUT,">$base/user_replaceables/files-to-include.txt.new") || die "Can't opne files-to-include.txt.new for writing.\n";
	$ch=0;
	while(defined($line=<IN>)) {
	    $ch +=($line =~s/ $fname//);
	    $ch +=($line =~s/$fname //);
	    $ch +=($line =~s/$fname//);
	    print OUT "$line";
	}
	close(IN);
	close(OUT);
	if ($ch>0) {
	    rename("$base/user_replaceables/files-to-include.txt.new","$base/user_replaceables/files-to-include.txt");
	} else {
	    unlink("$base/user_replaceables/files-to-include.txt.new");
	}
	    
	print "**********************************************************************\n";
	print "File deleted: $file\n";
	print "File is also deleted from $base/user_replaceables/files-to-include.txt\n" if $ch>0;
	print "You now need to run configure.\n";
	print "**********************************************************************\n";
    }
}

############################################################
### Subroutines
############################################################

sub create_replaceable{
    my ($srcfile) = shift @_;
    my ($usfile)  = shift @_;
    open(IN,"<$srcfile") || die "Can't open $srcfile for reading.";
    open(OUT,">$usfile") || die "Can't open $usfile for writing.";
    print OUT<<END;
**********************************************************************
*** This file is automatically generated from the file 
*** $srcfile
*** with the script scr/make_replaceable.pl on $date.
*** The file is copied as is, but to be of any use you should of
*** course modify this file to your liking. The way the default linking
*** is set up, this file will be linked to before the corresponding
*** DarkSUSY file, meaning that this is the file that will be used,
*** i.e. it will replace the default DarkSUSY one.
*** A few lines of code are added to the executable section below
*** to remind you that you need to change this file. Delete those lines
*** (and these comments) and modify the file to your liking and you
*** are ready to go!
**********************************************************************

END

    $printed=0; # set to 1 one code part is added
    while(defined($line=<IN>)) {
	if (line_is_code($line) and $printed == 0) {
	    print OUT "\n**********************************************************************\n";
            print OUT "*** This comment and the following few lines should be deleted when you\n";
	    print OUT "*** modified this routine to your liking.\n";
	    print OUT "      write(*,*) '*****************************************'\n";
	    print OUT "      write(*,*) 'ERROR in $usfile'\n";
            print OUT "      write(*,*) 'You have linked to this user replaceable file that'\n";
            print OUT "      write(*,*) 'is just a dummy file created by make_replaceable.pl.'\n";
	    print OUT "      write(*,*) 'You need to modify this user replaceable file'\n";
	    print OUT "      write(*,*) 'to your liking.'\n";
	    print OUT "      write(*,*) 'Note that you can also modify'\n";
            print OUT "      write(*,*) '$base/user_replaceables/files-to-include.txt'\n";
	    print OUT "      write(*,*) 'to determine if a file should be compiled and linked or not.'\n";
	    print OUT "      write(*,*) 'Program stopping.'\n";
	    print OUT "      write(*,*) '*****************************************'\n";	    print OUT "      stop\n";
	    print OUT "**********************************************************************\n\n";
	    $printed=1;
	}
	print OUT $line;
    }
    close(IN);
    close(OUT);
}

# sub line_is_code determines if a row from a Fortran source file is part
# of the header and declaration or code
sub line_is_code{
    my ($row) = shift @_;
    return 0 if $row =~ /^\S/;
    return 0 if $row =~ /^\s*real/i;
    return 0 if $row =~ /^\s*integer/i;
    return 0 if $row =~ /^\s*double/i;
    return 0 if $row =~ /^\s*common/i;
    return 0 if $row =~ /^\s*character/i;
    return 0 if $row =~ /^\s*complex/i;
    return 0 if $row =~ /^\s*logical/i;
    return 0 if $row =~ /^\s*implicit/i;
    return 0 if $row =~ /^\s*include/i;
    return 0 if $row =~ /^\s*parameter/i;
    return 0 if $row =~ /^\s*data/i;
    return 0 if $row =~ /^\s*external/i;
    return 0 if $row =~ /^\s*save/i;
    return 0 if $row =~ /^\s*function/i;
    return 0 if $row =~ /^\s*subroutine/i;
    return 0 if $row =~ /^\s*equivalence/i;
    return 0 if $row =~ /^     \&/;
    return 0 if $row =~ /^\s*$/;
    return 1;
}

 
    
