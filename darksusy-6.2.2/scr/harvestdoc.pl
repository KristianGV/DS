#!/usr/bin/perl -w
#
# Script harvestdoc goes through the DarkSUSY folders and sets together
# the documentation manual as one large latex file. It harvests documentation
# files from different places and also extracts documentation from source
# file headers, either for all of the routines or for the most important ones.
# Author: Joakim Edsjo, edsjo@fysik.su.se
# Date: March 8, 2006 (original headers2tex version)
# Modified: December 18, 2017

# This is the list of subdirectories of src/ from which headers and
# documentation will be read.

$texfile="docs/Manual.tex";
$texlinelimit=2; # total number of lines in tex files to require to include
                 # this directory in the manual
#$srcfile="docs/src-dirs-to-include.txt.auto"; # no longer used

$include_headers=0; # don't include all routine headers by default
$debug=0;
foreach $arg (@ARGV) {
    if ($arg =~ /\-\-long/) {
	$include_headers=1;
	$texfile =~ s/\.tex/-long.tex/;
	$texlinelimit=1; # no need for non-debug version here
    } 
    if ($arg =~ /\-\-debug/) {
	$debug=1;
	$texfile =~ s/\.tex/-debug.tex/;
	$texlinelimit=1; # Require only one line (name file)
    }
}
#if (($include_headers==1) && ($debug ==1)) {
#    print "Option --debug does not have any meaning when option --long is used.\n";
#    print "Will use only option --long.\n";
#    $debug=0;
#    $texfile =~ s/-debug//;
#}

# Now figure out which directories to include documentation from
@dirlist=();
$number_of_modules=0; # number of modules found
$srcm=0;
foreach $dir (<src/* src_models/*{}/* examples examples/test examples/aux>) {
    next if (not(-d $dir));
    next if ($dir =~ 'user_replaceable');
    # check how many tex lines we have (i.e. if we have anything to include)
    $texlines=0;
    foreach $texfile (<$dir/*.tex>) {
	$line=`wc -l $texfile`;
	$line =~ s/^\s+//;
	($nlines) = split(/\s+/,$line);
	$texlines += $nlines;	
    }
    # Determine number of source files
    $nsrcfiles=0;
    foreach $file (<$dir/*.f $dir/*.F $dir/*.f90 $dir/*.F90>) {
	$nsrcfiles++;
    }
    next if (($include_headers==0) && ($texlines < $texlinelimit));
    next if (($include_headers==1) && ($texlines < $texlinelimit) && ($nsrcfiles==0));
    $srcm++ if ($dir =~ /src\_models/);
    # Check if first directory in src_models, then add common first
    # do we need this? How do we make sure 
#    if ($srcm==1) {
#	push(@dirlist,"src_models/common/docs");
#    }
    $modname=$dir;
    $modname =~ s/^.*src_models\///;
    $modname =~ s/\/.*$//;

# Check if first instance of a particle physics module, then add docs first
    if ($dir =~ 'src_models') {
	$modno{$modname}++;
	if ($modno{$modname}==1) { # New module found
	    $number_of_modules++;
	    $modno2tag{$number_of_modules}=$modname;
	    $modtag2no{$modname}=$number_of_modules;
	    $modtag2fullname{$modname}=mod_full_name($modname);
	    push(@dirlist,"src_models/$modname/docs");
	}
    }

    
# docs directories taken care of manually, don't include them from file glob
    next if ($dir =~ /docs/);

    push(@dirlist,$dir);
}

# Now add chapter names for all chapters
foreach $dir(@dirlist) {
    $filefound=1;
    open(TI,"<$dir/DOC-name.tex") || ($filefound=0);
    if ($filefound==1) {
	$name=<TI>;
	chomp($name);
	close(TI);
    } else { # Use module path as chapter name
	$name = $dir;
	$name =~ s/_/\\_/g;
    }
    $ch_names{$dir}=$name;
    print "$dir $name\n";
}

## The code below is the old obsolete code that read the above from a file
## instead of creating it on the fly
# Now read in which directories to include and their names
#open(IN,$srcfile) || die "Can't open $srcfile for reading.\n";
#while(defined($line=<IN>)) {
#    next if $line =~ /^#/;
#    ($dir)=split(/\s/,$line);
#    $name = $line;
#    $name =~ s/^\S+\s//;
#    chomp($name);
#    push(@dirlist,$dir);
#    $ch_names{$dir}=$name;
#}


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

print_texheader($dsver,$date,$include_headers,$debug);

print_part("Prelude");

foreach $file (<docs/I*tex>) {
    include_texfile($file,$debug);
}


# Start going through the directories and files

$root=`pwd`;
chomp($root);

# Include chapter definition for src and general text for src structure

print_part("Main DarkSUSY routines in src/");

#include_texfile("docs/Src-Intro.tex"); Now above in general I-file instead.

$model="";
foreach $dir (@dirlist) {    
    chdir($root) || die "Can't cd to root directory $root\n";
    chdir($dir) || die "Can't cd to $dir\n";
    print "Directory: $dir\n";
    $newtex=1;
    $newfortran=1;
    if ($dir =~ /src_models\/common\/docs/) { # Now we reach the modules in src_models
	    print_part("Particle physics modules in src\\_models");
    } 
    if ($dir =~ /^examples$/) { # Now we reach the examples
	    print_part("Example programs");
    } 
    if ($dir =~ /docs$/) {
	$model = $dir;
	$model =~ s#src_models/##;
	$model =~ s#/docs##;
	print "Documentation directory, using tex files directly.\n";
    } else {
	if (($dir =~ /src_models/)) { # create sections, otherwise chapters
           print_texdir_sec($dir,$ch_names{$dir});
	} else {
           print_texdir_chap($dir,$ch_names{$dir});
	}
    }
    if ($include_headers){
	@files=<*.tex *.f *.F *.f90 *.F90>;
    } else {
	@files=<*.tex>;
    }
    foreach $file (@files) {
	next if ($file =~ /DOC\-name\.tex/); # this is already taken care of
        print "   File: $file\n";
        if ($newtex and $file =~ /\.tex/) {
#            print TEX "\\section{Overview (from tex files)}\n";
            $newtex=0;
	}
        if ($newfortran and $file =~ /\.(f|F|f90|F90)/) {
	    if (($dir =~ /src_models/)) {
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
    include_texfile($file,$debug);
}

# Now insert the rest of the files at the end
foreach $file (<docs/E*tex>) {
    include_texfile($file,$debug);
}


print_texend();

# Copy figures to the docs/fig/ directory (not needed anymore)
# system ('mkdir docs/fig') unless -d 'docs/fig';
# system ('cp -p docs/fig/* docs/fig/');

print "Done. Output written to $texfile\n";
exit;

##########
sub print_texheader {
    my $ver=$_[0];   # DarkSUSY version
    my $date=$_[1];
    my $inc_head=$_[2];
    my $debug=$_[3];

    $template="./docs/Template.tex";
    open(IN,$template) or die "Can't open $template for reading.\n";
    if ($inc_head) {
	$title="Manual and description of routines";
	$subtitle="with full routine headers";
    } else {
	$title="Manual and description of routines";
	$subtitle="";
    }
    while(defined($line=<IN>)) {
        $line=~ s/\[Date\]/$date/;
        $line=~ s/\[DarkSUSYVersion\]/$ver/;
        $line=~ s/\[Title\]/$title/;
        $line=~ s/\[SubTitle\]/$subtitle/;
        if ($line=~ /\[Macros\]/) {
            foreach $file (<docs/headers/*.tex>) {
                include_texfile($file,$debug);
            }
	    next;
	}
        if ($line=~ /\[HeaderFiles\]/) {
            foreach $file (<docs/H*.tex>) {
                include_texfile($file,$debug);
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
    $printdir =~ s/src_models\///;
    $printdir =~ s/_/\\_/g;

print TEX <<END;

% TB commented out to avoid too many almost empty pages
% \\newpage
\\section\[$printdir: $name\]{\\codeb{$printdir}:\\\\ $name}
\\label{sec:$dir}

END
}

##########
sub print_texdir_chap{
    my $dir=$_[0];
    my $name=$_[1];

    $printdir = $dir;
    $printdir =~ s/src\///;
    $printdir =~ s/_/\\_/g;

print TEX <<END;

\\chapter\[$printdir: $name\]{\\codeb{$printdir}:\\\\ $name}
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
    my $name;
    $tfile =~ s#\_#\\\_#g;
    $name=$tfile;
    $name =~ s/^.*\///;
    $name =~ s/\.\S+$//;

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
	    # Check for tags to replace
	    if ($line =~ /\[IFLIST\]/) {
		print_iflist();
		$line="";
	    }
	    if ($line =~ /\[IFHEADERS\]/) {
		print_ifheaders();
		$line="";
	    }
	    if ($line =~ /\[CULIST\]/) {
		print_culist();
		$line="";
	    }
	    if ($line =~ /\[CUHEADERS\]/) {
		print_cuheaders();
		$line="";
	    }
	    if ($line =~ /\[MODULESLIST\]/) {
		print_moduleslist();
		$line="";
	    }
	    if ($line =~ /\[SRCLIST\]/) {
		print_srclist();
		$line="";
	    }
	    print TEX $line;
	}
        return;
    }


print TEX <<END;

%%%%% routine $file %%%%%
\\begin{routine}{$tfile}
\\index[routines]{$name}
END
    if ($verb) {
	print TEX "{\\footnotesize\\begin{verbatim}\n";
    }

# Now go through the file

    open(IN,"$file") || die "Can't open $file for reading.\n";
    $i=0;
    $tex=0;
    $tex=1 if $type eq "latex";
    $clean=0; # 0=print line as is, 1=take away leading comment signs and spaces
    $defcont=0; # continuation line for subroutine / function found or not
    while(defined($line=<IN>)) {
	$line =~ s/^\t/      /; # replace tabs with spaces in beginning of line
        if (substr($line,0,1) ne " " or $type eq "text") {
            if ($line =~ /BeginTex/i) {
		print TEX "\\end{verbatim}}\n" if ($verb);
                $tex=1;
                next;
	    }
            if ($line =~ /EndTex/i) {
		print TEX "{\\footnotesize\\begin{verbatim}\n" if ($verb);
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
	print TEX "\\end{verbatim}}\n \\end{routine}\n";
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
    my $debug=$_[1];
    print "Including tex file: $file\n";
    print TEX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    print TEX "%%% Here comes $file %%%\n";
    print TEX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    open (IN2,"$file") || die "Can't open $file for reading.\n";
    while (defined($line=<IN2>)) {
        # Check for tags to replace
	if ($line =~ /\[IFLIST\]/) {
	    print_iflist();
	    $line="";
	}
	if ($line =~ /\[IFHEADERS\]/) {
	    print_ifheaders();
	    $line="";
	}
	if ($line =~ /\[CULIST\]/) {
	    print_culist();
	    $line="";
	}
	if ($line =~ /\[CUHEADERS\]/) {
	    print_cuheaders();
	    $line="";
	}
	if ($line =~ /\[MODULESLIST\]/) {
	    print_moduleslist();
	    $line="";
	}
	if ($line =~ /\[SRCLIST\]/) {
	    print_srclist();
	    $line="";
	}
	if (($line =~ /newcommand\{\\comment\}/) && ($debug == 0)) {
	    $line=<IN2>; # read one more line
	    $line="\\newcommand{\\comment}[1]{}\n";
	}
	if (($line =~ /newcommand\{\\mcomment\}/) && ($debug == 0)) {
	    $line=<IN2>; # read one more line
	    $line="\\newcommand{\\mcomment}[1]{}\n";
	}
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

##########
# print_ifheaders will look through the empty model and extract all interface
# functions from it and print their headers
sub print_ifheaders{
    my $iline="";
    foreach $ifile (<src_models/empty/*/*.{f,F,f90,F90,c}>) {
	if (interface($ifile))  {
	    print_texfile($ifile);
	}
    }
}

##########
# print_cuheaders will look through src and extract all commonly used
# functions from it and print their headers
sub print_cuheaders{
    my $iline="";
    foreach $ifile (<src/*/*.{f,F,f90,F90,c}>) {
	if (commonly_used($ifile))  {
	    print_texfile($ifile);
	}
    }
}

##########
# print_iflist will look through the empty model and extract all interface
# functions from it
sub print_iflist{
    my $iline="";
    my $desc;
    my $tfile;
    my $mname;
    my $ii;
    my $modappear;
    my $modsin;
    my $whused;
    my @lines;
    my $tline;
    $iline = "\\begin{tabular}{llp{3cm}p{6cm}}\n";
    $iline .= " & {\\bfseries Appear in} & {\\bfseries } & {\\bfseries } \\\\\n";
    $iline .= "{\\bfseries Routine} & {\\bfseries modules} & {\\bfseries Used by} & {\\bfseries Description} \\\\\n";
    $iline .= "\\hline\n";
    @lines=(); # fill temporary array to sort it later
    foreach $ifile (<src_models/empty/*/*.{f,F,f90,F90,c}>) {
	if (interface($ifile))  {
	    $name = $ifile;
	    $name =~ s/^.*\///;
	    $name =~ s/\_/\\\_/g;
	    $name =~ s/\.f//;
	    $desc=get_desc($ifile);
	    # Now find which modules it exists in
	    $ii=0;
	    %modappear=(); # clear hash
	    $modsin=""; # clear list
	    foreach $tfile (<src_models/*/*/$name.{f,F,f90,F90,c,cc}>) {
		next if $tfile=~ /user_replaceables/;
		$mname=$tfile;
		$mname =~ s/src_models\///;
		$mname =~ s/\/.*//;
		$modappear{$mname}++;
		next if ($modappear{$mname}>1); # catch case-insensitive system
		$ii++;
		if ($ii==1) {
		    $modsin = $modtag2no{$mname};
		} else {
		    $modsin .= ", $modtag2no{$mname}";
		}
	    }
	    # Now find where it is used in src
	    $whused=srcdir_where_used($name);
	    if ($ii>0) {
		push (@lines,"\\ftb{$name} & $modsin & $whused & $desc \\\\\[0.5ex\]\n");
	    }
#	    $iline .= "\\ftb{$name} & $modsin & $whused & $desc \\\\\[0.5ex\]\n";
	}
    }
    foreach $tline (sort @lines) {
	$iline .= "$tline";
    }
    $iline .= "\\hline\n";
    $iline .= "\\end{tabular}\n \\smallskip\n";
    print TEX $iline;
}

##########
# print_iflist will look through the empty model and extract all interface
# functions from it
sub print_iflist_old{
    my $iline="";
    my $desc="Routine description";
    $iline = "\\begin\{brief-subs\}\n";
    foreach $ifile (<src_models/empty/*/*.{f,F,f90,F90,c}>) {
	if (interface($ifile))  {
	    $name = $ifile;
	    $name =~ s/^.*\///;
	    $name =~ s/\_/\\\_/g;
	    $name =~ s/\.f//;
	    $iline .= "\\bsub\{$name\} \\index\{$name\}\n";
# Get description of routine from file.
	    $desc=get_desc($ifile);
	    $iline .= "$desc\n";
	}
    }
    $iline .= "\\end\{brief-subs\}\n";
    print TEX $iline;
}

##########
# print_culist will look through src and extract all commonly used
# functions from it
sub print_culist_old{
    my $iline="";
    my $desc="Routine description";
    $iline = "\\begin\{brief-subs\}\n";
    foreach $ifile (<src/*/*.{f,F,f90,F90,c}>) {
	if (commonly_used($ifile))  {
	    $dname = $ifile;
            $dname =~ s/\/ds.*$//;
	    $name = $ifile;
	    $name =~ s/^.*\///;
	    $name =~ s/\_/\\\_/g;
	    $name =~ s/\.f//;
	    $iline .= "\\bsub\{$name\} \\index\{$name\}\n";
# Get description of routine from file.
	    $desc=get_desc($ifile);
            $desc .= " (see Chapter \\ref{ch:$dname})";		    
	    $iline .= "$desc\n";
	}
    }
    $iline .= "\\end\{brief-subs\}\n";
    print TEX $iline;
}

##########
# print_culist will look through src and extract all commonly used
# functions from it
sub print_culist{
    my $iline="";
    my $desc="Routine description";
    $iline = "\\begin{tabular}{lp{8.5cm}l}\n";
    $iline .= "{\\bfseries Routine} & {\\bfseries Description} & {\\bfseries Chapter} \\\\\n";
    $iline .= "\\hline\n";
    foreach $ifile (<src/*/*.{f,F,f90,F90,c}>) {
	if (commonly_used($ifile))  {
	    $dname = $ifile;
            $dname =~ s/\/ds.*$//;
	    $name = $ifile;
	    $name =~ s/^.*\///;
	    $name =~ s/\_/\\\_/g;
	    $name =~ s/\.f//;
	    $iline .= "\\ftb{$name} \\index[routines]{$name} & ";	    
# Get description of routine from file.
	    $desc=get_desc($ifile);
	    $iline .= "$desc & \\ref{ch:$dname} \\\\\n";	    
	}
    }
    $iline .= "\\hline\n";
    $iline .= "\\end{tabular}\n";
    print TEX $iline;
}


# sub interface
# Returns 1 if interface function, otherwise 0
sub interface{
    my $file = shift @_;
    open(INI,"<$file") || die "Can't open $file for reading.\n";
    $found=0;
    while (defined($line=<INI>)) {
	$found=1 if ($line =~ /^\S+\s*type\s*:\s*interface/i);
    }
    close(INI);
    return $found;
}

# sub commonly_used
# Returns 1 if commonly used function, otherwise 0
sub commonly_used{
    my $file = shift @_;
    open(INI,"<$file") || die "Can't open $file for reading.\n";
    $found=0;
    while (defined($line=<INI>)) {
	$found=1 if ($line =~ /^\S+\s*type\s*:\s*common/i);
    }
    close(INI);
    return $found;
}

# sub get_desc
# get description of routine from header, will search for desc : descriptor
sub get_desc{
    my $dfile=shift @_;
    my $ddesc="";
    my $tmp;
    open(IND,"<$dfile") || die "Can't open $dfile for reading.\n";
    while (defined($line=<IND>)) {
	if ($line =~ /^\S+\s*desc\s*:(.*)$/i) {
	    $tmp = $1;
	    $tmp =~ s/\s+\*+$//;
	    $ddesc .= $tmp;
	}
    }
    close(IND);
    return $ddesc;
}

# sub mod_full_name
# this subroutine take the module short name (directory name) and extracts
# the title from the chapter definition in the documentation file
sub mod_full_name{
    my $modtag=shift;
    my $fullname="";
    my $ffile;
    my $fline;
    foreach $ffile (<src_models/$modtag/docs/I01*.tex>) {
	open(INF,"<$ffile") || die "Can't open $ffile for reading.\n";
	while ($fline=<INF>) {
	    if ($fline=~ /\\chapter\{(.*)\}/) {
		$fullname=$1;
	    }
	}
	close(INF);
    }
    return $fullname;
}

# sub src_full_name
# this subroutine take the src (directory path) and extracts
# the title from the chapter definition in the documentation file
sub src_full_name{
    my $srcdir=shift;
    my $fullname="";
    my $ffile;
    my $fline="";
    my $srcdirprint;
    $srcdirprint=$srcdir;
    $srcdirprint =~ s/\_/\\\_/g;
    $ffile="$srcdir/DOC-name.tex";
    $found=1;
    open(INF,"<$ffile") || ($found=0);
    if ($found == 1 ) {
	while ($fline=<INF>) {
	    $fullname .= $fline;
	}
	close(INF);
    } else {
	$fullname = "\\comment{Missing file: $srcdirprint/DOC-name.tex}";
    }
    return $fullname;
}

##########
# print_moduleslist will print a table with the modules currently available
sub print_moduleslist{
    my $iline="";
    my $ii=0;
    my $tag;
    my $mname;
    $iline = "\\begin{tabular}{llp{8cm}}\n";
    $iline .= "\\multicolumn{3}{l}{\\bfseries Particle physics modules} \\\\\n";
    $iline .= "\\cline{1-3}\n";
    $iline .= "{\\bfseries Module} & {\\bfseries Short} & {\\bfseries } \\\\\n";
    $iline .= "{\\bfseries No.} & {\\bfseries name} & {\\bfseries Description} \\\\\n";
    $iline .= "\\hline\n";
    for ($ii=1; $ii<=$number_of_modules; $ii++) {
	$tag=$modno2tag{$ii};
	$mname=$modtag2fullname{$tag};
        $tag =~ s/\_/\\\_/g; # for LaTeX
        $iline .= "$ii & \\ftb{$tag} & $mname \\\\\[0.5ex\]\n";
    }
    $iline .= "\\hline\n";
    $iline .= "\\end{tabular}\n";
    print TEX $iline;
}


##########
# print_srclist will print a table with the folders currently available in src
sub print_srclist{
    my $iline="";
    my $dname;
    my $dtag;
    my $dtitle;
    $iline = "\\begin{tabular}{lp{10cm}}\n";
    $iline .= "{\\bfseries Subdirectory} & {\\bfseries } \\\\\n";
    $iline .= "{\\bfseries name} & {\\bfseries Description} \\\\\n";
    $iline .= "\\hline\n";
    foreach $dname (<src/*>) {
	next if $dname =~ /include/;
	next if (not(-d $dname));
	next if $dname =~ /user_replaceables/;
	next if $dname =~ /templates/;
	$dtag=$dname;
	$dtag =~ s/src\///;
	$dtag =~ s/\_/\\\_/g;  # for LaTeX
	$dtitle=src_full_name($dname);
        $iline .= "\\ft{$dtag} & $dtitle \\\\\[0.5ex\]\n";
    }
    $iline .= "\\hline\n";
    $iline .= "\\end{tabular}\n";
    print TEX $iline;
}

# sub srcdir_where_used
# check where an interface function is used and return the directory names
# where used
sub srcdir_where_used{
    my $ifname=shift;
    my $res;
    my $srcs="";
    my $srcno;
    my $lline;
    my $ldir;
    my $ldirtex;
    my $ii;
    my $tfile;
    my $found;
    $ii=0;
    foreach $tfile (<src/*/*.{f,F,f90,F90,c,cc}>) {
	open(INS,"<$tfile") || die "Can't open $tfile for reading.\n";
	$found=0;
	while(defined($lline=<INS>)) {
	    $found=1 if $lline=~ /$ifname\(/;
	    $found=1 if $lline=~ /$ifname\,/;
	    $found=1 if $lline=~ /call\s+$ifname/;
	}
	close(INS);
	next if ($found == 0);
	$lline=$tfile;
	$lline =~ /src\/(\w+)\//;
        $ldir=$1;
	next if $ldir =~ /user_replaceables/;
	next if $ldir =~ /include/;
	$ii++ if $found++;
	$ldirtex=$ldir;
	$ldirtex =~ s/\_/\\\_/;
	$srcno{$ldir}++;
	if ($ii==1) {
	    $srcs=$ldirtex;
	    next;
	}
	if ($srcno{$ldir}==1) {
	    $srcs .= ", $ldirtex";
	}
    }

    return $srcs;
}
