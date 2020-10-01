README for user_replacables in src_models/mssm/
===============================================

In principle, all functions/subroutines in DarkSUSY can be replaced by the
user. The way to do this is to create a new version of the function/subroutine
and make sure to link to it before linking to the default DarkSUSY version.

In this directorly, we show how this can be done when replacing routines
in one of the particle physics modules in src_models.
Replacing routines in src/ is handled in the same way.

To replace routines, do the following:

1. Add the routines you want to replace to this directory.

2. As you might not want to replace all the files all the time, also add
   the filenames which include your replaced routines to the file
   'files-to-include.txt'

3. In the DarkSUSY root, configure and make as usual.

4. DarkSUSY will now have created proper makefiles that put your replaced
   routines in a library libds_mssm_user.a and made sure that the main programs
   are created by linking to these files first (see e.g. makefile in
   examples/test). 

To help you out, we provide a script scr/make_replaceable.pl that you run on
the file(s) you wish to create user-replaceable functions of. The script will
perfom step 1-2 above. It also has an option to delete a user-replaceable
function that will undo steps 1-2 above.

