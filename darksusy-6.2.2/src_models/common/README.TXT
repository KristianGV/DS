This directory src_models/common is a set of routines that *can* be used by all
particle physics modules. It does not contain its own makefile in this
directory as it is instead built by all the modules and included in their
libraries. (The subdirectories here contain makefiles as ususal though.)