*                         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsempty.h                             ***
***         this piece of code is needed as a separate file          ***
***              the rest of the code 'includes' dsempty.h           ***
c----------------------------------------------------------------------c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 2014

* For every model, we use the same structure to represent the particle
* code. HOW this is implemented (i.e. which particle codes are assigned)
* is up to the model. 
      include 'dsparticles.h'

      parameter (numpartspecies=0)  ! # particles in this model



***                                                                  ***
************************ end of dsempty.h ******************************
