*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                          dsidnumber.h                            ***
*** This common block just stores a unique number identifying each   ***
*** particle physics model.                                          ***      
***                                                                  ***
*** This number should be increased by one whenever a new model is   ***
*** setup, allowing to store model-specific once-and-for-all         ***
*** calculations in a very efficient way. Needed also by some        ***
*** routines in the core library.                                    ***
c----------------------------------------------------------------------c
c  author: Joakim Edsjo, edsjo@fysik.su.se, December 2014
c  mod: torsten.bringmann@fys.uio.no, 10/10/2017


* id number
      integer dsidno, checksum
      common/dsidno/ dsidno, checksum
c save common blocks
      save /dsidno/
***                                                                 ***
********************** end of dsidnumber.h ****************************
