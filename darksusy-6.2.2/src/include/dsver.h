*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsver.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsver.h              ***
c----------------------------------------------------------------------c
*** This file is created by config2.pl on Wed Sep 30 14:37:20 CEST 2020

      character*14 dsver
      parameter(dsver='darksusy-6.2.2')
      character*60 dsversion
      common /dsv/dsversion
      save /dsv/
***                                                                  ***
************************** end of dsver.h ***************************
