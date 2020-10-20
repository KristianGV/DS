*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsver.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsver.h              ***
c----------------------------------------------------------------------c
*** This file is created by config2.pl on Tue Oct 20 16:51:25 CEST 2020

      character*65 dsver
      parameter(dsver='darksusy-6.2.2 (branch master from Thu Oct 8 15'
     & //':39:04 2020 +0200)')
      character*110 dsversion
      common /dsv/dsversion
      save /dsv/
***                                                                  ***
************************** end of dsver.h ***************************
