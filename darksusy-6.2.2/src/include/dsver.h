*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsver.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsver.h              ***
c----------------------------------------------------------------------c
*** This file is created by config2.pl on Mon Oct 26 14:19:49 CET 2020

      character*66 dsver
      parameter(dsver='darksusy-6.2.2 (branch master from Thu Oct 22 1'
     & //'7:35:47 2020 +0200)')
      character*110 dsversion
      common /dsv/dsversion
      save /dsv/
***                                                                  ***
************************** end of dsver.h ***************************
