*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dsficom.h                              ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsficom.h            ***
c----------------------------------------------------------------------c
c by Kristian Gjestad Vangsnes (kristgva@uio.no), 2020-09-16

* Global parameters
      integer stat
      common /global_fi/ stat
      save /global_fi/

* K1 parameters
      real*8 x1_k1,x2_k1,x3_k1,eta1_k1,eta2_k1,eta3_k1
      common /dsfik1var/ x1_k1,x2_k1,x3_k1,eta1_k1,eta2_k1,eta3_k1
      save /dsfik1var/

* Decay parameters
      real*8 M_dec,eta_dec
      common /dsfidecvar/ M_dec,eta_dec
      save /dsfidecvar/

* Two to two parameters
      real*8 T_22,x1_22,x2_22,xX_22,eta1_22,eta2_22,etaX_22,g1_22,
     &g2_22,c12_22,m1_22,m2_22,mdm
      integer ichannel_22
      common /dsfi2to2_var/ T_22,x1_22,x2_22,xX_22,eta1_22,eta2_22,
     &etaX_22,g1_22,g2_22,c12_22,m1_22,m2_22,mdm,ichannel_22
      save /dsfi2to2_var/