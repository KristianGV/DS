*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dskdcom.h                              ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dskdcom.h            ***
c----------------------------------------------------------------------c
c by Torsten Bringmann (troms@physto.se), 2010-01-23
c updates: 2013-06-11 (removed model-dependence)

* MH switches
      real*8 mheps
      common /KDswitches/ mheps
      save /KDswitches/

* needed for the integration of the Boltzmann equation
      real*8 resk(9),Tint, m0 ! Tinr is the temperature of the scattering partners! 
      integer nKDres
      common /KDresonances/ resk,Tint,m0, nKDres
      save /KDresonances/

* non-SM scattering partners
      integer nBSMscatt,nBSMscattlight
      logical BSMscattfermion(100) ! scattering partner is fermion
      common/KDBSM/ nBSMscatt,nBSMscattlight,BSMscattfermion
      save /KDBSM/

* degrees of freedom
      integer kkhi,kklo
      common /KDdof/ kkhi,kklo
      save /KDdof/

* QCD phase transition
      real*8 tqcdmin,tqcdmax
      parameter (tqcdmin=0.154,tqcdmax=4*tqcdmin) ! in GeV
      integer quark_how
      common /KDqcd/ quark_how
      save /KDqcd/
