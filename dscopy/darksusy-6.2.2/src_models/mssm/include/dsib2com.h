*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dsib2com.h                              ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsib2com.h            ***
c----------------------------------------------------------------------c


c... IB2 switches
      integer ib2limit
      parameter (ib2limit=300)
      integer IB2flag(1:612), ib2svhow, IB2sfgen
      character*10 ib2dnhow
      real*8 IB2acc, intres
      common /IB2switches/ IB2acc, intres, IB2flag, ib2svhow, IB2sfgen, ib2dnhow

c... external particle masses
      real*8 mx, MB, Mf, Mff
      common /masses/ mx, MB, Mf, Mff  

c... BB system kinematical quantities defined in dsib2kinematics
      real*8  EvJ, E1J, E2J, kv, kJ,
     &     CW, SW, Epl, Emi, ppl, pmi
      common /kin/ EvJ, E1J, E2J, kv, kJ,
     &     CW, SW,Epl, Emi, ppl, pmi

c...  necessary kinematical data transferred via common block
c...  use only for integration routines
      integer c_ftype, c_fbartype, c_Btype, c_yieldk 
      real*8 c_EB, c_E1, c_E2, c_xstable
      character*5 c_pfinal    ! particle for second/outer energy integration 
      common /int/ c_EB, c_E1,c_E2,c_xstable,
     &     c_ftype, c_fbartype, c_Btype, c_yieldk, c_pfinal

c... helicity amplitudes
      complex*16 amp(-1:1, 0:3)
      common /ampl/ amp
