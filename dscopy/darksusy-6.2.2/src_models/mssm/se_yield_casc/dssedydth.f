*****************************************************************************
*** function dssedydth is the differential yield dyield/dcostheta which
*** should be integrated (by the routine gadap e.g.).
*** cth is cos(theta).
*** units: 1.0e-15 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dssedydth(cth,m0,m1,m2,e0,eth,
     &  thm,chi,wwh,fk,fv)
      implicit none
      include 'dsanyieldmodelcom.h'

c------------------------ variables ------------------------------------

      real*8 cth
      real*8 e1cm,e1
      integer istat
      character*2 wh
      real*8 m0,m1,m2,e0,eth,thm
      integer chi,wwh,fk,fv   ! phiwh=1 - sun, 2 - earth

c------------------------ functions ------------------------------------
      real*8 dsseyield_sim
c-----------------------------------------------------------------------

      if (wwh.eq.1) then
        wh='su'
      else
        wh='ea'
      endif

c...calculate the energy of particle 1 in the lab system

      e1cm=(m0**2+m1**2-m2**2)/(2.0*m0)
      e1=e0/m0*(e1cm + sqrt(e0**2-m0**2)/e0*
     +  sqrt(e1cm**2-m1**2)*cth)
      dssedydth=1.0d15*1.0d0/4.0d0*
     &     dsseyield_sim(e1,eth,thm,sechi2pdg(chi),'0',
     &      wh,fk,fv,istat) 

      end











