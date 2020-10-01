***********************************************************************
*** auxiliary function for the quick determination of Tkd
*** (integrand for the thermal average). Called by dskdTkd.
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
*** updates: 2014-05-09 (removed model-dependence)
***********************************************************************

      real*8 function dskdcint2(y)
      implicit none

      include 'dskdcom.h'

      real*8 mf,x,y,tmpres, dsmass
      real*8 cn
      integer SMtype,n

      SMtype = 1
      tmpres = 0d0
      dskdcint2=0d0    

      mf=0.0d0
 10   if (SMtype.le.3) mf=0d0       ! neutrinos
      if (SMtype.eq.4) mf=dsmass(11) ! PDG code for electrons
      if (SMtype.eq.5) mf=dsmass(13) ! PDG code for muons
      if (SMtype.eq.6) mf=dsmass(15) ! PDG code for taus
      if (y.lt.mf/Tint/10.) goto 20
      x = sqrt(y**2+(mf/Tint)**2)

      call dskdm2simp(SMtype,cn,n)
      tmpres = tmpres+cn*y**2*x**n/(1.+exp(x))

 20   SMtype = SMtype+1
      if (SMtype.le.6) goto 10


      dskdcint2 = tmpres
      
      return
      end
