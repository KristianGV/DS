***********************************************************************
*** The function dskdmcut returns the typical mass scale of the smallest
*** gravitationally bound objects [in units of M_sun]
***
***  type : commonly used
***  desc : Cutoff mass in linear power spectrum
***
***  input: m0  - DM mass [in GeV]
***         tkd - kinetic decoupling temperature [in MeV]
***         how = 1,2,3
***               1 - maximum of 2,3 [default]
***               2 - horizon mass at T_kd
***               3 - cutoff scale associated to free-streaming
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
*** updates: 2013-06-12 (removed model-dependence)
***********************************************************************

      real*8 function dskdmcut(m0,tkd,how)
      implicit none

      integer how
      real*8  m0,tkd

      real*8 tmpres1,tmpres2,geff

      tmpres1=0d0
      tmpres2=0d0

c... this implements (13,14) in 0903.0189

      call dskdgeff(tkd/1.d3,geff)
      
      tmpres1=3.38D-6*(50.d0/tkd/geff**0.25d0)**3      ! DAO
      
      tmpres2=2.896d-6*(1.d0+log(geff**0.25*tkd/50.d0)/19.075d0)**3  !FS
     -                /(m0/100.d0*geff**0.5d0*tkd/50.d0)**1.5d0

      if (how.eq.2) then
        dskdmcut=tmpres1
      elseif (how.eq.3) then
        dskdmcut=tmpres2
      else    ! default
        if (tmpres1.ge.tmpres2) dskdmcut=tmpres1
        if (tmpres1.lt.tmpres2) dskdmcut=tmpres2
      endif

      return
      end





        

