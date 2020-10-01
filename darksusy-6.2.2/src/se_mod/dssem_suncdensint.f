**********************************************************************
*** dssem_suncdensint gives the column density in the Sun from the
*** centre out the tha given radius r (in meters)
*** if type = 'N', the total column density (up to that r) is calculated
***         = 'p', the column density on protons is calculated
***         = 'n', the column density on neutrons is calculated
*** in this routine, the actual integration is performed. for speed,
*** use dssem_suncdens instead which uses a tabulation of this result.
*** Author: joakim edsjo (edsjo@fysik.su.se)
*** Date: November 24, 2005
**********************************************************************

      real*8 function dssem_suncdensint(r,type)
      implicit none
      include 'dssem_sun.h'

      real*8 r,dssem_suncdfunc,dsf_int
      character*1 type
      external dssem_suncdfunc

      cdt=type ! transfer to common block
c...integrate sun density, the last factor of 100.0 comes from m->cm

      dssem_suncdensint=
     &  dsf_int(dssem_suncdfunc,0.0d0,min(r,r_sun),1.d-3)*100.0d0

      return
      end
