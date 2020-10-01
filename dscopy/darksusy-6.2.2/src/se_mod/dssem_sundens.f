***********************************************************************
*** dssem_sundens gives the density in the Sun as a function of radius
*** the radius should be given in m and the density is returned in
*** g/cm^3
*** Density and element mass fractions up to O16 are from the standard
*** solar model BP2000 of Bahcall, Pinsonneault and Basu,
*** ApJ 555 (2001) 990.
*** The mass fractions for heavier elements are from N. Grevesse and
*** A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
*** their total mass fractions matches that of the heavier elements in 
*** the BP2000 model.
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2003-11-25
***********************************************************************

      real*8 function dssem_sundens(r)
      implicit none

      include 'dssem_sun.h'
      real*8 r,rpl
      integer j

      j=0                       ! initialize
      
c...Check if data file is loaded
      call dssem_sunread

      if (r.ge.r_sun) then
        dssem_sundens=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dssem_sundens=sdrho(1)
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20
      
      dssem_sundens=0.0d0
      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dssem_sundens=sdrho(j)*(1.0d0-rpl)+sdrho(j+1)*rpl

      return

      end
