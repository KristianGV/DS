***********************************************************************
*** dssem_sunmass gives the mass of the Sun as a function of radius
*** the radius should be given in m and the mass is given in kg
*** up to the specified radius.
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

      real*8 function dssem_sunmass(r)
      implicit none

      include 'dssem_sun.h'
      real*8 r,rpl
      integer j

      j=0                       ! initialize
      
c...Check if data file is loaded
      call dssem_sunread

      if (r.gt.r_sun) then
        dssem_sunmass=m_sun
        return
      endif

      if (r.ge.0.0d0.and.r.lt.sdr(2)*r_sun) then ! do a better interpolation
        dssem_sunmass=sdm(2)*m_sun*r**3/(sdr(2)*r_sun)**3
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      dssem_sunmass=m_sun
      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dssem_sunmass=(sdm(j)*(1.0d0-rpl)+sdm(j+1)*rpl)*m_sun

      return

      end
