***********************************************************************
*** dssem_sunmfrac gives the mass fraction of an element with atomic 
*** number z and with isotope number iso (see dssem_sunread.f
*** for definition of iso) as a function of the solar radius r.
*** the radius should be given in m and returned is the mass fraction.
***
*** Element mass fractions are from the solar models read in in
*** dssem_sunread. Different choices are available in dssem_sunset.
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2003-11-26
***********************************************************************

      real*8 function dssem_sunmfrac(r,z,iso)
      implicit none

      include 'dssem_sun.h'
      real*8 r,rpl
      integer j,z,iso

      j=0                       ! initialize
      
c...Check if data file is loaded
      call dssem_sunread

      if (r.gt.r_sun) then
        dssem_sunmfrac=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dssem_sunmfrac=sdmfr(z,iso,1)
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      j=sdn-1

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dssem_sunmfrac=sdmfr(z,iso,j)*(1.0d0-rpl)
     &  +sdmfr(z,iso,j+1)*rpl

      return

      end
