***********************************************************************
*** This routine uses a derived column density from the given solar model
*** The data in sdcens() is calculated by dssem_sunread.f. 
***********************************************************************

***********************************************************************
*** dssem_suncdens gives the column density in the Sun from the
*** centre out the tha given radius r (in meters).
*** The radius should be given in m and the column density is returned in
*** g/cm^2
*** if type = 'N', the total column density (up to that r) is calculated
***         = 'p', the column density on protons is calculated
***         = 'n', the column density on neutrons is calculated
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2005-11-25
***********************************************************************

      real*8 function dssem_suncdens(r,type)
      implicit none

      include 'dssem_sun.h'
      real*8 r,rpl
      character*1 type
      integer j,ti

      j=0                       ! initialize
      
c...Check if data file is loaded
      call dssem_sunread

      if (type.eq.'N') then 
        ti=0
      elseif (type.eq.'p') then
        ti=1
      elseif (type.eq.'n') then
        ti=2
      else
        write(*,*) 'ERROR in dssem_sundens: wrong type: ',type
        stop
      endif

      if (r.ge.r_sun) then
        dssem_suncdens=cd_sun(ti) ! total column density
        return
       endif

      if (r.le.0.0d0) then
        dssem_suncdens=0.0d0 ! no column density
        return
      endif

c...Interpolate in table

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      dssem_suncdens=sdcdens(sdn,ti)

      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dssem_suncdens=sdcdens(j,ti)*(1.0d0-rpl)+sdcdens(j+1,ti)*rpl

      return

      end
