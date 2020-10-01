****************************************************************
*** Ultracompact primordial minihalo density profile,        ***  
*** without adiabatic contraction.                           ***
***                                                          ***
*** radialdist = distance from centre of minihalo in kpc     ***
***                                                          ***
*** Input: radial distance in kpc                            ***
***                                                          ***
*** Output: density in gev/cm**3                             ***
***                                                          ***
*** Author: Pat Scott (pat@fysik.su.se)                      ***      
*** Date: 2009-08-08                                         ***
****************************************************************

      real*8 function dshmucmhrho(radialdist)
      implicit none
      real*8 radialdist, r, Rh_local
      real*8 dsucmh_numerical_dens
c      include 'dshmcom.h'
      include 'dsucmh.h'
      include 'dsmpconst.h'

      if (radialdist .gt. Rh_ucmh) then
        dshmucmhrho = 0.d0
        return
      endif

      if (ucmhprofile .eq. 'plain') then
      
        if (rhcut .gt. Rh_ucmh) then
          r = rhcut * 1.d3 * cmperpc
        else
          r = radialdist * 1.d3 * cmperpc     !r in cm
        endif
        Rh_local = Rh_ucmh * 1.d3 * cmperpc       !Rh in cm
        dshmucmhrho = 3.d0 / (16.d0 * pi) * Mh_ucmh * Rh_local **(-0.75d0) * r**(-2.25d0)
        dshmucmhrho = dshmucmhrho * dmfrac * GeVperSolarMass
      
      else
     
        if (rhcut .gt. Rh_ucmh) then
          dshmucmhrho = maxdens
        else
          dshmucmhrho = dsucmh_numerical_dens(radialdist)
        endif    
 
      endif
     
      return
      end

