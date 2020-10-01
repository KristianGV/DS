***********************************************************************
*** dssem_suncdfunc returns the density of protons, neutrons or the total
*** density depending on the common block variable cdt. If
***   cdt='N': the total density is returned
***   cdt='p': the density in protons is returned
***   cdt='n': the density in neutrons is returned
*** the radius should be given in m and the density is returned in
*** g/cm^3.
*** This routine is used by dssem_suncdensint to calculate the column
*** density in the Sun.
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2005-11-24
***********************************************************************

      real*8 function dssem_suncdfunc(r)
      implicit none
      include 'dssem_sun.h'

      real*8 dssem_sunmfrac,pfr,dssem_sundens,r

c...Check if data file is loaded
      call dssem_sunread

      if (cdt.eq.'N') then
        dssem_suncdfunc=dssem_sundens(r)
      elseif (cdt.eq.'p') then
        pfr=dssem_sunmfrac(r,1,1)
        dssem_suncdfunc=dssem_sundens(r)*(pfr + 0.5d0*(1.0d0-pfr))
      elseif (cdt.eq.'n') then
        pfr=dssem_sunmfrac(r,1,1)
        dssem_suncdfunc=dssem_sundens(r)*(0.5d0*(1.0d0-pfr))
      else
        write(*,*) 'ERROR in dssem_suncdfunc: wrong type: ',cdt
        stop
      endif

      return

      end
