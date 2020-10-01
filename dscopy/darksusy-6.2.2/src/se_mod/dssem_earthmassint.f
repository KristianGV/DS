

**********************************************************************
*** dssem_earthmassint gives the mass of the earth in units of kg within
*** a sphere with of radius r meters.
*** in this routine, the actual integration is performed. for speed,
*** use dssem_earthmass instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@fysik.su.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dssem_earthmassint(r)
      implicit none

      real*8 r,dssem_edfunc,dsf_int
      external dssem_edfunc

c...integrate earth density

      dssem_earthmassint=dsf_int(dssem_edfunc,0.0d0,r,1.d-2)

      return
      end
