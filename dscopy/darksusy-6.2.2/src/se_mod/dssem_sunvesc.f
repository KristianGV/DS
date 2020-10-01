

**********************************************************************
*** dssem_sunvesc gives the escape velocity in km/s as a function of
*** the radius r (in meters) from the sun's core.
*** author: joakim edsjo (edsjo@fysik.su.se)
*** input: radius in m
*** output escape velocity in km/s
*** date: 2003-11-26
**********************************************************************

      real*8 function dssem_sunvesc(r)
      implicit none

      real*8 r,dssem_sunpot,vescsurf,phisurf
      parameter(vescsurf=617.57d0)  ! km/s
      parameter(phisurf=-1.9069d11) ! m^2 s^-2

      dssem_sunvesc=vescsurf*sqrt(abs(dssem_sunpot(r)/phisurf))

      return
      end
