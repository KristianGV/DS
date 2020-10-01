**********************************************************************
*** dssem_earthpotint gives the gravitational potential inside and outside
*** of the earth as a function of the radius r (in meters).
*** in this routine, the actual integration is performed. for speed,
*** use dssem_earthpot instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@fysik.su.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dssem_earthpotint(r)
      implicit none
      include 'dsmpconst.h'
      real*8 r,dssem_epfunc,dsf_int,dssem_earthmass
      external dssem_epfunc

c...integrate earth density

      if (r.lt.6378.140d3) then
        dssem_earthpotint=
     &    dsf_int(dssem_epfunc,100.0d0,max(r,110.0d0),1.d-2)-
     &    1.123782d8
      else
        dssem_earthpotint=-dssem_earthmass(6378.140d3)*
     &     1.d-3*GNewton/(max(r,100.0d0))
      endif

      return
      end
