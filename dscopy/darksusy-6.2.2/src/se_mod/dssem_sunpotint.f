**********************************************************************
*** dssem_sunpotint gives the gravitational potential inside and outside
*** of the sun as a function of the radius r (in meters).
*** in this routine, the actual integration is performed. for speed,
*** use dssem_sunpot instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@fysik.su.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dssem_sunpotint(r)
      implicit none
      include 'dssem_sun.h'
      include 'dsmpconst.h'

      real*8 r,dssem_spfunc,dsf_int
      external dssem_spfunc

c...integrate sun density

      if (r.lt.r_sun) then
        dssem_sunpotint=
     &    -1.d-3*GNewton * m_sun/r_sun  ! surface potential, -1.9069d11 m^2 s^-2
     &    -dsf_int(dssem_spfunc,max(r,100.0d0),r_sun,1.d-3)
      else
        dssem_sunpotint=-m_sun*1.d-3*GNewton/(max(r,100.0d0))
      endif

      return
      end
