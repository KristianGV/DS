***********************************************************************
*** note. this routine assumes a maxwell-boltzmann velocity
*** distribution and uses approximations in the jkg review,
*** Jungman, Kamionkowski and Griest, Phys. Rep. 267 (1996) 195.
*** In particular, it is assumed that the Sun's velocity is
*** sqrt(2/3)*vd_3d.
*** for more accurate results, use dssenu_capearthfull instead.
*** for an arbitrary velocity distribution, use dssenu_capearthnumi instead.
***********************************************************************

      real*8 function dssenu_capearth(mx,rho,sigsi)
c----------------------------------------------------------------------
c         capture rate in the earth
c         based on jungman, kamionkowski, griest review
c       mx: WIMP mass
c       sigsi: spin independent cross section in units of cm^2
c       vobs: average halo velocity
c       lars bergstrom 1995-12-14
c----------------------------------------------------------------------
      implicit none
      include 'dsmpconst.h'
      include 'dshmcom.h' ! needed for vd_3d
      real*8 mx,rho,sigsi,dssenu_litlf_e

c...Eq. (9-27) in jkg with the assumption m_chi >> m_proton
      dssenu_capearth=sigsi*1.0d40*2.57d-13*pi/4./m_p**2 ! from sigsi to f_p
      dssenu_capearth=dssenu_capearth*dssenu_litlf_e(mx,vd_3d)
     &  *2.4d28*(rho/0.3)
      return
      end
