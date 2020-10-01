       real*8 function dssenu_capsun(mx,rho,sigsi,sigsd)
c----------------------------------------------------------------------
c         capture rate in the sun
c         based on jungman, kamionkowski, griest review
c       mx: WIMP mass
c       sigsi: spin independent cross section in units of cm^2
c       sigsd: spin dependent cross section in units of cm^2
c       vobs: average halo velocity
c       output:
c         capture rate in s^-1
c       lars bergstrom 1995-12-12
c----------------------------------------------------------------------
       implicit none
       include 'dsmpconst.h'
       include 'dshmcom.h' ! needed for vd_3d
       real*8 mx,rho,sigsi,sigsd,dssenu_litlf_s,dssenu_ss

       dssenu_capsun=sigsi*1.0d40*2.57d-13*pi/4./m_p**2 ! from sigsi to f_p
       dssenu_capsun=dssenu_capsun*dssenu_litlf_s(mx,vd_3d)*2.4d37
c   add spin dependent piece:
       dssenu_capsun=dssenu_capsun+
     &   1.3d25*sigsd*1.0d40*dssenu_ss(mx/m_p,vd_3d)/(mx*vd_3d/270.d0)
       dssenu_capsun=dssenu_capsun*(rho/0.3)
       return
       end








