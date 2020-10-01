      real*8 function dsdbtdps_td(t,td)
************************************************************************
*** antideuteron "confinement time per annihilation volume" for a 
*** non-static point source within which WIMP pair annihilations take 
*** place.
***
*** inputs:
***   it assumes that the observer is located at x,y,z=0       
***     t - time of observation (10^15 s), having assumed t=0 for the
***         time at which the source is the closest to the observer
***     td - antideuteron kinetic energy per nucleon (gev)
***
*** output: in kpc^-3 10^15 s
***
*** the antideuteron flux from wimp pair annihilation is obtained by
*** multiplying the result of this function by the factor:
***   units * 1/(4 pi) * vd * rho2vol * (sigma v)/(2*mwimp**2) * dn/dtd
*** where:
***   vd = vd(td) = antideuteron velocity for given antideuteron kinetic
***      energy td in unit of c 
***   rho2vol = int dr^3_s \rho_s^2(r_s) = integral over the source
***      volume of the wimp density squared in kpc^3 GeV^2 cm^-6
***   sigma v = dssigmav0tot in cm^3 s^-1
***   mwimp = WIMP mass in GeV
***   dn/dtd = dn/dtd(td) = antideuteron spectrum in GeV^-1
***   units = factor taking into account units =
***    1.d15 s sr^-1 2.99792d10 cm s^-1 GeV^2 cm^-6 cm^3 s^-1 GeV^-3
***      = 2.99792d25 cm^-2 s^-1 GeV^-1 sr^-1
*** 
*** for dark matter particle decays, the antideuteron flux is obtained 
*** by multiplying the result of this function by the factor:
***   units * 1/(4 pi) * vd * rhovol * dec-rate/mwimp * dn/dtd
*** having defined (as opposed to the previous case):
***   rho2vol = int dr^3_s \rho_s(r_s) = integral over the source
***      volume of the source density in kpc^3 GeV cm^-3
***   dec-rate = dsdecratewimp = dark matter particle decay rate in s^-1
***   mwimp = decaying particle mass in GeV
***   units = factor taking into account units =
***    1.d15 s sr^-1 2.99792d10 cm s^-1 GeV cm^-3 s^-1 GeV^-2
***      = 2.99792d25 cm^-2 s^-1 GeV^-1 sr^-1
***
************************************************************************
      implicit none
      include 'dsmpconst.h' 
      real*8 t,td
      real*8 dspbtdpsc_td
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
ccc set the antiproton mass and atomic number
ccc
      mnuc=m_d
      anuc=2
      znuc=1
ccc
      dsdbtdps_td=dspbtdpsc_td(t,td)
      return
      end
