***********************************************************************
*** dssenu_capearthnum calculates the capture rate at present.
*** Intead of using the assumptions of Gould (i.e. capture as in
*** free space), a tabulted velocity distribution based on detailed
*** numerical simulations of Johan Lundberg is used.
*** A numerical intregration has to be performed instead of the
*** convenient expressions in jkg.
*** Input: mx = WIMP mass in GeV
***        rho = local halo density in GeV/cm^3      
***        sigsi = spin-independent scattering cross section in cm^2
*** author: joakim edsjo (edsjo@fysik.su.se)
*** date: July 10, 2003
***********************************************************************
      real*8 function dssenu_capearthnum(mx,rho,sigsi)
      implicit none
      include 'dssecom.h'

      real*8 mx,rho,sigsi
      real*8 dssenu_capearthnumi
      real*8 dssenu_foveruearth

      external dssenu_foveruearth

c...integrate over the earth according to gould apj 321 (1987) 571.
c...set up things for radial and velocity integration

c...Transfer to common blocks for integration

      semx=mx                   ! so that internal routines knows the mass
      serho=rho

c...perform integration
      dssenu_capearthnum=dssenu_capearthnumi(mx,sigsi,
     &     dssenu_foveruearth,1)

      return
      end
