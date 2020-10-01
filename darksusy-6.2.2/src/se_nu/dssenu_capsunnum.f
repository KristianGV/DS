***********************************************************************
*** dssenu_capsunnum calculates the capture rate at present.
*** Instead of using the approximations in jkg, i.e. a gaussian
*** velocity distribution and approximating all elements as being 
*** at their typical radius, we here integrate numerically over
*** the actual velocity distribution and over the Sun's radius.
*** The velocity distribution used is the one set up by the
*** option veldf in dshmcom.h (see src/hm/dshmudf.f for details)
*** Input: mx = WIMP mass in GeV
***        rho = local halo density in GeV/cm^3      
***        sigsi = spin-independent scattering cross section (cm^2)
***        sgisd = spin-dependent scattering cross section on protons (cm^2)
*** author: joakim edsjo (edsjo@fysik.su.se)
*** date: 2003-11-26
***********************************************************************
      real*8 function dssenu_capsunnum(mx,rho,sigsi,sigsd)
      implicit none
      include 'dssecom.h'

      real*8 mx,rho,sigsi,sigsd
      real*8 dssenu_capsunnumi,dssenu_foveru
      real*8 dssenu_veoutjupiter

      external dssenu_foveru

c...integrate over the sun according to gould apj 321 (1987) 571.
c...set up things for radial and velocity integration

      semx=mx                   ! so that internal nt routines know the mass
      serho=rho

c...Include Jupiter effects
c...The veout below adds the the possibility to only consider scatterings
c...that will reach out to a distance from the Sun with the veout escape
c...velocity. Putting this to the escape velocity of Jupiter means that 
c...WIMPs that after the first scatter would reach further out than Jupiter
c...are considered as not captured. The rationale for this is work by
c...Peter and Tremaine showing that WIMPs that reach out to Jupiter
c...most likely will be disturbed by Jupiter before they have a chance
c...to scatter again and sink to the Sun's core. 
c...Putting this to zero gives back the usual Gould capture formulae. 
c...sejup=2 means that we set the maximal distance after first scatter based
c...on the scattering cross section.
      if (sejup.eq.0) then ! no Jupiter effect
         veout=0.d0
      elseif (sejup.eq.1) then ! approximate Jupiter effect
         veout=sqrt(2.0d0)*13.07d0 ! escape velocity at Jupiter
      elseif (sejup.eq.2) then ! more accurate Jupiter approach
         veout=dssenu_veoutjupiter(mx,sigsi,sigsd)
      else
         write(*,*) 
     &  'WARNING in dssenu_capsunnum: uknown option sejup = ',sejup
         veout=0.0d0
      endif

c...perform integration
      dssenu_capsunnum=dssenu_capsunnumi(mx,sigsi,sigsd,dssenu_foveru)

      return
      end
