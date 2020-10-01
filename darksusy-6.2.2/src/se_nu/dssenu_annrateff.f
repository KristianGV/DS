      subroutine dssenu_annrateff(mx,rho,gps,gns,gpa,gna,sigma_v,arateea,
     &  aratesu)
c_______________________________________________________________________
c
c     wimp annihilation rate in the sun and in the earth
c     in units of 10^24 annihilations per year
c
c     also gives the capture rate and the annih/capt equilibration time
c
c    november, 1995
c    uses routines by p. gondolo and j. edsjo
c    modified by l. bergstrom and j. edsjo and p. gondolo
c    capture rate routines are written by l. bergstrom
c    input:  mx      - wimp mass in GeV
c            rho - local halo density in GeV/cm^3      
c            sigsip  - spin-indep wimp-proton cross section in cm^2
c            sigsdp  - spin-dep wimp-proton cross section in cm^2
c            sigma_v  - wimp self-annihilation cross section in cm^3/s
c            rescale - rescale factor for local density
c    output: arateea  - 10^24 annihilations per year, earth
c            aratesu  - 10^24 annihilations per year, sun
c    slightly modified by j. edsjo.
c    modified by j. edsjo 97-05-15 to match new inv. rate convention
c    modified by j. edsjo 97-12-03 to match muflux3.21 routines.
c    modified by p. gondolo 98-03-04 to detach it from susy routines.
c    modified by J. Edsjo, 2010-05-28 to allow for numerical
c      form factor integration
c
c=======================================================================
      implicit none
      include 'dssecom.h'
      real*8 mx,rho,sigsip,sigsdp,gps,gns,gpa,gna,
     &     sigma_v,arateea,aratesu,
     &     ca,tt_sun,
     &     tt_earth,fluxs,fluxe,
     &     cap_sun,cap_earth
      real*8 dssenu_capearthnum, dssenu_capearthtab,
     &     dssenu_capsunnumff,dssenu_capsuntabff
      integer ierr
      complex*16 gg(27,2)
      real*8 sigij(27,27)
      real*8 spinx,dsdmspin

c --------------------------------------------- start calculations

c...Calculate cross sections (needed for the Earth)
      gg(1,1)=dcmplx(gps,0.d0)
      gg(1,2)=dcmplx(gns,0.d0)
      gg(4,1)=dcmplx(gpa,0.d0)
      gg(4,2)=dcmplx(gna,0.d0)
      spinx=dsdmspin()
      call dsddg2sigma(mx,spinx,0.d0,0.d0,1,1,gg,sigij,ierr)
      sigsip=sigij(1,1)
      sigsdp=sigij(4,4)

c...If setab=0, we do numerical integrations. If sejup=2 we also have to
c...do it numerically
      if (setab.eq.0.or.sejup.eq.2) then
        cap_earth=dssenu_capearthnum(mx,rho,sigsip) ! JE FIX: FF Earth
        cap_sun=dssenu_capsunnumff(mx,rho,gps,gns,gpa,gna) 
      else
        cap_earth=dssenu_capearthtab(mx,rho,sigsip)
        cap_sun=dssenu_capsuntabff(mx,rho,gps,gns,gpa,gna) 
      endif
        
      csu=cap_sun
      cea=cap_earth
c **************************************************************
      ca=sigma_v/6.6d28*(mx/20.d0)**(3./2.)
      tausu=1.0d0/dsqrt(cap_sun*ca)
      tt_sun=1.5d17*dsqrt(cap_sun*ca)
      ca=sigma_v/2.3d25*(mx/20.d0)**(3./2.)
      tauea=1.0d0/dsqrt(cap_earth*ca)
      tt_earth=1.5d17*dsqrt(cap_earth*ca)
      fluxs=cap_sun*0.5d0*tanh(tt_sun)**2
      fluxe=cap_earth*0.5d0*tanh(tt_earth)**2
      arateea = fluxe*1.d-24*3.15d7  ! 10^24 ann. per year
      aratesu = fluxs*1.d-24*3.15d7  ! 10^24 ann. per year

      return

      end
