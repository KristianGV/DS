      subroutine dssenu_rates(egev,theta,kind,rtype,ptype,rho,
     &  rateea,ratesu,istat)
c_______________________________________________________________________
c
c    This routine calculates
c       - the capture rate of WIMPs in the Sun/Earth
c       - the annihilation rate of WIMPs in the Sun/Earth
c       - the yield of neutrinos or muons from these annihilations
c       - the total flux of these particles at the detector at Earth
c
c    type : commonly used
c    desc : Rates of neutrinos and neutrino-induced leptons and hadronic
c    desc : showers from WIMP annihilations in the Sun/Earth
c
c    november, 1995
c    uses routines by p. gondolo and j. edsjo
c    modified by l. bergstrom and j. edsjo
c    capture rate routines are written by l. bergstrom
c    input:
c       egev  - energy in GeV
c       theta  - angular cut in degrees
c       kind    - 1 = integrated fluxes above energy egev and up to angle theta
c               - 2 = differential fluxes at egev and theta
c               - 3 = mixed, integrated up to theta, differential in energy
c       rtype   - 1 = flux of neutrinos (nu_mu and/or nu_mu-bar) km^-2 yr^-1
c               - 2 = contained mu- and/or mu+ events km^-3 yr^-1
c               - 3 = through-going mu- and/or mu+ events km^-2 yr^-1
c       ptype   - 1 = particles only (nu_mu or mu-)
c               - 2 = anti-particles only (nu_mu-bar or mu+)
c               - 3 = summed rates (both particles and anti-particles)
c       rho - local dark matter halo density [GeV/cm^3]
c    hidden input (set via dssenu_set, see that routine for details)
c       secalcmet - 1 = use jkg approximations
c                   2 = use jkg for sun, full gould for earth
c                   3 = use jkg for sun, full gould+dk for earth (deprecated)
c                   4 = use full numerical calculations for Sun, Earth,
c                       but analytical form factors (exponential)
c                   5 = use full numerical calculaitons for Sun, Earth,
c                       and numerical form factor integrations [default]      
c    output: rateea  - events from earth ann. per km^2(3) per yr
c            ratesu  - events from sun ann. per km^2(3) per yr
c            the km^2 is for ptype=1,3 and km^3 for ptype=2      
c    For kind=1, the units are as above
c            =2, the units have an additional GeV^-1 degree^-1
c            =3, the units have an additional GeV^-1      
c    slightly modified by j. edsjo.
c    modified by j. edsjo 97-05-15 to match new inv. rate convention
c    modified by j. edsjo 97-12-03 to match muflux3.21 routines.
c    modified by p. gondolo 98-03-04 to detach dssenu_annrate from susy
c    routines.
c    modified by j. edsjo 98-09-07, corrected istat handling.
c    modified by j. edsjo 98-09-23 to use damour-krauss distributions
c      and full earth formulas.
c    modified by j. edsjo 99-03-17 to include better damour-krauss
c      velocity distributions and numerical capture rate integrations
c      for these non-gaussian distributions
c    modified by p. scott 11-04-23 to allow rates to be calculated
c      for only particles or only antiparticles
c    modified by j. edsjo 2014-12-10 to have both integrated, differential
c      and mixed fluxes in the same routine
c    modified by j. edsjo 2015-06-11 to clean up and avoid excessive use
c      of common block variables and options      
c
c=======================================================================
      implicit none
      include 'dshmcom.h'
      include 'dssecom.h'
      include 'dsseyieldcom.h'

      real*8 egev,theta,rateea,ratesu,arateea,aratesu,
     &     yield,mwimp,sigv0,gps,gns,gpa,gna,
     &     sigsip,sigsin,sigsdp,sigsdn
      complex*16 gg(27,2)
      real*8 sigij(27,27)
      integer ierrgg
      real*8 dssenu_summedyields,dsmwimp,dssigmav0tot
      integer istat,itmp,kind,rtype,ptype,ierr
      real*8 rho
c      real*8,optional,intent(out) :: hej

c ----------------------------------------- zero common block data

      tausu=0.0d0
      csu=0.0d0
      tauea=0.0d0
      cea=0.0d0
      searateea=0.0d0
      searatesu=0.0d0

c... initialize local variables      
      aratesu=0.0d0
      arateea=0.0d0

      rateea=0.0d0
      ratesu=0.0d0

c --------------------------------------------- start calculations

      if (rho.eq.0.0d0) then
        rateea=0.0d0
        ratesu=0.0d0
        return
      endif

c ------------------ get quantities from particle physics module

      mwimp=dsmwimp()
      sigv0=dssigmav0tot()
c...In principle we want to call dsddsigma directly in all capture routines
c...but to save time and make tabulation possible we use the four usual
c...four-fermion couplings here, via a call to dsddgpgn. Note though
c...that this direct call to dsddgpgn will most likely be phased out
c...at some point.      
      call dsddgpgn(gg,ierrgg)
      gps=real(gg(1,1))
      gns=real(gg(1,2))
      gpa=real(gg(4,1))
      gna=real(gg(4,2))
c...For those approximate calculations that use cross sections,
c...calculate cross sections (these are only used for
c...certain values of secalcmet below)      
      call dsddsigma(0.d0,0.d0,1,1,sigij,ierr)
      sigsip=sigij(1,1)
      sigsdp=sigij(4,4)
      call dsddsigma(0.d0,0.d0,1,0,sigij,ierr)
      sigsin=sigij(1,1)
      sigsdn=sigij(4,4)

c **************************************************************

      if (mwimp.gt.egev) then

        if (secalcmet.eq.1.or.secalcmet.eq.2
     &    .or.secalcmet.eq.4) then
c...      jkg and/or gould w/ Gauss or full dist.
          call dssenu_annrate(mwimp,rho,sigsip,sigsdp,sigv0,secalcmet,
     &       arateea,aratesu)
       elseif (secalcmet.eq.5) then ! full dist + numerical FF
          if (ierrgg.eq.0) then
            call dssenu_annrateff(mwimp,rho,gps,gns,gpa,gna,sigv0,
     &            arateea,aratesu)
          else                   ! we don't have the g's, revert to method 4
            call dssenu_annrate(mwimp,rho,sigsip,sigsdp,sigv0,4,
     &           arateea,aratesu)
          endif
        else
          write(*,*) 'DS ERROR in dssenu_rates: invalid option,',
     &      ' secalcmet = ',secalcmet
          return
        endif
c...arateea and aratesu in units of 10^24 yr^-1

        yield=dssenu_summedyields(egev,theta,'su',kind,rtype,ptype,istat)
c...yield is in units of 10^-30 m^-2(3).
        ratesu=yield*aratesu
c...we now have units km^-2 yr^-1 if rtype=3.
c...one more m^-1 -> km^-1 still to go for rtype=2
        if (rtype.eq.2) then
          ratesu=ratesu*1.0d3
        endif
        itmp=istat
        yield=dssenu_summedyields(egev,theta,'ea',kind,rtype,ptype,istat)
c...yield is in units of 10^-30 m^-2(3).
        rateea=yield*arateea
c...we now have units km^-2 yr^-1 if rtype=3.
c...one more m^-1 -> km^-1 still to go for rtype=2
        if (rtype.eq.2) then
          rateea=rateea*1.0d3
        endif
        istat=or(istat,itmp)
      else
        rateea=0.d0
        ratesu=0.d0
      endif

      searateea=arateea
      searatesu=aratesu

      return

      end




