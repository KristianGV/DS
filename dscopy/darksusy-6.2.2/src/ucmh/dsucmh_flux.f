!Calculates the gamma-ray flux expected from a single UCMH
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: mostly late 2009, in this form Oct 2014
!      
!Input:   mucmh            UCMH mass (M_Sun)
!         rucmh            UCMH radius (kpc)
!         rcore            Core radius of the UCMH (kpc)
!         dist             distance at which UCMH is situated, for flux calculation (kpc)
!         ethreshold       Instrumental energy threshold (GeV)
!
!Output:  flux             flux from WIMP annihilation in a single UCMH (cm^-2 s^-1)
!         Eflux            'Energy flux' from WIMP annihilation in a single UCMH (GeV cm^-2 s^-1)

      subroutine dsucmh_flux(mucmh, rucmh, rcore, dist, ethreshold, flux, Eflux)

      implicit none

      include 'dsucmh.h'
      include 'dslocal.h'
      include 'dsmpconst.h'

      double precision, intent(IN) :: mucmh, rucmh, rcore, dist, ethreshold
      double precision, intent(OUT) :: flux, Eflux
      double precision, parameter :: Eflux_maxE = 10.0 !Integration limit for diffuse 'energy flux' in GeV
      double precision :: totalJ, dshmucmhrho, dsucmhEfluxyield, dsucmh_jpntsrc !, dshrgacsusy
c      double precision :: delta, cospsi0, dshmj
      integer :: istat=0
      double precision :: dscrsource
     
      !Correct details of UCMH profile.
      !The radius of the core is set to the greater of rcore and the core radius calculated at the time the density profile was initialised.
      !Likewise for the maximum density -- if the implied maximum density at rcore is less than the existing maximum, adopt the density at rcore.
      call dshmucmhset(mucmh, rucmh, max(rhcut,rcore), dist)
      maxdens = min(maxdens, dshmucmhrho(rcore))

      ! Calculate the J factor for the UCMH, assuming a constant density core.
      totalJ = dsucmh_jpntsrc(max(rhcut,rcore),Rh_ucmh)
         
      ! Calculate the fluxes
      flux = dscrsource(ethreshold,0,22,2,0.d0,istat)*totalJ
c      flux = 0.5d0*dsanyield(ethreshold,52,istat)*sv/mwimp/mwimp*totalJ
      if (ethreshold .lt. Eflux_maxE) then
c        Eflux = 0.5d0*dsucmhEfluxyield(ethreshold,Eflux_maxE)*sv/mwimp/mwimp*totalJ
        Eflux = dsucmhEfluxyield(ethreshold,Eflux_maxE)*totalJ           
      else
        Eflux = 0.d0
      endif

      !Change this switch to roughly check if your UCMH is a pt source to Fermi at dist or not; the final flux should come out lower than 
      !the Fermi pt source sensitivity for the assumption that it appears as a point source to be approximately correct.
c      if (.false.) then
c        cospsi0=1.d0 !dcos(viewing angle / 180.d0 * pi)
c        delta = 0.01d0 / 180.d0 / 180.d0 * pi * pi
c        totalJ = dshmj(cospsi0,delta)*delta
c        flux = dshrgacsusy(ethreshold,istat)*totalJ
c      endif
c      if (istat .ne. 0) write(*,*)
c     &  'Error from dshrgacsusy in dsucmh_flux: ',istat     

      end subroutine dsucmh_flux


      double precision function dsucmhEfluxyield(ethreshold,Etop)

      implicit none

      double precision, intent(IN) :: ethreshold, Etop
      double precision :: Efluxintegrand, dsf_int
      external Efluxintegrand

      dsucmhEfluxyield = dsf_int(Efluxintegrand,ethreshold,Etop,1.d-8)

      end function dsucmhEfluxyield


      double precision function Efluxintegrand(E)

      implicit none

      double precision, intent(IN) :: E
c      double precision :: dsanyield
      double precision :: dscrsource
      integer :: istat

c      Efluxintegrand = E*dsanyield(E,152,istat)
       Efluxintegrand = E*dscrsource(E,1,22,2,0.d0,istat)

      end function Efluxintegrand

