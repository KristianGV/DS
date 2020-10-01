!Limit on cosmological fraction of UCMHs from high-lattitude diffuse gamma-ray background
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
! 
!Input:   Eflux           energy-flux from WIMP annihilation in a single UCMH at distance dist (GeV cm^-2 s^-1)
!         dist            reference distance at which UCMH is situated (kpc)
!         skyfrac         Fraction of the sky accessible to the probe in question (dimensionless)
!         CL_fmax         The confidence level at which fmax should be computed (dimensionless)
!         quiet           Suppress warnings
!
!Output:  fmax            maximum cosmological fraction of UCMHs today (dimensionless)

      double precision function dsucmh_fmax_gamma_diffuse(Eflux, dist, skyfrac, CL_fmax, quiet)

      implicit none

      include 'dsmpconst.h'
      include 'dsucmh.h'
      include 'dslocal.h'

      logical, intent(IN) :: quiet
      double precision, intent(IN) :: Eflux, dist, skyfrac, CL_fmax
      double precision, parameter :: effAreaSigmaPercentage = 0.1, fermiObsEflux_raw = 1.d-5 !(dimensionless, GeV cm^-2 s^-1 sr^-1)
      double precision :: sigmaReq, fermiObsEflux, erfinv, dmax, dsf_int, Jfactor, totalMaxEflux
      external Jfactor

      !Choose upper observational limit on observed flux according to desired CL on f
      sigmaReq = erfinv(CL_fmax, 1.d0-CL_fmax)*sqrt(2.d0)
      fermiObsEflux = fermiObsEflux_raw * (1.d0 + sigmaReq*effAreaSigmaPercentage)

      !Work out how much diffuse flux would be expected if all DM in the Galaxy were in UCMHs
      dmax = sqrt(r_virial*r_virial + r_sun_gal*r_sun_gal)
      totalMaxEflux = 1.d0 / (dmfrac * Mh_ucmh) * Eflux * dist*dist * dsf_int(Jfactor,0.d0,dmax,1.d-8) *pc*pc*1d6 / Msun
      !Rescale to find fmax from the diffuse flux observed by Fermi
      dsucmh_fmax_gamma_diffuse = fermiObsEflux / totalMaxEflux
     
      !Rescale fmax according to the fraction of the sky being probed
      if (skyfrac .gt. 0.d0) dsucmh_fmax_gamma_diffuse = dsucmh_fmax_gamma_diffuse / skyfrac

      !Make sure that fmax is less than or equal to 1, as UCMHs must be less than 100% of the Milky Way. 
      !Otherwise, the limit is invalid for this point, so throw an error.
      if (dsucmh_fmax_gamma_diffuse .gt. 1.d0 .and. .not. quiet) then 
        write(*,*) 'Warning: fmax in dsucmh_fmax_gamma_diffuse = ', dsucmh_fmax_gamma_diffuse, '> 1, so limits are invalid.'
      endif

      end function dsucmh_fmax_gamma_diffuse


      double precision function Jfactor(d)
      !input    d  distance from Sun (m)
      !output      (kg m^-3)

      implicit none

      include 'dslocal.h'

      double precision, intent(IN) :: d
      double precision :: nfwrho

      Jfactor = nfwrho(sqrt(d*d + r_sun_gal*r_sun_gal))

      end function Jfactor


      double precision function nfwrho(r)
      !input     d  distance from centre of Galaxy (m)
      !output       NFW DM density  (kg m^-3)    

      implicit none

      include 'dslocal.h'

      double precision, intent(IN) :: r
      double precision :: r_s

      r_s = r_s3**(1.d0/3.d0)
      nfwrho = rhobg*del/(r/r_s*(1.d0+r/r_s)**2)
      
      end function nfwrho

