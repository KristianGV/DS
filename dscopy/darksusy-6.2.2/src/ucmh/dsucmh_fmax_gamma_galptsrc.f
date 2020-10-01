!Limit on cosmological fraction of UCMHs from galactic point sources.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
! 
!Input:   flux                 flux from WIMP annihilation in a single UCMH at distance dist (must have same units as gamma_pt_sensitivity)
!         dist                 reference distance at which UCMH is situated (kpc)
!         skyfrac              Fraction of the sky accessible to the probe in question (dimensionless)
!         gamma_pt_sensitivity Point source sensitivity of gamma-ray instrument (must have same units as flux)
!         CL_ptsrc             The confidence level with which gamma_pt_sensitivity is stated (dimensionless)
!         CL_fmax              The confidence level at which fmax should be computed (dimensionless)
!         quiet                Suppress warnings
!
!Output:  fmax                 maximum cosmological fraction of UCMHs today (dimensionless)

      double precision function dsucmh_fmax_gamma_galptsrc(flux, dist, skyfrac, gamma_pt_sensitivity, CL_ptsrc, CL_fmax, quiet)

      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'
      
      logical, intent(IN) :: quiet
      double precision, intent(IN) :: flux, dist, skyfrac, gamma_pt_sensitivity, CL_ptsrc, CL_fmax
      double precision :: nMWmin, dmin, massindmin, massfrac, safelog1m, dsucmh_localmass

      !Work out the minimum number of UCMHs in the Milky Way on which to base limits, on the basis of Poisson statistics 
      !nMWmin = float(int(0.5d0-safelog1m(CL_fmax/CL_ptsrc)))          !old version (arXiv:1110.2484)
      nMWmin = float(int(0.5d0-safelog1m(CL_fmax)/CL_ptsrc))           !new version (arXiv:1211.7361)
      
      !Work out the minimum distance at which UCMHs are not observable at CL CL_ptsrc 
      dmin = dist * sqrt(flux/gamma_pt_sensitivity)

      !Work out how much Milky Way dark matter is contained in dmin, and the fraction of the Milky Way's dark mass in UCMHs
      massindmin = dsucmh_localmass(dmin)
      massfrac = massindmin/M_MW
      if (massfrac .gt. 1.d0 .or. massfrac .lt. 0.d0) then
        write(*,*) 'Problem calculating massfrac in dsucmh_fmax_gamma_galptsrc: ',massfrac
        call exit(0)
      endif

      !Work out the corresponding maximum fraction of DM contained in UCMHs
      !dsucmh_fmax_gamma_galptsrc = dmfrac * Mh_ucmh / M_MW * safelog1m(CL_fmax/CL_ptsrc) / safelog1m(massfrac) !old version (arXiv:1110.2484)
      dsucmh_fmax_gamma_galptsrc = dmfrac * Mh_ucmh / M_MW * safelog1m(CL_fmax) / safelog1m(CL_ptsrc*massfrac)  !new version (arXiv:1211.7361)
      !If the limit is < nMWmin UCMHs in whole MW, reset it to nMWmin UCMHs
      if (dsucmh_fmax_gamma_galptsrc .lt. nMWmin * dmfrac * Mh_ucmh / M_MW) then
        dsucmh_fmax_gamma_galptsrc = nMWmin * dmfrac * Mh_ucmh / M_MW
      endif

      !Rescale fmax according to the fraction of the sky being probed
      if (skyfrac .gt. 0.d0) dsucmh_fmax_gamma_galptsrc = dsucmh_fmax_gamma_galptsrc / skyfrac

      !Make sure that fmax is less than or equal to 1, as UCMHs must be less than 100% of the Milky Way. 
      !Otherwise, the limit is invalid for this point, so throw an error.
      if (dsucmh_fmax_gamma_galptsrc .gt. 1.d0 .and. .not. quiet) then 
        write(*,*) 'Warning: fmax in dsucmh_fmax_gamma_galptsrc = ', dsucmh_fmax_gamma_galptsrc, '> 1, so limits are invalid.'
      endif

      end function dsucmh_fmax_gamma_galptsrc
