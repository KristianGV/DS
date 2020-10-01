!Limit on cosmological fraction of UCMHs from extragalactic point sources.
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

      double precision function dsucmh_fmax_gamma_exgalptsrc(flux, dist, skyfrac, gamma_pt_sensitivity, CL_ptsrc, CL_fmax, quiet)

      implicit none

      include 'dsmpconst.h'
      include 'dsucmh.h'
      include 'dslocal.h'
      
      logical, intent(IN) :: quiet
      double precision, intent(IN) :: flux, dist, skyfrac, gamma_pt_sensitivity, CL_ptsrc, CL_fmax
      double precision :: dmin, massindmin, safelog1m

      !Work out the minimum distance at which UCMHs are not observable at CL CL_ptsrc 
      dmin = dist * sqrt(flux/gamma_pt_sensitivity)

      !If UCMHs are visible beyond the confines of the Milky Way, comput limits from point extragalactic source searches.
      if (r_virial/pc*1.d-3 .lt. dmin) then
        !Work out how much cosmological + Milky Way dark matter is contained in dmin
        massindmin = (4.d0/3.d0*pi*((dmin*pc)**3*1.d9 - r_virial*r_virial*r_virial)*rhobg + M_200)/Msun
        !Work out the corresponding maximum fraction of DM contained in UCMHs
        !dsucmh_fmax_gamma_exgalptsrc = -1.d0 * dmfrac * Mh_ucmh / massindmin * safelog1m(CL_fmax/CL_ptsrc)  !old version (arXiv:1110.2484)
        dsucmh_fmax_gamma_exgalptsrc = -1.d0 / CL_ptsrc * dmfrac * Mh_ucmh / massindmin * safelog1m(CL_fmax) !new version (arXiv:1211.7361)
      else
        dsucmh_fmax_gamma_exgalptsrc = 10.d0
      endif

      !Rescale fmax according to the fraction of the sky being probed
      if (skyfrac .gt. 0.d0) dsucmh_fmax_gamma_exgalptsrc = dsucmh_fmax_gamma_exgalptsrc / skyfrac

      !Make sure that fmax is less than or equal to 1, as UCMHs must be less than 100% of the Milky Way. 
      !Otherwise, the limit is invalid for this point, so raise a warning.
      if (dsucmh_fmax_gamma_exgalptsrc .gt. 1.d0 .and. .not. quiet) then 
        write(*,*) 'Warning: fmax in dsucmh_fmax_gamma_exgalptsrc = ', dsucmh_fmax_gamma_exgalptsrc, '> 1, so limits are invalid.'
      endif

      end function dsucmh_fmax_gamma_exgalptsrc
