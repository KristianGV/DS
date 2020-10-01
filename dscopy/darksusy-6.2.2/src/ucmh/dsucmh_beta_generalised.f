!Works out UCMH/PBH fraction at equality for generalised power spectra.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
! 
!Input:  log10P_R   log10(power in curvature perturbations) (dimensionless)
!        delta_min  lower limit of density contrasts that could form UCMHs (dimensionless)
!        delta_max  upper limit of density contrasts that could form UCMHs (dimensionless)
!
!Output: beta       fraction of perturbations leading to UCMHs/PBHs (dimensionless)

      double precision function dsucmh_beta_generalised(log10P_R,delta_min,delta_max)

      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'

      double precision, intent(IN) :: log10P_R, delta_min, delta_max
      double precision :: dsucmh_mvar_generalised, dsf_int, betaintegrand

      external betaintegrand

      !Work out the mass variance
      twoucmhsigma_sq = 2.d0 * dsucmh_mvar_generalised(log10P_R)
      !Work out beta by integrating over the spectrum of perturbations, then normalising.
      dsucmh_beta_generalised = dsf_int(betaintegrand,delta_min,delta_max,1.d-8)
      dsucmh_beta_generalised = dsucmh_beta_generalised / sqrt(pi*twoucmhsigma_sq)

      end function dsucmh_beta_generalised

