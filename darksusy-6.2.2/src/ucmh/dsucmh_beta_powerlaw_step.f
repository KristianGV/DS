!Works out UCMH/PBH fraction at equality for running power-law power spectra with a step.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
! 
!Input:  n          spectral index (dimensionless)
!        alpha      running of spectral index (dimensionless)
!        p          size of step in amplitude of power spectrum (dimensionless)
!        k_s        wavenumber at which the step occurs (Mpc^-1) 
!        k          wavenumber corresponding to scale of perturbations (Mpc^-1)
!        delta_min  lower limit of density contrasts that could form UCMHs (dimensionless)
!        delta_max  upper limit of density contrasts that could form UCMHs (dimensionless)
!
!Output: beta       fraction of perturbations leading to UCMHs/PBHs (dimensionless)

      double precision function dsucmh_beta_powerlaw_step(n,alpha,p,k_s,k,delta_min,delta_max)

      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'

      double precision, intent(IN) :: n, alpha, p, k_s, k, delta_min, delta_max
      double precision :: lnp2, dsf_int, betaintegrand, dsucmh_mvar_powerlaw_step

      external betaintegrand

      !Check the size of p to prevent infinities
      lnp2 = merge(2.d0*log(p), -1.d2, p .gt. epsilon(p))
      !Work out the mass variance
      twoucmhsigma_sq = 2.d0 * dsucmh_mvar_powerlaw_step(n,alpha,lnp2,k_s,k)
      !Work out beta by integrating over the spectrum of perturbations, then normalising.
      dsucmh_beta_powerlaw_step = dsf_int(betaintegrand,delta_min,delta_max,1.d-8)
      dsucmh_beta_powerlaw_step = dsucmh_beta_powerlaw_step / sqrt(pi*twoucmhsigma_sq)

      end function dsucmh_beta_powerlaw_step

