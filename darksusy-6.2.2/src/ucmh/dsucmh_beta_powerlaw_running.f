!Works out UCMH/PBH fraction at equality for running power-law power spectra.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
! 
!Input:  n          spectral index (dimensionless)
!        alpha      running of spectral index (dimensionless)
!        k          wavenumber corresponding to scale of perturbations (Mpc^-1)
!        delta_min  lower limit of density contrasts that could form UCMHs (dimensionless)
!        delta_max  upper limit of density contrasts that could form UCMHs (dimensionless)
!        mH         horizon mass (M_solar)
!        rough      use the rough Green-Liddle formula
!
!Output: beta       fraction of perturbations leading to UCMHs/PBHs (dimensionless)

      double precision function dsucmh_beta_powerlaw_running(n,alpha,k,delta_min,delta_max,mH,rough)

      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'

      logical, intent(IN) :: rough
      double precision, intent(IN) :: n, alpha, k, delta_min, delta_max, mH 
      double precision :: dsf_int, betaintegrand, dsucmh_mvar_powerlaw_running

      external betaintegrand

      if (rough) then
        !Work out the mass variance using the simple Green & Liddle normalisation
        twoucmhsigma_sq = 2.d0 * 9.025d-9 * (mH * 1.98892d-23)**(0.5d0-0.5d0*n)
      else
        !Work out the mass variance properly
        twoucmhsigma_sq = 2.d0 * dsucmh_mvar_powerlaw_running(n,alpha,k)
      endif
      !Work out beta by integrating over the spectrum of perturbations, then normalising.
      dsucmh_beta_powerlaw_running = dsf_int(betaintegrand,delta_min,delta_max,1.d-8)
      dsucmh_beta_powerlaw_running = dsucmh_beta_powerlaw_running / sqrt(pi*twoucmhsigma_sq)

      end function dsucmh_beta_powerlaw_running

