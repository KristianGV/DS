!Power spectrum normalisation delta_H^2
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
! 
!Input:  n          spectral index (dimensionless)
!        alpha      running of spectral index (dimensionless)
!        k          wavenumber corresponding to scale of perturbations (Mpc^-1)
!
!Output: beta       power spectrum normalisation (dimensionless)
!        kratio     ratio of k to reference scale (dimensionless)

      double precision function dsucmh_deltahsq(n,alpha,k,kratio)

      implicit none

      double precision, intent(IN) :: n, alpha, k
      double precision, intent(OUT) :: kratio
      !double precision, parameter :: a = -0.95d0, b = -0.169  !Old normalisation, without gravitational waves
      !double precision, parameter :: a = 1.d0, b = 1.97d0     !Old normalisation, with gravitational waves from inflation
      double precision, parameter :: r = 0.d0                  !New normalisation; scalar to tensor ratio
      double precision :: ntilde

      ntilde = n-1.d0       !Modified spectral index
      !kratio = (k*c_0/H_0) !Old COBE normalisation (arXiv:1110.2484v1)
      !deltaHSq = 5.13d-9 * exp(2.d0*(a*ntilde + b*ntilde*ntilde)) * kratio**(ntilde + alpha*log(kratio))
      kratio = (k/0.05d0)   !WMAP3 normalisation (arXiv:1110.2484v2, from Liddle et al. 2006 PRD 74:083512 * (10/9)^2)
      dsucmh_deltahsq = 4.5844d-10 / (0.53*r+1.d0) * exp(2.d0*ntilde*(1.24d0-1.04d0*r)) * kratio**(ntilde + alpha*log(kratio))

      end function dsucmh_deltahsq
