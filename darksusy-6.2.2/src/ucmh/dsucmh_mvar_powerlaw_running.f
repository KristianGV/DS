!Mass variance calculation for running power-law power spectra.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input:   n               spectral index of primordial spectrum of perturbations (dimensionless)
!         alpha           slope of running spectral index (dimensionless) 
!         k               wavenumber of perturbation (Mpc^-1; 0.5d3/lambda)
!
!Output: mass variance (dimensionless)

      double precision function dsucmh_mvar_powerlaw_running(n,alpha,k)

      implicit none

      double precision, intent(IN) :: n, alpha, k
      double precision :: alphaSq, deltaHSq, dsucmh_transFunc_rad, midpnt
      double precision :: dsucmh_deltahsq, alphaIntegrand, qromo, kratio, midinf
      double precision, parameter :: IntCut = 2.d0
      external midpnt, midinf, alphaIntegrand

      !Get the normalisation of the power spectrum
      deltaHSq = dsucmh_deltahsq(n,alpha,k,kratio)

      !Get alpha^2
      alphaSq = dsucmh_transFunc_rad(1.d0)
      alphaSq = 81.d0/16.d0/(alphaSq*alphaSq)
      alphaSq = alphaSq*(qromo(alphaIntegrand,n,alpha,kratio,0.d0,IntCut,midpnt)+
     &          qromo(alphaIntegrand,n,alpha,kratio,IntCut,1.d300,midinf))

      !Apply the normalisation to alpha^2
      dsucmh_mvar_powerlaw_running = alphaSq*deltaHSq

      end function dsucmh_mvar_powerlaw_running

