!Mass variance calculation for running power-law power spectra with a step.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input:   n               spectral index of primordial spectrum of perturbations (dimensionless)
!         alpha           slope of running spectral index (dimensionless) 
!         logp2           ln(square of step size in primordial spectrum of density perturbations) (dimensionless)
!         k_s             wavenumber at which spectrum has a step discontinuity (Mpc^-1)
!         k               wavenumber of perturbation (Mpc^-1; 0.5d3/lambda)
!
!Output: mass variance (dimensionless)

      double precision function dsucmh_mvar_powerlaw_step(n,alpha,logp2,k_s,k)

      implicit none

      double precision, intent(IN) :: n, alpha, logp2, k_s, k 
      double precision :: xcut, alphaSq, deltaHSq, dsucmh_transFunc_rad, midpnt
      double precision :: dsucmh_deltahsq, alphaIntegrand, qromo, kratio, midinf
      double precision, parameter :: IntCut = 2.d0
      external midpnt, midinf, alphaIntegrand

      !Dimensionless step location
      xcut = k_s/k

      !Get the normalisation of the power spectrum
      deltaHSq = dsucmh_deltahsq(n,alpha,k,kratio)

      !Get the first part of alpha^2
      alphaSq = dsucmh_transFunc_rad(1.d0)
      alphaSq = 81.d0/16.d0/(alphaSq*alphaSq)

      !Get the second part of alpha^2
      if (xcut .lt. IntCut) then
        alphaSq = alphaSq * ( qromo(alphaIntegrand,n,alpha,kratio,0.d0,xcut,midpnt) +
     &             exp(logp2) * (qromo(alphaIntegrand,n,alpha,kratio,xcut,IntCut,midpnt) + 
     &             qromo(alphaIntegrand,n,alpha,kratio,IntCut,1.d300,midinf)) )

      else
        alphaSq = alphaSq * ( qromo(alphaIntegrand,n,alpha,kratio,0.d0,IntCut,midpnt) +
     &             qromo(alphaIntegrand,n,alpha,kratio,IntCut,xcut,midinf) + 
     &             exp(logp2)*qromo(alphaIntegrand,n,alpha,kratio,xcut,1.d300,midinf) )
      endif

      !Apply the normalisation to alpha^2
      dsucmh_mvar_powerlaw_step = alphaSq*deltaHSq

      end function dsucmh_mvar_powerlaw_step

