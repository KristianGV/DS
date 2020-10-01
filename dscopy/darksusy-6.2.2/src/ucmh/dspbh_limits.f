!Returns limit on beta in PBHs from arXiv:0912.5297
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: Oct 2014
!
!Input:	k 	    wavenumber of perturbation (Mpc^-1)
!Output:    betamax   maximum allowed value of beta

      double precision function dspbh_limits(k)
      
      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'
 
      double precision :: k, logM
      double precision, parameter :: const = 1.d0/(3.**3*2.) 
      integer :: lowerindx = 0

      logM = log10(sqrt(const) * (keq/k)**2 * (106.75/3.36)**0.5 * (106.75/3.91)**(-2./3.) * Mhzeq * SolarMass)
      if (logM .lt. minval(PBH_masses) .or. logM .gt. maxval(PBH_masses)) then
        dspbh_limits = 1.d0
      else
        call dshunt(PBH_masses,PBHfileLen,logM,lowerindx)
        dspbh_limits = 10.d0**(((logM-PBH_masses(lowerindx))*PBH_betamaxes(lowerindx+1) +
     &                 (PBH_masses(lowerindx+1)-logM)*PBH_betamaxes(lowerindx)) / 
     &                 (PBH_masses(lowerindx+1)-PBH_masses(lowerindx)))
      endif

      end function dspbh_limits

