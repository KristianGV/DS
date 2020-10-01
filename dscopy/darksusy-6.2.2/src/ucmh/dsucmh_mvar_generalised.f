!Mass variance calculation for generalised power spectra.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input: log10(power in curvature perturbations) (dimensionless)
!Output: mass variance (dimensionless)

      double precision function dsucmh_mvar_generalised(log10P_R)

      implicit none

      double precision, intent(IN) :: log10P_R
      double precision, parameter :: ln10 = log(10.d0)        

      dsucmh_mvar_generalised = 0.907d0 * exp(ln10*log10P_R) !Here the prefactor assumes n=1 locally -- it should be recalculated for n ne 1

      end function dsucmh_mvar_generalised


