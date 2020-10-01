!Integrand for dsucmh_beta
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
! 
!Input:  delta density contrast (dimensionless)
!Output: integrand  (dimensionless)

      double precision function betaintegrand(delta)


      implicit none

      include 'dsucmh.h'

      double precision, intent(IN) :: delta
 
      betaintegrand = exp(-delta*delta/twoucmhsigma_sq)

      end function betaintegrand      


