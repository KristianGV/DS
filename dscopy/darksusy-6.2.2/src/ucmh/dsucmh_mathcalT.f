!Returns fitting function \mathcal{T} from Weinberg's cosmology book
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input:   k             wavenumber (Mpc^-1)
!Output   \mathcal{T}   fitting func (dimensionless)
      double precision function dsucmh_mathcalT(k)

      implicit none
      include 'dsmpconst.h'

      double precision, intent(IN) :: k
      double precision, parameter :: root2 = sqrt(2.d0)
      double precision :: kap, transkap, intermed1, intermed2

      kap = root2*k/keq
      transkap = (0.124d0*kap)**2

      intermed1 = (1.d0 + (1.257d0*kap)**2 + (0.4452*kap)**4 + (0.2197*kap)**6)
      intermed2 = (1.d0 + (1.606d0*kap)**2 + (0.8568*kap)**4 + (0.3927*kap)**6)
      dsucmh_mathcalT = log(1.d0+transkap) / transkap * sqrt(intermed1/intermed2)

      end function dsucmh_mathcalT     

