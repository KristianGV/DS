!N=17 beta calculator for feeder scaling. 
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2012
! 
!Input:  M     M_3 non-Gaussianty parameter (dimensionless)
!        nu    ratio delta/sigma of density contrast to mass variance (dimensionless)
!Output: beta  (dimensionless)

      double precision function dsucmh_BetaF17(M,v)

      implicit none
      include 'dsmpconst.h'

      double precision, intent(IN) :: M,v
      double precision :: dsucmh_BetaF16
      external dsucmh_BetaF16

      dsucmh_BetaF17 = dsucmh_BetaF16(M,v) + 
     &  (46141773101.*M**5.666666666666667*(2027025. - 16216200.*v**2 + 18918900.*v**4 - 7567560.*v**6 + 1351350.*v**8 - 
     &       120120.*v**10 + 5460.*v**12 - 120.*v**14 + v**16))/(3.705077376e11*exp(v**2/2.)*Sqrt(2.*pi))

      end function dsucmh_BetaF17

