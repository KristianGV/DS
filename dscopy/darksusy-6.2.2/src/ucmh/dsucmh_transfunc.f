!Dark matter transfer function.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input:    x    wavenumber * horizon scale (dimensionless)
!Output:   T    transfer function (dimensionless)

      double precision function dsucmh_transFunc(x)

      implicit none

      double precision, intent(IN) :: x
      double precision :: theta, theta_sq, finalterm, mci, msi, sinc
      double precision, parameter :: root3 = sqrt(1.d0/3.d0)

      theta = x*root3
      theta_sq = theta*theta

      if (theta .lt. 0.1d0) then

        finalterm = -0.5d0 + theta_sq*(2.5d-2 - theta_sq/1680.d0)
        
      else

        finalterm = 3.d0/theta_sq*(sinc(theta) - 1.d0)

      endif

      call ModCosInt(theta,mci,msi)
      dsucmh_transFunc = -6.d0*mci + finalterm

      end function dsucmh_transFunc



