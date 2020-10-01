!Radiation transfer function.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input:    x    wavenumber * horizon scale (dimensionless)
!Output:   T    transfer function (dimensionless)

      double precision function dsucmh_transFunc_rad(x)

      implicit none

      double precision, intent(IN) :: x
      double precision :: theta, theta_sq, sinc
      double precision, parameter :: root3 = sqrt(1.d0/3.d0)

      theta = x*root3
      theta_sq = theta*theta

      if (theta .lt. 0.1d0) then

        dsucmh_transFunc_rad = 1.d0 - theta_sq*(0.1d0 - theta_sq/280.d0)

      else

        dsucmh_transFunc_rad = 3.d0/theta_sq*(sinc(theta) - cos(theta))

      endif      

      end function dsucmh_transFunc_rad

