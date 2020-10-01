!=======================================================================
! Calculates sinc(x) = sin(x)/x, to double precision.
! 
! Author: Pat Scott (p.scott@imperial.ac.uk)
! Date: Some time in 2010/11.
! Added to DarkSUSY: Oct 13 2014
! 
!=======================================================================


      double precision function sinc(x)
 
      implicit none

      double precision, intent(IN) :: x

      if (x .lt. 1.d-3) then
        sinc = 1.d0 - x*x*0.16666666666666666666d0
      else
        sinc = sin(x)/x
      endif

      end function sinc

