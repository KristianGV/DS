      !Calculates ln(1-x) with a Taylor expansion if x < 1e-3, otherwise uses the F90 native ln function     

      double precision function safelog1m(x)

      double precision, intent(IN) :: x

      if (x .gt. 1.d-3) then
        safelog1m = log(1.d0-x)
      else
        !ln(1.d0-x) = (-x-0.5*x^2) by Taylor expansion of ln(z) around z = 0, in cases where this term causes numerical issues
        safelog1m = (-1.d0*x-0.5d0*x*x)
      endif
      
      end function safelog1m
