	integer function dsidnumber()
************************************************************************
*** This function returns an integer number that is required to be
*** unique for each new particle physics model. This routine is used
*** by other routines both in src_models/ and in src/ to be able to optimze
*** the code by only recalculting quantities when needed (e.g. for a new
*** model).
***
***    Type: interface function
***	
*** The subroutine dsnewidnumber creates a new one and is called for each
*** new particle physics model
***
***  NB: You can first check with dsidnumberset whether this function can safely
***      be used!
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: December 10, 2014	
************************************************************************
      implicit none
      include 'dsidnumber.h'

      dsidnumber=dsidno

      return
      end
