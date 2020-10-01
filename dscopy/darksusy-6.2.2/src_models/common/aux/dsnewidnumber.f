	subroutine dsnewidnumber()
************************************************************************
*** This subroutine creates a new unique id number and should be called
*** for each new particle physics model. The Function dsidnumber returns
*** the number and can be used by different functions for optimzation
*** purposes.	
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: December 10, 2014	
*** mod: torsten.bringmann@fys.uio.no, 10/10/2017 (added checksum)
************************************************************************
       implicit none
       include 'dsidnumber.h'
       logical first
       data first/.true./
       save first

       if (first) then
          dsidno=0
          checksum=31415     
          first=.false.
       endif

       dsidno=dsidno+1

       return
       end
