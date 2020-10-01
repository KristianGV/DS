**********************************************************************
*** function dssigmav0tot returns the *total* annihilation cross section
*** sigma v at p=0 for WIMP-WIMP annihilation.
*** This is obtained by summing over all implemented 2-body channels
*** plus contributions from final states with more particles (currently none).
***                                                                         
***  type : INTERFACE                                                        
***                                                                         
*** Units of returned cross section: cm^3 s^-1
***
*** author: Paolo Gondolo 2016
**********************************************************************

      real*8 function dssigmav0tot()
      implicit none
      include 'dssilveira_zee.h'

      real*8 sv,dssigmavpartial
      integer ich
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dssigmav0tot')

c      write (*,*) 'PG dssigmav0tot> entered'
      sv=0.d0
      do ich=1,18
        sv=sv+dssigmavpartial(ich,0.d0)
      enddo
      dssigmav0tot=sv
c      write (*,*) 'PG dssigmav0tot> dssigmav0tot=',dssigmav0tot
c      write (*,*) 'PG dssigmav0tot> exited'
      
      return
      end
