**********************************************************************
*** function dssigmav returns the *total* annihilation cross section
*** sigma v at p=0 for WIMP-WIMP annihilation.
*** Here v is defined to be the relative velocity of one particle in
*** the frame of the other.
***
*** Notice that in 1306.4710, their sigma v_rel is related to the invariant rate W via
***  W = s sigma v_rel, from comparison with the invariant rate in their formula for <sigma v>
*** From the general relation
***  W = 4 p1.p2 sigma v = 2 (s-m1^2-m2^2) sigma v = 4 (p^2+sqrt(p^2+m1^2)*sqrt(p^2+m2^2)) sigma v
*** one obtains
***  sigma v = s/4/(p^2+sqrt(p^2+m1^2)*sqrt(p^2+m2^2)) sigma v_rel
*** and for m1=m2==mx
***  sigma v = W/4/(mx^2+2*p^2)
***          = (mx^2+p^2)/(mx^2+2*p^2) sigma v_rel
***                                                                         
***  type : INTERFACE                                                        
***                                                                         
*** Units of returned cross section: cm^3 s^-1
***
*** author: Paolo Gondolo 2016
**********************************************************************

      function dssigmav(p)
      implicit none
      real*8 dssigmav,p
      real*8 ans,dssigmavpartial
      integer ich
      include 'dssilveira_zee.h'

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dssigmav')

c      write (*,*) 'PG dssigmav> entered with p=',p
      ans=0.d0
      do ich=1,18
        ans=ans+dssigmavpartial(ich,p)
      enddo
      dssigmav=ans
c      write (*,*) 'PG dssigmav> dssigmav=',dssigmav
c      write (*,*) 'PG dssigmav> exited'
      
      return
      end
