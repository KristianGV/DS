***********************************************************************
*** This is just an empty dummy function. If we don't replace any
*** functions, link to this one so that the archive is not empty
*** which would throw an error.      
***      
*** type : interface (not interface in the usual sense, but needed)      
***      
***********************************************************************

      real*8 function empty_dummy()
      implicit none

      empty_dummy=0.d0
      
      return
      end