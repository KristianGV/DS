**********************************************************************
*** function dssigmav0tot returns the *total* annihilation cross section
*** sigma v at p=0 for neutralino-neutralino annihilation.
*** This is obtained by summing over all implemented 2-body channels
*** (as returned by dssigmav) plus contributions from final states with
*** more particles.
***                                                                         
***  type : interface                                                       
***
***  desc : Total annihilation cross section at $v=0$
***                                                                         
*** Units of returned cross section: cm^3 s^-1
***
*** author: Torsten.Bringmann@fys.uio.no
*** date: 2014-11-14
**********************************************************************

      real*8 function dssigmav0tot()
      implicit none
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dssigmav0tot')

      dssigmav0tot=0.0d0

      return
      end
