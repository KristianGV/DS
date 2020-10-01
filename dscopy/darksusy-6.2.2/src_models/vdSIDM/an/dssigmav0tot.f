**********************************************************************
*** function dssigmav0tot returns the *total* annihilation cross section
*** sigma v at p=0 for WIMP-WIMP annihilation.
*** This is obtained by summing over all implemented 2-body channels
*** (as returned by dssigmav) plus contributions from final states with
*** more particles.
***                                                                         
***  type : interface                                                       
***                                                                         
*** Units of returned cross section: cm^3 s^-1
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2018-02-18
**********************************************************************

      real*8 function dssigmav0tot()
      implicit none
      include 'dsvdSIDM.h'

      real*8 dsanwx, pcms

     
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dssigmav0tot')

      pcms=0.0d0
      dssigmav0tot=dsanwx(pcms)/(4.*mass(kdm)**2)

      return
      end
