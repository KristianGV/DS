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
*** date: 2014-11-14
**********************************************************************

      real*8 function dssigmav0tot()
      implicit none
      include 'dsgeneric_wimp.h'

      real*8 dsmwimp, dsmass

     
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dssigmav0tot')

      dssigmav0tot=0.d0
      if (dsmass(svch).gt.dsmwimp()) return ! final state can't be heavier
                                            ! than DM particle for v=0 !

      dssigmav0tot=sva
    
      return
      end
