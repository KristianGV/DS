*******************************************************************************
*** Function dsmwimp returns the WIMP mass                                  ***                            
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
*** author: torsten.bringmann@fys.uio.no, 2016-02-06                        ***
*******************************************************************************

      real*8 function dsmwimp()
      implicit none

      include 'dsparticles.h'

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_decayingDM','dsmwimp')

      dsmwimp=mass(kdm)

      return
      end
