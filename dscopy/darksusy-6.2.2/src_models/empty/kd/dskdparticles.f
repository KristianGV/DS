*******************************************************************************
*** Prepares integration of Boltzmann equation and has to be called         ***
*** before dskdboltz (called by dskdtkd)                                    ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : Initalization of kinetic decoupling for module                  ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date  : 2015-06-11                                                      ***
*******************************************************************************
      subroutine dskdparticles
      implicit none

      include 'dskdcom.h'

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dskdparticles')


c... are there any BSM scattering partners?
      nBSMscatt=0
      nBSMscattlight=0

c... set up potential resonances and order them
      nKDres=0

      return

      end





        

