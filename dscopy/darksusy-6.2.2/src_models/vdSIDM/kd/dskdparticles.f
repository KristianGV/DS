*******************************************************************************
*** Prepares integration of Boltzmann equation and has to be called         ***
*** before dskdboltz (called by dskdtkd)                                    ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date  : 2015-06-11                                                      ***
*******************************************************************************
      subroutine dskdparticles
      implicit none

      include 'dskdcom.h'

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dskdparticles')

c... are there any BSM scattering partners?
      nBSMscatt = 2  ! scattering with DR + mediators
      nBSMscattlight = 1 ! only DR
      BSMscattfermion(1) = .true. ! DR is a fermion
      BSMscattfermion(2) = .false. ! The mediator is not a fermion      

c... set up potential resonances and order them
      nKDres=0

      return

      end





        

