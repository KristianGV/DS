*******************************************************************************
***  This is a dummy version of dsinit_module.                              ***
***  This file exists so that you can link to DarkSUSY without linking to   ***
***  a particle physics module in case you actually do not need anything    ***
***  from the particle physics module. In the default linking where         ***
***  the particle physics modules are linked before the core, this dummy    ***
***  routine will be replaced with the proper one from the particle physics ***
***  module.                                                                ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : Intialzation of module                                          ***
***                                                                         ***
***  Author: Torsten Bringmann (torsten.bringmann@fys.uio.no)               ***
***  Date: 04/2014                                                          ***
*******************************************************************************
      subroutine dsinit_module
      implicit none

      include 'dsidtag.h'
      
c... naming the module is important for internal consistency checks!
      moduletag='none'

      write(*,*) 'Initialization of particle physics module none complete.'

      return
      end


