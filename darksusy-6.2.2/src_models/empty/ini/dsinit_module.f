*******************************************************************************
***  This is the initialization subroutine for the "empty" module which     ***
***  essentially does not contain any actual particle model information at  ***
***  all.                                                                   ***
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
c... the following should be included by all dsinit_module versions: 
      include 'dssmparam.h'
      include 'dsmpconst.h'
      include 'dsidtag.h'
c... this is the header file for particle-specific common blocks
c... (always contains dsparticles.h) 
      include 'dsempty.h'
      
      integer i

c... Initialize common block variables in dsparticles.h
c... Unlike for any "real" model, we don't need to assign anything here
      do i=0,maxnumpartspecies
        pdgcode(i)=0
        pname(i)='no_name_assigned'
        mass(i)=0.d0
        width(i)=0.d0
        spin(i)=0.d0
        kdof(i)=1
      enddo

c... naming the module is important for internal consistency checks!
      moduletag='empty'

c... but we always want to set the electron mass that appears in dsmpconst.h 
c... (needed in propagation routines)
      m_e=mass_e_def

      write(*,*) 'Initialization of particle physics module empty complete.'

      return
      end


