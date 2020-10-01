*******************************************************************************
***  This is the initialization subroutine for the "vdSIDM" module    ***
***  which just contains a simple WIMP model with a mass, annihilation      ***
***  branching fractions, cross sections etc.                               ***
***                                                                         ***
***  Author: Torsten Bringmann (torsten.bringmann@fys.uio.no)               ***
***  Date: 21/06/2016                                                       ***
*******************************************************************************
      subroutine dsinit_module
      implicit none

c... the following should be included by all dsinit_module versions: 
      include 'dssmparam.h'
      include 'dsmpconst.h'
      include 'dsidtag.h'
c... this is the header file for particle-specific common blocks
c... (always contains dsparticles.h) 
      include 'dsvdSIDM.h'
      
c... naming the module is important for internal consistency checks!
      moduletag='vdSIDM'


c... We now initialize the SM physics in the "standard" way by filling
c... in particular the first 17 entries of the common block variables
c... in dsparticles.h. Note that this is a convenient, but not required way
c... of including the SM (it could be done in a completely user-defined way).
c... Additionally, any particular setting can be overridden by simply 
c... re-assigning the corresponding values after this call
      call dsinit_sm

c... set electron mass in dsmpconst.h -- this is always needed by the CR routines
      m_e=mass_e_def

c... DM identifer. Already defined in dsparticles.h, and hence cannot be
c... redefined in dsvdSIDM.h
      kdm = 18 

c... The other internal particle (integer) codes are defined in dsvdSIDM.h

      pname(kdm) = 'DM'
      pname(kdr) = 'DR'
      pname(kmed) = 'dark_mediator'

c... We assume that the Dark Sector and the SM have been in thermal contact
c... at a very high temperature ('infinity')
      DSTdec = 1.0d5 ! 100 TeV


      write(*,*) 'Initialization of particle physics module vdSIDM complete.'

      return
      end


