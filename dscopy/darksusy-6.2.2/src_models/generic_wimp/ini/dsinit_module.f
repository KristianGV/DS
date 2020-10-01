*******************************************************************************
***  This is the initialization subroutine for the "generic_wimp" module    ***
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
      include 'dsgeneric_wimp.h'
      
c... naming the module is important for internal consistency checks!
      moduletag='generic_wimp'


c... We now initialize the SM physics in the "standard" way by filling
c... in particular the first 17 entries of the common block variables
c... in dsparticles.h. Note that this is a convenient, but not required way
c... of including the SM (it could be done in a completely user-defined way).
c... Additionally, any particular setting can be overridden by simply 
c... re-assigning the corresponding values after this call
      call dsinit_sm

c... set electron mass in dsmpconst.h -- this is always needed by the CR routines
      m_e=mass_e_def


c... DM identifer. For convenience, we simply choose the largest particle number
c... accessible -- though in general any number > 17 (SM dof) will do
      kdm = numpartspecies ! = 18, set in dsgeneric_wimp.h

      selfconj = -100 ! make sure that it is not accidentally decided whether
                      ! DM is self-conjugate or not (but only after setting up
                      ! a model).

      write(*,*) 'Initialization of particle physics module generic_wimp complete.'

      return
      end


