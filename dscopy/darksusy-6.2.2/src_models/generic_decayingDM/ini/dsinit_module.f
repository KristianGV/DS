*******************************************************************************
***  This is the initialization subroutine for the "generic_decayingDM"     ***
***  module which just contains a simple DM model with a mass,              ***
***  toatl decay rate and branching fractions.                              ***
***                                                                         ***
***  Author: Torsten Bringmann (torsten.bringmann@fys.uio.no)               ***
***  Date: 2016-02-06                                                       ***
*******************************************************************************
      subroutine dsinit_module
      implicit none

c... the following should be included by all dsinit_module versions: 
      include 'dssmparam.h'
      include 'dsmpconst.h'
      include 'dsidtag.h'
c... this is the header file for particle-specific common blocks
c... (always contains dsparticles.h) 
      include 'dsgeneric_decayingDM.h'
      
      integer i

c... naming the module is important for internal consistency checks!
      moduletag='generic_decayingDM'

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
      kdm = numpartspecies ! = 18, set in dsggeneric_decayingDM.h


c... full list of annihilation channels to be taken into account in yield routines
c... (note that the ordering of the channels does not matter!)
       dec_2body(1,1)=23  ! Z Z
       dec_2body(1,2)=23
       dec_2body(2,1)=24  ! W+ W-
       dec_2body(2,2)=-24
       dec_2body(3,1)=12  ! nue nuebar
       dec_2body(3,2)=-12
       dec_2body(4,1)=11  ! e- e+
       dec_2body(4,2)=-11
       dec_2body(5,1)=14  ! numu numubar
       dec_2body(5,2)=-14
       dec_2body(6,1)=13  ! mu- mu+
       dec_2body(6,2)=-13
       dec_2body(7,1)=16  ! nutau nutaubar
       dec_2body(7,2)=-16
       dec_2body(8,1)=15  ! tau- tau+
       dec_2body(8,2)=-15
       dec_2body(9,1)=2   ! u ubar
       dec_2body(9,2)=-2
       dec_2body(10,1)=1   ! d dbar
       dec_2body(10,2)=-1
       dec_2body(11,1)=4   ! c cbar
       dec_2body(11,2)=-4
       dec_2body(12,1)=3   ! s sbar
       dec_2body(12,2)=-3
       dec_2body(13,1)=6   ! t tbar
       dec_2body(13,2)=-6
       dec_2body(14,1)=5   ! b bbar
       dec_2body(14,2)=-5
       dec_2body(15,1)=21   ! gluon gluon
       dec_2body(15,2)=21
       dec_2body(16,1)=22   ! gamma gamma
       dec_2body(16,2)=22
       dec_2body(17,1)=22   ! gamma Z
       dec_2body(17,2)=23 

       numdecch2b=numdecch2bmax

c... now add angular momentum information about channels
       do i=1,numdecch2b  
         dec_2body(i,3)=0  ! 2*J
         dec_2body(i,4)=-1  ! CP
c... TB FIXME: the following is true for fermions and scalars, 
c...           but needs to be updated for everything else  !!!       
         dec_2body(i,5)=0  ! 2*S
         dec_2body(i,6)=0  ! 2*L
         if (i.eq.2) then
           dec_2body(i,5)=2  ! 2*S
           dec_2body(i,6)=2  ! 2*L
         endif    
       enddo      


c... full list of annihilation channels to be taken into account for monochromatic
c... yield. Convention: PDG code for CR particle in question first!

       yieldchannels_line(1,1)=22   ! gamma gamma
       yieldchannels_line(1,2)=22
       yieldchannels_line(2,1)=22   ! gamma Z
       yieldchannels_line(2,2)=23
       yieldchannels_line(3,1)=-11  ! e+ e-
       yieldchannels_line(3,2)=11
       yieldchannels_line(4,1)=12   ! nue nuebar
       yieldchannels_line(4,2)=-12
       yieldchannels_line(5,1)=14   ! numu numubar
       yieldchannels_line(5,2)=-14
       yieldchannels_line(6,1)=16   ! nutau nutaubar
       yieldchannels_line(6,2)=-16

       numyieldch_line=numyieldch_linemax

      write(*,*) 'Initialization of particle physics module generic_decayingDM complete.'

      return
      end


