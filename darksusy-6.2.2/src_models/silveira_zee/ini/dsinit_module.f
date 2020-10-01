*******************************************************************************
***  This is the initialization subroutine for the "silveira_zee" module    ***
***                                                                         ***
***  Author: Paolo Gondolo, Torsten Bringmann                               ***
***  Date: 05/19/2016                                                       ***
*******************************************************************************
      subroutine dsinit_module
      implicit none

c... this is the header file for particle-specific common blocks
      include 'dssilveira_zee.h'
      include 'dsio.h'
      
      integer i

c... naming the module is important for internal consistency checks!
      moduletag='silveira_zee'


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

c... Particles in addition to SM-no-Higgs
      pname(kdm)='singlet-higgs'
      mass(kdm)=0.d0
      width(kdm)=0.d0
      spin(kdm)=0.d0
      kdof(kdm)=1
      

c... full list of annihilation channels to be taken into account in yield routines
c... (note that the ordering of the channels in principle does not matter -- but it
c...  is presently hard-coded as "ichannel" in several routines, e.g. dssigmavpartial!)
c... FIXME: eventually, the "ichannel" should go...


       ann_2body(1,1)=12  ! nue nuebar
       ann_2body(1,2)=-12
       ann_2body(2,1)=11  ! e- e+
       ann_2body(2,2)=-11
       ann_2body(3,1)=14  ! numu numubar
       ann_2body(3,2)=-14
       ann_2body(4,1)=13  ! mu- mu+
       ann_2body(4,2)=-13
       ann_2body(5,1)=16  ! nutau nutaubar
       ann_2body(5,2)=-16
       ann_2body(6,1)=15  ! tau- tau+
       ann_2body(6,2)=-15
       ann_2body(7,1)=2   ! u ubar
       ann_2body(7,2)=-2
       ann_2body(8,1)=1   ! d dbar
       ann_2body(8,2)=-1
       ann_2body(9,1)=4   ! c cbar
       ann_2body(9,2)=-4
       ann_2body(10,1)=3   ! s sbar
       ann_2body(10,2)=-3
       ann_2body(11,1)=6   ! t tbar
       ann_2body(11,2)=-6
       ann_2body(12,1)=5   ! b bbar
       ann_2body(12,2)=-5
       ann_2body(13,1)=22   ! gamma gamma
       ann_2body(13,2)=22
       ann_2body(14,1)=24  ! W+ W-
       ann_2body(14,2)=-24
       ann_2body(15,1)=23  ! Z Z
       ann_2body(15,2)=23
       ann_2body(16,1)=21   ! gluon gluon
       ann_2body(16,2)=21
       ann_2body(17,1)=25   ! H H
       ann_2body(17,2)=25 
       ann_2body(18,1)=22   ! gamma Z
       ann_2body(18,2)=23 

       numannch2b=numannch2bmax

c... now add angular momentum information about channels
       do i=1,numannch2b  
         ann_2body(i,3)=0  ! 2*J
         ann_2body(i,4)=1  ! CP
c... TB FIXME: the following needs to be updated!!!       
         ann_2body(i,5)=0  ! 2*S
         ann_2body(i,6)=0  ! 2*L
         if (i.eq.2) then
           ann_2body(i,5)=2  ! 2*S
           ann_2body(i,6)=2  ! 2*L
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


c... There are various options for how to calculate the partial Higgs decay widths
       gammahhow = 1 ! Hdecay tables [DEFAULT]
c       gammahhow = 2 ! Dittmair tables [as used in Cline+]           

      write(*,*) 'Initialization of particle physics module silveira_zee ' //
     &           'complete.'

      return
      end


