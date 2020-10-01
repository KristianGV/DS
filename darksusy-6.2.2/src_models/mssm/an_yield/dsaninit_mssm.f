*****************************************************************************
***   subroutine dsaninitmodel initializes model dependent parts of 
***   annihilation routines. Specifically channel numbers, scalar decays
***   etc. Needs to be called once (typically from dsinit_module)
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: May, 2014
*** modified: 2014-11-11 Torsten Bringmann 
***           (added explicit annihilation channel list with PDG codes)
***           2016-02-06 Torsten Bringmann 
***           (added contributions from electroweak internal bremsstrahlung)
*****************************************************************************

      subroutine dsaninit_mssm
      implicit none
      include 'dsio.h'
      include 'dsanyieldmodelcom.h'

c------------------------ variables ------------------------------------

      integer i


c...to go from new channel numbers to compressed channels
c...these compressed channel numbers are the indices for the arrays
c...New channel numbers as of July 4, 2013
      data chcomp/
     &  11*0,9,8,3*0,10,0,11,2,1,4,3,6,5,7,3*0/
      data chi2pdg/1,2,3,4,5,6,21,24,23,13,15/
      data sechi2pdg/1,2,3,4,5,6,21,24,23,13,15,12,14,16/
c...JE FIXME. chi2pdg is just temporary. The fundamental channel numbers should
c...not be used here at all. The same for sechi2pdg

c----------------------------------------- set-up common block variables


c...Minimal branching fraction for scalar decay to include in an routines
c         ansbrmin=0.d0 ! include all
       ansbrmin=1.d-3 ! very good approximation
       anexhi=0 ! include decaying scalars

c... full list of annihilation channels to be taken into account in yield routines
c... (note that the ordering of the channels does not matter!)
       anch_2body(1,1)=25  ! h h
       anch_2body(1,2)=25
       anch_2body(2,1)=25  ! h H
       anch_2body(2,2)=35
       anch_2body(3,1)=35  ! H H
       anch_2body(3,2)=35
       anch_2body(4,1)=36  ! A A
       anch_2body(4,2)=36
       anch_2body(5,1)=25  ! h A
       anch_2body(5,2)=36
       anch_2body(6,1)=35  ! H A
       anch_2body(6,2)=36
       anch_2body(7,1)=37  ! H+ H-
       anch_2body(7,2)=-37
       anch_2body(8,1)=23  ! Z h
       anch_2body(8,2)=25
       anch_2body(9,1)=23  ! Z H
       anch_2body(9,2)=35
       anch_2body(10,1)=23  ! Z A
       anch_2body(10,2)=36
       anch_2body(11,1)=24  ! W+ H-
       anch_2body(11,2)=-37
       anch_2body(12,1)=23  ! Z Z
       anch_2body(12,2)=23
       anch_2body(13,1)=24  ! W+ W-
       anch_2body(13,2)=-24
       anch_2body(14,1)=12  ! nue nuebar
       anch_2body(14,2)=-12
       anch_2body(15,1)=11  ! e- e+
       anch_2body(15,2)=-11
       anch_2body(16,1)=14  ! numu numubar
       anch_2body(16,2)=-14
       anch_2body(17,1)=13  ! mu- mu+
       anch_2body(17,2)=-13
       anch_2body(18,1)=16  ! nutau nutaubar
       anch_2body(18,2)=-16
       anch_2body(19,1)=15  ! tau- tau+
       anch_2body(19,2)=-15
       anch_2body(20,1)=2   ! u ubar
       anch_2body(20,2)=-2
       anch_2body(21,1)=1   ! d dbar
       anch_2body(21,2)=-1
       anch_2body(22,1)=4   ! c cbar
       anch_2body(22,2)=-4
       anch_2body(23,1)=3   ! s sbar
       anch_2body(23,2)=-3
       anch_2body(24,1)=6   ! t tbar
       anch_2body(24,2)=-6
       anch_2body(25,1)=5   ! b bbar
       anch_2body(25,2)=-5
       anch_2body(26,1)=21   ! gluon gluon
       anch_2body(26,2)=21
       anch_2body(27,1)=10000   ! (not used)
       anch_2body(27,2)=10000
       anch_2body(28,1)=22   ! gamma gamma
       anch_2body(28,2)=22
       anch_2body(29,1)=22   ! gamma Z
       anch_2body(29,2)=23 

c... now add angular momentum information about channels
       do i=1,numanch2b  
         anch_2body(i,3)=0  ! 2*J
         anch_2body(i,4)=-1  ! CP
c... TB FIXME: the following is true for fermions and scalars, 
c...           but needs to be updated for everything else  !!!       
         anch_2body(i,5)=0  ! 2*S
         anch_2body(i,6)=0  ! 2*L
         if (i.eq.13) then
           anch_2body(i,5)=2  ! 2*L
           anch_2body(i,6)=2  ! 2*S
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

c... contributions from internal bremsstrahlung of any bosons in the theory:
      call dsibset('dynamic')  ! do include photon IB as default

      call dsib2set('off')     ! do *not* include SU(2) IB as default (too time consuming)
         
      return

      end
