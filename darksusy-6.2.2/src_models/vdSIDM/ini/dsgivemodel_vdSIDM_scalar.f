*******************************************************************************
***  subroutine dsgivemodel_vdSIDM_scalar reads in parameters to descibe    *** 
***  a simple model where a Dirac DM particle coupless to massive scalar    ***
***  mediators, which also couples (with the same strength) to massless     ***
***  dark radiation.                                                        *** 
***  It then transfers these parameters to common blocks                    ***
***                                                                         ***
***  input:                                                                 ***
***                                                                         ***
***    mDM  - DM mass (in GeV)                                              ***
***    mmed - scalar mediator mass (in GeV)                                 ***
***    g    - DM-mediator coupling                                          ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-05-22                                                         ***
*******************************************************************************
      subroutine dsgivemodel_vdSIDM_scalar(mDM,mmed,g)
      implicit none
      include 'dsvdSIDM.h'

      real*8 mDM,mmed,g
c-----------------------------------------------------------------------

c... transfer to common blocks
      mass(kdm)=mDM
      mass(kmed)=mmed
      gDM=g

c... set properties of particles
      DMtype = 1 ! fermion 
      DRtype = 1 ! fermion 
      Mediatortype = 3 ! scalar 

      kdof(kdm) = 2 ! Dirac Fermion 
      spin(kdm) = 0.5d0
      kdof(kdr) = 2 ! Dirac Fermion 
      spin(kdr) = 0.5d0
      kdof(kmed) = 1 ! real scalar 
      spin(kdm) = 0.0d0
    
c... assume massless DR coupling with the same strength
      mass(kdr)=0.0d0
      gDR=gDM


      selfconj=.false. ! for Dirac DM
      ! note on degrees of freedom and genselfconj :
      ! for the relic density of non-self-conjugate particles, one doubles
      ! kdof(kdm); for the cosmic ray routines, one requires selfconj=false
      if (.not.selfconj) kdof(kdm)=2*kdof(kdm)



      end


