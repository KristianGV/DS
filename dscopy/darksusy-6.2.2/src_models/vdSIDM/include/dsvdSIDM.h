*                         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                         dsvdSIDM.h                         ***
***         this piece of code is needed as a separate file          ***
***         the rest of the code 'includes' dsvdSIDM.h         ***
c----------------------------------------------------------------------c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 2018

* For every model, we use the same structure to represent the particle
* code. HOW this is implemented (i.e. which particle codes are assigned)
* is up to the model. 
      include 'dsparticles.h'
      include 'dssm.h'

      parameter (numpartspecies=20)  ! # particles in this model (including 17 from SM)

c...particle codes for DM, DR, and mediator 
      integer kdr,kmed  
      parameter (kdr=19,kmed=20)

c... self-conjugate particle [true=yes, false=no]
      logical selfconj
      common /cselfconj/ selfconj

c... decoupling temperature from SM heat bath [GeV]
      real*8 DSTdec
      common /DecTemp/ DSTdec

c... particle content
      real*8 gDM, gDR  ! DM and DR couplings
      integer DMtype, DRtype, Mediatortype  ! 1=fermion,2=vector,3=scalar 
      common /darksector/ gDM, gDR, DMtype, DRtype, Mediatortype
     
***                                                                  ***
******************* end of dsvdSIDM.h ****************************
