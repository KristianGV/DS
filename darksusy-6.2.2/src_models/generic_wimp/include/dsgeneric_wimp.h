!*                         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                         dsgeneric_wimp.h                         ***
!***         this piece of code is needed as a separate file          ***
!***         the rest of the code 'includes' dsgeneric_wimp.h         ***
!c----------------------------------------------------------------------c
!c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 2015

!* For every model, we use the same structure to represent the particle
!* code. HOW this is implemented (i.e. which particle codes are assigned)
!* is up to the model. 
      include 'dsparticles.h'
      include 'dssm.h'

      parameter (numpartspecies=18)  ! # particles in this model (including 17 from SM)

      integer kwimp
      parameter (kwimp=18)

!c... self-conjugate DM particle [1=yes, 2=no]
      integer selfconj 
      common /cselfconj/ selfconj

!c... describe annihilation cross section as a + b*v^2
      real*8 sva, svb
      integer svch
      common /svgen/ sva, svb, svch
     
!c... parameters to describe direct detection     
      real*8 sigsip,sigsin,sigsdp,sigsdn
      common /dd/ sigsip,sigsin,sigsdp,sigsdn 
   
      save /cselfconj/, /svgen/, /dd/ 

!***                                                                  ***
!******************* end of dsgeneric_wimp.h ****************************
