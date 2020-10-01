!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                         dsparticles.h                            ***
!***         this piece of code is needed as a separate file          ***
!***             the rest of the code 'includes' this file            ***
!c----------------------------------------------------------------------c
!c  author: paolo gondolo (paolo.gondolo@utah.edu) 2013
!c  modified: Torsten Bringmann 04/2014
!c  modified: Paolo Gondolo 11/2016 added particle spin
!c  modified: TB 04/2019 added further quantum numbers

!* maximum number of particle species
      integer maxnumpartspecies, numpartspecies   ! numpartspecies is set as 
      parameter (maxnumpartspecies=255)           ! param in include/[module].h
!* PDG particle codes
      integer pdgcode(0:maxnumpartspecies)
      common /pdgcodes/pdgcode
!* particle names
      character*20 pname(0:maxnumpartspecies)
      common /pnames/pname
!* particle pole masses
      real*8 mass(0:maxnumpartspecies)
      common /pmasses/ mass
!* particle decay widths
      real*8 width(0:maxnumpartspecies)
      common /pwidths/ width
!* particle spins
      real*8 spin(0:maxnumpartspecies)
      common /pspins/ spin
!* number of spin+color states
      integer kdof(0:maxnumpartspecies)
      common /pdofs/ kdof
!* other quantum numbers
      real*8 ncolor(0:maxnumpartspecies),wiso3(0:maxnumpartspecies)
      real*8 echarg(0:maxnumpartspecies)
      common /qnum/ ncolor,wiso3,echarg
!* dark matter particle
      integer kdm
      common /pdarkmatter/ kdm
!c save common blocks
      save /pdgcodes/,/pnames/,/pmasses/,/pwidths/,/pspins/,/pdofs/, /qnum/
      save /pdarkmatter/
!***                                                                 ***
!********************* end of dsparticles.h ****************************
