!*         -*- mode: fortran -*-
!*######################################################################*
!*                       i n c l u d e     f i l e                      *
!*######################################################################*

!************************************************************************
!***                             dssm.h                               ***
!***  this piece of code contains a *possible* way of representing    ***
!***  various SM parameters that are not already contained in         ***
!***  dsparticles.h, and is provided as a convenience for particle    ***
!***  modules that choose to use this structure (see also the header  ***
!***  of dsinit_sm)                                                   ***
!c----------------------------------------------------------------------c
!c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 11/2016


!c...Particle codes
!* IMPORTANT: the first 12 particle species must be the SM leptons and
!* quarks and must be listed in the order below (because they are 
!* sometimes referenced directly)
      integer knue,ke,knumu,kmu,knutau,ktau,ku,kd,kc
      integer ks,kt,kb,kgamma,kw,kz,kgluon,khsm
      parameter(knue=1)
      parameter(ke=2)
      parameter(knumu=3)
      parameter(kmu=4)
      parameter(knutau=5)
      parameter(ktau=6)
      parameter(ku=7)
      parameter(kd=8)
      parameter(kc=9)
      parameter(ks=10)
      parameter(kt=11)
      parameter(kb=12)
      parameter(kgamma=13)
      parameter(kw=14)
      parameter(kz=15)
      parameter(kgluon=16)
      parameter(khsm=17) ! warning: don't use this in models that provide
                ! more Higgs bosons (like the mssm module: there kh2=18
                ! is the most SM like Higgs!)


!c... Useful derived quantities
      real*8 v0 ! Higgs vev
      common /dssmderived/ v0

!* quark mixings
      complex*16 ckm(3,3)
      real*8 ckms12,ckms23,ckms13,ckmdelta
      common /sckm/ ckm, ckms12,ckms23,ckms13,ckmdelta

!* couplings constants & weak mixings
      real*8 alph3mz,GFermi,alphem
      real*8 s2thw, sinthw, costhw ! tree-level expressions
      real*8 s2wmz,swmz,cwmz       ! full expressions at MS-bar value of mZ
      real*8 g2weak, g2wmz, g3stro, gyweak, gywmz
      common /smcoupling/ alph3mz,GFermi,alphem,s2thw, sinthw, costhw,
     &               s2wmz,swmz,cwmz, g2weak, g2wmz, g3stro, gyweak, gywmz 

!* flag indicating if first or later call of dsralph34loop
      logical first_dsralph34loop
      common/smfirsts/first_dsralph34loop

!* quark mass spectrum
      real*8 mu2gev,md2gev,ms2gev,mcmc,mbmb,mtmt
      character*5 roption    ! option for how to treat running of quarks
      common /smquarkmasses/ mu2gev,md2gev,ms2gev,mcmc,mbmb,mtmt, roption
      

!* save common blocks!
      save /dssmderived/, /sckm/, /smcoupling/ 
      save /smfirsts/, /smquarkmasses/

!***                                                                 ***
!************************** end of dssm.h ******************************
