*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dsib3com.h                             ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsib3com.h           ***
c----------------------------------------------------------------------c
c by torsten.bringmann@fys.uio.no, 2015-06-01



* dNdE tables
      integer ntype, nmass, nenergy
      parameter (ntype=16,nmass=110,nenergy=1100)
      integer ehi(7:18,ntype),elo(7:18,ntype), mhi(7:18,ntype), mlo(7:18,ntype)
      real*8 dNdEdat(7:18,ntype,nmass,nenergy), xdat(7:18,ntype,nenergy)
      real*8 fitmass(7:18,ntype,nmass)
      common /dNdE3com/ dNdEdat,xdat,fitmass,ehi,elo, mhi, mlo

* parameters for interpolation of 3-body spectra, as in 1510.02473
      logical fitexists(7:12,ntype)
      real*8 ci(1:3,7:12,ntype), ni(1:3,7:12,ntype), Rprime(7:12,ntype,nmass)
      common /fit3bcom/ rprime, ci, ni, fitexists
