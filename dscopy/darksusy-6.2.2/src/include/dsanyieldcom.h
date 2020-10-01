*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                         dsanyieldcom.h                           ***
***         this piece of code is needed as a separate file          ***
***         the rest of the code 'includes' dsanyieldcom.h           ***
c----------------------------------------------------------------------c
c  Author: Joakim Edsjo (edsjo@fysik.su.se)
c  Date: April, 2014
c....simres - simulation result tables
      real*8 lb,ub,mi,zindex,dz,ndec
      real phiint,phidiff
      integer yieldtype,milow,zn,ntype,nmass,nch,dpn,dbntype,zndb
      integer dbnkind
      character*128 andir
      character anftype
      parameter(zn=200,       ! Number of bins in z (except dbar)
     &          zndb=100,     ! Number of bins in z for dbar
     &          ndec=10.0d0,  ! Number of decades tabulated
     &          ntype=23,     ! Code for highest yield type
     &          nmass=30,     ! Number of mass points
     &          nch=11,       ! Number of fundamental channels
     &          dpn=30,       ! Number of delta-p bins for dbar sims
     &          dbntype=3,    ! Number of yield types (MC gens) for dbar sims
     &          dbnkind=2)    ! Number of yield kinds (yield and error)
      common/ansim/lb(nch),ub(nch),mi(nmass),
     &  zindex(-1:zn,2),dz(-1:zn),
     &  yieldtype(4,ntype),milow(nch),
     &  andir,anftype
      common/ansim2/phiint(0:zn,nmass,nch,ntype),
     &  phidiff(-1:zn,nmass,nch,ntype)

      integer dbct,dbflxk
      real*8 dbp0fit(11:10+dbntype),dbp0low(11:10+dbntype),
     &  dbp0high(11:10+dbntype)
      real*8 dbp0bar,dbdpindex(-1:dpn)
      common/dbsim/dbp0fit,dbp0low,dbp0high,dbp0bar,dbdpindex,dbct,
     &  dbflxk

      real phiintdb(0:zndb,-1:dpn-1,nmass,nch,11:10+dbntype,1:dbnkind),
     &  phidiffdb(-1:zndb,-1:dpn-1,nmass,nch,11:10+dbntype,1:dbnkind)

      common/dbsim2/phiintdb,phidiffdb

      real*8 dbzindex(-1:zndb,2),dbdz(-1:zndb)
      common/dbsim3/dbzindex,dbdz

c....aninfo - option info switches for the program
      integer ansmooth
      common /aninfo/ansmooth

c...logical switches
      logical dsanyieldinitcalled
      common /anlogical/dsanyieldinitcalled

c...ansim3 - Information about simulation channels
      real*8 msim(nch)
      common /ansim3/msim

c save common block
      save /ansim/,/aninfo/,/ansim2/,
     & /ansim3/,/anlogical/

***                                                                 ***
********************** end of dsanyieldcom.h **************************





