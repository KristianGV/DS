*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                         dsseyieldcom.h                           ***
***         this piece of code is needed as a separate file          ***
***           the rest of the code 'includes' dsseyieldcom.h         ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se) 96-03-22
c  modified: 00-08-16, 08-04-01
c....sesim - simulation result tables
c...semax is the number of tables to maximally have in memory at the
c...same time. For most users, the default below is OK, but if you
c...want to use all tables simultaneously, semax has to be increased
c...to 52.
      integer zn2,thn
      parameter(zn2=200)        ! WimpAnn energy bins
      parameter(thn=90)         ! WimpAnn theta bins
      integer semax,selast(2),senm
      parameter(semax=6) ! number of tables to load simultaneously
      parameter(senm=28) ! number of masses for tabulation
      real*8 lb(14),ub(14),mi(senm),thindex(-1:thn,2),
     &   zindex(-1:zn2,2),dth(-1:thn),dz(-1:zn2)
      real phiint(-1:thn,0:zn2,senm,13,2,semax),
     &     phidiff(-1:thn,-1:zn2,senm,13,2,semax),
     &     phimixed(-1:thn,-1:zn2,senm,13,2,semax)
      integer milow(14),yload(2,26)
      integer kind2ki(3)
      character*128 sedir
      character seftype
      character*40 sebase
      common/sesim/lb,ub,mi,thindex,zindex,dth,dz,
     &  yload,selast,milow,kind2ki,
     &  sedir,seftype,sebase
      common/sesim2/ phiint,phidiff,phimixed


c...segen - general stuff
      integer ch2chi(29),chi2chii(14),chii2chi(13),chi2ch(14)
      common /segen/ch2chi,chi2chii,chii2chi,chi2ch

c...seinfo - tag etc.
      integer seerr,seistat
      common/seinfo/seerr,seistat

c...sesim3 - Information about simulation channels
      real*8 msim(14)
      common /sesim3/msim

c save common block
      save /sesim/,/seinfo/,/sesim2/,/sesim3/,/segen/
      
***                                                                 ***
************************** end of muoncom.h ***************************






