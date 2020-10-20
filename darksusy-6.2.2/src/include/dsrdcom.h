*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsrdcom.h                               ***
***         this piece of code is needed as a separate file          ***
***            the rest of the code 'includes' dsrdcom.h               ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@cfpa.berkeley.edu) 96-03-22
c  modified: 98-03-03
c  modified: 07-12-29 pg switch for choosing dof in early universe
c....rdmgev - relic and coannihilating masses in gevs
      integer tharsi
      parameter (tharsi=1000)  
      real*8 mco(tharsi),mdof(tharsi),rgev(tharsi),rwid(tharsi)
      integer nco,nres
      common /rdmgev/ mco,mdof,rgev,rwid,nco,nres
c....rdpth - momentum threshold to take into account when integrating
      integer nth,incth(0:tharsi+1)
      real*8 pth(0:tharsi+1)
      common /rdpth/ pth,incth,nth
c....rddof - degrees of freedom (5*300 real, 4 integer)
      real*8 tgev(300),fh(300),fg(300),fe(300),fp(300)
      integer nf,khi,klo,dofcode
      common /rddof/ tgev,fh,fg,fe,fp,nf,khi,klo,dofcode
c....rderrors - error flags (2 integer)
      integer rderr,rdwar,rdinit
      common /rderrors/ rderr,rdwar,rdinit
c....rdlun - logical units (2 integer)
      integer rdlulog,rdluerr
      common /rdlun/ rdlulog,rdluerr
c....rdrate - tabulated invariant rate (3*nrmax real,3*nrmax+3)
      integer nrmax,ispl
      parameter (nrmax=3000)
      real*8 pp(0:nrmax),yy(0:nrmax),yy2(0:nrmax)
      integer indx(0:nrmax),nlo,nhi,nr
      common /rdrate/ pp,yy,yy2,indx,ispl(0:tharsi+1),nlo,nhi,nr
c...rdinfo - info
      character*20 rdtag
      common /rdinfo/rdtag
c...rdpars - modifiable parameters
      real*8 cosmin,waccd,dpminr,dpthr,wdiffr,wdifft,
     &  hstep,hmin,compeps,xinit,xfinal,umax,cfr,pmax
      common /rdpars/cosmin,waccd,dpminr,dpthr,wdiffr,wdifft,
     &  hstep,hmin,compeps,xinit,xfinal,umax,cfr,pmax
c...rdpadd - how points are added
      real*8 pdivr,dpres
      integer nlow,nhigh,npres,nthup,cthtest,spltest
      common /rdpadd/pdivr,dpres,nlow,nhigh,npres,nthup,
     &  cthtest,spltest
c...rdlims
      real*8 plow(1:2*tharsi),phigh(1:2*tharsi)
      integer nlim
      common /rdlims/plow,phigh,nlim
c...rdswitch - switches
      integer thavint,rdprt
      common /rdswitch/thavint,rdprt
c...rdtime - time of calculation
      real*8 rdt_max,rdt_start,rdt_end
      common /rdtime/rdt_max,rdt_start,rdt_end
c...rdtabulation - tabulation of dsanwx in dsrdwx (separate from rdrate)
      real*8 ppp(0:nrmax),rdwx(0:nrmax)
      integer nrd
      common /rdtabulation/ppp,rdwx,nrd
c...rdresopt - optimize how resoances are handdled
c......used by fast=20 option
      logical resinc(tharsi),resfit(tharsi)
      real*8 resalpha(tharsi),resgamma(tharsi),resnorm(tharsi)
      real*8 resconst(tharsi)
      common /rdresopt/ resalpha,resgamma,resnorm,resconst,
     &   resinc,resfit
c....save common blocks
      save /rdmgev/,/rddof/,/rderrors/,/rdlun/,/rdrate/,/rdpth/,
     &   /rdinfo/,/rdpars/,/rdpadd/,/rdlims/,/rdswitch/,/rdtime/,
     &   /rdtabulation/,/rdresopt/
     
***                                                                 ***
************************** end of dsrdcom.h *****************************







