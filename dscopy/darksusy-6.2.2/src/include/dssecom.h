*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dssecom.h                               ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@fysik.su.se), 2000-08-16
c  modified: Joakim Edsjo (edsjo@fysik.su.se), 2015-06-11

c...Intermediate and final results (accessible if needed, main results
c...are returned by main routines).      
      real*8 tausu,csu,tauea,cea,searateea,searatesu
      common /seres/tausu,csu,tauea,cea,searateea,searatesu

      real*8 semx,serho ! needed for foveru where we don't want extra arguments
      common /seaux/semx,serho
      
      real*8 veout
      integer secalcmet,selambda,setab,sesunacc,sejup
      common /separa/veout,secalcmet,selambda,setab,sesunacc,sejup

      save /seres/,/separa/,/seaux/

c...Internal common blocks      
      integer zx,ix,ax,sctype
      common /seinternal2/zx,ix,ax,sctype
      save /seinternal2/

      real*8 gpsx,gnsx,gpax,gnax,wwx
      common /seinternal3/gpsx,gnsx,gpax,gnax,wwx
      save /seinternal3/

*************************** end of dssecom.h *****************************






