*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dssvcom.h                              ***
***         this piece of code is needed as a separate file          ***
***              the rest of the code 'includes' dssvcom.h           ***
c----------------------------------------------------------------------c
c  author: paolo gondolo 130812

c.... parameters for annihilation cross section and branching ratios 
c     at zero velocity
      real*8 mx,sigmav,sigv,wtot
      common/svcom/mx,sigmav,sigv(29),wtot

c save common blocks
      save /svcom/

***                                                                 ***
************************ end of dssvcom.h ******************************



