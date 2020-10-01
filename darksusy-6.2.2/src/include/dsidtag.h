*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            idtag.h                               ***
***         this piece of code is needed as a separate file          ***
***               the rest of the code 'includes' idtag.h            ***
***            where reference to a model id tag is needed           ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@teorfys.uu.se) 96-12-12
c          torsten bringmann (torsten.bringmann@fys.uio.no) 2014-04
c          - added particle module tag

* model id tag
      character*20 idtag, moduletag,halomoduletag
      common /dstag/ idtag, moduletag,halomoduletag

      save /dstag/
***************************** end **************************************
