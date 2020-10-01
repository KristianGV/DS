*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                         dssem_earth.h                            ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2006-04-21

      integer ean
      parameter(ean=42) ! number of entries in Earth density table
c TB: moved to mpconst.h
c      real*8 r_earth
c      parameter(r_earth=6378.14d3) ! Earth radius in m
      real*8 eadens(ean),eadepth(ean)
      common/dsearth/eadens,eadepth
      save /dsearth/

************************** end of dsearth.h ****************************






