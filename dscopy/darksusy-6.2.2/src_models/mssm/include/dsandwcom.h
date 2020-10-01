*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                          dsandwcom.h                             ***
***         this piece of code is needed as a separate file          ***
***          the rest of the code 'includes' dsandwcom.h             ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se) 97-09-09
c  modified: 01-09-12
c....dwcom - common block needed for the dwdcos optimization
      real*8 parts(6,6,54)
      logical tur(6,6,54)
      common /dwcom1/ parts,tur

      integer slcode
      logical incglue,incgaga,incgaz
      common /dwcom2/ incglue,incgaga,incgaz,slcode

      real*8 prtial(114)
      common /partials/ prtial

c...Store coannihilating particle info (needed by dsandwdcos)
      integer comax
      parameter(comax=50)
      real*8 mco(comax),mdof(comax)
      integer kco(comax),nco
      common /dsancoann/mco,mdof,kco,nco
      
      save /dwcom1/,/dwcom2/,/partials/,/dsancoann/

*$OMP THREADPRIVATE (/dwcom1/,/partials/)


***                                                                 ***
************************** end of dwcom.h *****************************







