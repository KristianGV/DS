*                         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                         dsgeneric_decayingDM.h                         ***
***         this piece of code is needed as a separate file          ***
***         the rest of the code 'includes' dsgeneric_decayingDM.h         ***
c----------------------------------------------------------------------c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 2016

* For every model, we use the same general structure to represent the particle
* code. HOW this is implemented (i.e. which particle codes are assigned)
* is up to the model. 
      include 'dsparticles.h'

      parameter (numpartspecies=18)  ! # particles in this model (including 17 from SM)

      integer kwimp
      parameter (kwimp=1)

      
c...decay rate and channels
      integer numdecch2b, numyieldch_line
      integer numdecch2bmax, numyieldch_linemax
      parameter (numdecch2bmax=17,numyieldch_linemax=6)
      real*8 Gammatot
      real*8 decBR(numdecch2bmax)
      integer dec_2body(numdecch2bmax,6),yieldchannels_line(numyieldch_linemax,2)
      common /decrates/ Gammatot, decBR, dec_2body, yieldchannels_line,
     &                  numdecch2b, numyieldch_line 
      
     

***                                                                  ***
******************* end of dsgeneric_decayingDM.h ****************************
