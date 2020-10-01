*                         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                         dssilveira_zee.h                         ***
***         this piece of code is needed as a separate file           ***
***         the rest of the code 'includes' dssilveira_zee.h         ***
c----------------------------------------------------------------------c
c  author: Paolo Gondolo 2016
c          Torsten Bringmann -- added yield parameters, outsourced SM part to dssm.h

* For every model, we use the same structure to represent the particle
* code. HOW this is implemented (i.e. which particle codes are assigned)
* is up to the model. 

      include 'dsparticles.h'
      include 'dssmparam.h'
      include 'dsmpconst.h'
      include 'dsidtag.h'
      include 'dssm.h'

      parameter (numpartspecies=18)  ! # particles in this model (including 17 from SM)


c...Particle codes beyond SM particles
      integer khs
      parameter (khs=18)

c...Basic model parameters
      real*8 lambda,mhs
      common /dssilveirazee/ lambda,mhs
      save /dssilveirazee/

c...Higgs width
      integer gammahhow
      common /widthscom/ gammahhow
      save /widthscom/
      
c...annihilation rate and channels
      integer numannch2b, numyieldch_line
      integer numannch2bmax, numyieldch_linemax
      parameter (numannch2bmax=18,numyieldch_linemax=6)
      integer ann_2body(numannch2bmax,6),yieldchannels_line(numyieldch_linemax,2)
      common /annrates/ ann_2body, yieldchannels_line,
     &                  numannch2b, numyieldch_line       

***                                                                  ***
******************* end of dssilveira_zee.h ****************************
