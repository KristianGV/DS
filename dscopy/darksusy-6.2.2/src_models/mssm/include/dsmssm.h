*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsmssm.h                              ***
***         this piece of code is needed as a separate file          ***
***              the rest of the code 'includes' dsmssm.h            ***
c----------------------------------------------------------------------c
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c    10-nov-95 complex vertex constants
c  modified by joakim edsjo 97-02-11, 13-06-11
c  derived from dssusy.h 2008-01-22 by paolo gondolo
c  modified by paolo gondolo (paolo@physics.utah.edu) 08-01-30eor
c  modified by pat scott 08-12-03, 09-10-20, 11-04-27
c   for including running masses
c  modified by paolo gondolo 13-10-02
c  modified by torsten bringmann 2015, 2016, 2019

      include 'dsparticles.h'
      include 'dssm.h'

      parameter (numpartspecies=50)  ! # particles in this model (including 17 from SM)

* mass spectrum
      real*8 runmass(0:numpartspecies)
      common /mspctm/ runmass
c save common blocks
      save /mspctm/

c
c     mssm particles
c

*   particle codes
*   (SM particle codes are decalred in dssm.h)

c  knu=(nue,numu,nutau)   kl=(e,mu,tau)    kqu=(u,c,t)    kqd=(d,s,b)
      integer knu(3),kl(3),kqu(3),kqd(3)
      integer kh1,kh2,kh3,khc,ksnu1,ksnu2,ksnu3,ksl1,ksl2,ksl3,ksl4,ksl5,ksl6,
     &        ksu1,ksu2,ksu3,ksu4,ksu5,ksu6,ksd1,ksd2,ksd3,ksd4,ksd5,ksd6,
     &        kn1,kn2,kn3,kn4,kcha1,kcha2,kgluin,kgold0,kgoldc 
      parameter (kh1=17,kh2=18,kh3=19,khc=20)
c old naming convention
*      parameter (ksnue=21,kse1=22,kse2=23,ksnumu=24,ksmu1=25,
*     &     ksmu2=26,ksnutau=27,kstau1=28,kstau2=29)
*      parameter (ksu1=30,ksu2=31,ksd1=32,ksd2=33,ksc1=34,ksc2=35,
*     &     kss1=36,kss2=37,kst1=38,kst2=39,ksb1=40,ksb2=41)
      parameter (ksnu1=21,ksl1=22,ksl2=23,ksnu2=24,ksl3=25,
     &     ksl4=26,ksnu3=27,ksl5=28,ksl6=29)
      parameter (ksu1=30,ksu2=31,ksd1=32,ksd2=33,ksu3=34,ksu4=35,
     &     ksd3=36,ksd4=37,ksu5=38,ksu6=39,ksd5=40,ksd6=41)
      parameter (kn1=42,kn2=43,kn3=44,kn4=45,kcha1=46,kcha2=47,
     &     kgluin=48)
      parameter (kgold0=49,kgoldc=50)
c old convention:      
c  ksqu=(~u1,~c1,~t1,~u2,~c2,~t2)    ksqd=(~d1,~s1,~b1,~d2,~s2,~b2)
c new convention: ksnu,ksl,ksqu,ksqd are *not* ordered
c ksquflav=((~u1,~c1,~t),(~u2,~c2,~t2)), ksquflav=((~d1,~s1,~b1),(~d2,~s2,~b2))
c      integer kse(2),ksmu(2),kstau(2),ksu(2),ksd(2),ksc(2),kss(2),
c     &     kst(2),ksb(2)   
      integer kn(4),kcha(2), ksnu(6),ksl(6),ksqu(6),ksqd(6),
     &        ksnu_flav(3,1),ksl_flav(3,2),ksqu_flav(3,2),ksqd_flav(3,2)

      common /pacodes_mssm/ knu,kl,kqu,kqd,ksnu,ksl,ksqu,ksqd,
     &                      kn,kcha,ksnu_flav,ksl_flav,ksqu_flav,ksqd_flav
     
     
* vertices
      complex*16 gl(numpartspecies,numpartspecies,numpartspecies),
     &     gr(numpartspecies,numpartspecies,numpartspecies)
      common /vrtxs_mssm/ gl,gr
c save common blocks
      save /pacodes_mssm/, /vrtxs_mssm/

c
c    standard model variables not included in dssm.h
c

* useful global variables + NLOoption
      real*8 delrho,delr
      character*7 NLOoption
      common /smruseful/ delrho,delr
      common /smcuseful/ NLOoption
* coupling constants
      real*8 alph3,yukawa(12)
      common /couplingconstants/ alph3,yukawa

c save common blocks
      save /smruseful/,/smcuseful/,
     &  /couplingconstants/

c
c    minimal supersymmetric model variables
c

* type of model (MSSM-7, mSUGRA, etc.) (same meaning as ModSel_Model in SLHA)
* 0=full MSSM, 1=mSUGRA
      integer modeltype
      common /mssmtype/ modeltype
* model parameters
      real*8 tanbe,mu,m2,m1,m3,ma,mass2u(3),mass2q(3),
     & mass2d(3),mass2l(3),mass2e(3),asoftu(3),asoftd(3),asofte(3)
      common /mssmpar/ tanbe,mu,m2,m1,m3,ma,mass2u,mass2q,mass2d,
     & mass2l,mass2e,asoftu,asoftd,asofte
* program switches
      integer higloop,neuloop,bsgqcd,higwid
      real*8 msquarks,msleptons
      common /mssmswitch/ msquarks,msleptons,
     &     higloop,neuloop,bsgqcd,higwid
* set sfermion masses by hand
      real*8 massup1(3),massup2(3),thetamixu(3),
     & massdn1(3),massdn2(3),thetamixd(3),
     & masssn(3),masssl1(3),masssl2(3),thetamixsl(3)
      common/sfermionmass/ massup1,massup2,thetamixu,
     & massdn1,massdn2,thetamixd,
     & masssn,masssl1,masssl2,thetamixsl
* useful global variables
      integer lsp,kln
      real*8 cosbe,sinbe,cos2be,sin2be,zg,lgzg
      common /mssmiuseful/ lsp,kln
      common /mssmruseful/ cosbe,sinbe,cos2be,sin2be,zg,lgzg
* coupling constants
      real*8 lam1,lam2,lam3,lam4,lam5,lam6,lam7
      common /mssmcouplingconstants/ lam1,lam2,lam3,lam4,lam5,lam6,lam7
* decay widths
      real*8 hdwidth(32,4)
      common /mssmwidths/ hdwidth
* mixings
      real*8 alpha,mix_stop,mix_sbot,mix_stau
      complex*16 neunmx(4,4),chaumx(2,2),chavmx(2,2),
     & slulmx(3,3),sldlmx(6,3),sldrmx(6,3),
     & squlmx(6,3),squrmx(6,3),sqdlmx(6,3),sqdrmx(6,3)
      common /mssmmixing/ neunmx,chaumx,chavmx,
     & slulmx,sldlmx,sldrmx,
     & squlmx,squrmx,sqdlmx,sqdrmx,alpha,mix_stop,mix_sbot,mix_stau
* msugra variables
      real*8 m0var,mhfvar,a0var,tgbetavar,sgnmuvar
      common/sugrainput/m0var,mhfvar,a0var,tgbetavar,sgnmuvar

c save common blocks
      save /mssmpar/,/mssmswitch/,/mssmiuseful/,/mssmruseful/,
     &  /mssmcouplingconstants/,
     &  /mssmmixing/,/sugrainput/


***                                                                 ***
************************ end of dsmssm.h ******************************
