*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                     dsanyieldmodelcom.h                          ***
***         this piece of code is needed as a separate file          ***
***       the rest of the code 'includes' dsanyieldmodelcom.h       ***
*** Note: this file only contains the model specific variables, like ***
*** scalar decays etc. All the common varialbes are contained in     ***
*** dsanyieldcom.h in src/                                           ***
c----------------------------------------------------------------------c
c  Author: Joakim Edsjo (edsjo@fysik.su.se)
c  Date: April, 2014

c....phi2par - variable passing to cascade routines
      real*8 phim0,phim1,phim2,phie0,phieth,phigap,phibep,
     &  phie1,phie2,phicthmin,phicthmax,phimp
      integer phich,phipkg,phifk,phitype,phihno
      common/anpar/phim0,phim1,phim2,phie0,phieth,
     &  phigap,phibep,phie1,phie2,phicthmin,phicthmax,phimp,
     &  phitype,phich,phipkg,phifk,phihno

c...separ - extra variables for se routines
      real*8 phithm
      integer phichi,phiwh,phifv ! phiwh=1 - sun, 2 - earth
      common /separ/phithm,phichi,phiwh,phifv


c...anch - channel conversion blocks
c...chcomp converts from new channel numbers to compressed channel numbers
c...Currently these are the old channel numbers, before Pythia runs
c...have been updated
      integer chcomp(29) 
      common /anch/chcomp

      integer chi2pdg(11),sechi2pdg(14)
      common/chistuff/chi2pdg,sechi2pdg
      
c...anbranch - annihilation branching rates and Scalar decay rates
c...Also, switches for IB (internal bremsstrahlung) are set here.
      integer numanch2b, numyieldch_line
      parameter (numanch2b=29,numyieldch_line=6)
      real*8 ans0br(numanch2b,3),anscbr(15),ans0m(3),
     &  anscm
      integer anch_2body(numanch2b,6),yieldchannels_line(numanch2b,2)
      common /anbranch/ans0br,anscbr,ans0m,anscm,
     &  anch_2body, yieldchannels_line

c...anopt - options
      integer anexhi
      real*8 ansbrmin
      common/anopt/ansbrmin,anexhi

c...Error/warning flags
      integer anerr,anistat
      common/anstat/anerr,anistat

c save common block
      save /anpar/,/anch/,/anbranch/,
     & /anopt/,/anstat/

***                                                                 ***
********************** end of dsanyieldcom.h **************************





