*****************************************************************************
*** function dsanyield_ch calculates the yield above threshold
*** or the differential flux for a given fluxtype 
***
***  input -- mwimp      : DM mass (in GeV)
***           egev       : energy of yield particle
***           pdg1, pdg2 : PDG codes of final state particles
***           yieldk     : fluxtype according to the following table:
***
***                       particle       integrated yield     differential yield
***                       --------       ----------------     ------------------
***                       positron                     51                    151
***                       cont. gamma                  52                    152
***                       nu_mu and nu_mu-bar          53                    153
***                       antiproton                   54                    154
***                       cont. gamma w/o pi0          55                    155
***                       nu_e and nu_e-bar            56                    156
***                       nu_tau and nu_tau-bar        57                    157
***                       pi0                          58                    158
***                       nu_mu and nu_mu-bar          71                    171 
***                                                    (same as 53/153)
***                       muons from nu at creation    72                    172
***                       muons from nu at detector    73                    173
***
***   output -- istat : zero if everything went OK
***
*** the units are (annihilation)**-1
*** for the differential yields, the units are the same plus gev**-1.
***
*** Note 1. The correct data files need to be loaded. This is handled by
*** a call to dsaninit. It is done automatically here upon first call.
***
*** Note 2. These routines do not contain internal bremsstrahlung (IB)
*** contributions (except those final state radiations (FSR) that are
*** included in the Pythia runs). The full IB contributions are added
*** in dsanyield.f.
***
*** author: joakim edsjo (edsjo@physto.se)
*** date: 98-01-26
*** modified: 08-01-15
*** modified: 08-11-27 pat scott 
*** modified: 09-10-20 pat scott 
*** modified: 2014-11-11 torsten bringmann (added pdg codes)
*****************************************************************************

      real*8 function dsanyield_ch(mwimp,egev,pdg1,pdg2,yieldk,istat)
      implicit none
      include 'dsanyieldmodelcom.h'
      include 'dsanyieldcom.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mwimp,egev,mp1,mp2,e1,e2,yield
      integer pdg1,pdg2,ch,istat,yieldk,chi
      integer yieldkk, yieldpdg, diff

c------------------------ functions ------------------------------------

      real*8 dsanyield_sim,dsanyields

c-----------------------------------------------------------------------
c... For internal use (=not visible to the DS core library), we simply convert 
c... pdg codes back to previous channel numbers (and only keep SUSY channels).

*** Ch No  Particles                 Old Ch No   Old chi   New chcomp
*** -----  ---------                 ---------   -------   ----------
***  1     S1 S1                      -          -         Not done yet
***  2     S1 S2                      -          -
***  3     S2 S2                      -          -
***  4     S3 S3                      -          -
***  5     S1 S3                      7          -
***  6     S2 S3                      11         -
***  7     S- S+                      -          -
***  8     S1 Z                       8          -
***  9     S2 Z                       9          -
*** 10     S3 Z	                      -          -
*** 11     W- S+ and W+ S-            10         -
*** 12     Z0 Z0                      6          6
*** 13     W+ W-                      5          5
*** 14     nu_e nu_e-bar              -          -
*** 15     e+ e-                      -          -
*** 16     nu_mu nu_mu-bar            -          -
*** 17     mu+ mu-                    13         7
*** 18     nu_tau nu_tau-bar	      -          -
*** 19     tau+ tau-	              4          4
*** 20     u u-bar                    -          -
*** 21     d d-bar                    -          -
*** 22     c c-bar                    1          1
*** 23     s s-bar                    -          -
*** 24     t t-bar                    3          3
*** 25     b b-bar                    2          2
*** 26     gluon gluon                12         8
*** 27     q q gluon (not implemented yet, put to zero)
*** 28     gamma gamma (1-loop)
*** 29     Z0 gamma (1-loop)


      if (pdg1.eq.25.and.pdg2.eq.25) then
        ch=3
      elseif (pdg1.eq.25.and.pdg2.eq.35.or.pdg1.eq.35.and.pdg2.eq.25) then
        ch=2
      elseif (pdg1.eq.35.and.pdg2.eq.35) then
        ch=1
      elseif (pdg1.eq.36.and.pdg2.eq.36) then
        ch=4
      elseif (pdg1.eq.25.and.pdg2.eq.36.or.pdg1.eq.36.and.pdg2.eq.25) then
        ch=6
      elseif (pdg1.eq.35.and.pdg2.eq.36) then
        ch=5
      elseif (pdg1.eq.37.and.pdg2.eq.-37.or.pdg1.eq.-37.and.pdg2.eq.37) then
        ch=7
      elseif (pdg1.eq.23.and.pdg2.eq.25.or.pdg1.eq.25.and.pdg2.eq.23) then
        ch=9
      elseif (pdg1.eq.23.and.pdg2.eq.35.or.pdg1.eq.35.and.pdg2.eq.23) then
        ch=8
      elseif (pdg1.eq.23.and.pdg2.eq.36.or.pdg1.eq.36.and.pdg2.eq.23) then
        ch=10
      elseif (pdg1.eq.24.and.pdg2.eq.-37.or.pdg1.eq.-37.and.pdg2.eq.24) then
        ch=11
      elseif (pdg1.eq.22.and.pdg2.eq.23.or.pdg1.eq.23.and.pdg2.eq.22) then
        ch=29 
      elseif (pdg1.eq.22.and.pdg2.eq.22) then
        ch=28  ! gamma gamma - not used
      elseif (pdg1.eq.10000.and.pdg2.eq.10000) then
         ch=-10000              ! dummy chanel -- not used!!!!
c         write(*,*) 'DS Error in dsanyield_ch:'
c         write(*,*) 'DS: called with pdg1=',pdg1,' pdg2=',pdg2
         dsanyield_ch=0.d0
         return
      else
         write(*,*)'DS WARNING -- ',
     &   'channel not implemented in dsanyield_ch: pdg1 = ',
     &   pdg1,' pdg2 = ',pdg2   
        dsanyield_ch=0d0
        return
      endif

c... The same needs to be done for yieldk: FIXME

      yieldkk=mod(yieldk,100)
      diff=yieldk/100
      if (diff.ne.0.and.diff.ne.1) then
         write(*,*) 'ERROR in dsanyield_ch:',
     &  ' unsupported argument yieldk =', yieldk
        dsanyield_ch=0.0d0
        return
      endif
      
      if (yieldkk.eq.51) then
         yieldpdg = -11 ! positron yields
      elseif (yieldkk.eq.52) then 
         yieldpdg = 22  ! cont. gammas     
      elseif (yieldkk.eq.53) then 
         yieldpdg = 14  ! muon neutrinos     
      elseif (yieldkk.eq.54) then 
         yieldpdg = -2212 ! antiproton yields     
      elseif (yieldkk.eq.61) then
         yieldpdg = -1000010020 ! anti-deuteron yields
         call dsanyield_dbset(61,-100.d0)
      elseif (yieldkk.eq.59) then 
         yieldpdg = -1000010020 ! anti-deuteron yields, old sph. coal mod
         call dsanyield_dbset(59,-100.d0)
      elseif (yieldkk.eq.71) then 
         yieldpdg = 14  ! neutrino yields (same as 53)    
      elseif (yieldkk.eq.72) then 
         yieldpdg = 130072 ! muon yields at creation     
      elseif (yieldkk.eq.73) then 
         yieldpdg = 130073 ! integrated muon yields in ice     
      else
         write(*,*) 'ERROR in dsanyield_ch:',
     &    ' unsupported argument yieldk =', yieldk
        dsanyield_ch=0.0d0
        return
      endif

c----------------------------------------------------------------------
      chi=chcomp(ch) ! convert from new to compact channel numbers

      istat=0
      anerr=0
      mp1=0.d0
      mp2=0.d0

c--------------------------------------- if first call, load tables

c      write(*,*) 'dsanyield_ch called with:'
c      write(*,*) '  mwimp = ',mwimp
c      write(*,*) '  egev = ',egev
c      write(*,*) '  ch = ',ch


c-----------------------------------------------------------------------

      if (chi.ge.1.and.chi.le.11) then ! "fundamental" channel

c        chok = .false.
c        if (mwimp.ge.msim(chi)) chok = .true.
c        if ((chi.eq.3.or.chi.eq.5.or.chi.eq.6).and.
c     &     mwimp.gt.(0.99*msim(chi))) chok = .true.
c        if (chok) then
c           dsanyield_ch=dsanyield_sim(mwimp,egev,chi2pdg(chi),'0',
c     &        yieldk,istat)
c        else
c          dsanyield_ch=0.0
c          istat=mod(istat,2)+2
c            write(6,5000) 'error in dsanyield_ch: channel ',ch,' is not',
c     &        ' energetically allowed.'
c            write(6,*) 'Mass of WIMP: ',mwimp
c            write(6,*) 'Mass of each final state particle: ',msim(chi)
c            write(6,*) 'Energy: ',egev
c            write(6,*) 'model: ',idtag
c
c        endif

        dsanyield_ch=dsanyield_sim(mwimp,egev,chi2pdg(chi),'0',
     &     yieldpdg,diff,istat)
        if (btest(istat,0)) then
           write(*,*) 'DS ERROR in dsanyield_ch:'
           write(*,*) '  Error raised by dsanyield_sim that simulation'
           write(*,*) '  tables were used outside of simulated regions.'
           write(*,*) '  Mass of WIMP: ',mwimp
           write(6,*) '  Final state channel: ',ch
           write(6,*) '  Energy: ',egev
           write(6,*) '  Model: ',idtag
        endif
 

      else                           ! "complex" channel

c...determine masses of the annihilation particles
        if (ch.eq.5) then
          mp1=ans0m(1)         ! S10 mass
          mp2=ans0m(3)         ! S30 mass
        elseif (ch.eq.8) then
          mp1=msim(9)        ! z0 mass
          mp2=ans0m(1)         ! S10 mass
        elseif (ch.eq.9) then
          mp1=msim(9)        ! z0 mass
          mp2=ans0m(2)         ! S20 mass
        elseif (ch.eq.11) then
          mp1=msim(8)        ! w+- mass
          mp2=anscm         ! S+- mass
        elseif (ch.eq.6) then
          mp1=ans0m(2)         ! S20 mass
          mp2=ans0m(3)         ! S30 mass
        elseif (ch.eq.29) then
          mp1=msim(9)        ! z0 mass
          mp2=0.0d0         ! gamma mass
        endif

c...if energetically allowed channel, go on...
        if (mwimp.ge.0.995d0*((mp1+mp2)/2.0d0)) then

c...calculate the energy of the annihilation particles
          e1=((2.0d0*mwimp)**2-mp2**2+mp1**2)/(4.0d0*mwimp)
          e2=2.0d0*mwimp-e1
          e1=max(e1,mp1+0.001d0)
          e2=max(e2,mp2+0.001d0)

c...check different annihilation channels
          yield=0.0d0

c---------- h10 h30 channel ----------
          if (ch.eq.5) then
            yield=yield+dsanyields(e1,egev,1,yieldk,istat)
            yield=yield+dsanyields(e2,egev,3,yieldk,istat)
          endif

c---------- z0 h10 channel ----------
          if (ch.eq.8) then
             yield=yield+0.5d0*dsanyield_sim(e1,egev,23,'0',
     &         yieldpdg,diff,istat)
            yield=yield+dsanyields(e2,egev,1,yieldk,istat)
          endif

c---------- z0 h20 channel ----------
          if (ch.eq.9) then
            yield=yield+0.5d0*dsanyield_sim(e1,egev,23,'0',yieldpdg,diff,istat)
            yield=yield+dsanyields(e2,egev,2,yieldk,istat)
          endif

c---------- w+h- w-h+ channel ----------
c...this calculation gives a mean of the two channels w+h- & w-h+
          if (ch.eq.11) then
            yield=yield+0.5d0*dsanyield_sim(e1,egev,24,'0',yieldpdg,diff,istat)
            yield=yield+dsanyields(e2,egev,4,yieldk,istat)
          endif

c---------- h20 h30 channel ----------
          if (ch.eq.6) then
            yield=yield+dsanyields(e1,egev,2,yieldk,istat)
            yield=yield+dsanyields(e2,egev,3,yieldk,istat)
          endif

c---------- z0 gamma channel ----------
          if (ch.eq.29) then
            yield=yield+0.5d0*dsanyield_sim(e1,egev,23,'0',yieldpdg,diff,istat)
          endif

          dsanyield_ch=yield

        else   ! not energetically allowed channel
          dsanyield_ch=0.0d0
          istat=mod(istat,2)+2
c          if (prtlevel.gt.0) then
            write(6,5000) 'error in dsanyield_ch: channel ',ch,' is not',
     +        ' energetically allowed.'
            write(6,*) 'Mass of WIMP: ',mwimp
            write(6,*) 'Mass of final state particles: ',mp1,mp2
            write(6,*) 'Energy: ',egev
            write(6,*) 'model: ',idtag
c          endif
        endif

      endif

      if (anerr.gt.0.and.prtlevel.gt.0) then
        write(*,*) 'warning in dsanyield_ch for model ',idtag,
     &    ', yield ',yieldk,' and channel ',ch
        write(*,*) '  the integration over higgs decay angles ran',
     &    'into numerical problems.'
        write(*,*) '  the results can only be trusted as a lower',
     &    'bound.'
        istat=ibset(istat,2)
      endif

 5000 format(' ',a,i2,a,a)

      end




