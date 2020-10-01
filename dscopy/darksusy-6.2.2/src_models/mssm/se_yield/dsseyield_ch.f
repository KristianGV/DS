      real*8 function dsseyield_ch(mneu,e,theta,pdg1,pdg2,
     &  wh,kind,type,istat)
*****************************************************************************
*** Function dsseyield_ch calculates the yield of neutrinos, leptons or
*** hadronic showers from WIMP annihilations in the Sun/Earth for a given
*** annihilation channel. Channels 1-xx are supported. Channel 10 (mu- mu+)
*** always gives zero (as muons are stopped in the Sun/Earth), but is
*** included to keep the channel numbering the same as for halo annihilation
*** routines.
*** Inputs:
***   mneu  - WIMP mass (GeV)
***   emu   - energy of neutrino, lepton or hadronic shower (GeV)
***   theta - angle from the Sun / centre of the Earth (degrees)
***   pdg1, pdg2 -- PDG codes of final state particles
***
***   wh    - 'su' for Sun, 'ea' for Earth
***   kind  - 1 = integrated
***           2 = differential
***           3 = mixed, integrated in theta, differential in energy
***   type  - Type of yield
***
***   type   Yield at detector
***   ----   -----------------
***   1      nu_e
***   2      nu_e-bar
***   3      nu_mu
***   4      nu_mu-bar
***   5      nu_tau
***   6      nu_tau-bar
***   7      e- at neutrino-nucleon vertex
***   8      e+ at neutrino-nucleon vertex
***   9      mu- at neutrino-nucleon vertex
***   10     mu+ at neutrino-nucleon vertex
***   11     tau- at neutrino-nucleon vertex
***   12     tau+ at neutrino-nucleon vertex
***   13     mu- at an imaginary plane in detector (i.e. after propagation)
***   14     mu+ at an imaginary plane in detector (i.e. after propagation)
***   15     hadronic shower from nu_e charged current (CC) interactions
***   16     hadronic shower from nu_e-bar charged current (CC) interactions
***   17     hadronic shower from nu_mu charged current (CC) interactions
***   18     hadronic shower from nu_mu-bar charged current (CC) interactions
***   19     hadronic shower from nu_tau charged current (CC) interactions
***   20     hadronic shower from nu_tau-bar charged current (CC) interactions
***   21     hadronic shower from nu_e neutral current (NC) interactions
***   22     hadronic shower from nu_e-bar neutral current (NC) interactions
***   23     hadronic shower from nu_mu neutral current (NC) interactions
***   24     hadronic shower from nu_mu-bar neutral current (NC) interactions
***   25     hadronic shower from nu_tau neutral current (NC) interactions
***   26     hadronic shower from nu_tau-bar neutral current (NC) interactions
***
*** Outputs:
***   yield in units of 
***      1e-30 m**-2 (annihilation)**-1 for types 1-6 and 13-14
***      1e-30 m**-3 (annihilation)**-1 for types 7-12, 15-26.
***      For the differential yields, the units are the same plus
***      GeV**-1 degree**-1 (kind 2) and GeV**-1 (kind 3)
***   istat - status flag (non-zero if something went wrong)
*** NOTE: Compared to previous versions of this routine, particles and
*** antiparticles are no longer summed, you thus need to call it for both
*** types and add them up.
*** Author: Joakim Edsjo (edsjo@fysik.su.se)
*** Date: 1995
*** Modified: April 3, 2008
*** Modified: May 6, 2011 Pat Scott (patscott@physics.mcgill.ca)
*** Modified: December, 2014 (edsjo)	
*****************************************************************************
      use omp_lib
      implicit none
	
      include 'dsanyieldmodelcom.h'
      include 'dsseyieldcom.h'	
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mneu,e,theta,flx
      real*8 mp1,mp2,e1,e2
      integer ch,chi,istat,i,kind,type,pdg1,pdg2,seerror
      character*2 wh

      logical first,firstwar
      data first /.true./
      data firstwar /.true./
      save first,firstwar


c------------------------ functions ------------------------------------

      real*8 dsseyield_sim,dsseyields

c-----------------------------------------------------------------------
c... For internal use (=not visible to the DS core library), we simply convert 
c... pdg codes back to previous channel numbers (and only keep SUSY channels).

*** ch     Particles                 WS chi  Internal chii
*** -----  ---------                 ------  ------------
***  1     S1 S1                     -       -
***  2     S1 S2                     -       -
***  3     S2 S2                     -       -
***  4     S3 S3                     -       -
***  5     S1 S3                     -       -
***  6     S2 S3                     -       -
***  7     S- S+                     -       -
***  8     Z S1                      -       -
***  9     Z S2                      -       -
*** 10     Z S3	                     -       -
*** 11     W- S+ and W+ S-           -       - 
*** 12     Z0 Z0 	             9       9 
*** 13     W+ W-                     8       8
*** 14     nu_e nu_e-bar             12      11
*** 15     e+ e-                     -       -
*** 16     nu_mu nu_mu-bar           13      12
*** 17     mu+ mu-                   10      -
*** 18     nu_tau nu_tau-bar	     14      13 
*** 19     tau+ tau-	             11      10
*** 20     u u-bar                   2       2
*** 21     d d-bar                   1       1
*** 22     c c-bar                   4       4
*** 23     s s-bar                   3       3
*** 24     t t-bar                   6       6
*** 25     b b-bar                   5       5
*** 26     gluon gluon               7       7
*** 27     q q gluon (not implemented yet, put to zero)
*** 28     gamma gamma (1-loop)
*** 29     Z0 gamma (1-loop)
***
*** ch is the channel number (including complex channel numbers) that
*** are used as input to this routine. These are the same as in dssigmav
*** in case of SUSY. chi (1-14) are the channel numbers used by WimpSim
*** chii are the internal channel numbers used as array indices (differs
*** from chi as mu- mu+ channels are not stored (they give zero yield)).

      istat=0
      seerror=0
      dsseyield_ch=0.d0
      
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
      elseif (pdg1.eq.35.and.pdg2.eq.36.or.pdg1.eq.36.and.pdg2.eq.35) then
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
c         write(*,*) 'DS Error in dsseyield_ch:'
c         write(*,*) 'DS: called with pdg1=',pdg1,' pdg2=',pdg2
         return
      else
         write(*,*)'DS WARNING -- ',
     &   'channel not implemented in dsanyield_ch: pdg1 = ',
     &   pdg1,' pdg2 = ',pdg2   
        return
      endif
c------------------------------------------------------------------
c      write(*,*) 'dsseyield_ch called with ',kind,type

      mp1=0.d0
      mp2=0.d0

      chi=ch2chi(ch) ! from full channel number to WimpSim channel number

c--------------------------------------- if first call, load yield tables


      if (first) then
        if (omp_get_thread_num() .eq. 0) then
           do i=1,26
              yload(1,i)=0
              yload(2,i)=0
           enddo
           selast(1)=0 ! last index for integrated yields stored in memory
           selast(2)=0 ! last index for differential yields stored in memory
           first=.false.
         endif
!$omp barrier
      endif

      if (yload(kind2ki(kind),type).eq.0) then
         if (omp_get_thread_num() .eq. 0) then
            call dsseinit(kind,type)
         endif
!$omp barrier
      endif

c-----------------------------------------------------------------------
      if (chi.gt.0) then              ! "fundamental" channel
        dsseyield_ch=dsseyield_sim(mneu,e,theta,chi,wh,
     &    kind,type,istat)
        if (btest(istat,0)) then
           write(*,*) 'DS ERROR in dsseyield_ch:'
           write(*,*) '  Error raised by dsseyield_sim that simulation'
           write(*,*) '  tables were used outside of simulated regions.'
           write(*,*) '  Mass of WIMP: ',mneu
           write(6,*) '  Final state channel: ',ch
           write(6,*) '  Energy: ',e
           write(6,*) '  Model: ',idtag
        endif

      else                           ! "complex" channel

c...determine masses of the annihilation particles
        if (ch.eq.1) then
          mp1=ans0m(1)         ! S10 mass
          mp2=ans0m(1)         ! S10 mass
        elseif (ch.eq.2) then
          mp1=ans0m(1)         ! S10 mass
          mp2=ans0m(2)         ! S20 mass
        elseif (ch.eq.3) then
          mp1=ans0m(2)         ! S20 mass
          mp2=ans0m(2)         ! S20 mass
        elseif (ch.eq.4) then
          mp1=ans0m(3)         ! S30 mass
          mp2=ans0m(3)         ! S30 mass
        elseif (ch.eq.5) then
          mp1=ans0m(1)         ! S10 mass
          mp2=ans0m(3)         ! S30 mass
        elseif (ch.eq.6) then
          mp1=ans0m(2)         ! S20 mass
          mp2=ans0m(3)         ! S30 mass
        elseif (ch.eq.7) then
          mp1=anscm            ! S+ mass
          mp2=anscm            ! S- mass
        elseif (ch.eq.8) then
          mp1=msim(9)           ! z0 mass
          mp2=ans0m(1)         ! h20 mass
        elseif (ch.eq.9) then
          mp1=msim(9)           ! z0 mass
          mp2=ans0m(2)         ! S20 mass
        elseif (ch.eq.10) then
          mp1=msim(9)           ! z0 mass
          mp2=ans0m(3)         ! S30 mass
        elseif (ch.eq.11) then
          mp1=msim(8)           ! w+- mass
          mp2=anscm            ! S+- mass
        elseif (ch.eq.29) then
          mp1=msim(9)           ! z0 mass
          mp2=0.0d0            ! gamma mass
        else                   ! not a supported channel
          dsseyield_ch=0.0d0
          return
        endif

c...if energetically allowed channel, go on...
        if (mneu.ge.0.995d0*((mp1+mp2)/2.0d0)) then

c...calculate the energy of the annihilation particles
          e1=((2.0d0*mneu)**2-mp2**2+mp1**2)/(4.0d0*mneu)
          e2=2.0d0*mneu-e1
          e1=max(e1,mp1+0.001d0)
          e2=max(e2,mp2+0.001d0)

c...check different annihilation channels
          flx=0.0d0

c---------- S10 S10 channel ----------
          if (ch.eq.1) then
             flx=flx+2.d0*dsseyields(e1,e,theta,1,wh,kind,type,istat,
     &          seerror)

c---------- S10 S20 channel ----------
          elseif (ch.eq.2) then
            flx=flx+dsseyields(e1,e,theta,1,wh,kind,type,istat,seerror)
            flx=flx+dsseyields(e2,e,theta,2,wh,kind,type,istat,seerror)

c---------- S20 S20 channel ----------
          elseif (ch.eq.3) then
             flx=flx+2.d0*dsseyields(e1,e,theta,2,wh,kind,type,istat,
     &          seerror)

c---------- S30 S30 channel ----------
          elseif (ch.eq.4) then
             flx=flx+2.d0*dsseyields(e1,e,theta,3,wh,kind,type,istat,
     &          seerror)

c---------- S10 S30 channel ----------
          elseif (ch.eq.5) then
            flx=flx+dsseyields(e1,e,theta,1,wh,kind,type,istat,seerror)
            flx=flx+dsseyields(e2,e,theta,3,wh,kind,type,istat,seerror)

c---------- S20 S30 channel ----------
          elseif (ch.eq.6) then
            flx=flx+dsseyields(e1,e,theta,2,wh,kind,type,istat,seerror)
            flx=flx+dsseyields(e2,e,theta,3,wh,kind,type,istat,seerror)

c---------- S+ S- channel ----------
          elseif (ch.eq.7) then
             flx=flx+2.d0*dsseyields(e1,e,theta,4,wh,kind,type,istat,
     &            seerror)

c---------- Z0 S10 channel ----------
          elseif (ch.eq.8) then
            flx=flx+0.5d0*dsseyield_sim(e1,e,theta,23,'0',wh,
     &        kind,type,istat)
            flx=flx+dsseyields(e2,e,theta,1,wh,kind,type,istat,seerror)

c---------- Z0 S20 channel ----------
          elseif (ch.eq.9) then
            flx=flx+0.5d0*dsseyield_sim(e1,e,theta,23,'0',wh,
     &        kind,type,istat)
            flx=flx+dsseyields(e2,e,theta,2,wh,kind,type,istat,seerror)

c---------- Z0 S30 channel ----------
          elseif (ch.eq.10) then
            flx=flx+0.5d0*dsseyield_sim(e1,e,theta,23,'0',wh,
     &        kind,type,istat)
            flx=flx+dsseyields(e2,e,theta,3,wh,kind,type,istat,seerror)

c---------- w+h- w-h+ channel ----------
c...this calculation gives a mean of the two channels w+h- & w-h+
          elseif (ch.eq.11) then
            flx=flx+0.5d0*dsseyield_sim(e1,e,theta,24,'0',wh,
     &        kind,type,istat)
            flx=flx+dsseyields(e2,e,theta,4,wh,kind,type,istat,seerror)

c---------- z0 gamma channel ----------
          elseif (ch.eq.29) then
            flx=flx+0.5d0*dsseyield_sim(e1,e,theta,23,'0',wh,
     &        kind,type,istat)
          endif

          dsseyield_ch=flx

        else   ! not energetically allowed channel
          dsseyield_ch=0.0d0
          istat=mod(istat,2)+2
          if (prtlevel.gt.0) then
c$omp critical (stdout)
            write(6,5000) 'WARNING in dsseyield_ch: channel ',ch,
     &        ' is not energetically allowed. '
            write(6,5000) '(NB: When running DS as a standalone this should not happen!). '     
            write(6,*) 'model: ',idtag
c$omp end critical (stdout)

          endif
        endif

      endif

      if (seerr.gt.0) then
        write(*,*) 'DS WARNING in dsseyield_ch for model ',idtag,
     &    ', yield type ',type,' and channel ',ch
        write(*,*) '  the integration over scalar decay angles ran',
     &    ' into numerical problems.'
        write(*,*) '  The results can only be trusted as a lower',
     &    'bound.'
        seistat=ibset(seistat,2)
        istat=ibset(istat,2)
      endif

 5000 format(' ',a,i2,a,a,a)

      end
