*****************************************************************************
*** function dsseyield_sim calculates the yield (including propagation
*** effects like oscillations and interactions) of particles from
*** WIMP annihilations in the Sun and the Earth.
*** This routine assumes that annihilation
*** takes place to two final state particles where the second is the
*** antiparticle of the first. For polarized states, the polarization
*** is assumed to be the same for both the particles.
***
*** Inputs:
***   - mwimp: mass of WIMP in GeV
***   - E:     kinetic energy of the neutrino, lepton or hadronic shower
***            where the yield is calculated (in GeV)
***   - theta: angle where the yield is calculated (in degrees)
***   - pdg:   PDG code of annihilation final state particle
***   Only the pdg code of the first particle is given, the second one
***   is assumed to have pdg code -pdg.
***   Only channels for which simulation data from Pythia simulations exist
***   are available here. More complex channels
***   like channels containing Higgs bosons etc that decay to standard model
***   particles are treated in each respective particle physics module
***   in src_model.
***   The currently implemented channels are
***   
***      pdg  Channel           Channel No chi  Array Ch No chii
***                             (only listed temporarily)
***     ----  -------           --------------  ----------------
***        1  d d-bar            1              1
***        2  u u-bar            2              2
***        3  s s-bar            3              3
***        4  c c-bar            4              4
***        5  b b-bar            5              5
***        6  t t-bar            6              6
***       12  nu_e nu_e-bar     12              11
***       13  mu- mu+           10               - not simulated
***       14  nu_mu nu_mu-bar   13              12
***       15  tau- tau+         11              10
***       16  nu_tau nu_tau-bar 14              13
***       24  W+ W-              8              8
***       23  Z0 Z0              9              9
***       21  gluon gluon        7              7
***
***   Note: If a channel that is not simulated is asked for, the yield
***   0 is returned and a warning is issues (istat bit 3 set). Note that
***   in the internal numbering here, mu- mu+ is not included.
***   However when WimpSim is run, it is included (as number 10) to have
***   a consistent numbering with HaloAnn.      
***
***   hel: helicity state of final state. Currently, only states where the
***        two final state particles are in the same helicity state (e.g.
***        both left, both right, both transverse, both longitudinal)
***        is included.
***        Note: right now, only unpolarized yields are given, will
***        change eventually
***
***   Possible values for hel:
***     'L': left-polarized fermion for fermions, longitudinal for vectors
***     'R': right-polarized fermion for fermions
***     'T': transverse for vectors
***     '0': unpolarized yields (unphysical, but included for comaprison)
***
***   - wh:    Flag determinging source body:
***            'su' for Sun and 'ea' for Earth
***   - kind:  The yields are of different kinds:
***      = 1:  integrated up to a given theta and above a given energy
***        2:  differetial in energy and angle
***        3:  differential in energy, but integrated in angle up to the
***            given angular cut
***   - type:  Type of yield:
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
*** If this routine is called outside of the kinematical regions
*** where tables exist, the following is done:
***   - for lower energies, than the lowest simulted ones,
***     extrapolations are used
***   - for masses below the lowest simulated, extrapolations
***     are used
***   - for energies above the highest simulated, the results
***     for the highest energy simulated are used
*** Output:
***   - istat (=0 if no warnings/errors are reported)
***     bit 0 is set if extrapolations below the lowest simulated mass
***       or above the highest simulated mass is needed
***     bit 3 is set if the requested channel is not simulated,
***       yield zero is returned
***     bit 4 is set if the requested polarization state is not available,
***       the yield of the simulated polarization yield (typically
***       unpolarized) is returned instead
***   - dsseyield_sim in units of
*** units: 1.0e-30 m**-2 (annihilation)**-1  integrated (kind 1)
*** units: 1.0e-30 m**-2 gev**-1 (degree)**-1 (annihilation)**-1 differential
***        (kind 2)
*** units: 1.0e-30 m**-2 gev**-1 (annihilation)**-1 mixed (kind 3)
*** types 7-12, 15-26 have an additional unit of m**-1.
*** author: joakim edsjo, edsjo@physics.berkeley.edu
*** date: 1995
*** modified: dec 03, 1997, April 3, 2008
*** Modified 2011-05-08 to included mixed yields (kind=3). Corrected a bug
*** for integrated yields below 0.2 degrees.
*** Modified: December 2014 to use PDG codes (edsjo)      
*****************************************************************************

      real*8 function dsseyield_sim(mwimp,e,theta,pdg,hel,wh,kind,
     &  type,istat)
      use omp_lib
      implicit none
      include 'dsseyieldcom.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mwimp,e,theta,phi1,phi2,mp1,mp2,zpl,thpl,mn,z,th,
     &  tmp
      integer chi,istat,zi,thi,m1i,m2i,whi,kind,type,yli,chii,i
      integer pdg,apdg
      character*2 wh
      character*1 hel
      logical wb,chok
      external dsseyieldcom ! set up common block variables

      logical first
      data first/.true./
      save first
c-----------------------------------------------------------------------

c--------------------------------------- if first call, load yield tables


      if (first) then
        do i=1,26
           yload(1,i)=0
           yload(2,i)=0
        enddo
        selast(1)=0 ! last index for integrated yields stored in memory
        selast(2)=0 ! last index for differential yields stored in memory
        first=.false.
      endif

      if (yload(kind2ki(kind),type).eq.0) then
        call dsseinit(kind,type)
      endif

c-----------------------------------------------------------------------

      
      yli=yload(kind2ki(kind),type)

      chi=0
      apdg=abs(pdg)

      if (apdg.eq.1) then      ! d d-bar
         if (hel.eq.'0') then
            chi=1 
         else ! eventually we will have sims here
            chi=1 
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.2) then      ! u u-bar
         if (hel.eq.'0') then
            chi=2 
         else ! eventually we will have sims here
            chi=2 
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.3) then      ! s s-bar
         if (hel.eq.'0') then
            chi=3
         else ! eventually we will have sims here
            chi=3 
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.4) then      ! c c-bar
         if (hel.eq.'0') then
            chi=4
         else ! eventually we will have sims here
            chi=4
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.5) then ! b b-bar
         if (hel.eq.'0') then
            chi=5
         else ! eventually we will have sims here
            chi=5
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.6) then ! t t-bar
         if (hel.eq.'0') then
            chi=6
         else ! eventually we will have sims here
            chi=6
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.11) then ! e- e+ ! not sim, but yield zer0
            chi=0 
      elseif (apdg.eq.12) then ! nu_e nu_e-bar
         if (hel.eq.'0') then
            chi=12
         else ! eventually we will have sims here
            chi=12
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.13) then ! mu- mu+ ! not sim, but yield zer0
            chi=0 
      elseif (apdg.eq.14) then ! nu_mu nu_mu-bar
         if (hel.eq.'0') then
            chi=13
         else ! eventually we will have sims here
            chi=13
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.15) then ! tau- tau+
         if (hel.eq.'0') then
            chi=11
         else ! eventually we will have sims here
            chi=11
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.16) then ! nu_tau nu_tau-bar
         if (hel.eq.'0') then
            chi=14
         else ! eventually we will have sims here
            chi=14
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.24) then ! W- W+
         if (hel.eq.'0') then
            chi=8
         else ! eventually we will have sims here
            chi=8
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.23) then ! Z0 Z0
         if (hel.eq.'0') then
            chi=9
         else ! eventually we will have sims here
            chi=9
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.21) then ! g g
         if (hel.eq.'0'.or.hel.eq.'T') then ! TB corr 20190607
            chi=7
         else ! 
            chi=7
            istat=ibset(istat,4)
         endif
      else
         istat=ibset(istat,3) 
         dsseyield_sim=0.d0
         return
      endif

      if (chi.eq.0) then        ! yield zero for these channels
         dsseyield_sim=0.d0
         return
      endif
      
c...Convert chi above (1-14) to array indices chi (1-13, muon channel missing)
      chii=chi2chii(chi)

c...Safecheck
      if (yli.eq.0) then
        if (omp_get_thread_num() .eq. 0) then
          call dsseinit(kind,type)
        endif
!$omp barrier
         yli=yload(kind,type)
      endif

      dsseyield_sim=0.0d0

      wb=.true.
      if (wh.eq.'su'.or.wh.eq.'SU') then
        whi=1
      else
        whi=2
      endif

      if (e.ge.mwimp) then
        dsseyield_sim=0.0d0
        return
      endif

      mn = mwimp

c...Check kinematics
      chok=.false.
      if (mn.ge.msim(chi)) chok=.true.
      if ((chi.eq.6).and.
     &     mn.ge.(0.95*msim(chi))) chok = .true. ! t t-bar
      if ((chi.eq.8.or.chi.eq.9).and.
     &     mn.ge.(0.97*msim(chi))) chok = .true. ! W+W-, Z0 Z0
      if (.not.chok) istat=ibset(istat,0)

c...take care of the case where mwimp is between the mass of the annihilation
c...product and the lower bound of the simulations (there might be a small
c...gap of less than a gev) for tt-bar,ww and zz.
      if (mn.lt.lb(chi)) then
         if (mn.lt.0.97d0*lb(chi)) istat=(istat/2)*2+1
         mn=lb(chi)
      endif

c      if (chi.eq.6.or.chi.eq.8.or.chi.eq.9) then
c        if (mn.lt.lb(chi).and.mn.gt.(0.97*lb(chi))) then
c          mn=lb(chi)
c        endif
c      endif

c...check if mwimp within correct bounds
c...This part is obsolete, we extrapolate instead
c      if (mn.lt.lb(chi).and.e.lt.lb(chi)) then
c        if (prtlevel.gt.0) then
c          write(6,*)
c          write(6,5000) 
c     &      'WARNING in dsseyield_sim: a WIMP mass of ',mn,
c     +      ' gev wants to be used,'
c          write(6,5010) 'while the lower bound for channel ',chi,
c     +      ' is ',lb(chi),' GeV.'
c          write(6,*) 'the yield is put to 0.0 for these too low masses.'
c          write(6,*) 'the results can thus only be trusted as a',
c     +      ' lower bound.'
c          write(6,*) 'model: ',idtag
c        endif
c        wb=.false.
c        dsseyield_sim=0.0d0
c        istat=(istat/2)*2+1
c      endif

      if (mn.gt.ub(chi)) then
        if (prtlevel.gt.0) then
          write(6,*)
          write(6,5000) 
     &      'WARNING in dsseyield_sim: a WIMP mass of ',mn,
     +      ' gev wants to be used,'
          write(6,5010) 'while the upper bound for channel ',chi,
     +      ' is ',ub(chi),' GeV.'
          write(6,5020) 'a WIMP mass of ',ub(chi),' GeV is used',
     +      ' instead for these too high masses.'
          write(6,*) 'the results can thus only be trusted as a',
     +      ' lower bound.'
          write(6,*) 'model: ',idtag
        endif
        istat=(istat/2)*2+1
      endif


c---------------------------------------------------- integrated yields
      if (wb.and.kind.eq.1) then
c...determine which entries in phiint to use and how

        z=e/mn
        call dsseifind(z,zindex(0,1),zpl,zi,0,zn2-1)

        if (zi.eq.-5.or.zi.ge.zn2) then
          dsseyield_sim=0.0d0
          return
        endif

        th=theta
        call dsseifind(th,thindex(-1,1),thpl,thi,-1,thn-1)

        if (theta.gt.30.0d0) then
          thi=thn-1
          thpl=1.0d0
        endif
        if (thi.eq.-5) then
          dsseyield_sim=0.0d0
          return
        endif

        call dsseifind(mn,mi(1),tmp,m1i,1,senm-1)
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(senm)) then
          m1i=senm
          m2i=senm
          mp1=mi(senm)
          mp2=mp1

          dsseyield_sim = (1.0-thpl)*((1.0-zpl)*
     &      dble(phiint(thi,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phiint(thi,zi+1,m1i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phiint(thi+1,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phiint(thi+1,zi+1,m1i,chii,whi,yli)))
        else
          phi1 = (1.0-thpl)*((1.0-zpl)*
     &      dble(phiint(thi,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phiint(thi,zi+1,m1i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phiint(thi+1,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phiint(thi+1,zi+1,m1i,chii,whi,yli)))
          phi2 = (1.0-thpl)*((1.0-zpl)*
     &      dble(phiint(thi,zi,m2i,chii,whi,yli))+
     &      zpl*dble(phiint(thi,zi+1,m2i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phiint(thi+1,zi,m2i,chii,whi,yli))+
     &      zpl*dble(phiint(thi+1,zi+1,m2i,chii,whi,yli)))
c          dsseyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
c     &      ((mp2-mp1)*(mp2+mp1)) ! old interpolation
c          dsseyield_sim = phi1*(mp2-mn)/(mp2-mp1)
c     &      +phi2*(mn-mp1)/(mp2-mp1) ! linear interpolation is better
          if (mn.lt.250.d0.or.phi1.eq.0.d0.or.phi2.eq.0.d0) then
            dsseyield_sim = phi1*(mp2-mn)/(mp2-mp1)
     &            +phi2*(mn-mp1)/(mp2-mp1) ! linear interpolation here
          else
            dsseyield_sim = 10**(log10(phi1)*(mp2-mn)/(mp2-mp1)
     &         +log10(phi2)*(mn-mp1)/(mp2-mp1)) ! Log-lin interpolation
          endif
c          write(*,*) 'mx=',mwimp,'  e=',e,'  zi=',zi,'  zpl=',zpl
c          write(*,*) 'phi1=',phi1
c          write(*,*) 'phi2=',phi2
c          write(*,*) 'phi1 parts:',
c     &      phiint(thi,zi,m1i,chii,whi,yli),
c     &      phiint(thi,zi+1,m1i,chii,whi,yli),
c     &      phiint(thi+1,zi,m1i,chii,whi,yli),
c     &      phiint(thi+1,zi+1,m1i,chii,whi,yli)
c          write(*,*) 'phi2 parts:',
c     &      phiint(thi,zi,m2i,chii,whi,yli),
c     &      phiint(thi,zi+1,m2i,chii,whi,yli),
c     &      phiint(thi+1,zi,m2i,chii,whi,yli),
c     &      phiint(thi+1,zi+1,m2i,chii,whi,yli)



        endif
        if (mn.ne.mwimp) then
          dsseyield_sim=mwimp**2/mn**2*dsseyield_sim
        endif
      endif

c-------------------------------------------------- differential yields
      if (wb.and.kind.eq.2) then
c...determine which entries in phidiff to use and how
        z=e/mn
        call dsseifind(z,zindex(-1,2),zpl,zi,-1,zn2-1)

        if (zi.eq.-5.or.zi.ge.zn2) then
          dsseyield_sim=0.0d0
          return
        endif

        if ((e.gt.mn).or.(e.le.0.0).or.zi.eq.-5) then
          dsseyield_sim=0.0d0
          return
        endif

        th=theta
        call dsseifind(th,thindex(-1,2),thpl,thi,-1,thn-1)

c...The following is to make sure that we get the last bin correct
c...It contains all the events with theta>30 degrees. We want it to
c...integrate correctly to the correct number of events so we place
c...it between 30.0 and 30.5 to make sure that happens.
c...We are also careful with how we interpolate the last bin before
c...that to avoid problems
        if (theta.gt.30.0d0.and.theta.lt.30.5d0) then
           thi=thn
           thpl=0.0d0
        endif

        if (theta.gt.29.75.and.theta.lt.30.0d0) then
           thpl=0.d0 ! do not include last theta>30 bin
        endif
           
        if (theta.gt.30.5d0.or.thi.eq.-5) then
          dsseyield_sim=0.0d0
          return
        endif

        call dsseifind(mn,mi(1),tmp,m1i,1,senm-1)
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(senm)) then
          m1i=senm
          m2i=senm
          mp1=mi(senm)
          mp2=mp1
          dsseyield_sim = (1.0-thpl)*((1.0-zpl)*
     &      dble(phidiff(thi,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phidiff(thi,zi+1,m1i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phidiff(thi+1,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phidiff(thi+1,zi+1,m1i,chii,whi,yli)))
c... convert from dyield / dz dtheta to dyield / de dtheta
          dsseyield_sim=dsseyield_sim/mwimp
        else
          phi1 = (1.0-thpl)*((1.0-zpl)*
     &      dble(phidiff(thi,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phidiff(thi,zi+1,m1i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phidiff(thi+1,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phidiff(thi+1,zi+1,m1i,chii,whi,yli)))
          phi2 = (1.0-thpl)*((1.0-zpl)*
     &      dble(phidiff(thi,zi,m2i,chii,whi,yli))+
     &      zpl*dble(phidiff(thi,zi+1,m2i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phidiff(thi+1,zi,m2i,chii,whi,yli))+
     &      zpl*dble(phidiff(thi+1,zi+1,m2i,chii,whi,yli)))
c... convert from dyield / dz dtheta to dyield / de dtheta
          phi1=phi1/mp1
          phi2=phi2/mp2
c          dsseyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
c     &      ((mp2-mp1)*(mp2+mp1)) ! lin-lin
          if (mn.lt.250.d0.or.phi1.eq.0.d0.or.phi2.eq.0.d0) then
             dsseyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
     &            ((mp2-mp1)*(mp2+mp1)) ! lin-lin
          else
            dsseyield_sim = 10**(log10(phi1)*(mp2-mn)/(mp2-mp1)
     &            +log10(phi2)*(mn-mp1)/(mp2-mp1)) ! Log-lin interpolation
	  endif
	  if (mn.gt.500.d0.and.mn.lt.750.d0) then ! JE TMP
	     write(49,*) mn,phi1,phi2,dsseyield_sim,dsseyield_sim/mn,
     &       (phi1/mp1)+(phi2/mp2-phi1/mp1)*(mn-mp1)*(mn+mp1)/
     &       ((mp2-mp1)*(mp2+mp1)) ! JE TMP	     
	  endif
        endif
        if (mn.ne.mwimp) then
          dsseyield_sim=mwimp**2/mn**2*dsseyield_sim
        endif
c... convert from dyield / dz dtheta to dyield / de dtheta
c        dsseyield_sim=dsseyield_sim/mwimp ! done before interpolation now
      endif

c-------------------------------------------------- mixed yields
      if (wb.and.kind.eq.3) then
c...determine which entries in phidiff to use and how
        z=e/mn
        call dsseifind(z,zindex(-1,2),zpl,zi,-1,zn2-1)

        if (zi.eq.-5.or.zi.ge.zn2) then
          dsseyield_sim=0.0d0
          return
        endif

        if ((e.gt.mn).or.(e.le.0.0).or.zi.eq.-5) then
          dsseyield_sim=0.0d0
          return
        endif

        th=theta
        call dsseifind(th,thindex(-1,1),thpl,thi,-1,thn-1)

        if (theta.gt.30.0d0) then
          thi=thn-1
          thpl=1.0d0
        endif
        if (thi.eq.-5) then
          dsseyield_sim=0.0d0
          return
        endif

        call dsseifind(mn,mi(1),tmp,m1i,1,senm-1)
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(senm)) then
          m1i=senm
          m2i=senm
          mp1=mi(senm)
          mp2=mp1
          dsseyield_sim = (1.0-thpl)*((1.0-zpl)*
     &      dble(phimixed(thi,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phimixed(thi,zi+1,m1i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phimixed(thi+1,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phimixed(thi+1,zi+1,m1i,chii,whi,yli)))
c... convert from dyield / dz to dyield / de
          dsseyield_sim=dsseyield_sim/mwimp
        else
          phi1 = (1.0-thpl)*((1.0-zpl)*
     &      dble(phimixed(thi,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phimixed(thi,zi+1,m1i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phimixed(thi+1,zi,m1i,chii,whi,yli))+
     &      zpl*dble(phimixed(thi+1,zi+1,m1i,chii,whi,yli)))
          phi2 = (1.0-thpl)*((1.0-zpl)*
     &      dble(phimixed(thi,zi,m2i,chii,whi,yli))+
     &      zpl*dble(phimixed(thi,zi+1,m2i,chii,whi,yli)))+
     &      thpl*((1.0-zpl)*dble(phimixed(thi+1,zi,m2i,chii,whi,yli))+
     &      zpl*dble(phimixed(thi+1,zi+1,m2i,chii,whi,yli)))
c... convert from dyield / dz to dyield / de
          phi1=phi1/mp1
          phi2=phi2/mp2
c          dsseyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
c     &      ((mp2-mp1)*(mp2+mp1))
          if (mn.lt.250.d0.or.phi1.eq.0.d0.or.phi2.eq.0.d0) then ! linear
             dsseyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
     &            ((mp2-mp1)*(mp2+mp1)) ! lin-lin
          else
             dsseyield_sim = 10**(log10(phi1)*(mp2-mn)/(mp2-mp1)
     &            +log10(phi2)*(mn-mp1)/(mp2-mp1)) ! Log-lin
          endif
        endif
        if (mn.ne.mwimp) then
          dsseyield_sim=mwimp**2/mn**2*dsseyield_sim
        endif
c... convert from dyield / dz to dyield / de
c        dsseyield_sim=dsseyield_sim/mwimp ! done above
      endif


 5000 format(' ',a,f8.2,a)
 5010 format(' ',a,i2,a,f8.2,a)
 5020 format(' ',a,f8.2,a,a)

      return

      end
