*****************************************************************************
*** function dsanyield_sim calculates the yield above threshold
*** (or differential at that energy) for the requested annihilation channel
*** and the fluxtype given by yieldpdg. This routine assumes that annihilation
*** takes place to two final state particles where the second is the
*** antiparticle of the first. For polarized states, the polarization
*** is assumed to be the same for both the particles.
***
*** Inputs:
***   - mwimp = WIMP mass in GeV
***   - e = kinetic energy where the yield is calculated (in GeV)
***         Note: e is kinetic energy per particle (i.e. not per nucleon for
***         anti-deuterons)
*** 
***   - pdg = pdg codes of annihilation final state particle.
***   Only the pdg code of the first particle is given, the second one
***   is assumed to have pdg code -pdg.
***   Only channels for which simulation data from Pythia simulations exist
***   are available here. More complex channels
***   like channels containing Higgs bosons etc that decay to standard model
***   particles are treated in each respective particle physics module
***   in src_model.
***   The currently implemented channels are
***   
***      pdg  Channel        Internal number (only listed temporarily)
***     ----  -------        ---------------
***        1  d d-bar        1
***        2  u u-bar        2
***        3  s s-bar        3
***        4  c c-bar        4
***        5  b b-bar        5
***        6  t t-bar        6
***       13  mu- mu+        10
***       15  tau- tau+      11
***       21  gluon gluon    7
***       23  Z0 Z0          9
***       24  W+ W-          8
***
***   Note: If a channel that is not simulated is asked for, the yield
***   0 is returned and a warning is issued (istat bit 3 set)
***
***   hel: helicity state of final state. Currently, only states where the
***        two final state particles are in the same helicity state (e.g.
***        both left, both right, both transverse, both longitudinal)
***        is included.
***   Possible values for hel:
***     'L': left-polarized fermion for fermions, longitudinal for vectors
***     'R': right-polarized fermion for fermions
***     'T': transverse for vectors
***     '0': unpolarized yields (unphysical, but included for comaprison)
***
***   - yieldpdg, PDG code for the yield type; currently the following is implemented:
***
***          yieldpg       yield type   
***          -------       ----------------  
***          22            cont. gamma rays
***          -11           positrons
***          -2212         antiprotons
***          -1000010020   anti-deuteron
***          111           pi0
***          12 or -12     nu_e and nu_e-bar
***          14 or -14     nu_mu and nu_mu-bar
***          16 or -16     nu_tau and nu_tau-bar
***          130072        muons from nu at creation
***          130073        muons from nu, as seen by a detector in ice
***                        (i.e. integrating 130072 over the mean muon path)
***
***    - diff: dictates whether differential source term at egev (diff=1) 
***            or integrated source term above egev (diff=0) is returned
***
***
*** Output:
***   - istat (=0 if no warnings/errors are reported)
***     bit 0 is set if extrapolations below the lowest simulated mass
***       or above the highest simulated mass is needed
***     bit 3 is set if the requested channel is not simulated,
***       yield zero is returned
***     bit 4 is set if the requested polarization state is not available,
***       the yield of the simulated polarization yield (typically
***       unpolarized) is returned instead
***   - dsanyield_sim, yield in units of
***   units: (annihilation)**-1  integrated
***   units: gev**-1 (annihilation)**-1 differential
***
*** If this routine is called outside of the kinematical regions
*** where tables exist, the following is done:
***   - for lower energies, than the lowest simulated ones,
***     extrapolations are used
***   - for masses below the lowest simulated, extrapolations
***     are used
***   - for energies above the highest simulated, the results
***     for the highest energy simulated are used
***
*** Note: at initialization of DarkSUSY, dsaninit should be called
*** to initialize these routines (done automatically in dsinit). This is
*** only needed once per run.
***
***
***   type : commonly used
***   desc : Simulated particle yields
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Modifications: 2010-11-05 (JE): better extrapolation below lowest
***   simulated energy (now dN/dz is taken as constant below)
*** mod tb -- added FSR from e+e-, linked light quarks to cc (prelim fix)
*** mod tb -- added pdg yield codes
*** mod je -- new simulation yields (including light quarks and anti-deuterons)
***	and internal channel numbers	
*****************************************************************************

      real*8 function dsanyield_sim(mwimp,e,pdg,hel,yieldpdg,diff,istat)
      implicit none
      include 'dsanyieldcom.h'
      include 'dsmpconst.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mwimp,e,xp,phi1,phi2,mp1,mp2,zpl,mn,z,
     &  tmp,lge,dsanyieldget,scalefactor,Euse
      integer ch,istat,zi,m1i,m2i,fltype,fi,yieldk, yieldpdg, diff,
     &    yieldkuse,fiuse
      real*8 pdb,td,nucleon,ed,ep,tp
      integer pdg,apdg
      integer i,j
      logical wb,chok
      parameter(lge=0.434294481903d0)
      character*1 hel

      real*8 pcoal0,sig_pcoal,pcoal
      parameter(pcoal0=0.079d0) ! Donato et al. arXiv:0803.2640, spherical coal
      parameter(sig_pcoal=0.015d0) ! Just a guess...
	
      logical first
      data first/.true./
      save first

c------------------------ functions ------------------------------------
      real*8 dsdilog

c-----------------------------------------------------------------------

      dsanyield_sim=0.d0

c...Determine yieldtype number from pdg codes

      yieldk=0
      if (yieldpdg.eq.22.or.yieldpdg.eq.-22) then      ! cont. gamma rays
         yieldk=52 
      elseif (yieldpdg.eq.-11) then                    ! positrons
         yieldk=51 
      elseif (yieldpdg.eq.-2212) then                  ! antiprotons
         yieldk=54
      elseif (yieldpdg.eq.-1000010020) then            ! anti-deuterons
         yieldk=dbflxk   ! (set from common, default in dsanyield_init)
      elseif (yieldpdg.eq.111) then                    ! neutral pions
         yieldk=58 
      elseif (yieldpdg.eq.12.or.yieldpdg.eq.-12) then  ! nu_e and nu_e-bar
         yieldk=56 
      elseif (yieldpdg.eq.14.or.yieldpdg.eq.-14) then  ! nu_mu and nu_mu-bar
         yieldk=53 
      elseif (yieldpdg.eq.16.or.yieldpdg.eq.-16) then  ! nu_tau and nu_tau-bar
         yieldk=57 
      elseif (yieldpdg.eq.130072) then                 ! mu from nu @ creation
         yieldk=72 
      elseif (yieldpdg.eq.130073) then                 ! mu from nu in ice
         yieldk=73    
      else
         istat=ibset(istat,3) 
         return
      endif

      if (diff.eq.1) then
        yieldk = yieldk +100
      else
        if (diff.ne.0) then
          istat=ibset(istat,3) 
          return        
        endif
      endif  

      call dsandec(yieldk,fltype,fi)

      if (yieldk.eq.59) then
        write(*,*) 'DS ERROR in dsanyield_sim:',
     &     ' old spherical coalescence model for anti-deuterons'
        write(*,*) 'only exists for differential yields. Stopping.'
        stop
      endif
	 
c--------------------------------------- if first call, load tables
      if (first) then
        do i=1,ntype
          do j=1,2
            yieldtype(j,i)=0
          enddo
        enddo

        first=.false.
        if (yieldk.eq.159) then
          call dsanreadfiles(159)
        else ! for all other cases, load corresponding file
           call dsanreadfiles(yieldk)
        endif
      endif

      if (yieldtype(fltype,fi).eq.0) then
        if (yieldk.eq.159) then
          call dsanreadfiles(159)
        else ! for all other cases, load corresponding file
           call dsanreadfiles(yieldk)
        endif
      endif


c...Set up energy and yield variables to use (the same as input, except for
c...anti-deuterons in old spherical coalescence model).

      euse=e
      yieldkuse=yieldk
      fiuse=fi
      pdb=0.0d0
      pcoal=0.0d0
	
      if (yieldk.eq.159) then ! dbar in old spherical coal. mod
        pcoal=pcoal0+dbp0bar*sig_pcoal
c...calculate kinematical variables, need momentum, e is kinetic energy of dbar
        pdb=sqrt(e**2+2.0d0*m_d*e)
c... shift yieldk to the value for pbar
        yieldkuse=154 ! pbar
        fiuse=4 ! pbar
        td=e/2.d0  ! kinetic energy per nucleon of dbar
        nucleon=2.d0
        ed=nucleon*td+m_d ! total energy of dbar
        ep=ed/2.d0 ! proton energy
        tp=ep-m_p  ! proton kinetic energy
        if(tp.lt.0.d0) then
          write(*,*) 'negative tp in dsanyield_sim'
          write(*,*) 'td, tp = ',td,tp
          tp=0.d0
        endif
        euse=tp ! kinetic energy of proton
      endif 
	
c...Determine channel number from pdg codes

      ch=0
      apdg=abs(pdg)
      if (apdg.eq.1) then      ! d d-bar
         if (hel.eq.'0') then
            ch=1 
         else ! eventually we will have sims here
            ch=1 
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.2) then      ! u u-bar
         if (hel.eq.'0') then
            ch=2 
         else ! eventually we will have sims here
            ch=2 
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.3) then      ! s s-bar
         if (hel.eq.'0') then
            ch=3 
         else ! eventually we will have sims here
            ch=3 
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.4) then      ! c c-bar
         if (hel.eq.'0') then
            ch=4
         else ! eventually we will have sims here
            ch=4
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.5) then ! b b-bar
         if (hel.eq.'0') then
            ch=5
         else ! eventually we will have sims here
            ch=5
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.6) then ! t t-bar
         if (hel.eq.'0') then
            ch=6
         else ! eventually we will have sims here
            ch=6
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.15) then ! tau- tau+
         if (hel.eq.'0') then
            ch=11
         else ! eventually we will have sims here
            ch=11
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.24) then ! W- W+
         if (hel.eq.'0') then
            ch=8
         else ! eventually we will have sims here
            ch=8
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.23) then ! Z0 Z0
         if (hel.eq.'0') then
            ch=9
         else ! eventually we will have sims here
            ch=9
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.11) then ! e- e+
         dsanyield_sim=0.d0
              ! NB: in the absence of sims, we take Weizsæcker-Williams 
              !     for photons here, zero yield for everything else
         xp=euse/mwimp

         if (yieldk.eq.152.and.xp.lt.(1-(m_e/mwimp)**2))  ! differential yield 
     %     dsanyield_sim=alpha_em/pi*(1+(1-xp)**2)/xp*
     &                   log(4.*mwimp**2*(1-xp)/m_e**2)
     
         if (yieldk.eq.52.and.xp.lt.(1-(m_e/mwimp)**2))   ! integrated yield
     &     dsanyield_sim=alpha_em/pi*(15 - 4*Pi**2 + 3*(-6 + xp)*xp - 
     &        6*(-1 + (-2 + xp)**2)*Log(4 - 4*xp) - 48*Log(xp)*Log(2*mwimp/m_e)
     &       + 12*(-3 + xp)*(-1 + xp)*Log(m_e/mwimp) + 24*dsdilog(xp))/12.
     
         if (hel.ne.'0') istat=ibset(istat,4)
         return
      elseif (apdg.eq.13) then ! mu- mu+
         if (hel.eq.'0') then
            ch=10
         else ! eventually we will have sims here
            ch=10
            istat=ibset(istat,4)
         endif
      elseif (apdg.eq.21) then ! g g
         if (hel.eq.'T') then
            ch=7
         else ! 
            ch=7
            istat=ibset(istat,4)
         endif
      else
         istat=ibset(istat,3) 
         dsanyield_sim=0.d0
         return
      endif
         


c      write(*,*) 'dsanyield_sim called with: mwimp = ',mwimp
c      write(*,*) '  e = ',e

      wb=.true.

      if (euse.ge.mwimp) then
        dsanyield_sim=0.0d0
        return
      endif

      mn = mwimp

c...Check kinematics
      chok=.false.
      if (mn.ge.msim(ch)) chok=.true.
      if ((ch.eq.6).and.
     &     mn.ge.(0.95*msim(ch))) chok = .true. ! t t-bar
      if ((ch.eq.8.or.ch.eq.9).and.
     &     mn.ge.(0.98*msim(ch))) chok = .true. ! W+W- or Z0 Z0 
      if (.not.chok) istat=ibset(istat,0)
      
c...take care of the case where mwimp is between the mass of the annihilation
c...product and the lower bound of the simulations (there might be a small
c...gap of less than a gev) for tt-bar, ww and zz.
      scalefactor=1.d0 ! use yields directly as they are
      if (mn.lt.lb(ch)) then
         if (mn.lt.0.98d0*lb(ch)) istat=ibset(istat,0)
         mn=lb(ch)
         scalefactor=mwimp/mn ! rescale energies, i.e. assume dN/dz is OK
      endif

      if (mn.gt.ub(ch)) then
        if (prtlevel.gt.0) then
          write(6,*)
          write(6,*) 'warning in dsanyield_sim for model: ',idtag
          write(6,5000) 'a WIMP mass of ',mn,
     +      ' GeV wants to be used,'
          write(6,5010) 'while the upper bound for channel ',ch,
     +      ' is ',ub(ch),' gev.'
          write(6,5020) 'a WIMP mass of ',ub(ch),' gev is used',
     +      ' instead for these too high masses.'
          write(6,*) 'the results can thus only be trusted as a',
     +      ' lower bound.'
        endif
        istat=ibset(istat,0)
      endif


c---------------------------------------------------- integrated yields
      if (wb.and.fltype.eq.1) then
c...determine which entries in phiint to use and how
        if (fi.ge.1.and.fi.le.20) then ! log tabulated yields
          z=(log10(euse/mwimp)+ndec)/ndec
        else
          z=euse/mwimp
        endif
        if (z.lt.0.0d0) then
          dsanyield_sim=0.0d0
          return
        endif

        if (fi.lt.9.or.fi.gt.13) then ! not dbar
           call dsanifind(z,zindex(0,1),zpl,zi,0,zn-1)
        else                    ! dbar
           call dsanifind(z,dbzindex(0,1),zpl,zi,0,zndb-1)
        endif
	
        if (zi.eq.-5.or.zi.ge.zn) then
          dsanyield_sim=0.0d0
          return
        endif

        call dsanifind(mn,mi(1),tmp,m1i,1,nmass-1) ! mn here
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(nmass)) then
          m1i=nmass
          m2i=nmass
          mp1=mi(nmass)
          mp2=mp1

          dsanyield_sim =
     &      (1.0-zpl)*dsanyieldget(zi,m1i,ch,fiuse,fltype,istat)+
     &      zpl*dsanyieldget(zi+1,m1i,ch,fiuse,fltype,istat)
        else
          phi1 =
     &      (1.0-zpl)*dsanyieldget(zi,m1i,ch,fiuse,fltype,istat)+
     &      zpl*dsanyieldget(zi+1,m1i,ch,fiuse,fltype,istat)
          phi2 =
     &      (1.0-zpl)*dsanyieldget(zi,m2i,ch,fiuse,fltype,istat)+
     &      zpl*dsanyieldget(zi+1,m2i,ch,fiuse,fltype,istat)
          dsanyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
     &      ((mp2-mp1)*(mp2+mp1))
        endif
c        if (mn.ne.mwimp) then
c          dsanyield_sim=mwimp**2/mn**2*dsanyield_sim
c        endif
      endif

c-------------------------------------------------- differential yields
      if (wb.and.fltype.eq.2) then
c...determine which entries in phidiff to use and how
        if (fi.ge.1.and.fi.le.20) then ! log tabulated yields
          z=(log10(euse/mwimp)+ndec)/ndec
        else
          z=euse/mwimp
        endif
        if (z.lt.0.0d0) then
          dsanyield_sim=0.0d0
          return
        endif

        if (fi.lt.9.or.fi.gt.13) then ! not dbar
           call dsanifind(z,zindex(-1,2),zpl,zi,-1,zn-1)
        else                    ! dbar
           call dsanifind(z,dbzindex(-1,2),zpl,zi,-1,zndb-1)
        endif

        if (zi.eq.-5.or.zi.ge.zn) then
          dsanyield_sim=0.0d0
          return
        endif

        if ((e.gt.mn).or.(e.le.0.0).or.zi.eq.-5) then
          dsanyield_sim=0.0d0
          return
        endif

        call dsanifind(mn,mi(1),tmp,m1i,1,nmass-1) ! here we use mn
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(nmass)) then
          m1i=nmass
          m2i=nmass
          mp1=mi(nmass)
          mp2=mp1
          dsanyield_sim =
     &      (1.0-zpl)*dsanyieldget(zi,m1i,ch,fiuse,fltype,istat)+
     &      zpl*dsanyieldget(zi+1,m1i,ch,fiuse,fltype,istat)
        else
          phi1 =
     &      (1.0-zpl)*dsanyieldget(zi,m1i,ch,fiuse,fltype,istat)+
     &      zpl*dsanyieldget(zi+1,m1i,ch,fiuse,fltype,istat)
          phi2 =
     &      (1.0-zpl)*dsanyieldget(zi,m2i,ch,fiuse,fltype,istat)+
     &      zpl*dsanyieldget(zi+1,m2i,ch,fiuse,fltype,istat)
          dsanyield_sim = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
     &      ((mp2-mp1)*(mp2+mp1))
        endif
c        if (mn.ne.mwimp) then
c          dsanyield_sim=mwimp**2/mn**2*dsanyield_sim
c        endif
c... convert from dyield/dz or dyield/dx to dyield/de
        if (fiuse.ge.1.and.fiuse.le.20) then ! log tabulated yields
          dsanyield_sim=dsanyield_sim*lge/(ndec*euse)
        else
          dsanyield_sim=dsanyield_sim/mwimp
        endif
        if (yieldk.eq.159) then
          dsanyield_sim=0.25d0*dsanyield_sim**2
          dsanyield_sim=(4.d0/3.d0*pcoal**3/pdb)*m_d/m_p**2 
     &    *dsanyield_sim*0.5d0
c...Final factor of 0.5d0 comes from that we here give yield per dbar kinetic
c...energy (not per kinetic energy per nucleon)

        endif

      endif


 5000 format(' ',a,f8.2,a)
 5010 format(' ',a,i2,a,f8.2,a)
 5020 format(' ',a,f8.2,a,a)

      end

