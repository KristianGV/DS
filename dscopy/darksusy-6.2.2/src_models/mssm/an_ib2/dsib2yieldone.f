************************************************************************
*** function dsIB2yieldone calculates the electroweak IB yield from ONE 
*** individual 3-body annihilation channel. NB: For the default case
*** (on-shell constribution already subtracted), only the *sum* over all
*** these channels gives a well-defined yield (as returned by dsib2yield)!!! 
***
***   Input: egev  - energy in GeV
***          yieldk  - which yield to calculate
***                    (see dsanyield_ch for an explanation)
***          IB2ch   - 1XX for ffZ
***                    2XX for WfF
***                    3XX for ffh (XX like for ffZ)
***                    4XX for ffH (XX like for ffZ)
***                    5XX for ffA (XX like for ffZ)
***                    6XX for ffH+- (XX like for fFW)
***
***                    XX   | \bar f f Z      |   \bar f f W
***                   ------+-----------------+-------------------
***                    01   | nu_e   nu_e   Z |   e+     nu_e   W-
***                    02   | e      e      Z |   nu_e   e-     W+
***                    03   | nu_mu  nu_mu  Z |   mu+    nu_mu  W-
***                    04   | mu     mu     Z |   nu_mu  mu-    W+
***                    05   | nu_tau nu_tau Z |   tau+   nu_tau W-
***                    06   | tau    tau    Z |   nu_tau tau-   W+
***                    07   | u      u      Z |   dbar   u      W-
***                    08   | d      d      Z |   ubar   d      W+
***                    09   | c      c      Z |   sbar   c      W-
***                    10   | s      s      Z |   cbar   s      W+
***                    11   | t      t      Z |   bbar   t      W-
***                    12   | b      b      Z |   tbar   b      W+
***                    (note that this numbering must be consistent 
***                     with dmssm.h !!! )
***
***  TODO: CHECK that 2nd column is indeed everywhere consistent!!!
***
***          onshell - 0 subtracts the on-shell contributions ("2-body final 
***                      states")in the narrow width approximation [default]
***                      (note that the result can be negative!)
***                    1 returns the result including on-shell contributions
***
***   Output: yield [total tree-level(!) annihilation]**-1, 
***                  differential also [GeV]**-1
***           -> use dssigmav() & dsib2sigmav() to manually adjust 
***              normalization to refer to individual channels instead.
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-12-16, update 2015-03-22 (streamlined handling of c_ftype),
***         update 2015-11-19 (added neutrino and positron spectra)
***         update 2016-02-05 updated to DS6 conventions
************************************************************************

      real*8 function dsIB2yieldone(egev,IB2ch,yieldk,onshell,istat)
      implicit none
      include 'dsanyieldmodelcom.h'
      include 'dsanyieldcom.h'
      include 'dsmssm.h'
      include 'dsib2com.h'
      include 'dsmpconst.h'
      include 'dsio.h'

c------------------------ functions ------------------------------------

      real*8 dssigmav0, dssigmav0tot,dsIB2convint,dsIB2intres, dsanyield_ch
      real*8 dsIB2dsde_aux, dsIB2BRtree, dsanyield_sim
      external dsIB2convint,dsIB2dsde_aux

c------------------------ variables ------------------------------------

      real*8 egev
      integer IB2ch,istat,yieldk, onshell
      
      real*8 result,tmpresult,tmpres2, mfinal1, mfinal2, mfinal3
      real*8 zmin,zmax,mstable, Tkinmax, xkinmax, sv
      real*8 nwapart, nwastable
      integer err,pdgyield, diff, ier, i, j
      real*8 ans0brsav(30,3), anscbrsav(21), Fwidth

      logical direct, indirect

c-----------------------------------------------------------------------

      result=0.d0
      tmpresult=0.d0
      tmpres2=0d0
      dsIB2yieldone=0.d0
      istat=0

      direct = .true.   ! keep contribution from stable final particles
      indirect = .true. ! keep contribution from unstable final particles


      pdgyield=0
      diff=yieldk/100   ! diff=0: integrated; diff=1: differential
c TB FIXME: positron yield needs fixing -> problem with propagation routines !?
c           preliminary set to zero 
c      if (mod(yieldk,100).eq.51) pdgyield=-11           ! positrons
      if (mod(yieldk,100).eq.52) pdgyield=22            ! cont. gamma rays
      if (mod(yieldk,100).eq.54) pdgyield=-2212         ! antiprotons
      if (mod(yieldk,100).eq.56) pdgyield=12            ! nu_e and nu_e-bar
      if (mod(yieldk,100).eq.71.or.mod(yieldk,100).eq.53) pdgyield=14      ! nu_mu and nu_mu-bar
      if (mod(yieldk,100).eq.57) pdgyield=16            ! nu_tau and nu_tau-bar
      if (mod(yieldk,100).eq.72) pdgyield=130072        ! muons from nu at creation
      if (mod(yieldk,100).eq.73) pdgyield=130073        ! muons from nu, as seen by a detector in ice
                                                        ! (i.e. integrating 130072 over the mean muon path)
      if (pdgyield.eq.0.or.diff.lt.0.or.diff.gt.1) then
        if (prtlevel.ge.2) then
          write(*,*) 'WARNING in dsib2yieldone: called with unsopported ',
     &               'parameter yieldk = ',yieldk
          write(*,*) 'Returning with zero yield...'
        endif
        return
      endif
          
 

c... set masses and internal channel specifications
      call dsib2chinit(IB2ch,err)
      if (err.eq.1) then
         write(*,*) 'ERROR in dsIB2sigmav: unknown channel', IB2ch
         stop
      endif  

      mstable=1.d15
      if (mod(yieldk,100).eq.51) mstable=mass(ke)  ! e+
      if (mod(yieldk,100).eq.52) mstable=0.d0      ! gamma
      if (mod(yieldk,100).eq.53) mstable=0.d0      ! nu_mu
      if (mod(yieldk,100).eq.54) mstable=m_p       ! pbar
      if (mod(yieldk,100).eq.56) mstable=0.d0      ! nu_e
      if (mod(yieldk,100).eq.57) mstable=0.d0      ! nu_tau
      Tkinmax=mx+mstable**2/mx/4.0d0-mstable ! CHECK: - why do we need this additional criterion?
                                             !        - do we use E <-> T consistently 
                                             !          (only relevant for pbars!)?
      if (egev.gt.Tkinmax) return

      if (2.d0.lt.1.001*(mB+Mf+Mff)) return
      c_xstable=egev/mx       ! rescaled energy of stable particle

cTB debug
c      goto 330

      if (.not.direct) goto 100 ! skip contributions from stable final states

c... let's first add the contribution from stable final states
c... A) positron or antineutrino yields from antifermion final states     
      tmpresult=0.d0
      xkinmax=Tkinmax/mx
      if ((mod(yieldk,100).eq.51.and.c_fbartype.eq.2).or.         ! positrons
     &    (mod(yieldk,100).eq.56.and.c_fbartype.eq.1).or.         ! nu_e_bar
     &    (mod(yieldk,100).eq.53.and.c_fbartype.eq.3).or.         ! nu_mu_bar
     &    (mod(yieldk,100).eq.57.and.c_fbartype.eq.5))  then      ! nu_tau_bar 
            c_pfinal='fbar'
c            write(*,*) 'TEST :', IB2ch,c_fbartype
           if (yieldk/100.eq.1) tmpresult=dsIB2dsde_aux(c_xstable)/mx    ! differential yield
           if (yieldk/100.eq.0)                                          ! integrated yield
     &        call dgadap(c_xstable,xkinmax,dsIB2dsde_aux,IB2acc,tmpresult)  
               
           result=result + tmpresult
      endif

c      write(*,*) 'mstable, result = ', mstable, c_xstable, result
      
c... B) neutrino yields from neutrino final states
      if ((mod(yieldk,100).eq.56.and.c_ftype.eq.1).or.
     &    (mod(yieldk,100).eq.53.and.c_ftype.eq.3).or.
     &    (mod(yieldk,100).eq.57.and.c_ftype.eq.5))  then        
            if (c_Btype.eq.2.or.c_Btype.eq.6) then  ! only need to calculate for charged Boson 
                                                    ! final state; for neutral bosons, we assume
                                                    ! gauge invariance and recycle the result 
                                                    ! from the last step
              c_pfinal='f'
              if (yieldk/100.eq.1) tmpresult= dsIB2dsde_aux(c_xstable)/mx 
              if (yieldk/100.eq.0) 
     &            call dgadap(c_xstable,xkinmax,dsIB2dsde_aux,IB2acc,tmpresult) 
            endif     
            result=result + tmpresult
      endif

 100  if (.not.indirect) goto 330 ! skip contributions from unstable final states 

c... now we add the contribution from unstable final states       
c... Here, we need the same preliminary fix as in IB routines: 
c... dsanyield _chuses different definitions of masses: adopt those for
ccc  kinematical integration limits

      mfinal1=1.d15
      if (c_Btype.eq.1) mfinal1 = msim(chcomp(12))/mx      ! 'mZ'
      if (c_Btype.eq.2) mfinal1 = msim(chcomp(13))/mx      ! 'mW'
      if (c_Btype.eq.3) mfinal1 = ans0m(2)/mx    !  SM Higgs
      if (c_Btype.eq.4) mfinal1 = ans0m(1)/mx    !  heavy Higgs
      if (c_Btype.eq.5) mfinal1 = ans0m(3)/mx    !  A0
      if (c_Btype.eq.6) mfinal1 = anscm/mx       ! 'mH+'
      mfinal2 = msim(chcomp(13+c_ftype))/mx       ! 'Mf'
      mfinal3 = msim(chcomp(13+c_fbartype))/mx    ! 'Mff'
c... temporary fix: map all light (= not implemented) quark channels to cc      
      if (c_ftype.ge.7.and.c_ftype.le.10) mfinal2 = msim(chcomp(22))/mx
      if (c_fbartype.ge.7.and.c_fbartype.le.10) mfinal3 = msim(chcomp(22))/mx

      if (2.d0.lt.1.001*(mfinal1+mfinal2+mfinal3)) return

      c_yieldk=yieldk   

c... now sum the three contributions to the convolution integral
c... start with contribution from the bosonic final state
      tmpresult=0.0d0
      c_pfinal='B'
      zmin=c_xstable
      if (c_xstable.lt.1.00d0*mfinal1) zmin=1.00d0*mfinal1
      zmax= 1.d0+(mfinal1**2-(mfinal2+mfinal3)**2)/4.d0
      if (zmax.le.zmin.or.c_xstable.gt.zmax) goto 200
      result=result + dsib2intres(dsIB2convint,c_pfinal,zmin,zmax,IB2acc,ier)
      if (ier.ne.0) write(*,*) 'dsib2yieldone: WARNING - error in ',
     &                          c_pfinal,' integration: ', ier, c_ftype


c... add contribution from *unstable* final state fermions
 200  continue
      if (c_ftype.eq.1.or.c_ftype.eq.2.or.c_ftype.eq.3.or.c_ftype.eq.5) goto 300

      c_pfinal='f'
      zmin=c_xstable
      if (c_xstable.lt.1.00*mfinal2) zmin=1.00d0*mfinal2
      zmax= 1.d0+(mfinal2**2-(mfinal1+mfinal3)**2)/4.d0
      if (zmax.le.zmin.or.c_xstable.gt.zmax) goto 300  
      tmpresult= dsib2intres(dsIB2convint,c_pfinal,zmin,zmax,IB2acc,ier)
      if (ier.ne.0) write(*,*) 'dsib2yieldone: WARNING - error in ',
     &                            c_pfinal,' integration: ', ier, c_ftype
      
      result=result + tmpresult


c... add contribution from final state ANTIfermions, 
c... only need to compute new if final state fermions not identical
 300  continue 
      if (c_fbartype.eq.1.or.c_fbartype.eq.2.or.
     &    c_fbartype.eq.4.or.c_fbartype.eq.6) goto 330
 
      if (c_Btype.eq.2.or.c_Btype.eq.6) then  ! for neutral Bosons, we will just recycle
                                              ! the fermion result (just before label 330)    
         c_pfinal='fbar'
         zmin=c_xstable
         if (c_xstable.lt.1.00*mfinal3) zmin=1.00d0*mfinal3
         zmax= 1.d0+(mfinal3**2-(mfinal2+mfinal1)**2)/4.d0
         if (zmax.le.zmin.or.c_xstable.gt.zmax) goto 330  

         tmpresult= dsib2intres(dsIB2convint,c_pfinal,zmin,zmax,IB2acc,ier)
         if (ier.ne.0) write(*,*) 'dsib2yieldone: WARNING - error in ',
     &                             c_pfinal,' integration: ', ier, c_ftype
      endif
      
      result=result+tmpresult

      
330   sv=dssigmav0tot()
      result=result*gev2cm3s/sv
      result=result/1.d20  !correct for factor introduced in dsib2Msqaux
      

c... default: subtract on-shell contributions in NWA approximation
c MG changed call to BR calculation (BR for each fermion species individually)
c MG added factor 0.5 in dssigmav0(24,-37)=WH=W+H- + W-H+, because single contribution is required 

      nwapart=0.0d0
      if (onshell.ne.1) then

c... added by MG, 09/16; adopted to DS6 by TB
c... for our NWA subtraction to work, we need to make sure that the yields from the
c... *almost* on-shell decaying Scalars (but not those from the on-shell scalars!)
c... are based on branching ratios that
c...  a) are determined at tree-level 
c...  b) are re-scaled such the contribution from decays to fermions adds up to 1

c... First we need to save the common block values for Higgs decay branching ratios
        do j=1,3
           do i=1,29
               ans0brsav(i,j) = ans0br(i,j) ! neutral Higgses
           enddo
        enddo
        do i=1,15
            anscbrsav(i) = anscbr(i) ! charged Higgs
        enddo

c         goto 499

c... NB: The following uses the same hardcoded channel numbering as in dsanyieldset
c... and dsfeynhiggs.F (as well as the hardcoded channel coding of IB2ch)!
        if ((IB2ch/100.eq.1).or.(IB2ch/100.eq.5)) then ! change H and h rates
           Fwidth = 0.0d0
           do i=1,12
              Fwidth = Fwidth + dsIB2BRtree('hff',i)
           enddo
           do i=1,29
              ans0br(i,2) = 0d0
              ans0br(i,1) = 0d0
           enddo 
           do i=1,12
              ans0br(i+13,2) = dsIB2BRtree('hff',i)/Fwidth
           enddo
           Fwidth = 0.0d0
           do i=1,12
              Fwidth = Fwidth + dsIB2BRtree('Hff',i)
           enddo
           do i=1,12 
              ans0br(i+13,1) = dsIB2BRtree('Hff',i)/Fwidth
           enddo
        endif
        if ((IB2ch/100.eq.4).or.(IB2ch/100.eq.3)) then ! change A rates
           Fwidth = 0.0d0
           do i=1,12
              Fwidth = Fwidth + dsIB2BRtree('Aff',i)
           enddo
           do i=1,29
              ans0br(i,3) = 0d0
           enddo 
           do i=1,12
              ans0br(i+13,3) = dsIB2BRtree('Aff',i)/Fwidth
           enddo
        endif
        if (IB2ch/100.eq.2) then ! change H+- rates
           Fwidth = 0.0d0
           do i=1,10,2
              Fwidth = Fwidth + dsIB2BRtree('HfF',i)
           enddo
           Fwidth = Fwidth + dsIB2BRtree('Hbt',12)
           do i=1,15
              anscbr(i) = 0d0
           enddo 
           anscbr(1) = dsIB2BRtree('HfF',7)/Fwidth  ! u d
           anscbr(5) = dsIB2BRtree('HfF',9)/Fwidth  ! c s
           anscbr(9) = dsIB2BRtree('Hbt',12)/Fwidth ! t b
           anscbr(10) = dsIB2BRtree('HfF',1)/Fwidth ! e nu
           anscbr(11) = dsIB2BRtree('HfF',3)/Fwidth ! mu nu
           anscbr(12) = dsIB2BRtree('HfF',5)/Fwidth ! tau nu
        endif

c 499    continue

      
c... this is the standard case      
        if (.not.indirect) goto 500

c... FIXME: replace yields with polarized yields
        if (IB2ch.ge.201.and.IB2ch.le.210) then  ! WfF other than Wtb
          if (dssigmav0(24,-24).gt.0.d0) nwapart = dssigmav0(24,-24)/sv 
     &        *dsIB2BRtree('WfF',c_ftype)*dsanyield_sim(mx,egev,24,0,pdgyield,diff,ier)
          if (dssigmav0(24,-37).gt.0.d0) nwapart=nwapart+0.5*dssigmav0(24,-37)/sv*
     &               dsIB2BRtree('HfF',c_ftype)*dsanyield_ch(mx,egev,24,-37,yieldk,ier)
        elseif (IB2ch.eq.211.or.IB2ch.eq.212) then      ! Wtb
          if (dssigmav0(24,-37).gt.0.d0) nwapart=0.5*dssigmav0(24,-37)/sv*
     &               dsIB2BRtree('Hbt',11)*dsanyield_ch(mx,egev,24,-37,yieldk,ier)
          if (dssigmav0(6,-6).gt.0.d0) nwapart = nwapart + dssigmav0(6,-6)/sv*
     &               dsIB2BRtree('tWb',12)*dsanyield_sim(mx,egev,6,0,pdgyield,diff,ier)
        elseif (IB2ch.ge.101.and.IB2ch.le.112) then  ! Zff
          if (dssigmav0(23,23).gt.0.d0) nwapart = 2.*dssigmav0(23,23)/sv*
     &               dsIB2BRtree('Zff',c_ftype)*dsanyield_sim(mx,egev,23,0,pdgyield,diff,ier)
          if (dssigmav0(23,35).gt.0.d0) nwapart = nwapart + dssigmav0(23,35)/sv*
     &               dsIB2BRtree('Hff',c_ftype)*dsanyield_ch(mx,egev,23,35,yieldk,ier)
          if (dssigmav0(23,25).gt.0.d0) nwapart = nwapart + dssigmav0(23,25)/sv*
     &               dsIB2BRtree('hff',c_ftype)*dsanyield_ch(mx,egev,23,25,yieldk,ier)
        elseif (IB2ch.ge.601.and.IB2ch.le.610) then  ! HfF
          if (dssigmav0(24,-37).gt.0.d0) nwapart=0.5*dssigmav0(24,-37)/sv*
     &               dsIB2BRtree('WfF',c_ftype)*dsanyield_ch(mx,egev,24,-37,yieldk,ier)
       elseif (IB2ch.ge.501.and.IB2ch.le.512) then  ! Aff
          if (dssigmav0(25,36).gt.0.d0) nwapart = dssigmav0(25,36)/sv*
     &               dsIB2BRtree('hff',c_ftype)*dsanyield_ch(mx,egev,25,36,yieldk,ier)
          if (dssigmav0(35,36).gt.0.d0) nwapart = nwapart + dssigmav0(35,36)/sv*
     &               dsIB2BRtree('Hff',c_ftype)*dsanyield_ch(mx,egev,35,36,yieldk,ier)
        elseif (IB2ch.ge.401.and.IB2ch.le.412) then  ! Hff
          if (dssigmav0(35,36).gt.0.d0) nwapart = dssigmav0(35,36)/sv*
     &               dsIB2BRtree('Aff',c_ftype)*dsanyield_ch(mx,egev,35,36,yieldk,ier)
          if (dssigmav0(23,35).gt.0.d0) nwapart = nwapart + dssigmav0(23,35)/sv*
     &               dsIB2BRtree('Zff',c_ftype)*dsanyield_ch(mx,egev,23,35,yieldk,ier)
        elseif (IB2ch.ge.301.and.IB2ch.le.312) then  ! hff
          if (dssigmav0(25,36).gt.0.d0) nwapart = dssigmav0(25,36)/sv*
     &               dsIB2BRtree('Aff',c_ftype)*dsanyield_ch(mx,egev,25,36,yieldk,ier)
          if (dssigmav0(23,25).gt.0.d0) nwapart = nwapart + dssigmav0(23,25)/sv*
     &               dsIB2BRtree('Zff',c_ftype)*dsanyield_ch(mx,egev,23,25,yieldk,ier)
        endif     

 500    continue
c... COMMENT OUT THE "GOTO 800" IF INTERESTED IN A BETTER ESTIMATE OF THE YIELD 
c... *FOR INDIVIDUAL CHANNELS*
        goto 800
 
c... if both final state fermions are stable, the approximation about the yield
c... from those (which went into the above expressions) is particular bad, and
c... we need to correct for that to get a better estimate 
        if ((c_fbartype.le.3.or.c_fbartype.eq.5).and.
     &      (c_ftype.le.3.or.c_ftype.eq.5)) then
c... if one of those stable particles is required yield, we have to add it
c... (note that expressions assume massless stable fermions, resulting from the
c...  decay of an unpolarized Boson)       
          nwastable=0.0d0
          if (.not.direct) goto 550
          if ((mod(yieldk,100).eq.51.and.c_fbartype.eq.2).or.     ! positrons
     &        (mod(yieldk,100).eq.56.and.c_fbartype.eq.1).or.     ! nu_e_bar
     &        (mod(yieldk,100).eq.53.and.c_fbartype.eq.3).or.     ! nu_mu_bar
     &        (mod(yieldk,100).eq.57.and.c_fbartype.eq.5).OR.     ! nu_tau_bar 
     &        (mod(yieldk,100).eq.56.and.c_ftype.eq.1).or.        ! nu_e
     &        (mod(yieldk,100).eq.53.and.c_ftype.eq.3).or.        ! nu_mu
     &        (mod(yieldk,100).eq.57.and.c_ftype.eq.5))  then     ! nu_tau 

            if (IB2ch.ge.201.and.IB2ch.le.210) then  ! WfF other than Wtb
              nwastable = dssigmav0(24,-24)/sv* dsIB2BRtree('WfF',c_ftype)
     &                    +0.5*dssigmav0(24,-37)/sv*dsIB2BRtree('HfF',c_ftype)
            elseif (IB2ch.eq.211.or.IB2ch.eq.212) then      ! Wtb
              nwastable = 0.5*dssigmav0(24,-37)/sv*dsIB2BRtree('Hbt',11)
     &                   + dssigmav0(6,-6)/sv*dsIB2BRtree('tWb',12)
            elseif (IB2ch.ge.101.and.IB2ch.le.112) then  ! Zff
              nwastable = 2.*dssigmav0(23,23)/sv*dsIB2BRtree('Zff',c_ftype)
     &                   + dssigmav0(23,35)/sv*dsIB2BRtree('Hff',c_ftype)
     &                   + dssigmav0(23,25)/sv*dsIB2BRtree('hff',c_ftype)
            elseif (IB2ch.ge.601.and.IB2ch.le.610) then  ! HfF
              nwastable=0.5*dssigmav0(24,-37)/sv*dsIB2BRtree('WfF',c_ftype)
            elseif (IB2ch.ge.501.and.IB2ch.le.512) then  ! Aff
              nwastable = dssigmav0(25,36)/sv*dsIB2BRtree('hff',c_ftype)
     &                   + dssigmav0(35,36)/sv*dsIB2BRtree('Hff',c_ftype)
            elseif (IB2ch.ge.401.and.IB2ch.le.412) then  ! Hff
              nwastable = dssigmav0(35,36)/sv*dsIB2BRtree('Aff',c_ftype)
     &                   + dssigmav0(23,35)/sv*dsIB2BRtree('Zff',c_ftype)
            elseif (IB2ch.ge.301.and.IB2ch.le.312) then  ! hff
              nwastable = dssigmav0(25,36)/sv*dsIB2BRtree('Aff',c_ftype)
     &                   + dssigmav0(23,25)/sv*dsIB2BRtree('Zff',c_ftype)
            endif     
 
            if (yieldk/100.eq.1) nwastable=nwastable/mx        ! diff. yield
            if (yieldk/100.eq.0) nwastable=nwastable*(1.-egev) ! integrated yield

c... for neutrino yields, *both* nu and nubar contribute
            if (c_ftype.eq.c_fbartype.and.(mod(c_ftype,2).eq.1))
     &        nwastable = 2.*nwastable       

          endif
          
 550      nwapart = 0.5*nwapart + nwastable       
        endif 

 800  continue

c added by MG 09/16
c... Restore common block values of Higgs decay braching fractions
      do j=1,3
         do i=1,29
               ans0br(i,j) = ans0brsav(i,j)
         enddo
      enddo
      do i=1,15
            anscbr(i) = anscbrsav(i)
      enddo



      endif  ! end of NWA substraction prescription

                         
      dsIB2yieldone=result - nwapart
      istat = ier

      return
      end


