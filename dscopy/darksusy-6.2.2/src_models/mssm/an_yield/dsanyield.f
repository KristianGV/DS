*****************************************************************************
***   function dsanyield gives the total yield of positrons, cont. gammas
***   or neutrinos coming from WIMP annihilation in the halo
***   the yields are given as number / annihilation. the energy egev
***   is the threshold for integrated yields and the energy for
***   differential yields. the yields are
***     yieldk =  51: integrated positron yields
***     yieldk =  52: integrated cont. gammas
***     yieldk =  53: integrated muon neutrinos
***     yieldk =  54: integrated antiproton yields
***     yieldk =  59: integrated anti-deuteron yields (old. sph. coal. model)
***                   obsolete, use 61 instead. Only works for differential
***                   yields      
***     yieldk =  61: integrated anti-deuteron yields (Pythia6 MC, default dbar)
***     yieldk =  71: integrated neutrino yields (same as 53)
***     yieldk =  72: integrated muon yields at creation
***     yieldk =  73: integrated muon yields in ice
***     yieldk = above+100: differential in energy
***
*** The scalar decay properties are setup in this routine before
*** calling the cascade routines in an_yield_casc
***
*** istat will be zero in case of no errors. Otherwise its bits are set as
***   bit  decimal   reason
***     0        1   some inaccesible parts the differential spectra
***                  has been wanted, and the returned yield should then
***                  be treated as a lower bound.
***     1        2   energetically forbidden annihilation channels have been
***                  wanted.
***     2        4   problems with dsIBf_intdxdy integration for IB yields
***                  (only used for MSSM)
***     3        8   problems with dsIBf_intdy integration for IB yields
***                  (only used for MSSM)

*** author: joakim edsjo  edsjo@fysik.su.se
*** date: 98-01-29
*** modified: 98-04-15
*** modified: 2007-05-01 Torsten Bringmann (IB contribution added)
*** modified: 2008-01-15 Joakim Edsjo, made more modular
*** modified: April, 2014. Joakim Edsjo, split simulated part into
*** src/an_yield and model-dependent cascade decays into 
*** src_models/mssm/an_yield_casc.
*** modified: 2011-11-11 Torsten Bringmann (changed to PDG code sum)
*** modified: 2016-02-06 Torsten Bringmann (added electroweak IB)
*** modified: 2016-04-06 Torsten Bringmann (added gluon IB)
*** modified: 2016-04-12 Joakim Edsjo (added dbars, both MC and old sph. coal)
*****************************************************************************

      real*8 function dsanyield(egev,yieldk,istat)
      implicit none
      include 'dsanyieldmodelcom.h'
      include 'dsanyieldcom.h'
      include 'dsidtag.h'
      include 'dsmpconst.h'

*** FIXME: use PDG input instead of yieldk here...
* once done: make sure to remove obsolete code that maps PDG codes to yieldk
* before calling dsanyield!


c------------------------ functions ------------------------------------

      real*8 dsibyield,dsib2yield,dsanyield_ch
      real*8 dsanyield_sim_ls, dssigmav0, dssigmav0tot, dsmwimp, dsib3yield

c------------------------ variables ------------------------------------

      real*8 egev,yield,tmp, mDM
      integer ch,istat,jstat,yieldk, yieldkk, yieldpdg, diff

      real*8 pcoal0,sig_pcoal,pcoal,pdb,tp,nucleon,ep,ed,td
      parameter(pcoal0=0.079d0) ! Donato et al. arXiv:0803.2640, spherical coal
      parameter(sig_pcoal=0.015d0) ! Just a guess...

c----------------------------------------------- set-up common variables

      anerr=0

c...loop through all channels that give  calculate the yield above threshold
c...for each channel.

c      write(*,*)
c      write(*,*) 'model: ',idtag,'  eth = ',egev

c... map to PDG code input required by dsanyield_sim

      yieldkk=mod(yieldk,100)
      diff=yieldk/100
      if (diff.ne.0.and.diff.ne.1) then
        write(*,*) 'ERROR in dsanyield: unspoorted argument yieldk =', yieldk
        dsanyield=0.0d0
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
         yieldpdg = -1000010020 ! anti-deuteron yields, old sph. coal. mod
         call dsanyield_dbset(59,-100.d0)
      elseif (yieldkk.eq.71) then 
         yieldpdg = 14  ! neutrino yields (same as 53)    
      elseif (yieldkk.eq.72) then 
         yieldpdg = 130072 ! muon yields at creation     
      elseif (yieldkk.eq.73) then 
         yieldpdg = 130073 ! integrated muon yields in ice     
      else
        write(*,*) 'ERROR in dsanyield: unspoorted argument yieldk =', yieldk
        dsanyield=0.0d0
        return
      endif

      anistat=0
      yield=0.0d0
      mDM=dsmwimp()
      do 100 ch=1,numanch2b
         jstat=0
         if (dssigmav0(anch_2body(ch,1),anch_2body(ch,2),jstat)
     &        .gt.0d0) then
            tmp= dsanyield_sim_ls(mDM,egev, 
     &           anch_2body(ch,1),anch_2body(ch,2),
     &           anch_2body(ch,3),anch_2body(ch,4),
     &           anch_2body(ch,5),anch_2body(ch,6),
     &          yieldpdg,diff,jstat)
        
            if (btest(jstat,3)) then ! channel not simulated!
               if (yieldk.eq.159) then ! old sph. coal. model
                  pcoal=pcoal0+dbp0bar*sig_pcoal
c...calculate kinematical variables, need momentum, e is kinetic energy of dbar
                  pdb=sqrt(egev**2+2.0d0*m_d*egev)
                  td=egev/2.d0  ! kinetic energy per nucleon of dbar
                  nucleon=2.d0
                  ed=nucleon*td+m_d ! total energy of dbar
                  ep=ed/2.d0 ! proton energy
                  tp=ep-m_p  ! proton kinetic energy
                  if(tp.lt.0.d0) then
                     write(*,*) 'negative tp in dsanyield'
                     write(*,*) 'td, tp = ',td,tp
                     tp=0.d0
                  endif
                  tmp=dsanyield_ch(mDM,tp,anch_2body(ch,1),
     &                 anch_2body(ch,2),154,istat)
                  tmp=0.25d0*tmp**2
     &                 *(4.d0/3.d0*pcoal**3/pdb)*m_d/m_p**2 *0.5d0
c...Final factor of 0.5d0 comes from that we here give yield per dbar kinetic
c...energy (not per kinetic energy per nucleon)
               else ! default for all other channels
                 tmp=dsanyield_ch(mDM,egev,anch_2body(ch,1),
     &                anch_2body(ch,2),yieldk,istat)
                 anistat=or(anistat,istat)
               endif
            endif

            yield=yield+
     &      tmp*dssigmav0(anch_2body(ch,1),anch_2body(ch,2),jstat)
         endif    
  100 continue

c... now normalize to total annihilation rate
      yield = yield / dssigmav0tot()

c...add photon bremsstrahlung (IB)
      yield=yield+dsibyield(egev,yieldk,istat)
      anistat=or(anistat,istat*4)

c...add electroweak bremsstrahlung (IB2)
c... NB: This contribution is switched off by default because it takes O(10s)
c... to compute. Call dsib2set('standard') before calling this routine (or any
c... gamma or cosmic ray routine) to fully include these contributions!
      yield=yield+dsib2yield(egev,yieldk,0,istat)
      anistat=or(anistat,istat*8)

      yield=yield+dsib3yield(egev,yieldk,istat)
      anistat=or(anistat,istat*16)


      dsanyield=yield
      istat=anistat

      end



















