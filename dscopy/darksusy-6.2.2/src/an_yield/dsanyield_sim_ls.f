*****************************************************************************
*** NOTE: This routine is not fully functional yet, but shows what we intend
*** to have eventually. The goal is to have a much more general structure
*** regarding which channels and polarization states that are available.
*** This routine should eventually be able to return the yield for different
*** final state particles and polarization states, expressed as the final
*** state's quantum numbers j,P,l and s.
***      
*** function dsanyield_sim_ls calculates the yield above threshold
*** (or differential at that energy) for the requested annihilation channel
*** and the fluxtype given by yieldpdg. This routine assumes that annihilation
*** takes place to two final state particles.
***
*** Inputs:
***   - mwimp = WIMP mass in GeV
***   - e = kinetic energy where the yield is calculated (in GeV)
*** 
***   - pdg1, pdg2 = pdg codes of annihilation final state particles.
***   Only channels for which simulation data from Pythia simulations exist
***   are available here. More complex channels
***   like channels containing Higgs bosons etc that decay to standard model
***   particles are treated in each respective particle physics module
***   in src_model.
***   The currently implemented channels are
***   
***     pdg1  pdg2  Channel        Internal number (only listed temporarily (old new))
***     ----  ----  -------        ---------------
***        1    -1  d d-bar        -   1
***        2    -2  u u-bar        -   2
***        3    -3  s s-bar        -   3
***        4    -4  c c-bar        1   4
***        5    -5  b b-bar        2   5
***        6    -6  t t-bar        3   6
***       15   -15  tau- tau+      4  11
***       24   -24  W+ W-          5   8
***       23    23  Z0 Z0          6   9
***       13   -13  mu- mu+        7  10
***       21    21  gluon gluon    8   7
***
***   Note: If a channel that is not simulated is asked for, the yield
***   0 is returned and a warning is issues (istat bit 3 set)
***
***   For the final state polarization, we need a few arguments to
***   describe it fully.
***     twoj: total angular momentum quantum number of final state particles
***         times 2.
***     p: parity quantum number: JE/TB FIXME: Use CP instead?
***     twol: orbital angular momentum quantum number of final state times 2
***     twos: spin quantum number of final state times 2
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
***   - dsanyield_sim_ls, yield in units of
***   units: (annihilation)**-1  integrated
***   units: gev**-1 (annihilation)**1 differential
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
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Modifications: 2010-11-05 (JE): better extrapolation below lowest
***   simulated energy (now dN/dz is taken as constant below)
*** mod tb -- added pdg yield codes
*****************************************************************************

      real*8 function dsanyield_sim_ls(mwimp,e,pdg1,pdg2,twoj,p,
     &  twol,twos,yieldpdg,diff,istat)
      implicit none
      include 'dsanyieldcom.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mwimp,e,dsanyield_sim
      integer istat,yieldpdg,diff
      integer pdg1,pdg2,apdg1,apdg2,twoj,twol,twos, p
      character*1 cL,cR,cT,c0

      logical first
      data first/.true./
      save first

      data cL,cR,cT,c0/'L','R','T','0'/

      save cL,cR,cT,c0
c-----------------------------------------------------------------------

      dsanyield_sim_ls=0.d0

c...Right now we don't have these simulated explicitly, so call
c...dsanyield_sim for the moment
c...JE FIXME: update with more simulations and better calls to dsanyield_sim
c...Need to work out the Clebsch-Gordan coefficients to use here
      apdg1=abs(pdg1)
      apdg2=abs(pdg2)
      if (apdg1.eq.apdg2) then
         dsanyield_sim_ls=dsanyield_sim(mwimp,E,apdg1,'0',yieldpdg,
     &     diff,istat)
      else
         istat=ibset(istat,3)
         dsanyield_sim_ls=0.d0
      endif
         
      return
      end

