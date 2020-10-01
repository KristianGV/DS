*****************************************************************************
*** NOTE: This routine is not fully functional yet, but shows what we intend
*** to have eventually. The goal is to have a much more general structure
*** regarding which channels and polarization states that are avaialble.
*** This routine should eventually be able to return the yield for different
*** final state particles and polarization states, expressed as the final
*** state's quantum numbers j,P,l and s. Right now it just calls the old
*** routine, but eventually we want to add more stuff here      
***      
*** function dsseyield_sim_ls calculates the yield above threshold
*** (or differential at that energy and angle) for the requested
*** annihilation channel and the kind and type of yield.
*** This routine assumes that annihilation takes place to two final
*** state particles.
***
*** Inputs:
*** Inputs:
***   - mwimp: mass of WIMP in GeV
***   - E:     kinetic energy of the neutrino, lepton or hadronic shower
***            where the yield is calculated (in GeV)
***   - theta: angle where the yield is calculated (in degrees)
*** 
***   - pdg1, pdg2 = pdg codes of annihilation final state particles.
***   Only channels for which simulation data from Pythia simulations exist
***   are available here. More complex channels
***   like channels containing Higgs bosons etc that decay to standard model
***   particles are treated in each respective particle physics module
***   in src_model.
***   The currently implemented channels are
***   
***     pdg1  pdg2  Channel          Channel No chi  Array Ch No chii
***                                  (only listed temporarily)
***     ----  ----  -------           --------------  ----------------
***        1    -1  d d-bar            1              1
***        2    -2  u u-bar            2              2
***        3    -3  s s-bar            3              3
***        4    -4  c c-bar            4              4
***        5    -5  b b-bar            5              5
***        6    -6  t t-bar            6              6
***       12   -12  nu_e nu_e-bar     12              11
***       13   -13  mu- mu+           10               - not simulated
***       14   -14  nu_mu nu_mu-bar   13              12
***       15   -15  tau- tau+         11              10
***       16   -16  nu_tau nu_tau-bar 14              13
***       24   -24  W+ W-              8              8
***       23    23  Z0 Z0              9              9
***       21    21  gluon gluon        7              7
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
***
*** Note: at initialization of DarkSUSY, dsseinit should be called
*** to initialize these routines (done automatically in dsinit). This is
*** only needed once per run.
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*****************************************************************************

      real*8 function dsseyield_sim_ls(mwimp,e,theta,pdg1,pdg2,twoj,p,
     &  twol,twos,wh,kind,type,istat)
      implicit none
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mwimp,e,theta,lge,dsseyield_sim
      integer istat
      integer pdg1,pdg2,apdg1,apdg2,twoj,twol,twos,p,kind,type
      parameter(lge=0.434294481903d0)
      character*1 cL,cR,cT,c0
      character*2 wh
      
      logical first
      data first/.true./
      save first

      data cL,cR,cT,c0/'L','R','T','0'/

      save cL,cR,cT,c0
c-----------------------------------------------------------------------

      dsseyield_sim_ls=0.d0

c...Right now we don't have these simulated explicitly, so call
c...dsanyield_sim for the moment
c...JE FIXME: update with more simulations and better calls to dsanyield_sim
c...Need to work out the Clebsch-Gordan coefficients to use here
      apdg1=abs(pdg1)
      apdg2=abs(pdg2)
      if (apdg1.eq.apdg2) then
         dsseyield_sim_ls=dsseyield_sim(mwimp,E,theta,apdg1,'0',
     &     wh,kind,type,istat)
      else
         istat=ibset(istat,3)
         dsseyield_sim_ls=0.d0
      endif
         
      return
      end

