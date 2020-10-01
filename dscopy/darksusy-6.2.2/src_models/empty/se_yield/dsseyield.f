*****************************************************************************
*** Function dsseyield calculates the yield of neutrinos, leptons or
*** hadronic showers from WIMP annihilations in the Sun/Earth summing
*** over all annihilation channels.
***                                                                         
***  type : interface                                                       
***
***  desc : Yields from annihilation in the Sun/Earth
***                                                                         
*** Inputs:
***   emu   - energy of neutrino, lepton or hadronic shower (GeV)
***   theta - angle from the Sun / centre of the Earth (degrees)
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
***      GeV**-1 degree**-1.
***   istat - status flag (non-zero if something went wrong)
*** NOTE: Compared to previous versions of this routine, particles and
*** antiparticles are no longer summed, you thus need to call it for both
*** types and add them up.
***
*** Author: Torsten.Bringmann@fys.uio.no
*** Date: 11/06/2015
*****************************************************************************

      real*8 function dsseyield(e,theta,wh,kind,type,istat)
      implicit none

c------------------------ variables ------------------------------------

      real*8 e,theta
      integer istat,kind,type
      character*2 wh

c------------------------ functions ------------------------------------


c-----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsseyield')

      dsseyield=0.0d0

      return
      end




