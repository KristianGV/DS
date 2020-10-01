*****************************************************************************
***   function dscrsource returns the source term for DM-induced cosmic rays 
***   (including neutral species), assuming that the corresponding flux scales
***   like a power of the DM density rho_DM. Note that monochromatic 
***   contributions are not included here (see dscrsource_line for this).
***                                                                         
***   type : interface                                                      
***                                                                         
***   input: power  - determines scaling as flux ~ (rho_DM)^power 
***          pdg    - PDG code of CR species
***          egev   - CR energy [in GeV]
***          diff   - dictates whether differential source term at egev (diff=1) 
***                   or integrated source term above egev (diff=0) is returned
***          v      - relative velocity of annihilating DM particles
***                   in units of c
***                   [only for power=2, otherwise ignored]
***
***   output: istat - will be zero in case of no errors.
***
***   unit of return value: #particles / (cm^3 s) / (GeV / cm^3)^power
***         -> for power=2, this is multiplied by 1/c
***         -> for diff=1, this is multiplied by 1/GeV
***
*** author: Torsten.Bringmann.fys.uio.no
*** date  : 2015-06-11
*** modified TB 2019-12-13: checked FSR contributions
*****************************************************************************

      real*8 function dscrsource(egev,diff,pdg,power,v,istat)
      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsmpconst.h' 
c      include 'dsanyieldcom.h' ! to have access to dbflxk

c------------------------ functions ------------------------------------

      real*8 dssigmav0tot, dsmwimp, dsmass, dsanyield_sim
      
c------------------------ variables ------------------------------------

      real*8 egev, v, yield, ml, mx, x
      integer diff, pdg, power, istat


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dscrsource')

      istat=0
      dscrsource=0d0

      if (power.ne.2) return

      mx = dsmwimp()

c... simulated yields
      yield = dsanyield_sim(dsmwimp(),egev,svch,0,pdg,diff,istat)
c... FSR *is* actually included...
c      if (pdg.eq.22.and.(svch.eq.11.or.svch.eq.13.or.svch.eq.15)) then ! add FSR photons
c        ml = dsmass(svch)
c        x = egev/mx
c        if (x.lt.0.999.and.diff.eq.1) yield = yield + 
c     &       (alphem/pi*(1.+(1.-x)**2)/x*log(4.*(1.-x)*mx**2/ml**2))/mx
c      endif


c... this is the standard way of representing the "particle physics factor" 
c... for self-conjugate WIMP annihilation
      dscrsource = dssigmav0tot()*yield/mx**2/2.d0

c... if non-self-conjugate, divide by an extra factor of 2,
c... since rho refers to the total density of particles + antiparticles
      if (selfconj.eq.2) dscrsource=0.5d0*dscrsource

      return
      end












