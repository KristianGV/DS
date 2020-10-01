*****************************************************************************
***   function dscrsource returns the source term for DM-induced cosmic rays 
***   (including neutral species), assuming that the corresponding flux scales
***   like a power of the DM density rho_DM. Note that monochromatic 
***   contributions are not included here (see dscrsource_line for this).
***
***   type : INTERFACE
***
***   input: power  - determines scaling as flux ~ (rho_DM)^power 
***          pdg    - PDG code of CR species
***          egev   - CR energy [in GeV]
***                   Energy is per CR particle (i.e. not per nucleon)      
***          diff   - dictates whether differential source term at egev (diff=1) 
***                   or integrated source term above egev (diff=0) is returned
***          v      - relative velocity of annihilating DM particles
***                   in units of c
***                   [only for power=2, otherwise ignored]
***
***   output: istat - will be zero in case of no errors.
***
***   unit of return value: #particles / (cm^3 s) / (GeV / cm^3)^power
***         -> for diff=1, this is multiplied by 1/GeV
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2014-11-13
*****************************************************************************

      real*8 function dscrsource(egev,diff,pdg,power,v,istat)
      implicit none
      include 'dsio.h'
      include 'dsanyieldmodelcom.h'
      include 'dsanyieldcom.h' ! to get dbflxk
 
c------------------------ functions ------------------------------------

      real*8 dsanyield, dssigmav0tot, dsmwimp
      
c------------------------ variables ------------------------------------

      real*8 egev, v, result
      integer diff, pdg, power, istat, yieldch

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dscrsource')

      istat=0
      dscrsource=0d0


c... for WIMPs, only annihilation contributes 
c... (but note that for multi-component DM, several values of power can lead
c...  to a non-zero return value)
      if (power.ne.2) return
      
c... assign channels for internal use in the MSSM module 
      if (pdg.eq.-11) then  
        yieldch = 51   ! positrons
      elseif (pdg.eq.22) then
        yieldch = 52    ! gamma-rays
      elseif (pdg.eq.14) then
        yieldch = 53    ! muon neutrinos
      elseif (pdg.eq.-2212) then
        yieldch = 54           ! antiprotons
      elseif (pdg.eq.-1000010020) then
        yieldch = dbflxk  ! anti-deuterons, get code form common i src
c... FIXME: electron and tau neutrinos still missing here!

c... the CR particles are not created at the place of annihilation/decay,
c... but after "propagation"
      elseif (pdg.eq.130072) then
        yieldch = 72    ! muons from neutrinos, at creation 
      elseif (pdg.eq.130073) then
        yieldch = 73    ! muons from neutrinos, as seen by a detector in ice
                        ! (i.e. integrating channel 72 over the mean muon path)
      else
        istat = 1  ! FIXME: use error codes consistently
        if (prtlevel.gt.0) write(*,*) 'WARNING in dscrsource: pdg = ',pdg,'not supported.',
     &             'returning 0...'
        return
      endif
   
      if (diff.eq.1) then
        yieldch=yieldch + 100
      elseif (diff.ne.0) then
        istat = 1  
        write(*,*) 'ERROR in dscrsource: diff = ',diff,
     &    'not a possible choice!'
        return      
      endif
        
c... this is the standard way of representing the "particle physics factor" 
c... for WIMP annihilation

c      if (v.lt.1d-4) then   ! v!=0 is not yet set up in MSSM and thus 
                            ! simply mapped to the v=0 case
        result= dssigmav0tot()*dsanyield(egev,yieldch,istat)
     &    /dsmwimp()**power        
       
c      endif

c... the symmetry factor depends on the WIMP type -- here for Majorana DM
      result = result/2.d0
      

      dscrsource=result
      return

      end



















