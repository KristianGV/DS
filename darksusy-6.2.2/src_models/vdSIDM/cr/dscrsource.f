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
*** date  : 2018-02-18
*****************************************************************************

      real*8 function dscrsource(egev,diff,pdg,power,v,istat)
      implicit none
      include 'dsvdSIDM.h'
      include 'dsanyieldcom.h' ! to have access to dbflxk
      
c------------------------ variables ------------------------------------

      real*8 egev, v
      integer diff, pdg, power, istat


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dscrsource')

c... currently, the coupling of dark sector to SM particles is assumed to be
c... negligible in the vdSIDM module, and indirect detection signals are therefore
c... set to zero. See the generic_WIMP module for a simple way of adding 
c... non-zero CR source terms.

      istat=0
      dscrsource=0d0


      return
      end



















