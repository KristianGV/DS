*****************************************************************************
***   function dscrsource returns the source term for DM-induced cosmic rays 
***   (including neutral species), assuming that the corresponding flux scales
***   like a power of the DM density rho_DM. Note that monochromatic 
***   contributions are not included here (see dscrsource_line for this).
***                                                                         
***   type : interface                                                      
***
***   desc : Source term for dark matter induced cosmic rays
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
*** author: Torsten.Bringmann@fys.uio.no
*** date  : 2015-06-11   
*****************************************************************************

      real*8 function dscrsource(egev,diff,pdg,power,v,istat)
      implicit none
 
c------------------------ functions ------------------------------------

      
c------------------------ variables ------------------------------------

      real*8 egev, v
      integer diff, pdg, power, istat


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dscrsource')

      istat=0
      dscrsource=0d0
 
      return
      end



















