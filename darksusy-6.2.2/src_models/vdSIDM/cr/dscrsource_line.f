*****************************************************************************
***   In analogy to dscrsource, subroutine dscrsource_line returns 
***   *monochromatic* cosmic ray contributions to the source term, assuming 
***   that the corresponding flux scales like a power of the DM density rho_DM.
***                                                                         
***  type : interface                                                       
***                                                                         
***   input: power  - determines scaling as flux ~ (rho_DM)^power 
***          pdg    - PDG code of monochromatic CR species
***          n      - returns the nth monochromatic signal for this pdg code
***                   [if called with n=0, the first line will be returned,
***                    and n be set to the total number of existing lines]
***          v      - relative velocity of annihilating DM particles
***                   in units of c
***                   [only for power=2, otherwise ignored]
***
***   output: egev      - CR energy of particle pdg [in GeV]
***           widthline - signal width
***           pdg2      - pdgcode of associated 2nd final state particle 
***           istat     - equals 0 if there are no errors, bit 1 is set if line
***                       n does not exist, higher non-zero bits specify model-
***                       specific errors
***
***   unit of return value: #particles / (cm^3 s) / (GeV / cm^3)^power
***         -> for power=2, this is multiplied by 1/c
***
*** author: Torsten.Bringmann.fys.uio.no
*** date  : 2018-02-18
*****************************************************************************

      real*8 function dscrsource_line(pdg,n,power,v,egev,widthline,pdg2,istat)
      implicit none      

c------------------------ variables ------------------------------------

      real*8 egev, v, widthline
      integer pdg, pdg2, n, power, istat

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dscrsource_line')

c... currently, the coupling of dark sector to SM particles is assumed to be
c... negligible in the vdSIDM module, and indirect detection signals are therefore
c... set to zero. See the generic_WIMP module for a simple way of adding 
c... non-zero CR source terms.

      istat=0
      pdg2=0
      dscrsource_line=0d0
      egev=0d0
      widthline=0d0
 
      return
      end


