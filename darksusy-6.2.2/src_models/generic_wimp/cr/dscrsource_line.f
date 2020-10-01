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
*** date  : 2015-06-21
*****************************************************************************

      real*8 function dscrsource_line(pdg,n,power,v,egev,widthline,pdg2,istat)
      implicit none      
      include 'dsgeneric_wimp.h'

c------------------------ functions ------------------------------------

       
c------------------------ variables ------------------------------------

      real*8 egev, v, widthline
      integer pdg, pdg2, n, power, istat

      real*8 result, dsmwimp, dsmass

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dscrsource_line')

      istat=0
      pdg2=0
      dscrsource_line=0d0
      egev=0d0
      widthline=0d0
      result=0d0
      

      if (n.gt.1.or.n.lt.0) then         ! only main annihilation channel is identified as line for generic WIMP
        istat=ibset(istat,1)
        return 
      endif

      if (power.ne.2) return ! only annihilation

c... NB: currently only identical particle masses supported for generic WIMP
      if (svch.eq.pdg.and.
     &   (pdg.eq.11.or.pdg.eq.12.or.pdg.eq.14.or.pdg.eq.16)) then
        if (n.eq.0) n=1
        pdg2=-pdg 
        egev=dsmwimp()*(1.-dsmass(pdg)**2/4./dsmwimp()**2)                    
      elseif (svch.eq.pdg.and.pdg.eq.22) then
        if (n.eq.0) n=1
        pdg2=pdg
        egev=dsmwimp()
      else
        return
      endif


c... this is the standard way of representing the "particle physics factor" 
c... for WIMP annihilation
      result= sva/dsmwimp()**power        

c... for identical particles, the monochromatic yield must be doubled
      if (pdg.eq.pdg2) result=result*2.

c... if non-self-conjugate, divide by an extra factor of 2,
c... since rho refers to the total density of particles + antiparticles
      if (selfconj.eq.2) result = result/2.d0
      

      dscrsource_line=result

      return
      end


