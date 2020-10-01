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
*** date  : 2016-02-06   
*****************************************************************************

      real*8 function dscrsource(egev,diff,pdg,power,v,istat)
      implicit none
      include 'dsgeneric_decayingDM.h'

c------------------------ functions ------------------------------------

      real*8 dsGammatot, dsmwimp, dsanyield_sim_ls
      
c------------------------ variables ------------------------------------

      real*8 egev, v, yieldtot, mDM
      integer diff, pdg, power, istat, decistat, ch


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_decayingDM','dscrsource')

      istat=0
      dscrsource=0d0

      if (power.ne.1) return ! allow only decay

c... now calculate the yield by summing over all implemented channels

      yieldtot=0.d0

      mDM=dsmwimp()
      
      do ch=1,numdecch2b
         decistat=0
         if (decBR(ch).gt.0d0) then
            yieldtot= yieldtot + decBR(ch)*
     &           dsanyield_sim_ls(mDM/2.d0,egev,dec_2body(ch,1),
     &           dec_2body(ch,2),
     &           dec_2body(ch,3),dec_2body(ch,4),
     &           dec_2body(ch,5),dec_2body(ch,6),
     &           pdg,diff,decistat)
        
            if (btest(decistat,3)) then  ! channel not simulated!       
                                         ! nothing happen -> channel is treated as invisible
            endif
         endif    
      enddo


c... this is the standard way of representing the "particle physics factor" 
c... for annihilating DM

      dscrsource = dsGammatot()/dsmwimp()*yieldtot
       
      return
      end



















