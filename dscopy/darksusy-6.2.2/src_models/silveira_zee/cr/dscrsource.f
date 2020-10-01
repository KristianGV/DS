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
*** date  : 2016-06-28 
*****************************************************************************

      real*8 function dscrsource(egev,diff,pdg,power,v,istat)
      implicit none
      include 'dssilveira_zee.h'
      include 'dsio.h'

c------------------------ functions ------------------------------------

      real*8 dssigmav0tot, dsmwimp, dssigmavpartial, dsmass
      real*8 dsanyield_sim_ls, dsanyield_sim
      
c------------------------ variables ------------------------------------

      real*8 egev, v, yieldtot, mDM, sv, eZ
      integer diff, pdg, power, istat, chstat, ch


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dscrsource')

      istat=0
      dscrsource=0.d0

      if (power.ne.2) return ! allow only annihilation

c... now calculate the yield by summing over all implemented channels

      yieldtot=0.d0
      mDM=dsmwimp()
      sv=dssigmav0tot()
      
      do ch=1,numannch2b
         chstat=0
         if (dssigmavpartial(ch,0.d0)/sv.gt.1.d-7) then
c...first check whether channel is simulated
            yieldtot= yieldtot + dssigmavpartial(ch,0.d0)/sv*
     &           dsanyield_sim_ls(mDM,egev,ann_2body(ch,1),
     &           ann_2body(ch,2),
     &           ann_2body(ch,3),ann_2body(ch,4),
     &           ann_2body(ch,5),ann_2body(ch,6),
     &           pdg,diff,chstat)

c...then include channels that are not simulated        
            if (btest(chstat,3)) then  
            
              if ((ann_2body(ch,1).eq.22.and.ann_2body(ch,2).eq.23).or.
     &           (ann_2body(ch,2).eq.22.and.ann_2body(ch,1).eq.23)) then  ! gamma Z
     
                 eZ=((2.0d0*mdm)**2+dsmass(23)**2)/(4.0d0*mdm)
                 eZ=max(eZ,dsmass(23)+0.001d0)
     
                 yieldtot=yieldtot+0.5*dssigmavpartial(ch,0.d0)*
     &                   dsanyield_sim(eZ,egev,23,'0',pdg,diff,chstat)    
              endif
             
              if (ann_2body(ch,1).eq.25.and.ann_2body(ch,2).eq.25) then  ! H H

c...TB FIXME: this is a very approximate fix until we have HH channels implemented!     
                 yieldtot=yieldtot+dssigmavpartial(ch,0.d0)*
     &                   dsanyield_sim(mdm,egev,23,'0',pdg,diff,chstat)
                 if (dssigmavpartial(ch,0.d0)/sv.gt.1.d-1) then
                   if (prtlevel.gt.1) then
                     write(*,*) 'WARNING in dscrsource: yield from Higgs-Higgs'
                     write(*,*) 'final states approximated as ZZ final state!'
                   endif
                 endif      
              endif
            
            endif
         endif    
      enddo


c... this is the standard way of representing the "particle physics factor" 
c... for annihilating (self-conjugate) DM

      dscrsource = sv/mDM**2*yieldtot/2.d0


      return
      end



















