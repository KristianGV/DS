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
*** date  : 2016-06-28
*** mod   : 2019-06-20 tb, added kinematics check
*****************************************************************************

      real*8 function dscrsource_line(pdg,n,power,v,egev,widthline,pdg2,istat)
      implicit none      
      include 'dssilveira_zee.h'

c------------------------ functions ------------------------------------

      real*8 dssigmavpartial, dsmwimp, dsmass
       
c------------------------ variables ------------------------------------

      real*8 egev, v, widthline
      integer pdg, pdg2, n, power, istat, ch, nn

      real*8 result

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dscrsource_line')

      istat=0
      pdg2=0
      dscrsource_line=0d0
      egev=0d0
      widthline=0d0
      result=0d0
      

      if (n.lt.0) then 
        istat=ibset(istat,1)
        return 
      endif

      if (power.ne.2) return ! allow only annihilation

c... identify correct line
          nn=0
          ch=0
          do while (ch.lt.numyieldch_line.and.(nn.lt.n.or.n.eq.0))
            ch=ch+1
            if (yieldchannels_line(ch,1).eq.pdg) then
              if (dsmass(pdg)+dsmass(yieldchannels_line(ch,2)).lt.dsmwimp()) then 
                               ! make sure the line is kinematically accessible
                nn=nn+1
                if (nn.eq.1) pdg2=yieldchannels_line(ch,2) ! keep track of first line
              endif
            endif
          enddo
          
          if (n.eq.0) then ! return first line, along with total number of lines
            if (nn.eq.0) return
            n = nn
          elseif (nn.eq.n.and.n.gt.0) then ! return line number n
            pdg2=yieldchannels_line(ch,2)
          else  ! line does not exist
            istat=ibset(istat,0)
            pdg2=0
            return 
          endif
                                
c... determine corresponding cross section, energy and width
          if (pdg.eq.22.and.pdg2.eq.22) then ! gamma gamma
            egev=dsmwimp()
            widthline=0d0
          elseif (pdg.eq.22.and.pdg2.eq.23) then ! gamma Z
            egev=dsmwimp()*(1.-dsmass(23)**2/4./dsmwimp()**2)
            widthline=width(15) ! internal code, see dsinit_module
          elseif (pdg.eq.-11.and.pdg2.eq.11) then  ! e+ e-
            egev=dsmwimp()
            widthline=0d0            
          elseif (pdg.eq.12.and.pdg2.eq.-12) then  ! nue nuebar
            egev=dsmwimp()
            widthline=0d0
         elseif (pdg.eq.14.and.pdg2.eq.-14) then ! numu numubar
            egev=dsmwimp()
            widthline=0d0
          elseif (pdg.eq.16.and.pdg2.eq.-16) then ! nutau nutaubar
            egev=dsmwimp()
            widthline=0d0
          else
            istat=ibset(istat,0)
            return
          endif

          do ch = 1, numannch2b 
           if ((pdg.eq.ann_2body(ch,1).and.pdg2.eq.ann_2body(ch,2)).or.
     &        (pdg.eq.ann_2body(ch,2).and.pdg2.eq.ann_2body(ch,1))) 
     &      result = dssigmavpartial(ch,0.d0)
          enddo 
           
c... this is the standard way of representing the "particle physics factor" 
c... for DM annihilation (cross section is already included)

      result= result/dsmwimp()**2/2.0      


c... for identical particles, the monochromatic yield must be doubled
      if (pdg.eq.pdg2) result=result*2.

      

      dscrsource_line=result


      return
      end


