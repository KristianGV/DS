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
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2014-11-13
*****************************************************************************

      real*8 function dscrsource_line(pdg,n,power,v,egev,widthline,pdg2,istat)
      implicit none      
      include 'dsmssm.h'
      include 'dsanyieldmodelcom.h'
 
c------------------------ functions ------------------------------------

      real*8 dssigmav0, dsmwimp
      
c------------------------ variables ------------------------------------

      real*8 egev, v, widthline, result
      integer pdg, pdg2, n, nn, ch, power, istat


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dscrsource_line')

      istat=0
      dscrsource_line=0d0
      egev=0d0
      widthline=0d0

      if (n.lt.0) return


c... for WIMPs, only annihilation contributes 
c... (but note that for multi-component DM, several values of power can lead
c...  to a non-zero return value)
      if (power.ne.2) return
      
      
c      if (v.eq.0d0) then   ! v!=0 is not yet set up in MSSM and thus 
                            ! simply mapped to the v=0 case
      
c... identify correct line
          nn=0
          ch=0
          do while (ch.lt.numyieldch_line.and.(nn.lt.n.or.n.eq.0))
            ch=ch+1
            if (yieldchannels_line(ch,1).eq.pdg) then
              nn=nn+1
              if (nn.eq.1) pdg2=yieldchannels_line(ch,2) ! keep track of first line
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
            egev=dsmwimp()*(1.-mass(kz)**2/4./dsmwimp()**2)
            widthline=width(kz)
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

          result = dssigmav0(pdg,pdg2)

c      endif

           
c... this is the standard way of representing the "particle physics factor" 
c... for WIMP annihilation
      result= result/dsmwimp()**power        



c... for identical particles, the monochromatic yield must be doubled
      if (pdg.eq.pdg2) result=result*2.

c... the symmetry factor depends on the WIMP type -- here for Majorana DM
      result = result/2.d0
      

      dscrsource_line=result
      
      return

      end


