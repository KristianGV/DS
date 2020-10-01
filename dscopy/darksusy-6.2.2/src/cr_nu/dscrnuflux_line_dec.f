      real*8 function dscrnuflux_line_dec(n,egev,widthline,jpsi_dec,xi,istat)
**********************************************************************
***   function dscrnuflux_line_dec gives the flux of *monoenergetic* 
***   muon neutrinos from WIMP annihilation in the halo, in the limit of 
***   zero relative velocity.
***   
***   input:   n - returns the nth monochromatic gamma-ray signal 
***                implemented for the particle physics mode that is linked to
***
***   output: egev      - neutrino energy [in GeV]
***           widthline - signal width
***           jpsi_dec  - value of line-of sight integration over a solid angle
***                       (\int d\Omega\int ds rho [kpc sr GeV cm^-3],
***                        obtained e.g. with a call to dsomlosisph)
***           xi        - factor by which the DM density should be rescaled 
***                       (obtained e.g. as ratio of calculated to measured 
***                        relic density)
***           istat     - equals 0 if there are no errors, bit 1 is set if line
***                       n does not exist, higher non-zero bits specify model-
***                       specific errors
***   
***   
***   unit of return value: cm^-2 s^-1
***   
***   author: Torsten.Bringmann@fys.uio.no 
***   date: 2016-02-09
***   mod: 2019-11-15 (added oscillations)
**********************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 egev,widthline,jpsi_dec,xi,dscrsource_line
      integer n,istat, pdg2

      real*8 dscroscillation_simp
      real*8 crsrc(3),crosc

      crsrc(1) = dscrsource_line(12,n,1,0.d0,egev,widthline,pdg2,istat) ! nue at source
      crsrc(2) = dscrsource_line(14,n,1,0.d0,egev,widthline,pdg2,istat) ! numu at source
      crsrc(3) = dscrsource_line(16,n,1,0.d0,egev,widthline,pdg2,istat) ! nutau at source
      
      crosc = dscroscillation_simp(2,crsrc) ! numu at earth
      
      dscrnuflux_line_dec = crosc*jpsi_dec*xi*kpc*1.d21/(4.d0*pi) 
      
      return
      end


