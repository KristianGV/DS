      real*8 function dscrnuflux_v0ann(egev,diff,jpsi,xi,istat)
**********************************************************************
***   function dscrnuflux_v0ann gives the flux of muon neutrinos from WIMP
***   annihilation in the halo, in the limit of zero relative velocity.
***   
***   input: egev   - neutrino energy [in GeV]
***          diff   - dictates whether differential source term at egev (diff=1) 
***                   or integrated source term above egev (diff=0) is returned
***          jpsi   - value of line-of sight integration over a solid angle
***                    (\int d\Omega\int ds rho^2 [kpc sr GeV^2 cm^-6],
***                     obtained e.g. with a call to dsomlosisph)
***           xi    - factor by which the DM density should be rescaled 
***                   (obtained e.g. as ratio of calculated to measured relic density)
***   
***   
***   unit of return value: cm^-2 s^-1 / GeV^diff
***   
***   author: Torsten.Bringmann@fys.uio.no 
***   date: 2016-02-09
***   mod: 2019-11-15 (added oscillations)
**********************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 egev,jpsi,xi,dscrsource
      integer diff,istat

      real*8 dscroscillation_simp
      real*8 crsrc(3),crosc

      crsrc(1) = dscrsource(egev,diff,12,2,0d0,istat) ! nue at source
      crsrc(2) = dscrsource(egev,diff,14,2,0d0,istat) ! numu at source
      crsrc(3) = dscrsource(egev,diff,16,1,0d0,istat) ! nutau at source

      crosc = dscroscillation_simp(2,crsrc) ! numu at earth

      dscrnuflux_v0ann = crosc*jpsi*xi**2*kpc*1.d21/(4.d0*pi) 

      return
      end


