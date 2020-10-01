      real*8 function dscrgaflux_v0ann(egev,diff,jpsi,xi,istat)
**********************************************************************
***   function dscrgaflux_v0ann gives the flux of gamma-rays from WIMP
***   annihilation in the halo, in the limit of zero relative velocity.
***
***   type : commonly used
***   desc : Flux of gamma-rays from annihilation in the halo, in the
***   desc : limit of zero relative velocity
***   
***   input: egev   - gamma-ray energy [in GeV]
***          diff   - dictates whether differential source term at egev (diff=1) 
***                   or integrated source term above egev (diff=0) is returned
***          jpsi   - value of line-of sight integration over a solid angle
***                    (\int d\Omega\int ds rho^2 [kpc sr GeV^2 cm^-6],
***                     obtained e.g. with a call to dsomlosisph)
***          xi     - factor by which the DM density should be rescaled 
***                   (obtained e.g. as ratio of calculated to measured relic density)
***   
***   
***   unit of return value: cm^-2 s^-1/ GeV^diff
***   
***   author: Torsten.Bringmann@fys.uio.no 
***   date: 2016-02-09
**********************************************************************
      implicit none
      include 'dsmpconst.h'
      include 'dsdmdcom.h'  ! psannihi
      real*8 egev,jpsi,xi,dscrsource
      integer diff,istat
      
      dscrgaflux_v0ann=dscrsource(egev,diff,22,psannihi,0d0,istat)
     &  *jpsi*xi**2*kpc*1.d21/(4.d0*pi) ! standard factor of 4 pi goes here

      return
      end


