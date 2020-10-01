      real*8 function dscrgaflux_dec(egev,diff,jpsi_dec,xi,istat)
**********************************************************************
***   function dscrgaflux_dec gives the flux of gamma-rays from WIMP
***   annihilation in the halo, in the limit of zero relative velocity.
***   
***   input: egev     - gamma-ray energy [in GeV]
***          diff     - dictates whether differential source term at egev (diff=1) 
***                     or integrated source term above egev (diff=0) is returned
***          jpsi_dec - value of line-of sight integration over a solid angle
***                     (\int d\Omega\int ds rho [kpc sr GeV cm^-3],
***                     obtained e.g. with a call to dsomlosisph)
***          xi       - factor by which the DM density should be rescaled 
***                     (obtained e.g. as ratio of calculated to measured relic density)
***   
***   
***   unit of return value: cm^-2 s^-1/ GeV^diff
***   
***   author: Torsten.Bringmann@fys.uio.no 
***   date: 2016-02-09
**********************************************************************
      implicit none
      include 'dsmpconst.h'
      include 'dsdmdcom.h'  ! psdecay
      real*8 egev,jpsi_dec,xi,dscrsource
      integer diff,istat
      
      dscrgaflux_dec=dscrsource(egev,diff,22,psdecay,0d0,istat)
     &   *jpsi_dec*xi*kpc*1.d21/(4.d0*pi) ! standard factor of 4 pi goes here

      return
      end


