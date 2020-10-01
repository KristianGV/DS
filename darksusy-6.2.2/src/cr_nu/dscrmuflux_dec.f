      real*8 function dscrmuflux_dec(egev,diff,jpsi_dec,xi,istat)
**********************************************************************
***   function dscrmuflux_dec gives the flux of neutrino-induced mouns
***   from WIMP annihilation in the halo, in the limit of zero relative 
***   velocity. The function returns the flux as measured by a detector 
***   placed in ice.
***   
***   input: egev     - muon threshold energy [in GeV]
***          diff     - dictates whether differential source term at egev (diff=1) 
***                     or integrated source term above egev (diff=0) is returned
***          jpsi_dec - value of line-of sight integration over a solid angle
***                     (\int d\Omega\int ds rho [kpc sr GeV cm^-3],
***                     obtained e.g. with a call to dsomlosisph)
***          xi       - factor by which the DM density should be rescaled 
***                     (obtained e.g. as ratio of calculated to measured relic density)
***   
***   
***   unit of return value: km^-2 yr^-1 / GeV^diff
***   
***   author: Torsten.Bringmann@fys.uio.no
***   date: 2016-02-09
**********************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 egev,jpsi_dec,xi,dscrsource
      integer diff,istat

c...FIXME: oscillations still need to be added

      dscrmuflux_dec=dscrsource(egev,diff,130073,1,0d0,istat)*jpsi_dec*xi
     &          *kpc*1.d21/(4.d0*pi) ! standard factor of 4 pi goes here
     &          * year * 1.d10   ! convert from cm^-2 s^-1

      return
      end


