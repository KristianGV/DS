      real*8 function dscrgaflux_line_dec(n,egev,widthline,jpsi_dec,xi,istat)
**********************************************************************
***   function dscrgaflux_line_dec gives the flux of *monoenergetic* 
***   gamma-rays from WIMP annihilation in the halo, in the limit of 
***   zero relative velocity.
***   
***   input:   n - returns the nth monochromatic gamma-ray signal 
***                implemented for the particle physics mode that is linked to
***                [if called with n=0, the first line will be returned,
***                 and n be set to the total number of existing lines]
***
***   output: egev      - gamma-ray energy [in GeV]
***           widthline - signal width
***           jpsi_dec  - value of line-of sight integration over a solid angle
***                       (\int d\Omega\int ds rho [kpc sr GeV cm^-3],
***                       obtained e.g. with a call to dsomlosisph)
***           xi        - factor by which the DM density should be rescaled 
***                       (obtained e.g. as ratio of calculated to measured 
***                       relic density)
***           istat     - equals 0 if there are no errors, bit 1 is set if line
***                       n does not exist, higher non-zero bits specify model-
***                       specific errors
***   
***   
***   unit of return value: cm^-2 s^-1
***   
***   author: Torsten.Bringmann@fys.uio.no 
***   date: 2016-02-09
**********************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 egev,widthline,jpsi_dec,xi,dscrsource_line
      integer n,istat, pdg2
      
      dscrgaflux_line_dec=
     &          dscrsource_line(22,n,1,0.d0,egev,widthline,pdg2,istat)
     &          *jpsi_dec*xi**kpc*1.d21/(4.d0*pi) ! standard factor of 4 pi goes here
            
      return
      end
