***********************************************************************
*** This routine calculates the capture rates in the Earth.
*** It does the same thing as dssenu_capearthnum (i.e. dssenu_capearthnumi),
*** except that it uses tabulted versions of the results intead
*** of performing a numerical integration every time. 
*** Inputs: mx - WIMP mass in GeV
***         rho - local halo density in GeV/cm^3      
***         sigsi - spin-independent capture rate in cm^2
***         type - type of distribution (same as in dssenu_capearthnum)
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      real*8 function dssenu_capearthtab(mx,rho,sigsi)

      implicit none
      real*8 mx,rho,sigsi,dssenu_ctabget

      dssenu_capearthtab=dssenu_ctabget('ea',1,mx)*(sigsi/1d-40)
     &  *(rho/0.3d0)

      return
      end
