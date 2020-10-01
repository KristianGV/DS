***********************************************************************
*** This routine calculates the capture rates in the Sun.
*** It does the same thing as dssenu_capsunnum (i.e. dssenu_capsunnumi),
*** except that it uses tabulted versions of the results intead
*** of performing a numerical integration every time. 
*** Inputs: mx - WIMP mass in GeV
***         rho - local halo density in GeV/cm^3      
***         sigsi - spin-independent capture rate in cm^2
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      real*8 function dssenu_capsuntab(mx,rho,sigsi,sigsd)

      implicit none
      real*8 mx,rho,sigsi,sigsd,dssenu_ctabget

      dssenu_capsuntab=dssenu_ctabget('su',1,mx)*(sigsi/1.d-40)
     &  *(rho/0.3d0)
     &  +dssenu_ctabget('su',2,mx)*(sigsd/1.d-40)*(rho/0.3d0)

      return
      end
