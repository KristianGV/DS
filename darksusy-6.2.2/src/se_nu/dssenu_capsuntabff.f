***********************************************************************
*** This routine calculates the capture rates in the Sun.
*** It does the same thing as dssenu_capsunnumff (i.e. dssenu_capsunnumffi),
*** except that it uses tabulted versions of the results intead
*** of performing a numerical integration every time. 
*** Inputs: mx - WIMP mass in GeV
***         rho - local halo density in GeV/cm^3      
***         gps,gns,gpa,gna - WIMP-nucleon four-fermion couplings
***           (obtained from e.g. dsddgpgn), unit: GeV^-4
*** Author: Joakim Edsjo
*** Date: 2015-06-12
***********************************************************************

      real*8 function dssenu_capsuntabff(mx,rho,gps,gns,gpa,gna)

      implicit none
      real*8 mx,rho,gps,gns,gpa,gna,dssenu_ctabffget,tmp

      tmp=dssenu_ctabffget('su',1,mx)*gps**2
     &   +dssenu_ctabffget('su',2,mx)*gns**2
     &   +dssenu_ctabffget('su',3,mx)*gps*gns
     &   +dssenu_ctabffget('su',4,mx)*gpa**2
     &   +dssenu_ctabffget('su',5,mx)*gna**2
     &   +dssenu_ctabffget('su',6,mx)*gpa*gna
      
      dssenu_capsuntabff=tmp*(rho/0.3d0)

      return
      end
