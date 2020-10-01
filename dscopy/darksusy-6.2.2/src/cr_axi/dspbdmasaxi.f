      real*8 function dspbdmasaxi(R,z,power)
************************************************************************
*** spatially dependent part of the antiproton and antideuteron source 
*** function. an axisymmetric profile is required here.
***
***  input: galactic coordinates R,z in kpc
***  power = integer selecting annihilations (=2) or decays (=1)      
***  output:  [GeV^2 cm^-6] for dark matter pair annihilations
***           [GeV cm^-3] for dark matter decays
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 R,z
      integer power
ccc
      real*8 R0,z0,rcheck,dsdmasaxi
ccc
ccc cut out inner spherical region:
ccc
      rcheck=dsqrt(R**2+z**2)
      if(rcheck.lt.diffrcpb) then
        z0=z
        R0=dabs(diffrcpb**2-z**2)
      else
        z0=z
        R0=R
      endif
      dspbdmasaxi=dsdmasaxi(R0,z0,power)
      return
      end



