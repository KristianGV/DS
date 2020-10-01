************************************************************************
***  spherically symmetric version of the astrophysical source function
***  for dark matter pair annihilation or decays, setting the spatial     
***  part of the local injection yield. this function just links to the 
***  replaceable dmd driver dsdmddriver
***
***  NOTE: for a smooth dark matter distribution dsdmddriver must set:
***    dsdmassph -> dsdmdsph**2  for dark matter pair annihilations 
***    dsdmassph -> dsdmdsph     for dark matter decays 
***  while including substructures may result in more involved forms 
***
***  input: r - spherical galactocentric radius [kpc]
***         power - selecting annihilations or decays      
***
***  output:  [GeV^2 cm^-6] for pair annihilations (power=2)
***           [GeV cm^-3] for decays (power=1)
************************************************************************
      real*8 function dsdmassph(r,power)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdvcom.h'
      real*8 r
      integer power
      real*8 out
      ksoupow=power
      dvsph(idvrsph)=r
      call dsdmddriver(idmdsousph,ntysph,dvsph,ndummy,chdummy,labdummy,
     &  out)
      dsdmassph=out
      end
