************************************************************************
***  cartesian coordinates version of the astrophysical source function 
***  for dark matter pair annihilation or decays, setting the spatial     
***  part of the local injection yield. this function just links to the  
***  replaceable dmd driver dsdmddriver
***
***  NOTE: for a smooth dark matter distribution dsdmddriver must set:
***    dsdmastriax -> dsdmdtriax**2  for dark matter pair annihilations 
***    dsdmastriax -> dsdmdtriax     for dark matter decays 
***  while including substructures may result in more involved forms 
***
***  input: x,y,z - carthesian coordinates [kpc]
***         power - selecting annihilations or decays      
***
***  output:  [GeV^2 cm^-6] for pair annihilations (power=2)
***           [GeV cm^-3] for decays (power=1)
************************************************************************
      real*8 function dsdmastri(x,y,z,power)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdvcom.h'
      real*8 x,y,z
      integer power
      real*8 out
      ksoupow=power
      dvtri(idvx)=x
      dvtri(idvy)=y
      dvtri(idvz)=z
      call dsdmddriver(idmdsoutri,ntytri,dvtri,ndummy,chdummy,labdummy,
     &  out)
      dsdmastri=out
      end
