************************************************************************
***  axisymmetric version of the astrophysical source function for dark
***  matter pair annihilation or decays, setting the spatial part of the    
***  local injection yield. this function just links to the
***  corresponding dark matter density profile to the appropriate power
***
***  NOTE: for a smooth dark matter distribution dsdmddriver must set:
***    dsdmasaxi -> dsdmdaxi**2  for dark matter pair annihilations 
***    dsdmasaxi -> dsdmdaxi     for dark matter decays 
***  while including substructures may result in more involved forms 
***
***  input: R,z - cylindrical galactic coordinates [kpc]
***         power - selecting annihilations or decays      
***
***  output:  [GeV^2 cm^-6] for pair annihilations (power=2)
***           [GeV cm^-3] for decays (power=1)
************************************************************************
      real*8 function dsdmasaxi(R,z,power)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdvcom.h'
      real*8 R,z
      integer power
      real*8 out
      ksoupow=power
      dvaxi(idvraxi)=R
      dvaxi(idvzaxi)=z
      call dsdmddriver(idmdsouaxi,ntyaxi,dvaxi,ndummy,chdummy,labdummy,
     &  out)
      dsdmasaxi=out
      end
