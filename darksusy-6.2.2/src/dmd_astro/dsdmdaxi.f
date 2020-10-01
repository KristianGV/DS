************************************************************************
***  axisymmetric version of the dark matter density profile
***  this function just links to the replaceable dmd driver dsdmddriver
***
***  input: R,z - cylindrical galactic coordinates [kpc]
***
***  output:  [GeV cm^-3]
************************************************************************
      real*8 function dsdmdaxi(R,z)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdvcom.h'
      real*8 R,z
      real*8 out
      dvaxi(idvraxi)=R
      dvaxi(idvzaxi)=z
      call dsdmddriver(idmddenaxi,ntyaxi,dvaxi,ndummy,chdummy,labdummy,
     &  out)
      dsdmdaxi=out
      end
