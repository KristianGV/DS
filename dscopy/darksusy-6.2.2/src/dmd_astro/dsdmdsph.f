************************************************************************
***  spherically symmetric version of the dark matter density profile
***  this function just links to the replaceable dmd driver dsdmddriver
***
***  input: r - spherical galactocentric radius [kpc]
***
***  output:  [GeV cm^-3]
************************************************************************
      real*8 function dsdmdsph(r)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdvcom.h'
      real*8 r
      real*8 out
      dvsph(idvrsph)=r
      call dsdmddriver(idmddensph,ntysph,dvsph,ndummy,chdummy,labdummy,
     &  out)
      dsdmdsph=out
      end
