************************************************************************
***  dark matter density profile in cartesian coordinates
***  this function just links to the replaceable dmd driver dsdmddriver
***
***  input: x,y,z - carthesian coordinates [kpc]
***
***  output:  [GeV cm^-3]
************************************************************************
      real*8 function dsdmdtri(x,y,z)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdvcom.h'
      real*8 x,y,z
      real*8 out
      dvtri(idvx)=x
      dvtri(idvy)=y
      dvtri(idvz)=z
      call dsdmddriver(idmddentri,ntytri,dvtri,ndummy,chdummy,labdummy,
     &  out)
      dsdmdtri=out
      end
