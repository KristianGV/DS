*         -*- mode: fortran -*-
c
c type of dependent variables/coordinates used in a given function:
c      
      integer itydvsph  ! type: spherically symmetric -> r
      integer itydvaxi  ! type: axisymmetric -> R,z
      integer itydvhom  ! type: homoeoid -> m,q 
      integer itydvtri  ! type: triaxial -> x,y,z
      parameter(itydvsph=1,itydvaxi=2,itydvhom=3,itydvtri=4)
c      
c for each type specify here once and for all the (arbitrary) entry number
c in the dependent variable/coordinate vector
c
      integer idvrsph,ntysph
      parameter(idvrsph=1,ntysph=1)
      integer idvraxi,idvzaxi,ntyaxi
      parameter(idvraxi=1,idvzaxi=2,ntyaxi=2)
      integer idvmhom,idvqhom,ntyhom
      parameter(idvmhom=1,idvqhom=2,ntyhom=2)
      integer idvx,idvy,idvz,ntytri
      parameter(idvx=1,idvy=2,idvz=3,ntytri=3)
c
c this is the dependent variable/coordinate vector (it is an internal
c changing quantity, so it does not need to be stored in common blocks):
c (TB: but not doing so results in compiler warnings when it is included
c      but not used. So I put it in a common block, but without saving it...)
c      
      real*8 dvsph(ntysph),dvaxi(ntyaxi),dvhom(ntyhom),dvtri(ntytri)
      common /dsdvdcomvars/ dvsph,dvaxi,dvhom,dvtri

