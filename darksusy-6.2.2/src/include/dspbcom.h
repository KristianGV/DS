*         -*- mode: fortran -*-
      integer pbnzero,pbpropmodel
      real*8 pbr0,pbrh,pbhh,pbhg,pbkh,pbkg,pbrig0,pbcvel,pbnh,pbng,
     &  smod,pbdelta
      character*10 pbc
      common /pbpar/pbr0,pbrh,pbhh,pbhg,pbkh,pbkg,pbrig0,pbcvel,
     &     pbdelta,pbnh,pbng,pbnzero,pbpropmodel,smod,pbc
      save /pbpar/

c...Clumpy things
      real*8 rcl,zcl,thetacl,deltathetacl
      common/pbclpos/rcl,zcl,thetacl,deltathetacl
      integer nbesselk,nzerojk
      parameter(nbesselk=1000,nzerojk=5000)

c...Inner integration region for cuspy profiles
      real*8 pbrcy,pbzcy,pbrho2int
      common/pbcycom/pbrcy,pbzcy,pbrho2int

      save /pbclpos/,/pbcycom/

