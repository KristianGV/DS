*         -*- mode: fortran -*-
c
c variables needed in l.o.s.i. for the spherically symmetric case
c
      real*8 maxradius,thetav(10)
      integer npsf,sphilosi
      common/dslosisphs1com/maxradius,thetav,npsf,sphilosi
c
      real*8 npairs_norm,sinobjdist2
      integer locilosi,sphilosii
      common/dslosisphs2com/npairs_norm,sinobjdist2,locilosi,sphilosii
c
      real*8 taucut,taureschoweq1,taumax
      common/dslosisphs3com/taucut,taureschoweq1,taumax
c
      real*8 objdistance,radinnertr,radoutertr
      common/dsdmscom/objdistance,radinnertr,radoutertr
c
c extra variables needed in l.o.s.i. for the axisymmetric case
c
      real*8 angleb,anglel,xsign
      integer axihow,axihowi
      common/dslosiaxicom/angleb,anglel,xsign,axihow,axihowi

      save /dslosisphs1com/,/dslosisphs2com/,/dslosisphs3com/,
     & /dslosiaxicom/

