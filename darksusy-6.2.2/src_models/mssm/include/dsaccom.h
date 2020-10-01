*         -*- mode: fortran -*-
      character*8 aclabel
      real*8 gzinv,bsg,gm2muon,delta0,bsmumu,bsmumu_untag,
     &  brbtaunu,brbdtaunu,rmu23,amusi,bsgsi
      common /dsaccom/gzinv,bsg,gm2muon,delta0,bsmumu,bsmumu_untag,
     &  brbtaunu,brbdtaunu,rmu23,amusi,bsgsi,
     &  aclabel
      save /dsaccom/
