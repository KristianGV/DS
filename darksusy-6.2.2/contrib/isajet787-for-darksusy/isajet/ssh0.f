#include "PILOT.inc"
      DOUBLE PRECISION FUNCTION SSH0(M1,M2)
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "ssinf.inc"
      DOUBLE PRECISION SSA0,SSB00,M1SQ,M2SQ
      COMPLEX*16 SSB0
      REAL M1,M2
      M1SQ=M1*M1
      M2SQ=M2*M2
C      SSH0=4.D0*(((SSA0(M1)+SSA0(M2))/2.D0+(M1SQ+M2SQ)
C     $*SSB00(M1,M2)+M1SQ+M2SQ)/6.D0)
C     $+(-M1SQ-M2SQ)*SSB00(M1,M2)-SSA0(M1)-SSA0(M2)
      SSH0=4.D0*(((SSA0(M1)+SSA0(M2))/2.D0+(M1SQ+M2SQ)
     $*SSB0(0.,M1,M2)+M1SQ+M2SQ)/6.D0)
     $+(-M1SQ-M2SQ)*SSB0(0.,M1,M2)-SSA0(M1)-SSA0(M2)
      RETURN
      END