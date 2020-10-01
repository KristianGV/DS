#include "PILOT.inc"
#ifdef NORANLUX_X
      SUBROUTINE RANFST(SEED)
C
C          Set seed for RANF() from real or double precision SEED
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#ifdef SINGLE_X
      REAL SEED
#elif defined(DOUBLE_X)
      DOUBLE PRECISION SEED
#endif
#ifdef RANFCALL_X
      CALL RANSET(SEED)
#elif defined(CRAY_X)
      INTEGER ISEED
      ISEED=SEED
      CALL RANSET(ISEED)
#endif
      RETURN
      END
#endif
