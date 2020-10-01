#include "PILOT.inc"
#ifdef NORANLUX_X
      SUBROUTINE RANFGT(SEED)
C
C          Get seed for RANF() in real or double precision SEED.
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
      CALL RANGET(SEED)
#elif defined(CRAY_X)
      INTEGER ISEED,RANGET,IDUMMY
      ISEED=RANGET(IDUMMY)
      SEED=ISEED
#endif
      RETURN
      END
#endif
