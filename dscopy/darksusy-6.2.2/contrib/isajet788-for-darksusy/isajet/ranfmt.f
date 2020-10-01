#include "PILOT.inc"
#ifdef NORANLUX_X
      SUBROUTINE RANFMT
C
C          Get RANF seed and translate it to a character variable
C          to ensure exactly the same seed with a formatted read.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "seed.inc"
#ifdef SINGLE_X
      REAL SEED
#elif defined(DOUBLE_X)
      DOUBLE PRECISION SEED
#endif
      CALL RANFGT(SEED)
      WRITE(XSEED,'(E24.15)') SEED
      READ(XSEED,'(E24.15)') SEED
      CALL RANFST(SEED)
      RETURN
      END
#endif
