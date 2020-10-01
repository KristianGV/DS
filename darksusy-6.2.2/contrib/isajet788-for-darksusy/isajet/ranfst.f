#include "PILOT.inc"
#ifdef NORANLUX_X
      SUBROUTINE RANFST(SEEDg)
C
C          Set seed for RANF() from real or double precision SEEDg
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#ifdef SINGLE_X
      REAL SEEDg
#elif defined(DOUBLE_X)
      DOUBLE PRECISION SEEDg
#endif
#ifdef RANFCALL_X
      CALL RANSET(SEEDg)
#elif defined(CRAY_X)
      INTEGER ISEEDg
      ISEEDg=SEEDg
      CALL RANSET(ISEEDg)
#endif
      RETURN
      END
#endif
