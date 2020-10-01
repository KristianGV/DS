#include "PILOT.inc"
#ifdef NORANLUX_X
      SUBROUTINE RANFGT(SEEDg)
C
C          Get seed for RANF() in real or double precision SEEDg.
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
      CALL RANGET(SEEDg)
#elif defined(CRAY_X)
      INTEGER ISEEDg,RANGET,IDUMMY
      ISEEDg=RANGET(IDUMMY)
      SEEDg=ISEEDg
#endif
      RETURN
      END
#endif
