#include "PILOT.inc"
      DOUBLE PRECISION FUNCTION SSA0(M)
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "ssinf.inc"
      REAL M
      DOUBLE PRECISION MSQ
      MSQ=M**2
      IF(M.NE.0.) THEN
        SSA0=MSQ*(1.D0-DLOG(MSQ)+XLAM)
      ELSE
        SSA0=0.D0
      ENDIF
      RETURN
      END
