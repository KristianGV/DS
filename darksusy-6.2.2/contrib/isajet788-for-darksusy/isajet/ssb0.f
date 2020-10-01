#include "PILOT.inc"
      COMPLEX*16 FUNCTION SSB0(QSQ,M1,M2)
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "ssinf.inc"
      COMPLEX*16 SSF0
      REAL QSQ,M1,M2
      SSB0=XLAM*(1.D0,0.D0)-SSF0(QSQ,M1,M2)
      RETURN
      END
