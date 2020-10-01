#include "PILOT.inc"
      LOGICAL FUNCTION LOGP(IERR)
C
C       SET AND CHECK LIMITS FOR JET MOMENTA
C
#include "itapes.inc"
#include "jetlim.inc"
#include "primar.inc"
#include "jetpar.inc"
#include "const.inc"
#include "dylim.inc"
#include "keys.inc"
#include "q1q2.inc"
      DATA PLOW/1.0/
      DATA UNDEF/-.9E9/
C
      LOGP=.TRUE.
      DO 10 I=1,NJET
      FIXP(I)=.FALSE.
      IF(PMIN(I).LT.UNDEF.AND.PMAX(I).LT.UNDEF)  PMAX(I)=HALFE
      IF(PMAX(I).GT.HALFE) PMAX(I)=HALFE
      IF(PMAX(I).LT.UNDEF) FIXP(I)=.TRUE.
      IF(PMIN(I).LT.UNDEF) PMIN(I)=PLOW
      IF(FIXP(I)) THEN
        PMAX(I)=PMIN(I)
        P(I)=PMIN(I)
      ENDIF
      IF(KEYS(3).AND.I.EQ.3.AND.QTMIN.GT.0) PMIN(I)=QTMIN
   10 CONTINUE
C
      RETURN
      END
