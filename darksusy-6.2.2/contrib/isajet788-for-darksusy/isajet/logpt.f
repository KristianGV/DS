#include "PILOT.inc"
      LOGICAL FUNCTION LOGPT(IERR)
C
C       SET AND CHECK LIMITS FOR JET PT
C
#include "itapes.inc"
#include "jetlim.inc"
#include "primar.inc"
#include "jetpar.inc"
#include "const.inc"
#include "dylim.inc"
#include "keys.inc"
#include "q1q2.inc"
      DATA PTLOW/1.0/
      DATA UNDEF/-.9E9/
C
      LOGPT=.TRUE.
      DO 20 I=1,NJET
      FIXPT(I)=.FALSE.
C
      IF(PTMIN(I).LT.UNDEF.AND.PTMAX(I).LT.UNDEF) THEN
        PTMAX(I)=PMAX(I)
        PTMIN(I)=PTLOW
        IF(KEYS(3).AND.I.EQ.3.AND.QTMIN.GT.0.) PTMIN(I)=QTMIN
        IF(PMIN(I).LT.PTMIN(I)) PMIN(I)=PTMIN(I)
      ELSEIF(PTMAX(I).LT.UNDEF) THEN
        FIXPT(I)=.TRUE.
        PTMAX(I)=PTMIN(I)
      ELSEIF(PTMIN(I).LT.UNDEF) THEN
        PTMIN(I)=PTLOW
        IF(KEYS(3).AND.I.EQ.3.AND.QTMIN.GT.0.) PTMIN(I)=QTMIN
      ENDIF
C
      IF(FIXPT(I)) PTMAX(I)=PTMIN(I)
      IF(FIXPT(I)) PT(I)=PTMIN(I)
      IF(PTMAX(I).GT.PMAX(I)) PTMAX(I)=PMAX(I)
      IF(PMIN(I).LT.PTMIN(I)) PMIN(I)=PTMIN(I)
C
   20 CONTINUE
C
      RETURN
      END
