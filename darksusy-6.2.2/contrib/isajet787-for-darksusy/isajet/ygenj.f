#include "PILOT.inc"
      LOGICAL FUNCTION YGENJ(I)
C
C            GENERATE Y FOR TWOJET
C
#include "itapes.inc"
#include "jetlim.inc"
#include "jetpar.inc"
#include "primar.inc"
#include "ptpar.inc"
#include "totals.inc"
      ACOSH(X)=ALOG(X+SQRT(X**2-1.0))
      YGENJ=.TRUE.
      YMAX=ACOSH(HALFE/PT(I))
      YMIN=-YMAX
      IF(YMAX.LT.YJMIN(I).OR.YMIN.GT.YJMAX(I)) GOTO 10
      YJ(I)=YJMIN(I)+(YJMAX(I)-YJMIN(I))*RANF()
      IF(YJ(I).LT.YMIN.OR.YJ(I).GT.YMAX) GOTO 10
      TH(I)=2.*ATAN(EXP(-YJ(I)))
      CTH(I)=COS(TH(I))
      STH(I)=SIN(TH(I))
      WT=WT*(YJMAX(I)-YJMIN(I))
      RETURN
   10 YGENJ=.FALSE.
      RETURN
      END
