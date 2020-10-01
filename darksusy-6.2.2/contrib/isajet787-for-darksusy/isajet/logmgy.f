#include "PILOT.inc"
      LOGICAL FUNCTION LOGMGY(IERR)
C
C       Set and check limits for dijet masses.
C
C       Note we use the convention that not setting an upper limit
C       gives a fixed variable, even though that currently is not
C       implemented in N-jet phase space.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "itapes.inc"
#include "jetlim.inc"
#include "primar.inc"
#include "jetpar.inc"
#include "const.inc"
#include "dylim.inc"
#include "keys.inc"
#include "q1q2.inc"
#include "mglims.inc"
C
      REAL UNDEF
      INTEGER IERR
      DATA UNDEF/-.9E9/
C
      LOGMGY=.TRUE.
C
C          Attempt to fix YHMG
C
      IF(YHMGMN.GT.UNDEF.AND.YHMGMX.LT.UNDEF) THEN
        LOGMGY=.FALSE.
        RETURN
      ENDIF
C
C          No limits
C
      IF(EHMGMN.LT.0.OR.EHMGMX.LT.0) THEN
        LOGMGY=.FALSE.
        RETURN
      ENDIF
      IF(YHMGMN.LT.UNDEF) THEN
        YHMGMX=LOG(ECM/EHMGMN)
        YHMGMN=-YHMGMX
      ENDIF
C
      RETURN
      END
