#include "PILOT.inc"
      LOGICAL FUNCTION LOGMGM(IERR)
C
C       Set and check limits for multijet mass
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
      LOGMGM=.TRUE.
C
      IF(EHMGMN.LT.UNDEF.OR.EHMGMX.LT.UNDEF) THEN
        LOGMGM=.FALSE.
      ENDIF
C
      RETURN
      END
