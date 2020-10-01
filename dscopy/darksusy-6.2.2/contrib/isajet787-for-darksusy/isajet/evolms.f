#include "PILOT.inc"
      FUNCTION EVOLMS(J,FUDGE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-     Set evolution mass scale for parton J
C-
C-   Returned value  : maximum mass
C-
C-   Inputs  : 
C-     J    = index to PJSET array
C-     FUDGE= fudge factor
C-
C-   Created  16-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
      REAL    EVOLMS,FUDGE
      INTEGER J
#include "limevl.inc"
#include "jetset.inc"
#include "jetpar.inc"
C----------------------------------------------------------------------
C
      IF ( USELIM ) THEN
        EVOLMS=SQRT(PJSET(1,J)**2+PJSET(2,J)**2)*CONCUT
      ELSE
        EVOLMS=FUDGE*SQRT(QSQ)
      ENDIF
  999 RETURN
      END
