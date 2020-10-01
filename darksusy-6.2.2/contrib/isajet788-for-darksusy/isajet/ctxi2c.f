#include "PILOT.inc"
      SUBROUTINE CTXI2C(IVAL,CVAL,NSIZE)
C-----------------------------------------------------------------------
C          Convert integer array IVAL to character variable CVAL
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
      CHARACTER*(*) CVAL
      INTEGER I,NSIZE
      INTEGER IVAL(NSIZE)
C
      DO 100 I=1,NSIZE
100   CVAL(I:I)=CHAR(IVAL(I))
C
      RETURN
      END
