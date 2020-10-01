#include "PILOT.inc"
      SUBROUTINE SETNXT
C
C            RESET LIMITS BEFORE NEXT SET
C
#include "itapes.inc"
#include "lstprt.inc"
#include "totals.inc"
#include "dylim.inc"
#include "jetlim.inc"
#include "primar.inc"
#include "jetset.inc"
#include "partcl.inc"
      DATA UNDEF/-1.E9/
      DO 1 I=1,36
      IF(SETLMJ(I)) BLIMS(I)=UNDEF
    1 CONTINUE
      DO 2 I=1,12
      IF(SETLMQ(I)) BLIM1(I)=UNDEF
    2 CONTINUE
C          RESET /TOTALS/
      NKINPT=0
      NWGEN=0
      NKEEP=0
      SUMWT=0.
C          RESET /LSTPRT/
      LSTPRT=0
C          RESET NJSET AND NPTCL
      NJSET=0
      NPTCL=0
      NPAIR=0
      RETURN
      END
