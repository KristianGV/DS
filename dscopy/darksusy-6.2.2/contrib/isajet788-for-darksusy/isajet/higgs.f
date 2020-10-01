#include "PILOT.inc"
      SUBROUTINE HIGGS
C
C          FINISH HIGGS GENERATION STARTED BY DRLLYN FOR DECAY
C          HIGGS --> W W.
C
C          VER 7.14: TEST BOTH JET1 AND JET2 FOR W,Z FOR SAFETY
C
#include "itapes.inc"
#include "qcdpar.inc"
#include "jetpar.inc"
#include "pjets.inc"
#include "primar.inc"
#include "q1q2.inc"
#include "jetsig.inc"
#include "qsave.inc"
#include "wcon.inc"
#include "const.inc"
#include "hcon.inc"
C
      DIMENSION X(2)
      EQUIVALENCE (X(1),X1)
C
      IDABS1=IABS(IDJETS(1))
      IDABS2=IABS(IDJETS(2))
      IF(IDABS1.NE.80.AND.IDABS1.NE.90.AND.
     $IDABS2.NE.80.AND.IDABS2.NE.90) THEN
        NPAIR=0
        DO 100 I=1,4
        IDPAIR(I)=0
        JPAIR(I)=0
        DO 110 K=1,5
110     PPAIR(K,I)=0.
100     CONTINUE
      ELSE
        CALL WPAIR
      ENDIF
      RETURN
      END
