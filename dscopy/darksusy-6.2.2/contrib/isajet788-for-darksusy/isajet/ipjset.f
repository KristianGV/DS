#include "PILOT.inc"
      SUBROUTINE IPJSET
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-       Initialize PJSET starting from PJETS
C-
C-   Created  14-AUG-1991   Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "primar.inc"
#include "pjets.inc"
#include "jetset.inc"
      INTEGER I,K
C----------------------------------------------------------------------
      DO 110 I=1,NJET
        NJSET=NJSET+1
        JORIG(NJSET)=JPACK*I
        JTYPE(NJSET)=IDJETS(I)
        JDCAY(NJSET)=0
        DO 115 K=1,5
115     PJSET(K,NJSET)=PJETS(K,I)
        IFRAME(I)=1
110   CONTINUE
  999 RETURN
      END
