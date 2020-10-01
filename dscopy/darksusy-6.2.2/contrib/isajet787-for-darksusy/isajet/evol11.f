#include "PILOT.inc"
      SUBROUTINE EVOL11
C----------------------------------------------------------------------
C-
C-   Purpose and Methods :
C-        Setup for process 11 (EXTRADIM)
C-        Lorentz frames and perform initial and final QCD jet
C-        evolution in leading-log approximation.
C-
C----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "primar.inc"
#include "jetpar.inc"
#include "pjets.inc"
#include "jetset.inc"
#include "jwork.inc"
#include "jwork2.inc"
#include "q1q2.inc"
#include "frame.inc"
#include "wcon.inc"
C
      INTEGER K,NJFINL,J
      REAL EVOLMS
C----------------------------------------------------------------------
C
C          Add recoil jet (jet 3)
      NJSET=NJSET+1
      N0JETS=NJSET
      JORIG(NJSET)=JPACK*3
      JTYPE(NJSET)=IDJETS(3)
      JDCAY(NJSET)=0
      DO 105 K=1,5
105   PJSET(K,NJSET)=PJETS(K,3)
      IFRAME(3)=1

C          Add W (=KKG)
      NJSET=NJSET+1
      N0W=NJSET
      JORIG(NJSET)=0
      JTYPE(NJSET)=IDENTW
      JDCAY(NJSET)=0
      DO 120 K=1,5
120   PJSET(K,NJSET)=QWJET(K)
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      CALL ISTRAD(1.0)
      IF(NJSET.LT.0) RETURN
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      CALL IFRAMS(N0JETS,N0W,1,.FALSE.)
C
C          Set maximum off-shell masses and JDCAY flags.
C
      NJFINL=N0JETS
      DO 321 J=N0JETS,N0W
        IF(IABS(JTYPE(J)).LT.10) THEN
          PJSET(5,J)=EVOLMS(J,WFUDGE)
          JDCAY(J)=-1
        ENDIF
321   CONTINUE
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END
