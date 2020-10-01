#include "PILOT.inc"
      SUBROUTINE EVOL02
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 2 (E+E-)
C-        Lorentz frames and perform initial and final QCD jet
C-        evolution in leading-log approximation.
C-
C-   Created  13-AUG-1991   Frank E. Paige,Serban D. Protopopescu
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
#include "keys.inc"
#include "frame.inc"
      REAL    EVOLMS
      INTEGER I,K,J,NJSAVE,NJFINL
C----------------------------------------------------------------------
C
C          Copy momenta from /PJETS/ to /JETSET/
      N0JETS=NJSET+1
      CALL IPJSET
      NJSAVE=NJSET
C
C          Final state evolution.
C          Define Lorentz frames and JMATCH pointers for jet evolution
C          and fragmentation.
C
      CALL IFRAMS(N0JETS,NJSAVE,1,.FALSE.)
C
C          Set maximum off-shell masses and JDCAY flags.
C
      NJFINL=N0JETS
      DO 310 J=N0JETS,NJSAVE
        IF((IABS(JTYPE(J)).LT.10).OR.
     $     (IABS(JTYPE(J)).GE.21.AND.IABS(JTYPE(J)).LE.29).OR.
     $     (IABS(JTYPE(J)).GE.41.AND.IABS(JTYPE(J)).LE.46))THEN
          PJSET(5,J)=EVOLMS(J,1.0)
          JDCAY(J)=-1
        ENDIF
310   CONTINUE
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END