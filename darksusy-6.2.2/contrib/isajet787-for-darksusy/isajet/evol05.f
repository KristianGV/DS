#include "PILOT.inc"
      SUBROUTINE EVOL05
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Setup for process 5 (SUPERSYM)
C-        Lorentz frames and perform initial and final QCD jet
C-        evolution in leading-log approximation.
C-
C-   Created  13-AUG-1991   Frank E. Paige,Serban D. Protopopescu
C-
C----------------------------------------------------------------------
#include "primar.inc"
#include "jetpar.inc"
#include "pjets.inc"
#include "jetset.inc"
#include "jwork.inc"
#include "jwork2.inc"
#include "frame.inc"
      REAL    EVOLMS
      INTEGER I,K,J,NJSAVE,NJFINL,JTABS
C----------------------------------------------------------------------
C
C          Copy momenta from /PJETS/ to /JETSET/
      N0JETS=NJSET+1
      CALL IPJSET
      NJSAVE=NJSET
C
C          Set flags and maximum off-shell masses and generate
C          initial QCD parton shower.
C
      CALL ISTRAD(1.0)
C
      IF(NJSET.LT.0) RETURN
C
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
        DO 325 J=N0JETS,NJSAVE
          JTABS=IABS(JTYPE(J))
          IF(JTABS.GT.20.AND.JTABS.LT.30) THEN
          PJSET(5,J)=EVOLMS(J,1.0)
          JDCAY(J)=-1
        ENDIF
325   CONTINUE
C
C          Produce final-state QCD parton cascade
C
      CALL QCDJET(NJFINL)
C
      RETURN
      END