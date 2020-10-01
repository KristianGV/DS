#include "PILOT.inc"
      SUBROUTINE EEBEG
C          INITIALIZE E-E+ EVENTS FOR DOLOG
#include "itapes.inc"
#include "primar.inc"
#include "jetlim.inc"
#include "jetpar.inc"
      DO 100 I=1,2
      PMIN(I)=HALFE
      PMAX(I)=-1.E9
100   CONTINUE
      QSQ=SCM
      IDIN(1)=12
      IDIN(2)=-12
      RETURN
      END
