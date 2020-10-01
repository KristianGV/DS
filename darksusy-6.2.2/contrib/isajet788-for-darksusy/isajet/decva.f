#include "PILOT.inc"
      LOGICAL FUNCTION DECVA(IP,NADD,IDABS,PREST)
C
C          Compute matrix element unpolarized for V-A.
C          Auxiliary routine for DECAY. 
C
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "itapes.inc"
#include "wcon.inc"
#include "partcl.inc"
#include "dkytab.inc"
#include "jetset.inc"
#include "jwork.inc"
#include "const.inc"
#include "keys.inc"
#include "pjets.inc"
#include "xmssm.inc"
#include "sspols.inc"
C
      REAL PREST(4,6)
      REAL DOT,RANF,WT
      INTEGER IP,NADD,IDABS(5),I,K,I1,I2,IDIPA
C
      DOT(I1,I2)=PREST(4,I1)*PREST(4,I2)-PREST(1,I1)*PREST(1,I2)
     $-PREST(2,I1)*PREST(2,I2)-PREST(3,I1)*PREST(3,I2)
C
      IDIPA=IABS(IDENT(IP))
      IF(IDENT(NPTCL+1).LT.0) THEN
        WT=PPTCL(5,IP)*PREST(4,1)*DOT(2,3)
      ELSEIF(IDENT(NPTCL+2).LT.0) THEN
        WT=PPTCL(5,IP)*PREST(4,2)*DOT(1,3)
      ELSE
        WT=PPTCL(5,IP)*PREST(4,3)*DOT(1,2)
      ENDIF
      IF(WT.LT.RANF()*PPTCL(5,IP)**4/16.) THEN
        DECVA=.FALSE.
      ELSE
        DECVA=.TRUE.
      ENDIF
      RETURN
      END
