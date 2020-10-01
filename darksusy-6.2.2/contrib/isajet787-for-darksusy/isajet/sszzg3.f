#include "PILOT.inc"
        REAL FUNCTION SSZZG3(XARG)
C-----------------------------------------------------------------------
C          SSZIBF: ziss -> zjss + gm
C          Corresponds to integral K_H in Tadas' notes
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
        REAL XARG
        DOUBLE PRECISION XX,F,MF,MSF,MZINI,MZFIN,TOP,BOT
C
        MF=TMP(1)
        MSF=TMP(2)
        MZINI=TMP(3)
        MZFIN=TMP(4)
C
        XX=XARG
        TOP=-MZINI**2*XX*(1.D0-XX)+MF**2*XX+MSF**2*(1.D0-XX)
        BOT=-MZFIN**2*XX*(1.D0-XX)+MF**2*XX+MSF**2*(1.D0-XX)
        F=-(1+TOP*DLOG(TOP/BOT)/XX/(1.D0-XX)/(MZINI**2-
     $MZFIN**2))/(MZINI**2-MZFIN**2)
        SSZZG3=F
        RETURN
        END
