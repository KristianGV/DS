#include "PILOT.inc"
        REAL FUNCTION SSZZG1(XARG)
C-----------------------------------------------------------------------
C          SSZIBF: ziss -> zjss + gm
C          Corresponds to integral I_H in Tadas' notes
C
C          Change variables: XARG -> y, 0 < y < 1,
C          y = sqrt(1-x), dx = -2 y dy
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
C          xx=1-xarg**2 with better precision
      XX=XARG
      XX=(1-XX)*(1+XX)
      TOP=-MZINI**2*XX*(1.D0-XX)+MF**2*XX+MSF**2*(1.D0-XX)
      BOT=-MZFIN**2*XX*(1.D0-XX)+MF**2*XX+MSF**2*(1.D0-XX)
      F=DLOG(TOP/BOT)/(MZINI**2-MZFIN**2)/(1.D0-XX)
      SSZZG1=2*XARG*F
      RETURN
      END
