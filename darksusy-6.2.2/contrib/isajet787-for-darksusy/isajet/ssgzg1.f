#include "PILOT.inc"
        REAL FUNCTION SSGZG1(XARG)
C-----------------------------------------------------------------------
C          SSGLBF: glss -> ziss + gl
C          Baer's FUNI- removed mass dependence to simplify
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
        REAL XARG
        DOUBLE PRECISION TOP,BOT,XX,F,MQ,MX,MS,MG
C
        MQ=TMP(1)
        MX=TMP(2)
        MS=TMP(3)
        MG=AMGLSS
C
        XX=XARG
        TOP=-MG**2*XX*(1.D0-XX)+MS**2*XX+MQ**2*(1.D0-XX)
        BOT=-MX**2*XX*(1.D0-XX)+MS**2*XX+MQ**2*(1.D0-XX)
        F=DLOG(TOP/BOT)/XX
        SSGZG1=F
        RETURN
        END