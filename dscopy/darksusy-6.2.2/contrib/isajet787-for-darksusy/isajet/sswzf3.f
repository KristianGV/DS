#include "PILOT.inc"
        REAL FUNCTION SSWZF3(QQ)
C-----------------------------------------------------------------------
C          SSWZBF: wiss -> zjss f fbar
C          Baer's PHIFUN
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
C
      REAL QQ
      DOUBLE PRECISION M1,M2,M3,T,B,XLOG,SQBKT,Q,PI
      DATA PI/3.14159265D0/
C
      Q=QQ
      M1=TMP(1)
      M2=TMP(2)
      M3=TMP(3)
C
      T=M2**2*(M1-2*Q)-M1*M3**2
      B=(M1-2*Q)*(M2**2-2*M1*Q-M3**2)
      XLOG=DLOG(T/B)
      SQBKT=-Q*(M1**2-M3**2-2*M1*Q)/M1/(M1-2*Q)-
     $(2*M1*Q-M2**2+M3**2)*XLOG/2./M1
      SSWZF3=.5*PI**2*M1*M3*SQBKT/(M1**2-M2**2-2*M1*Q)
      RETURN
      END
