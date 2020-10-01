#include "PILOT.inc"
#if defined(VAX_X)
      SUBROUTINE DATIME(ID,IT)
C          CALL VAX DATE AND TIME.
#include "itapes.inc"
      CHARACTER*8 BUF
      CALL IDATE(IMON,IDAY,IYR)
      CALL TIME(BUF)
      ID=10000*IYR+100*IMON+IDAY
      READ(BUF,'(I2,1X,I2,1X,I2)') K1,K2,K3
      IT=10000*K1+100*K2+K3
      RETURN
      END
#endif
