#include "PILOT.inc"
      SUBROUTINE SSPRT(ID)
C-----------------------------------------------------------------------
C
C     Print decay modes for ID. Note these need not be contiguous,
C     so the loop is over all modes in /SSMODE/.
C
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "ssmode.inc"
C
      INTEGER ID,I,K,NOUT
      CHARACTER*5 SSID,LBLIN,LBLOUT(4)
C
      NOUT=0
      DO 100 I=1,NSSMOD
        IF(ISSMOD(I).NE.ID) GO TO 100
        NOUT=NOUT+1
        LBLIN=SSID(ISSMOD(I))
        DO 110 K=1,4
110     LBLOUT(K)=SSID(JSSMOD(K,I))
        WRITE(LOUT,1000) LBLIN,(LBLOUT(K),K=1,4),GSSMOD(I),BSSMOD(I)
1000    FORMAT(1X,A5,'  -->  ',4(A5,2X),2E15.5)
100   CONTINUE
C
      IF(NOUT.GT.0) WRITE(LOUT,*) ' '
C
      RETURN
      END
