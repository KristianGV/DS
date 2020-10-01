#include "PILOT.inc"
      SUBROUTINE ISAWND
C
C          WRITE END RECORD, TYPE 300
C          CONTAINS CROSS SECTIONS AND LUMINOSITY
C
#include "itapes.inc"
#include "final.inc"
#include "totals.inc"
#include "zevel.inc"
#include "jetlim.inc"
#include "primar.inc"
#include "dylim.inc"
#include "keys.inc"
#include "sslun.inc"
C
C       keep entry point WREND for backward compatibility
      ENTRY WREND
      ITA=IABS(ITEVT)
      IZEVEL(1)=300
      IZEVEL(2)=1
      IZEVEL(3)=NKINF
      ZEVEL(4)=SIGF
      ZEVEL(5)=ALUM
      ZEVEL(6)=ACCEPT
      IZEVEL(7)=NRECS
      IL=7
      CALL BUFOUT(IL)
      WRITE(ITLIS,1010) NRECS,ITA
1010  FORMAT(////' THIS RUN WROTE',I10,
     1' PHYSICAL RECORDS ON TAPE',I3)
      IF (WRTLHE) THEN
        WRITE(LHEOUT,2001)
        CLOSE(UNIT=LHEOUT)
      END IF
2001  FORMAT('</LesHouchesEvents>')
      RETURN
      END
