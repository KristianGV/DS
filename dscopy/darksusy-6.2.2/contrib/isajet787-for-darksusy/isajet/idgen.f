#include "PILOT.inc"
      SUBROUTINE IDGEN
C
C          Call system date and time routines (non-standard) to set up
C          run identification:
C          IDVER=100*VERSN     (integer ISAJET version number)
C          IDG(1)=YYMMDD       (integer year-month-day)
C          IDG(2)=HHMMSS       (integer hour-minute-second)
C
#include "itapes.inc"
#include "idrun.inc"
      CHARACTER*10 CHAR,DATE,TIME
      DIMENSION ISUN(3)
C
C          Default run id is zero.
      IYMD=0.
      IHMS=0.
#ifdef CDC_X
C          Call CDC date and time and convert to integer.
      CHAR=DATE()
      READ(CHAR,'(1X,I2,1X,I2,1X,I2,1X)') IA,IB,IC
      IYMD=10000*IC+100*IA+IB
      CHAR=TIME()
      READ(CHAR,'(1X,I2,1X,I2,1X,I2,1X)') IA,IB,IC
      IHMS=10000*IA+100*IB+IC
#elif defined(ETA_X)
C          Call ETA date and time and convert to integer.
      CHAR=DATE()
      READ(CHAR,'(I2,1X,I2,1X,I2)') IA,IB,IC
      IYMD=10000*IC+100*IA+IB
      CHAR=TIME()
      READ(CHAR,'(I2,1X,I2,1X,I2)') IA,IB,IC
      IHMS=10000*IA+100*IB+IC
#elif defined(IDATE_X)
C          Call date and time for SUN, SGI, and g77/gfortran
      CALL IDATE(ISUN)
      IYMD=10000*ISUN(3)+100*ISUN(2)+ISUN(1)
      CALL ITIME(ISUN)
      IHMS=10000*ISUN(1)+100*ISUN(2)+ISUN(3)
#elif defined(CERN_X)
C          Call DATIME for date and time. (In Cern library)
      CALL DATIME(IYMD,IHMS)
#endif
      IDG(1)=IYMD
      IDG(2)=IHMS
      RETURN
      END
