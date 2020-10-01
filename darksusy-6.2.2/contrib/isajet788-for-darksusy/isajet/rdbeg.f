#include "PILOT.inc"
      SUBROUTINE RDBEG
C
C          Read first record (type 200)
C          Inverse of WRBEG
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "itapes.inc"
#include "dylim.inc"
#include "frgpar.inc"
#include "idrun.inc"
#include "jetlim.inc"
#include "keys.inc"
#include "primar.inc"
#include "qcdpar.inc"
#include "qlmass.inc"
#include "q1q2.inc"
#include "types.inc"
#include "xmssm.inc"
C
#include "zevel.inc"
C
      INTEGER NL,IDSAVE,IL
C
      IL=3
      IDSAVE=IDVER
      CALL MOVLEI(IZEVEL(IL),IDVER,4)
      IF(IDVER.NE.IDSAVE) WRITE(ITLIS,1777) IDVER,IDSAVE
1777  FORMAT(///,
     $' WARNING:  DATA WERE GENERATED WITH VERSION',I5,/,
     $'           DATA ARE BEING READ WITH VERSION',I5,/,
     $' RESULTS CANNOT BE PREDICTED.'///)
      IL=IL+4
      CALL MOVLEI(IZEVEL(IL),NJET,7)
      IL=IL+7
      NL=NJET*MXGOQ
      IF(NJET.NE.0) CALL MOVLEL(LZEVEL(IL),GOQ(1,1),NL)
      IL=14+NL
      CALL MOVLEL(LZEVEL(IL),KEYS(1),10)
      IL=IL+10
      CALL MOVLEV(ZEVEL(IL),PMIN(1),36)
      IL=IL+36
      IF(.NOT.KEYS(3)) GO TO 11
      CALL MOVLEV(ZEVEL(IL),QMIN,12)
      IL=IL+12
11    CONTINUE
      CALL MOVLEL(LZEVEL(IL),GODY(1),5)
      IL=IL+5
      CALL MOVLEV(ZEVEL(IL),PUD,22)
      IL=IL+22
      CALL MOVLEV(ZEVEL(IL),ALAM,4)
      IL=IL+4
      CALL MOVLEV(ZEVEL(IL),AMLEP(6),3)
      IL=IL+3
      CALL MOVLEI(IZEVEL(IL),LOC(1),100)
      IL=IL+100
      CALL MOVLEL(LZEVEL(IL),GOMSSM,1)
      IL=IL+1
      CALL MOVLEV(ZEVEL(IL),XGLSS,11)
      IL=IL+11
      CALL MOVLEL(LZEVEL(IL),GOSUG,1)
      IL=IL+1
      CALL MOVLEV(ZEVEL(IL),XM0SU,5)
      IL=IL+5
C
      RETURN
      END
