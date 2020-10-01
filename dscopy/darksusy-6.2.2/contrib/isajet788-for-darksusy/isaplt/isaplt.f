#include "PILOT.inc"
      SUBROUTINE ISAPLT(JTIN,JTLIS)
C
C          Skeleton of a job to analyse ISAJET events from tape JTIN
C          with a trivial calorimeter simulation and a simple UA1
C          jet-finding algorithm.
C
C          Histogram results with CERN HBOOK 3 or 4. Each histogram
C          should be booked in SETUP with 1 < ID < 100. A second 
C          histogram is then automatically booked with ID + 100. The
C          first histogram is filled in USER. At the end of a run, the
C          cross section SIGF is used to fill the second histogram with
C          the proper weight, and the first histogram is cleared.
C
C          Print 1-dim histograms with contents and errors in E format.
C          Print 2-dim histograms with automatic scaling.
C          All printing is on tape JTLIS.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
C          ISAJET common blocks
#include "itapes.inc"
#include "rectp.inc"
#include "final.inc"
#include "primar.inc"
C
C          ISAPLT common blocks
#include "myhist.inc"
C
C          HBOOK and other variables. NWMEMO = number of words for 
C          histograms; it may have to be increased.
      INTEGER NWMEMO
      PARAMETER (NWMEMO=10000)
#ifdef HBOOK4_X
      COMMON/PAWC/HMEMOR(NWMEMO)
      REAL HMEMOR
#elif defined(HBOOK3_X)
      COMMON//HMEMOR(NWMEMO)
      REAL HMEMOR
#endif
      LOGICAL HEXIST
      CHARACTER*40 V,VISAJE
      INTEGER JTIN,JTLIS,ITAPE,IFL,I
      REAL SIGWT
C
C          Initialize
      V=VISAJE()
      CALL HLIMIT(NWMEMO)
      ITLIS=IABS(JTLIS)
      ITAPE=IABS(JTIN)
      CALL HOUTPU(ITLIS)
      CALL HERMES(ITLIS)
C
C          Set up histograms
      CALL SETUP
C
C          Read next record until eof encountered
      REWIND ITAPE
1     CONTINUE
      IFL=0
      CALL RDTAPE(ITAPE,IFL)
      IF(IFL.EQ.-1) GO TO 3
      IF(IRECTP.EQ.200) GO TO 1
      IF(IRECTP.EQ.300) GO TO 2
C
C          Event record
      CALL USER
      GO TO 1
C
C          End of run record. Weight histograms.
2     CONTINUE
      SIGWT=SIGF/FLOAT(NEVENT)
      IF(SIGWT.EQ.0.) SIGWT=1.
      DO 100 I=1,MXHIST
        IF(.NOT.HEXIST(I)) GO TO 100
        CALL HOPERA(I,'+',I+NHSHFT,I+NHSHFT,SIGWT,1.)
#ifdef HBOOK4_X
        CALL HRESET(I,' ')
#elif defined(HBOOK3_X)
        CALL HRESET(I)
#endif
100   CONTINUE
      GO TO 1
C
C          End of job. Output histograms.
3     CALL EHIST
      RETURN
      END
