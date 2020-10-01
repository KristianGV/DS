#include "PILOT.inc"
      SUBROUTINE SETUP
C
C          User routine to book histograms for HBOOK 3 or 4.
C
C          Sets up two histograms for each plot.
C          Histogram I, 1<I<MXHIST, is unweighted and never printed.
C          Histogram I+NHSHFT is weighted by SIGF/NEVENT unless SIGF=0.
C          Scatter plots are automatically scaled.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "itapes.inc"
#include "myhist.inc"
C
      LOGICAL HEXIST
      INTEGER I
C
C          MXHIST=maximum unweighted histogram ID.
      MXHIST=100
      NHSHFT=100
C
C          Book unweighted histograms with HBOOK1 and HBOOK2 below.
C          ID should be less than MXHIST.
C          2-DIM histograms *MUST* be booked with 1 word/channel.
C-----------------------------------------------------------------------
C-
C-
C-
C-
C-
C-----------------------------------------------------------------------
C
C          Provide automatic scaling for scatter plots.
      CALL HSCALE(0,0.)
C
C          Create copies of histograms with ID2=ID1+NHSHFT
C          Set error bar mode to display errors.
      DO 100 I=1,MXHIST
        IF(.NOT.HEXIST(I)) GO TO 100
#ifdef HBOOK4_X
        CALL HCOPY(I,I+NHSHFT,' ')
#elif defined(HBOOK3_X)
        CALL HCOPY(I,I+NHSHFT)
#endif
        CALL HBARX(I)
        CALL HBARX(I+NHSHFT)
100   CONTINUE
      RETURN
      END
