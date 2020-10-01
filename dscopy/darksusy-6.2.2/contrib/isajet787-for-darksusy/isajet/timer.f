#include "PILOT.inc"
      SUBROUTINE TIMER(IT)
C
C          CALL SYSTEM CPU CLOCK -- MACHINE DEPENDENT.
C          IT=1 FOR RUN START TIME.
C          IT=2 FOR RUN STOP TIME.
C
#include "itapes.inc"
#include "times.inc"
      DIMENSION TIMES(2)
      EQUIVALENCE (TIMES(1),TIME1)
      DIMENSION TTT(2)
#ifdef VAX_X
      INTEGER CPUTIM(2),ITMLST(4),NHSEC
      EXTERNAL JPI$_CPUTIM
#endif
C
C          DEFAULT IS TO RETURN ZERO.
      TNOW=0.
#ifdef CDC_X
C          SECOND GIVES CPU TIME ON CDC.
      CALL SECOND(TNOW)
#elif defined(ETA_X)
C          SECOND GIVES CPU TIME ON ETA.
      TNOW=SECOND()
#elif defined(IBMRT_X)
C          MCLOCK GIVES CPU TIME ON IBM RS/6000.
      TNOW=FLOAT(MCLOCK())/60.
#elif defined(ETIME_X)
C          ETIME gives CPU time on SUN, SGI, and g77/gfortran
      TNOW=ETIME(TTT)
#elif defined(VAX_X)
C          VAX HAS NO FORTRAN FUNCTION FOR CPU TIME.
C          FOLLOWING PROVIDED BY T. KILLIAN
      ITMLST(1)=ISHFT(%LOC(JPI$_CPUTIM),16)+4
      ITMLST(2)=%LOC(NHSEC)
      ITMLST(3)=0
      ITMLST(4)=0
      CALL SYS$GETJPI(,,,ITMLST,,,)
      TNOW=.01*NHSEC
#elif defined(CERN_X)
      CALL TIMEST(1.E7)
      CALL TIMEX(TNOW)
#endif
      TIMES(IT)=TNOW
      RETURN
      END
