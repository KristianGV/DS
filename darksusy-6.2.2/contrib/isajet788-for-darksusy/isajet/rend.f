#include "PILOT.inc"
      SUBROUTINE REND
C
C          INVERSE OF WREND
C          READ END RECORD (TYPE 300)
C
#include "itapes.inc"
#include "final.inc"
#include "zevel.inc"
      NKINF=IZEVEL(3)
      SIGF=ZEVEL(4)
      ALUM=ZEVEL(5)
      ACCEPT=ZEVEL(6)
      NRECS=IZEVEL(7)
      RETURN
      END
