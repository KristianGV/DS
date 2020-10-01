#include "PILOT.inc"
      SUBROUTINE SSSVME(ME)
C
C          Set MSSMOD flag for last mode in /SSMODE/
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "ssmode.inc"
C
      INTEGER ME
      MSSMOD(NSSMOD)=ME
      RETURN
      END
