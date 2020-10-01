#include "PILOT.inc"
      FUNCTION SSZHX(X)
C----------------------------------------------------------------
C          Auxiliary function for Z -> HL Z*. Called by SSTEST.
C----------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "ssmode.inc"
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
#include "sstype.inc"
C
      REAL X,SSZHX
      REAL MHL,MZ,GZ,R,DEN,TERM
C
      MHL=TMP(1)
      MZ=AMZ
      GZ=GAMZ
      R=MHL/MZ
      TERM=(1.-X+X**2/12.+2*R**2/3.)*SQRT(X**2-4*R**2)
      DEN=(X-R**2)**2+(GZ/MZ)**2
      SSZHX=TERM/DEN
      RETURN
      END
