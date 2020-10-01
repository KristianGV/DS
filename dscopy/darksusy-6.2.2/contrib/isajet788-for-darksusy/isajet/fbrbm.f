#include "PILOT.inc"
      REAL FUNCTION FBRBM(X)
C
C     Integrand for convolution of 
C     bremsstrahlung with beamstrahlung spectra
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "brembm.inc"
C
      REAL EBEAM,ESTRUC,X
C
      FBRBM=EBEAM(X,EB)*ESTRUC(XMIN/X,QSQBM)/X
      RETURN
      END
