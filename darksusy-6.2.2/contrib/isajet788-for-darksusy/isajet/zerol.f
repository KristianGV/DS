#include "PILOT.inc"
      SUBROUTINE ZEROL(Z,N)
C          SET N VALUES OF Z IN LEVEL2 TO ZERO
#include "itapes.inc"
      DIMENSION Z(N)
#ifdef LEVEL2_X
      LEVEL2,Z
#endif
      DO 1 I=1,N
    1 Z(I)=0
      RETURN
      END
