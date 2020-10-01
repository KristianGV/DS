#include "PILOT.inc"
      FUNCTION WZSS(T,U,T1,U1,T3,U3,P1,P2)
C          DECAY DISTRIBUTION FOR W- Z0 PAIRS FROM SCHOONSCHIP(1980).
C          SQUARE OF S GRAPH.
#include "itapes.inc"
#include "wwpar.inc"
      DIMENSION P1(4),P2(4)
#ifdef DOUBLE_X
      DOUBLE PRECISION WZSS
      DOUBLE PRECISION T,U,T1,U1,T3,U3,P1,P2
      DOUBLE PRECISION WM4,ZM4,WZM2,CSXCS
#endif
      WM4=WM2**2
      ZM4=ZM2**2
      WZM2=WM2*ZM2
      CSXCS=CS**2
      WZSS=
     1 +CSXCS*CV3*(-32.*WM2*ZM2*WM4-32.*WM2*ZM2*ZM4-128.*WM2*T1*T3**2-1
     1 28.*WM2*T1*ZM4+128.*WM2*U1*T3*U3-64.*WM2*T3*ZM4-64.*WM2*S13*ZM4+
     1 64.*ZM2*T1*U1*T3+64.*ZM2*T1*U1*U3-64.*ZM2*T1*WM4-64.*ZM2*T1**2*T
     1 3-64.*ZM2*U1**2*U3-128.*ZM2*T3*WM4-64.*ZM2*S13*WM4+128.*T1*U1*T3
     1 *U3-192.*T1*T3*WZM2-64.*T1*S13*WZM2-64.*T1**2*T3**2-32.*T1**2*WZ
     1 M2-32.*T1**2*ZM4+64.*U1*U3*WZM2-64.*U1*S13*WZM2-64.*U1**2*U3**2-
     1 32.*U1**2*WZM2-32.*U1**2*ZM4-128.*T3*S13*WZM2-64.*T3**2*WZM2-64.
     1 *T3**2*WM4-64.*S13**2*WZM2-96.*WM4*ZM4)
      WZSS=WZSS
     1 +CSXCS*CV3*T*(64.*WM2*T1*T3-32.*WM2*U1*T3-32.*WM2*U1*U3-64.*WM2*
     1 T3*S13+64.*WM2*T3**2+96.*WM2*ZM4+64.*ZM2*T1*T3-32.*ZM2*T1*S13+32
     1 .*ZM2*T1**2-32.*ZM2*U1*T3-32.*ZM2*U1*U3+32.*ZM2*U1*S13+32.*ZM2*U
     1 1**2+96.*ZM2*WM4-32.*T1*U1*T3-32.*T1*U1*U3-64.*T1*T3*S13+64.*T1*
     1 T3**2+128.*T1*WZM2+32.*T1*ZM4+32.*T1**2*T3-64.*U1*T3*U3+64.*U1*U
     1 3*S13+32.*U1**2*U3+128.*T3*WZM2+32.*T3*WM4+32.*S13*WZM2)

      WZSS=WZSS
     1 +CSXCS*CV3*T*U*(-32.*WM2*T3+32.*WM2*S13-32.*ZM2*T1+32.*ZM2*S13-3
     1 2.*T1*T3+32.*T1*S13+64.*U1*T3+32.*U1*U3+32.*U1*S13+64.*T3*S13+64
     1 .*S13**2-32.*WZM2)
     1 +CSXCS*CV3*T**2*(-32.*WM2*T3+32.*WM2*S13-32.*ZM2*T1+32.*ZM2*S13-
     1 32.*T1*T3+32.*T1*S13+32.*U1*T3+32.*U1*U3+64.*T3*S13-64.*WZM2)
     1 +CSXCS*CV3*T**2*U*(-32.*S13)
     1 +CSXCS*CV3*T**3*(-32.*S13)
      WZSS=WZSS
     1 +CSXCS*CV3*U*(64.*WM2*T1*T3-32.*WM2*U1*T3-32.*WM2*U1*U3+64.*WM2*
     1 T3*S13+64.*WM2*T3**2+32.*WM2*ZM4+64.*ZM2*T1*T3+32.*ZM2*T1*S13+32
     1 .*ZM2*T1**2-32.*ZM2*U1*T3-32.*ZM2*U1*U3-32.*ZM2*U1*S13+32.*ZM2*U
     1 1**2+32.*ZM2*WM4-32.*T1*U1*T3-32.*T1*U1*U3+64.*T1*T3*S13+64.*T1*
     1 T3**2+64.*T1*WZM2+32.*T1*ZM4+32.*T1**2*T3-64.*U1*T3*U3-64.*U1*U3
     1 *S13+32.*U1**2*U3+64.*T3*WZM2+32.*T3*WM4+32.*S13*WZM2)
     1 +CSXCS*CV3*U**2*(32.*U1*T3+32.*U1*S13)
      WZSS=WZSS
     1 +CSXCS*CA3*(32.*WM2*ZM2*WM4+32.*WM2*ZM2*ZM4+128.*WM2*T1*ZM4+64.*
     1 WM2*T3*ZM4+64.*WM2*S13*ZM4+64.*ZM2*T1*U1*T3-64.*ZM2*T1*U1*U3+64.
     1 *ZM2*T1*WM4+64.*ZM2*T1**2*T3-64.*ZM2*U1**2*U3+128.*ZM2*T3*WM4+64
     1 .*ZM2*S13*WM4+192.*T1*T3*WZM2+64.*T1*S13*WZM2+32.*T1**2*WZM2+32.
     1 *T1**2*ZM4-64.*U1*U3*WZM2-64.*U1*S13*WZM2-32.*U1**2*WZM2-32.*U1*
     1 *2*ZM4+96.*WM4*ZM4)
      WZSS=WZSS
     1 +CSXCS*CA3*T*(-64.*WM2*T1*T3-32.*WM2*U1*T3+32.*WM2*U1*U3-96.*WM2
     1 *ZM4-64.*ZM2*T1*T3+32.*ZM2*T1*S13-32.*ZM2*T1**2-32.*ZM2*U1*T3+32
     1 .*ZM2*U1*U3+32.*ZM2*U1*S13+32.*ZM2*U1**2-96.*ZM2*WM4-32.*T1*U1*T
     1 3+32.*T1*U1*U3-128.*T1*WZM2-32.*T1*ZM4-32.*T1**2*T3+32.*U1**2*U3
     1 -128.*T3*WZM2-32.*T3*WM4-32.*S13*WZM2)
     1 +CSXCS*CA3*T*U*(32.*WM2*T3-32.*WM2*S13+32.*ZM2*T1-32.*ZM2*S13+32
     1 .*T1*T3-32.*T1*S13+64.*U1*T3-32.*U1*U3+32.*U1*S13+32.*WZM2)
      WZSS=WZSS
     1 +CSXCS*CA3*T**2*(32.*WM2*T3-32.*WM2*S13+32.*ZM2*T1-32.*ZM2*S13+3
     1 2.*T1*T3-32.*T1*S13+32.*U1*T3-32.*U1*U3+64.*WZM2)
     1 +CSXCS*CA3*T**2*U*(32.*S13)
     1 +CSXCS*CA3*T**3*(32.*S13)
      WZSS=WZSS
     1 +CSXCS*CA3*U*(-64.*WM2*T1*T3-32.*WM2*U1*T3+32.*WM2*U1*U3-32.*WM2
     1 *ZM4-64.*ZM2*T1*T3-32.*ZM2*T1*S13-32.*ZM2*T1**2-32.*ZM2*U1*T3+32
     1 .*ZM2*U1*U3-32.*ZM2*U1*S13+32.*ZM2*U1**2-32.*ZM2*WM4-32.*T1*U1*T
     1 3+32.*T1*U1*U3-64.*T1*WZM2-32.*T1*ZM4-32.*T1**2*T3+32.*U1**2*U3-
     1 64.*T3*WZM2-32.*T3*WM4-32.*S13*WZM2)
     1 +CSXCS*CA3*U**2*(32.*U1*T3+32.*U1*S13)
      RETURN
      END