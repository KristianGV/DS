#include "PILOT.inc"
      SUBROUTINE NAMSB(M32,G,INUHM)
C...  Calculate soft terms for natural anomaly mediation
c...  Formulas are from WSS
#include "sugxin.inc"
      REAL M32
      REAL*8 G(31)
      REAL*8 MS,DPI,B1,B2,B3,BLHAT,BBHAT,BTHAT
      REAL*8 M012,M03,MHUB,MHDB,A0
      INTEGER INUHM
      DPI=4.D0*DATAN(1.D0)
      MS=DBLE(M32)/16.D0/DPI/DPI        
      M012=XAMIN(10)
      M03=XAMIN(11)
      MHUB=XAMIN(2)
      MHDB=XAMIN(3)
      A0=XAMIN(1)
C-----set nAMSB boundary conditions -----------------------------
      B1=33.D0/5.D0
      B2=1.D0
      B3=-3.D0
      BLHAT=G(4)*(-9*G(1)**2/5.D0-3*G(2)**2+3*G(5)**2+4*G(4)**2)
      BBHAT=G(5)*(-7*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
      BTHAT=G(6)*(-13*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             6*G(6)**2+G(5)**2)
      G(7)=MS*B1*G(1)**2
      G(8)=MS*B2*G(2)**2
      G(9)=MS*B3*G(3)**2
      G(10)=-BLHAT*MS/G(4)+A0
      G(11)=-BBHAT*MS/G(5)+A0
      G(12)=-BTHAT*MS/G(6)+A0
      IF (INUHM.EQ.0) THEN
        G(13)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(5)*BBHAT+
     ,G(4)*BLHAT)*MS**2+MHDB**2
        G(14)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(6)*BTHAT)*
     ,MS**2+MHUB**2
      END IF
      G(15)=(-198*G(1)**4/25.D0)*MS**2+M012**2
      G(16)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0)*MS**2+M012**2
      G(17)=(-22*G(1)**4/25.D0+8*G(3)**4)*MS**2+M012**2
      G(18)=(-88*G(1)**4/25.D0+8*G(3)**4)*MS**2+M012**2
      G(19)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4)*
     ,MS**2+M012**2
      G(20)=(-198*G(1)**4/25.D0+2*G(4)*BLHAT)*MS**2+M03**2
      G(21)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+G(4)*BLHAT)*MS**2+
     ,M03**2
      G(22)=(-22*G(1)**4/25.D0+8*G(3)**4+2*G(5)*BBHAT)*MS**2+M03**2
      G(23)=(-88*G(1)**4/25.D0+8*G(3)**4+2*G(6)*BTHAT)*MS**2+M03**2
      G(24)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4+G(5)*BBHAT+
     ,G(6)*BTHAT)*MS**2+M03**2
C
c        WRITE(6,*) 'g1-6=',g(1),g(2),g(3),g(4),g(5),g(6),inuhm
c        WRITE(6,*) 'g7-12=',g(7),g(8),g(9),g(10),g(11),g(12)
c        WRITE(6,*) 'g13-24=',sign(1.D0,g(13))*sqrt(abs(g(13))),
c     $sign(1.D0,g(14))*sqrt(abs(g(14))),
c     $sign(1.D0,g(15))*sqrt(abs(g(15))),
c     $sign(1.D0,g(16))*sqrt(abs(g(16))),
c     $sign(1.D0,g(17))*sqrt(abs(g(17))),
c     $sign(1.D0,g(18))*sqrt(abs(g(18))),
c     $sign(1.D0,g(19))*sqrt(abs(g(19))),
c     $sign(1.D0,g(20))*sqrt(abs(g(20))),
c     $sign(1.D0,g(21))*sqrt(abs(g(21))),
c     $sign(1.D0,g(22))*sqrt(abs(g(22))),
c     $sign(1.D0,g(23))*sqrt(abs(g(23))),
c     $sign(1.D0,g(24))*sqrt(abs(g(24)))  
      RETURN
      END
