#include "PILOT.inc"
      SUBROUTINE GENMIR(ALF,M32,G,INUHM)
C...  Calculate soft terms for generalized mirage mediation
c...  Formulas are from arXiv:1610.06205
#include "sugxin.inc"
      REAL ALF,M32
      REAL*8 G(31)
      REAL*8 MS,DPI,AL,B1,B2,B3,BYL,BYB,BYT,GAMHU,GAMHD,GAMHUD,GAMHDD
      REAL*8 CM,CM3,CHU,CHD,A3
      REAL*8 GAMQ1,GAMU1,GAMD1,GAML1,GAME1
      REAL*8 GAMQ3,GAMU3,GAMD3,GAML3,GAME3
      REAL*8 GAMQ1D,GAMU1D,GAMD1D,GAML1D,GAME1D
      REAL*8 GAMQ3D,GAMU3D,GAMD3D,GAML3D,GAME3D
      REAL*8 XIQ1,XIU1,XID1,XIL1,XIE1,XIQ3,XIU3,XID3,XIL3,XIE3,XIHU,XIHD
      INTEGER INUHM
      DPI=4.D0*DATAN(1.D0)
      AL=ALF
      MS=DBLE(M32)/16.D0/DPI/DPI        
      CM=XAMIN(1)
      CHU=XAMIN(2)
      CHD=XAMIN(3)
      A3=XAMIN(4)
      CM3=XAMIN(5)
      B1=33.D0/5.D0
      B2=1.D0
      B3=-3.D0
      BYL=3*G(5)**2+4*G(4)**2-3*G(2)**2-9*G(1)**2/5.D0
      BYB=G(6)**2+6*G(5)**2+G(4)**2-16*G(3)**2/3.D0-3*G(2)**2-7*G(1)**2/
     $15.D0
      BYT=6*G(6)**2+G(5)**2-16*G(3)**2/3.D0-3*G(2)**2-13*G(1)**2/15.D0
      G(7)=MS*(AL+B1*G(1)**2)
      G(8)=MS*(AL+B2*G(2)**2)
      G(9)=MS*(AL+B3*G(3)**2)
C-------anomalous dimensions
        GAMQ1=8*G(3)**2/3.D0+3*G(2)**2/2.D0+G(1)**2/30.D0
        GAMU1=8*G(3)**2/3.D0+8*G(1)**2/15.D0
        GAMD1=8*G(3)**2/3.D0+2*G(1)**2/15.D0
        GAML1=3*G(2)**2/2.D0+3*G(1)**2/10.D0
        GAME1=6*G(1)**2/5.D0
        GAMHU=3*G(2)**2/2.D0+3*G(1)**2/10.D0-3*G(6)**2
        GAMHD=3*G(2)**2/2.D0+3*G(1)**2/10.D0-3*G(5)**2-G(4)**2
        GAMQ3=GAMQ1-G(6)**2-G(5)**2
        GAMU3=GAMU1-2*G(6)**2
        GAMD3=GAMD1-2*G(5)**2
        GAML3=GAML1-G(4)**2
        GAME3=GAME1-2*G(4)**2
C-------dotted ADs
        GAMQ1D=-8*G(3)**4+3*G(2)**4/2.D0+11*G(1)**4/50.D0
        GAMU1D=-8*G(3)**4+88*G(1)**4/25.D0
        GAMD1D=-8*G(3)**4+22*G(1)**4/25.D0
        GAML1D=3*G(2)**4/2.D0+99*G(1)**4/50.D0
        GAME1D=198*G(1)**4/25.D0
        GAMHUD=3*G(2)**4/2.D0+99*G(1)**4/50.D0-3*G(6)**2*BYT
        GAMHDD=3*G(2)**4/2.D0+99*G(1)**4/50.D0-3*G(5)**2*BYB-G(4)**2*BYL
        GAMQ3D=GAMQ1D-G(6)**2*BYT-G(5)**2*BYB
        GAMU3D=GAMU1D-2*G(6)**2*BYT
        GAMD3D=GAMD1D-2*G(5)**2*BYB
        GAML3D=GAML1D-G(4)**2*BYL
        GAME3D=GAME1D-2*G(4)**2*BYL
C-------XI FACTORS
        XIQ1=(-8*G(3)**2/3.D0-3*G(2)**2/2.D0-G(1)**2/30.D0)/2.D0
        XIU1=(-8*G(3)**2/3.D0-8*G(1)**2/15.D0)/2.D0
        XID1=(-8*G(3)**2/3.D0-2*G(1)**2/15.D0)/2.D0
        XIL1=(-3*G(2)**2/2.D0-3*G(1)**2/10.D0)/2.D0
        XIE1=(-6*G(1)**2/5.D0)/2.D0
        XIHU=(-3*G(2)**2/2.D0-3*G(1)**2/10.D0+9*G(6)**2)/2.D0
        XIHD=(-3*G(2)**2/2.D0-3*G(1)**2/10.D0+9*G(5)**2+3*G(4)**2)/2.D0
        XIQ3=XIQ1+3*(G(6)**2+G(5)**2)/2.D0
        XIU3=XIU1+6*G(6)**2/2.D0
        XID3=XID1+6*G(5)**2/2.D0
        XIL3=XIL1+3*G(4)**2/2.D0
        XIE3=XIE1+6*G(4)**2/2.D0
      G(10)=MS*(-A3*AL+GAML3+GAMHD+GAME3)
      G(11)=MS*(-A3*AL+GAMQ3+GAMHD+GAMD3)
      G(12)=MS*(-A3*AL+GAMQ3+GAMHU+GAMU3)
      IF (INUHM.EQ.0) THEN
        G(13)=MS**2*(CHD*AL**2+4*AL*XIHD-GAMHDD)
        G(14)=MS**2*(CHU*AL**2+4*AL*XIHU-GAMHUD)
      END IF
      G(15)=MS**2*(CM*AL**2+4*AL*XIE1-GAME1D)
      G(16)=MS**2*(CM*AL**2+4*AL*XIL1-GAML1D)
      G(17)=MS**2*(CM*AL**2+4*AL*XID1-GAMD1D)
      G(18)=MS**2*(CM*AL**2+4*AL*XIU1-GAMU1D)
      G(19)=MS**2*(CM*AL**2+4*AL*XIQ1-GAMQ1D)
      G(20)=MS**2*(CM3*AL**2+4*AL*XIE3-GAME3D)
      G(21)=MS**2*(CM3*AL**2+4*AL*XIL3-GAML3D)
      G(22)=MS**2*(CM3*AL**2+4*AL*XID3-GAMD3D)
      G(23)=MS**2*(CM3*AL**2+4*AL*XIU3-GAMU3D)
      G(24)=MS**2*(CM3*AL**2+4*AL*XIQ3-GAMQ3D)
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
