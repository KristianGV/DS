#include "PILOT.inc"
        REAL FUNCTION SIG0L(P2,I,J,G1,G2,CHW)
C-----------------------------------------------------------------------
C          SIG0L: Neutralino mass matrix correction
C     Taken from Damien M. Pierce, Jonathan A. Bagger, Konstantin T. Matchev,
C     Ren-jie Zhang, Nucl.Phys.B491:3-67,1997, hep-ph/9606211
C     Programmed by Tadas Krupovnickas 1/2003
C          P2 = 4-momentum squared
C          CHW = Cos(theta_W) in DR bar scheme
C     Ordering: u=1,s=2,t=3,d=4,c=5,b=6,e=7,mu=8,tau=9,nue=10,num=11,nut=12
C     I and J are indexes in Pierce's base. To convert to ISAJET's basis
C     make the transformation 1->4, 2->3, 3->2, 4->1 for both inexes
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sugmg.inc"
#include "ssinf.inc"
      COMPLEX*16 SSB1,AP0C0Z(4,4),AP0C0H(4,4,4)
      DOUBLE PRECISION TEMP,AP0FF(4,12,2),AP0CPW(4,2),AP0CPH(4,2,2)
      REAL AP0R(4,12),AP0L(4,12),AP0PPW(4,2),AP0P0Z(4,4)
     $,AP0PPH(4,2,2),AP0P0H(4,4,4),G1,G2
      REAL P2,CHW,COST,SINT,COSB,SINB,COSL,SINL,GG,GGP
      REAL YUR,YUL,YDR,YDL,YER,YEL,YNR,YNL,I3UL,I3DL,I3EL,I3NL
      REAL THX,THY,VMAT(2,2),UMAT(2,2),BE,FB,FT,FMIX(12,2,2)
      INTEGER THIW1,THIW2
      INTEGER I,J,II,III
      COMPLEX IMAG,N(4,4)
      PARAMETER (IMAG=(0.,1.))
      REAL PI,SR2
      PI=4*ATAN(1.)
      SR2=SQRT(2.)
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      GG=G2
      GGP=SQRT(3./5.)*G1
      BE=ATAN(VUQ/VDQ)
      FB=MBQ/VDQ
      FT=MTQ/VUQ
      THX=SIGN(1.,TAN(GAMMAL))
      THY=SIGN(1.,TAN(GAMMAR))
      IF (SIGN(1.,AMW1SS).EQ.1.) THEN
         THIW1=0
      ELSE
         THIW1=1
      END IF
      IF (SIGN(1.,AMW2SS).EQ.1.) THEN
         THIW2=0
      ELSE
         THIW2=1
      END IF
      UMAT(1,1)=-(-1)**THIW1*SIN(GAMMAR)
      UMAT(1,2)=(-1)**THIW1*COS(GAMMAR)
      UMAT(2,1)=-(-1)**THIW2*THY*COS(GAMMAR)
      UMAT(2,2)=-(-1)**THIW2*THY*SIN(GAMMAR)
      VMAT(1,1)=-SIN(GAMMAL)
      VMAT(1,2)=COS(GAMMAL)
      VMAT(2,1)=-THX*COS(GAMMAL)
      VMAT(2,2)=-THX*SIN(GAMMAL)
      DO II=1,4
        IF (SIGN(1.,AMZISS(II)).EQ.1.) THEN
          III=0
        ELSE
          III=1
        END IF
        N(II,1)=-IMAG**III*ZMIXSS(4,II)
        N(II,2)=-IMAG**III*ZMIXSS(3,II)
        N(II,3)=IMAG**III*ZMIXSS(2,II)
        N(II,4)=IMAG**III*ZMIXSS(1,II)
      ENDDO
      YUR=-4./3.
      YUL=1./3.
      YDR=2./3.
      YDL=1./3.
      YER=2.
      YEL=-1.
      YNR=0.
      YNL=-1.
      I3UL=1./2.
      I3DL=-1./2.
      I3EL=-1./2.
      I3NL=1./2.
      DO II=1,4
        DO III=1,12
          AP0R(II,III)=0.
          AP0L(II,III)=0.
        ENDDO
      ENDDO
      AP0R(1,1)=GGP/SR2*YUR
      AP0R(1,2)=GGP/SR2*YUR
      AP0R(1,3)=GGP/SR2*YUR
      AP0R(1,4)=GGP/SR2*YDR
      AP0R(1,5)=GGP/SR2*YDR
      AP0R(1,6)=GGP/SR2*YDR
      AP0R(1,7)=GGP/SR2*YER
      AP0R(1,8)=GGP/SR2*YER
      AP0R(1,9)=GGP/SR2*YER
      AP0R(1,10)=GGP/SR2*YNR
      AP0R(1,11)=GGP/SR2*YNR
      AP0R(1,12)=GGP/SR2*YNR
      AP0L(3,6)=FB
      AP0L(4,3)=FT
      DO II=1,12
        FMIX(II,1,1)=1.
        FMIX(II,1,2)=0.
        FMIX(II,2,1)=0.
        FMIX(II,2,2)=1.
      ENDDO
      FMIX(3,1,1)=COST
      FMIX(3,1,2)=-SINT
      FMIX(3,2,1)=SINT
      FMIX(3,2,2)=COST
      FMIX(6,1,1)=COSB
      FMIX(6,1,2)=-SINB
      FMIX(6,2,1)=SINB
      FMIX(6,2,2)=COSB
      FMIX(9,1,1)=COSL
      FMIX(9,1,2)=-SINL
      FMIX(9,2,1)=SINL
      FMIX(9,2,2)=COSL
      DO II=1,4
        DO III=1,12
          AP0FF(II,III,1)=FMIX(III,1,1)*AP0L(II,III)
     $+FMIX(III,1,2)*AP0R(II,III)
          AP0FF(II,III,2)=FMIX(III,2,1)*AP0L(II,III)
     $+FMIX(III,2,2)*AP0R(II,III)
        ENDDO
      ENDDO
      DO II=1,4
        DO III=1,2
          AP0PPW(II,III)=0.
        ENDDO
      ENDDO
      AP0PPW(2,1)=-GG
      AP0PPW(4,2)=GG/SR2
      DO II=1,4
        DO III=1,2
          AP0CPW(II,III)=VMAT(III,1)*AP0PPW(II,1)
     $+VMAT(III,2)*AP0PPW(II,2)
        ENDDO
      ENDDO
      DO II=1,4
        DO III=1,4
          AP0P0Z(II,III)=0.
        ENDDO
      ENDDO
      AP0P0Z(3,3)=GG/2./CHW
      AP0P0Z(4,4)=-GG/2./CHW
      DO II=1,4
        DO III=1,4
          AP0C0Z(II,III)=N(III,3)*AP0P0Z(II,3)
     $+N(III,4)*AP0P0Z(II,4)
        ENDDO
      ENDDO
      DO II=1,4
        DO III=1,2
          AP0PPH(II,III,1)=0.
          AP0PPH(II,III,2)=0.
        ENDDO
      ENDDO
      AP0PPH(1,2,1)=GGP/SR2
      AP0PPH(2,2,1)=GG/SR2
      AP0PPH(3,1,1)=-GG
      DO II=1,4
        DO III=1,2
          AP0CPH(II,III,1)=COS(BE)*(UMAT(III,1)*AP0PPH(II,1,1)
     $+UMAT(III,2)*AP0PPH(II,2,1))+SIN(BE)*(UMAT(III,1)*AP0PPH(II,1,2)
     $+UMAT(III,2)*AP0PPH(II,2,2))
          AP0CPH(II,III,2)=-SIN(BE)*(UMAT(III,1)*AP0PPH(II,1,1)
     $+UMAT(III,2)*AP0PPH(II,2,1))+COS(BE)*(UMAT(III,1)*AP0PPH(II,1,2)
     $+UMAT(III,2)*AP0PPH(II,2,2))
        ENDDO
      ENDDO
      DO II=1,4
        DO III=1,4
          AP0P0H(II,III,1)=0.
          AP0P0H(II,III,2)=0.
          AP0P0H(II,III,3)=0.
          AP0P0H(II,III,4)=0.
        ENDDO
      ENDDO
      AP0P0H(1,3,1)=-GGP/2.
      AP0P0H(1,4,2)=GGP/2.
      AP0P0H(2,3,1)=GG/2.
      AP0P0H(2,4,2)=-GG/2.
      AP0P0H(1,3,3)=-GGP/2.
      AP0P0H(1,4,4)=-GGP/2.
      AP0P0H(2,3,3)=GG/2.
      AP0P0H(2,4,4)=GG/2.
      AP0P0H(3,1,1)=-GGP/2.
      AP0P0H(4,1,2)=GGP/2.
      AP0P0H(3,2,1)=GG/2.
      AP0P0H(4,2,2)=-GG/2.
      AP0P0H(3,1,3)=-GGP/2.
      AP0P0H(4,1,4)=-GGP/2.
      AP0P0H(3,2,3)=GG/2.
      AP0P0H(4,2,4)=GG/2.
      DO II=1,4
        DO III=1,4
          AP0C0H(II,III,1)=(CONJG(N(III,1))*AP0P0H(II,1,1)
     $+CONJG(N(III,2))*AP0P0H(II,2,1)+CONJG(N(III,3))*AP0P0H(II,3,1)
     $+CONJG(N(III,4))*AP0P0H(II,4,1))*COS(ALFAH)
     $-(CONJG(N(III,1))*AP0P0H(II,1,2)+CONJG(N(III,2))*AP0P0H(II,2,2)
     $+CONJG(N(III,3))*AP0P0H(II,3,2)+CONJG(N(III,2))
     $*AP0P0H(II,4,2))*SIN(ALFAH)
          AP0C0H(II,III,2)=(CONJG(N(III,1))*AP0P0H(II,1,1)
     $+CONJG(N(III,2))*AP0P0H(II,2,1)+CONJG(N(III,3))*AP0P0H(II,3,1)
     $+CONJG(N(III,4))*AP0P0H(II,4,1))*SIN(ALFAH)
     $+(CONJG(N(III,1))*AP0P0H(II,1,2)+CONJG(N(III,2))*AP0P0H(II,2,2)
     $+CONJG(N(III,3))*AP0P0H(II,3,2)+CONJG(N(III,2))
     $*AP0P0H(II,4,2))*COS(ALFAH)
          AP0C0H(II,III,3)=(CONJG(N(III,1))*AP0P0H(II,1,3)
     $+CONJG(N(III,2))*AP0P0H(II,2,3)+CONJG(N(III,3))*AP0P0H(II,3,3)
     $+CONJG(N(III,4))*AP0P0H(II,4,3))*COS(BE)
     $+(CONJG(N(III,1))*AP0P0H(II,1,4)+CONJG(N(III,2))*AP0P0H(II,2,4)
     $+CONJG(N(III,3))*AP0P0H(II,3,4)+CONJG(N(III,2))
     $*AP0P0H(II,4,4))*SIN(BE)
          AP0C0H(II,III,4)=-(CONJG(N(III,1))*AP0P0H(II,1,3)
     $+CONJG(N(III,2))*AP0P0H(II,2,3)+CONJG(N(III,3))*AP0P0H(II,3,3)
     $+CONJG(N(III,4))*AP0P0H(II,4,3))*SIN(BE)
     $+(CONJG(N(III,1))*AP0P0H(II,1,4)+CONJG(N(III,2))*AP0P0H(II,2,4)
     $+CONJG(N(III,3))*AP0P0H(II,3,4)+CONJG(N(III,2))
     $*AP0P0H(II,4,4))*COS(BE)
        ENDDO
      ENDDO
      TEMP=DBLE(3.*(AP0FF(I,1,1)*AP0FF(J,1,1)
     $*SSB1(P2,AMUP,ABS(MSS(2)))
     $+AP0FF(I,1,2)*AP0FF(J,1,2)*SSB1(P2,AMUP,ABS(MSS(3)))
     $+AP0FF(I,2,1)*AP0FF(J,2,1)*SSB1(P2,AMST,ABS(MSS(6)))
     $+AP0FF(I,2,2)*AP0FF(J,2,2)*SSB1(P2,AMST,ABS(MSS(7)))
     $+AP0FF(I,3,1)*AP0FF(J,3,1)*SSB1(P2,AMTP,ABS(MSS(13)))
     $+AP0FF(I,3,2)*AP0FF(J,3,2)*SSB1(P2,AMTP,ABS(MSS(12)))
     $+AP0FF(I,4,1)*AP0FF(J,4,1)*SSB1(P2,AMDN,ABS(MSS(4)))
     $+AP0FF(I,4,2)*AP0FF(J,4,2)*SSB1(P2,AMDN,ABS(MSS(5)))
     $+AP0FF(I,5,1)*AP0FF(J,5,1)*SSB1(P2,AMCH,ABS(MSS(8)))
     $+AP0FF(I,5,2)*AP0FF(J,5,2)*SSB1(P2,AMCH,ABS(MSS(9)))
     $+AP0FF(I,6,1)*AP0FF(J,6,1)*SSB1(P2,AMBT,ABS(MSS(10)))
     $+AP0FF(I,6,2)*AP0FF(J,6,2)*SSB1(P2,AMBT,ABS(MSS(11))))
     $+AP0FF(I,7,1)*AP0FF(J,7,1)*SSB1(P2,AME,ABS(MSS(17)))
     $+AP0FF(I,7,2)*AP0FF(J,7,2)*SSB1(P2,AME,ABS(MSS(18)))
     $+AP0FF(I,8,1)*AP0FF(J,8,1)*SSB1(P2,AMMU,ABS(MSS(19)))
     $+AP0FF(I,8,2)*AP0FF(J,8,2)*SSB1(P2,AMMU,ABS(MSS(20)))
     $+AP0FF(I,9,1)*AP0FF(J,9,1)*SSB1(P2,AMTAU,ABS(MSS(21)))
     $+AP0FF(I,9,2)*AP0FF(J,9,2)*SSB1(P2,AMTAU,ABS(MSS(22)))
     $+2.*(AP0CPW(I,1)*AP0CPW(J,1)*SSB1(P2,ABS(MSS(27)),AMW)
     $+AP0CPW(I,2)*AP0CPW(J,2)*SSB1(P2,ABS(MSS(28)),AMW))
     $+CONJG(AP0C0Z(I,1))*AP0C0Z(J,1)*SSB1(P2,ABS(MSS(23)),AMZ)
     $+CONJG(AP0C0Z(I,2))*AP0C0Z(J,2)*SSB1(P2,ABS(MSS(24)),AMZ)
     $+CONJG(AP0C0Z(I,3))*AP0C0Z(J,3)*SSB1(P2,ABS(MSS(25)),AMZ)
     $+CONJG(AP0C0Z(I,4))*AP0C0Z(J,4)*SSB1(P2,ABS(MSS(26)),AMZ)
     $+AP0CPH(I,1,1)*AP0CPH(J,1,1)*SSB1(P2,ABS(MSS(27)),AMW)
     $+AP0CPH(I,1,2)*AP0CPH(J,1,2)*SSB1(P2,ABS(MSS(27)),ABS(MSS(32)))
     $+AP0CPH(I,2,1)*AP0CPH(J,2,1)*SSB1(P2,ABS(MSS(28)),AMW)
     $+AP0CPH(I,2,2)*AP0CPH(J,2,2)*SSB1(P2,ABS(MSS(28)),ABS(MSS(32)))
     $+(CONJG(AP0C0H(I,1,1))*AP0C0H(J,1,1)
     $*SSB1(P2,ABS(MSS(23)),ABS(MSS(30)))
     $+CONJG(AP0C0H(I,1,2))*AP0C0H(J,1,2)
     $*SSB1(P2,ABS(MSS(23)),ABS(MSS(29)))
     $+CONJG(AP0C0H(I,1,3))*AP0C0H(J,1,3)*SSB1(P2,ABS(MSS(23)),AMZ)
     $+CONJG(AP0C0H(I,1,4))*AP0C0H(J,1,4)
     $*SSB1(P2,ABS(MSS(23)),ABS(MSS(31)))
     $+CONJG(AP0C0H(I,2,1))*AP0C0H(J,2,1)
     $*SSB1(P2,ABS(MSS(24)),ABS(MSS(30)))
     $+CONJG(AP0C0H(I,2,2))*AP0C0H(J,2,2)
     $*SSB1(P2,ABS(MSS(24)),ABS(MSS(29)))
     $+CONJG(AP0C0H(I,2,3))*AP0C0H(J,2,3)*SSB1(P2,ABS(MSS(24)),AMZ)
     $+CONJG(AP0C0H(I,2,4))*AP0C0H(J,2,4)
     $*SSB1(P2,ABS(MSS(24)),ABS(MSS(31)))
     $+CONJG(AP0C0H(I,3,1))*AP0C0H(J,3,1)
     $*SSB1(P2,ABS(MSS(25)),ABS(MSS(30)))
     $+CONJG(AP0C0H(I,3,2))*AP0C0H(J,3,2)
     $*SSB1(P2,ABS(MSS(25)),ABS(MSS(29)))
     $+CONJG(AP0C0H(I,3,3))*AP0C0H(J,3,3)*SSB1(P2,ABS(MSS(25)),AMZ)
     $+CONJG(AP0C0H(I,3,4))*AP0C0H(J,3,4)
     $*SSB1(P2,ABS(MSS(25)),ABS(MSS(31)))
     $+CONJG(AP0C0H(I,4,1))*AP0C0H(J,4,1)
     $*SSB1(P2,ABS(MSS(26)),ABS(MSS(30)))
     $+CONJG(AP0C0H(I,4,2))*AP0C0H(J,4,2)
     $*SSB1(P2,ABS(MSS(26)),ABS(MSS(29)))
     $+CONJG(AP0C0H(I,4,3))*AP0C0H(J,4,3)*SSB1(P2,ABS(MSS(26)),AMZ)
     $+CONJG(AP0C0H(I,4,4))*AP0C0H(J,4,4)
     $*SSB1(P2,ABS(MSS(26)),ABS(MSS(31)))
     $)/2.)/16./PI**2
      SIG0L=TEMP
      RETURN
      END
