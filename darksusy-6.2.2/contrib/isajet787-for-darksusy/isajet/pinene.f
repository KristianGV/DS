#include "PILOT.inc"
        REAL FUNCTION PINENE(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PINENE: electron sneutrino self-energy
C     Taken from Damien M. Pierce, Jonathan A. Bagger, Konstantin T. Matchev,
C     Ren-jie Zhang, Nucl.Phys.B491:3-67,1997, hep-ph/9606211
C          P2 = 4-momentum squared
C          CTHW = Cos(theta_W) in DR bar scheme
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sugmg.inc"
#include "ssinf.inc"
      COMPLEX*16 SSB0,SSG,SSF
      DOUBLE PRECISION SSA0
      COMPLEX TEMP,N(4,4),AC0TTL(4),BC0TTL(4),FTTLL(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,C(4),G1,G2,G3
     $,MH0(4),MHP(2),LHTLTI(4,2),LHTLBI(2,2),FBTLL(4),VMAT(2,2)
     $,UMAT(2,2),GBTLL(4),AP0TTL(4),BP0TTL(4)
     $,APPBTL(2),BPPBTL(2),ACPBTL(2),BCPBTL(2)
     $,LSTLTL(2),LSTLTR(2),LSTLTI(2,2)
     $,LHTLBL(2),LHTLBR(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2,I3UL,I3DL,I3EL,I3NL,YUR,YUL
     $,YDR,YDL,YER,YEL,YNL,SWW2,GNL,SINA,COSA,SINBE,COSBE
      INTEGER I,II,THIW1,THIW2
      COMPLEX IMAG
      PARAMETER (IMAG=(0.,1.))
      REAL PI,SR2
      PI=4*ATAN(1.)
      SR2=SQRT(2.)
      P=SQRT(P2)
      CHWW2=CTHW**2
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      COST2=COST**2
      SINT2=1-COST2
      COSB2=COSB**2
      SINB2=1-COSB2
      COSL2=COSL**2
      SINL2=1-COSL2
      GG=G2
      GGP=SQRT(3./5.)*G1
      BE=ATAN(VUQ/VDQ)
      SINA2=SIN(ALFAH)**2
      COSA2=COS(ALFAH)**2
      SINBE2=SIN(BE)**2
      COSBE2=COS(BE)**2
      SINA=SIN(ALFAH)
      COSA=COS(ALFAH)
      SINBE=SIN(BE)
      COSBE=COS(BE)
      C(1)=-(COSA2-SINA2)
      C(2)=COSA2-SINA2
      C(3)=-(COSBE2-SINBE2)
      C(4)=COSBE2-SINBE2
      MH0(1)=MSS(30)
      MH0(2)=MSS(29)
      MH0(3)=AMZ
      MH0(4)=MSS(31)
      MHP(2)=AMW
      MHP(1)=MSS(32)
      I3UL=1./2.
      I3DL=-1./2.
      I3EL=-1./2.
      I3NL=1./2.
      YUR=-4./3.
      YUL=1./3.
      YDR=2./3.
      YDL=1./3.
      YER=2.
      YEL=-1.
      YNL=-1.
      SWW2=1-(AMW/AMZ)**2
      GNL=I3NL
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
          I=0
        ELSE
          I=1
        END IF
        N(II,1)=-IMAG**I*ZMIXSS(4,II)
        N(II,2)=-IMAG**I*ZMIXSS(3,II)
        N(II,3)=IMAG**I*ZMIXSS(2,II)
        N(II,4)=IMAG**I*ZMIXSS(1,II)
      ENDDO
      DO I=1,4
        AP0TTL(I)=0.
        BP0TTL(I)=0.
        AC0TTL(I)=0.
        BC0TTL(I)=0.
      ENDDO
      BP0TTL(1)=GGP/SR2*YNL
      BP0TTL(2)=SR2*GG*I3NL
      DO I=1,4
        DO II=1,4
          AC0TTL(I)=AC0TTL(I)+CONJG(N(I,II))*AP0TTL(II)
          BC0TTL(I)=BC0TTL(I)+N(I,II)*BC0TTL(II)
        ENDDO
      ENDDO
      DO I=1,4
        FTTLL(I)=CONJG(AC0TTL(I))*AC0TTL(I)
     $+CONJG(BC0TTL(I))*BC0TTL(I)
      ENDDO
      DO I=1,2
        APPBTL(I)=0.
        BPPBTL(I)=0.
        ACPBTL(I)=0.
        BCPBTL(I)=0.
      ENDDO
      APPBTL(1)=GG
      DO I=1,2
        DO II=1,2
          ACPBTL(I)=ACPBTL(I)+VMAT(I,II)*APPBTL(II)
          BCPBTL(I)=BCPBTL(I)+UMAT(I,II)*BPPBTL(II)
        ENDDO
      ENDDO
      DO I=1,2
        FBTLL(I)=ACPBTL(I)*ACPBTL(I)+BCPBTL(I)*BCPBTL(I)
        GBTLL(I)=BCPBTL(I)*ACPBTL(I)+ACPBTL(I)*BCPBTL(I)
      ENDDO
      LSTLTL(1)=GG*AMZ/CTHW*GNL*COSBE
      LSTLTR(1)=0.
      LSTLTL(2)=-GG*AMZ/CTHW*GNL*SINBE
      LSTLTR(2)=0.
C     First index is for Higgs sector, second is for sfermion
      LSTLTI(1,1)=LSTLTL(1)
      LSTLTI(1,2)=LSTLTR(1)
      LSTLTI(2,1)=LSTLTL(2)
      LSTLTI(2,2)=LSTLTR(2)
      LHTLTI(1,1)=COSA*LSTLTI(1,1)-SINA*LSTLTI(2,1)
      LHTLTI(1,2)=COSA*LSTLTI(1,2)-SINA*LSTLTI(2,2)
      LHTLTI(2,1)=SINA*LSTLTI(1,1)+COSA*LSTLTI(2,1)
      LHTLTI(2,2)=SINA*LSTLTI(1,2)+COSA*LSTLTI(2,2)
      LHTLTI(3,1)=0.
      LHTLTI(3,2)=0.
      LHTLTI(4,1)=0.
      LHTLTI(4,2)=0.
      LHTLBL(1)=-GG*AMW/SR2*(COSBE2-SINBE2)
      LHTLBR(1)=0.
      LHTLBL(2)=GG*AMW/SR2*2.*COSBE*SINBE
      LHTLBR(2)=0.
      LHTLBI(1,1)=LHTLBL(1)
      LHTLBI(1,2)=LHTLBR(1)
      LHTLBI(2,1)=LHTLBL(2)
      LHTLBI(2,2)=LHTLBR(2)
      TEMP=4.*GG**2/CHWW2*GNL**2*SSA0(AMZ)+2.*GG**2*SSA0(AMW)
     $+GG**2/CHWW2*GNL**2*SSF(P,MSS(14),AMZ)
     $+GG**2/2.*SSF(P,MSS(17),AMW)
     $+GG**2/4.*(SSA0(MSS(14))+2.*SSA0(MSS(17)))
     $+GG**2*I3NL*(3.*(I3UL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+I3DL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11))))
     $+I3EL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+I3NL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
     $+GGP**2/4.*YNL**2*SSA0(MSS(14))
     $+GGP**2/4.*YNL*(3.*(YUL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+YUR*(SSA0(MSS(3))+SSA0(MSS(9))+SINT2*SSA0(MSS(13))
     $+COST2*SSA0(MSS(12)))+YDL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11)))
     $+YDR*(SSA0(MSS(5))+SSA0(MSS(7))
     $+SINB2*SSA0(MSS(10))+COSB2*SSA0(MSS(11))))
     $+YEL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+YER*(SSA0(MSS(18))+SSA0(MSS(20))
     $+SINL2*SSA0(MSS(21))+COSL2*SSA0(MSS(22)))
     $+YNL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
      DO I=1,4
        TEMP=TEMP-GG**2*GNL/2./CHWW2*C(I)*SSA0(MH0(I))/2.
      ENDDO
      DO I=3,4
        TEMP=TEMP+GG**2*(GNL/2./CHWW2-I3NL)*C(I)*SSA0(MHP(I-2))
      ENDDO
      DO I=1,4
          TEMP=TEMP+LHTLTI(I,1)**2*SSB0(P2,MH0(I),MSS(14))
      ENDDO
      DO I=1,2
        DO II=1,2
          TEMP=TEMP+(LHTLBI(I,II))**2*SSB0(P2,MSS(16+II),MHP(I))
        ENDDO
      ENDDO
      DO I=1,4
        TEMP=TEMP+FTTLL(I)*SSG(P,ABS(MSS(22+I)),0.)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FBTLL(I)*SSG(P,ABS(MSS(26+I)),AME)
     $-2.*GBTLL(I)*ABS(MSS(26+I))*AME*SSB0(P2,ABS(MSS(26+I)),AME)
      ENDDO
      PINENE=REAL(TEMP)/16./PI**2
      RETURN
      END
