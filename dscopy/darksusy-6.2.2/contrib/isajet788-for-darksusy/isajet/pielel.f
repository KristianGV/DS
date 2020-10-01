#include "PILOT.inc"
        REAL FUNCTION PIELEL(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PIELEL: left selectron self-energy
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
      COMPLEX TEMP,N(4,4),AC0BBL(4),BC0BBL(4),FLLLL(4),GLLLL(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,C(4),G1,G2,G3
     $,MH0(4),MHP(2),LHLLLI(4,2),LHLLN(2),FNLLL(4),VMAT(2,2)
     $,UMAT(2,2),AP0BBL(4),BP0BBL(4)
     $,APPTBL(2),BPPTBL(2),ACPTBL(2),BCPTBL(2)
     $,LSBLBL(2),LSBLBR(2),LSBLBI(2,2)
     $,LHBLTL(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2,I3UL,I3DL,I3EL,I3NL,YUR,YUL
     $,YDR,YDL,YER,YEL,YNL,EEL,SWW2,GEL,SINA,COSA,SINBE,COSBE
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
      EEL=-1.
      SWW2=1-(AMW/AMZ)**2
      GEL=I3EL-EEL*SWW2
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
        AP0BBL(I)=0.
        BP0BBL(I)=0.
        AC0BBL(I)=0.
        BC0BBL(I)=0.
      ENDDO
      BP0BBL(1)=GGP/SR2*YEL
      BP0BBL(2)=SR2*GG*I3EL
      DO I=1,4
        DO II=1,4
          AC0BBL(I)=AC0BBL(I)+CONJG(N(I,II))*AP0BBL(II)
          BC0BBL(I)=BC0BBL(I)+N(I,II)*BC0BBL(II)
        ENDDO
      ENDDO
      DO I=1,4
        FLLLL(I)=CONJG(AC0BBL(I))*AC0BBL(I)
     $+CONJG(BC0BBL(I))*BC0BBL(I)
        GLLLL(I)=CONJG(BC0BBL(I))*AC0BBL(I)
     $+CONJG(AC0BBL(I))*BC0BBL(I)
      ENDDO
      DO I=1,2
        APPTBL(I)=0.
        BPPTBL(I)=0.
        ACPTBL(I)=0.
        BCPTBL(I)=0.
      ENDDO
      BPPTBL(1)=GG
      DO I=1,2
        DO II=1,2
          ACPTBL(I)=ACPTBL(I)+VMAT(I,II)*APPTBL(II)
          BCPTBL(I)=BCPTBL(I)+UMAT(I,II)*BPPTBL(II)
        ENDDO
      ENDDO
      DO I=1,2
        FNLLL(I)=ACPTBL(I)*ACPTBL(I)+BCPTBL(I)*BCPTBL(I)
      ENDDO
      LSBLBL(1)=GG*AMZ/CTHW*GEL*COSBE
      LSBLBR(1)=0.
      LSBLBL(2)=-GG*AMZ/CTHW*GEL*SINBE
      LSBLBR(2)=0.
C     First index is for Higgs sector, second is for sfermion
      LSBLBI(1,1)=LSBLBL(1)
      LSBLBI(1,2)=LSBLBR(1)
      LSBLBI(2,1)=LSBLBL(2)
      LSBLBI(2,2)=LSBLBR(2)
      LHLLLI(1,1)=COSA*LSBLBI(1,1)-SINA*LSBLBI(2,1)
      LHLLLI(1,2)=COSA*LSBLBI(1,2)-SINA*LSBLBI(2,2)
      LHLLLI(2,1)=SINA*LSBLBI(1,1)+COSA*LSBLBI(2,1)
      LHLLLI(2,2)=SINA*LSBLBI(1,2)+COSA*LSBLBI(2,2)
      LHLLLI(3,1)=0.
      LHLLLI(3,2)=0.
      LHLLLI(4,1)=0.
      LHLLLI(4,2)=0.
      LHBLTL(1)=-GG*AMW/SR2*(COSBE2-SINBE2)
      LHBLTL(2)=GG*AMW/SR2*2.*COSBE*SINBE
      LHLLN(1)=LHBLTL(1)
      LHLLN(2)=LHBLTL(2)
      TEMP=4.*GG**2/CHWW2*GEL**2*SSA0(AMZ)+2.*GG**2*SSA0(AMW)
     $+4.*PI/137.036*SSF(P,MSS(17),0.)+GG**2/CHWW2*GEL**2
     $*SSF(P,MSS(17),AMZ)+GG**2/2.*SSF(P,MSS(14),AMW)
     $+GG**2/4.*(SSA0(MSS(17))+2.*SSA0(MSS(14)))
     $+GG**2*I3EL*(3.*(I3UL*(SSA0(MSS(2))+SSA0(MSS(8))
     $+COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12)))
     $+I3DL*(SSA0(MSS(4))+SSA0(MSS(6))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11))))
     $+I3EL*(SSA0(MSS(17))+SSA0(MSS(19))
     $+COSL2*SSA0(MSS(21))+SINL2*SSA0(MSS(22)))
     $+I3NL*(SSA0(MSS(14))+SSA0(MSS(15))+SSA0(MSS(16))))
     $+GGP**2/4.*YEL**2*SSA0(MSS(17))
     $+GGP**2/4.*YEL*(3.*(YUL*(SSA0(MSS(2))+SSA0(MSS(8))
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
        TEMP=TEMP-GG**2*GEL/2./CHWW2*C(I)
     $*SSA0(MH0(I))/2.
      ENDDO
      DO I=3,4
        TEMP=TEMP+GG**2*(GEL/2./CHWW2-I3EL)*C(I)
     $*SSA0(MHP(I-2))
      ENDDO
      DO I=1,4
          TEMP=TEMP+LHLLLI(I,1)**2*SSB0(P2,MH0(I),MSS(17))
     $+LHLLLI(I,2)**2*SSB0(P2,MH0(I),MSS(18))
      ENDDO
      DO I=1,2
          TEMP=TEMP+LHLLN(I)**2*SSB0(P2,MSS(14),MHP(I))
      ENDDO
      DO I=1,4
        TEMP=TEMP+FLLLL(I)*SSG(P,ABS(MSS(22+I)),AME)
     $-2.*GLLLL(I)*ABS(MSS(22+I))*AME*SSB0(P2,ABS(MSS(22+I)),AME)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FNLLL(I)*SSG(P,ABS(MSS(26+I)),0.)
      ENDDO
      PIELEL=REAL(TEMP)/16./PI**2
      RETURN
      END
