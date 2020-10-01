#include "PILOT.inc"
        REAL FUNCTION PITRTR(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PITRTR: t_R squark self-energy
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
      COMPLEX TEMP,N(4,4),AC0TTR(4),BC0TTR(4),FTTRR(4),GTTRR(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,FT,FB,C(4),DU(4),DD(4),G1,G2,G3
     $,MH0(4),MHP(2),LHTRTI(4,2),LHTRBI(2,2),FBTRR(4),VMAT(2,2)
     $,UMAT(2,2),GBTRR(4),AP0TTR(4),BP0TTR(4)
     $,APPBTR(2),BPPBTR(2),ACPBTR(2),BCPBTR(2)
     $,LSTRTR(2),LSTLTR(2),LHTRTR(2),LHTLTR(2),LSTRTI(2,2)
     $,LHTRBR(2),LHTRBL(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2,YUR,YUL
     $,YDR,YDL,YER,YEL,YNL,EUR,SWW2,GUR,SINA,COSA,SINBE,COSBE
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
      FT=MTQ/VUQ
      FB=MBQ/VDQ
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
      DU(1)=SINA2
      DU(2)=COSA2
      DU(3)=SINBE2
      DU(4)=COSBE2
      DD(1)=COSA2
      DD(2)=SINA2
      DD(3)=COSBE2
      DD(4)=SINBE2
      MH0(1)=MSS(30)
      MH0(2)=MSS(29)
      MH0(3)=AMZ
      MH0(4)=MSS(31)
      MHP(2)=AMW
      MHP(1)=MSS(32)
      YUR=-4./3.
      YUL=1./3.
      YDR=2./3.
      YDL=1./3.
      YER=2.
      YEL=-1.
      YNL=-1.
      EUR=-2./3.
      SWW2=1-(AMW/AMZ)**2
      GUR=-EUR*SWW2
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
        AP0TTR(I)=0.
        BP0TTR(I)=0.
        AC0TTR(I)=0.
        BC0TTR(I)=0.
      ENDDO
      AP0TTR(1)=GGP/SR2*YUR
      BP0TTR(4)=FT
      DO I=1,4
        DO II=1,4
          AC0TTR(I)=AC0TTR(I)+CONJG(N(I,II))*AP0TTR(II)
          BC0TTR(I)=BC0TTR(I)+N(I,II)*BC0TTR(II)
        ENDDO
      ENDDO
      DO I=1,4
        FTTRR(I)=CONJG(AC0TTR(I))*AC0TTR(I)
     $+CONJG(BC0TTR(I))*BC0TTR(I)
        GTTRR(I)=CONJG(BC0TTR(I))*AC0TTR(I)
     $+CONJG(AC0TTR(I))*BC0TTR(I)
      ENDDO
      DO I=1,2
        APPBTR(I)=0.
        BPPBTR(I)=0.
        ACPBTR(I)=0.
        BCPBTR(I)=0.
      ENDDO
      APPBTR(2)=-FT
      DO I=1,2
        DO II=1,2
          ACPBTR(I)=ACPBTR(I)+VMAT(I,II)*APPBTR(II)
          BCPBTR(I)=BCPBTR(I)+UMAT(I,II)*BPPBTR(II)
        ENDDO
      ENDDO
      DO I=1,2
        FBTRR(I)=ACPBTR(I)*ACPBTR(I)+BCPBTR(I)*BCPBTR(I)
        GBTRR(I)=BCPBTR(I)*ACPBTR(I)+ACPBTR(I)*BCPBTR(I)
      ENDDO
      LSTRTR(1)=GG*AMZ/CTHW*GUR*COSBE
      LSTLTR(1)=-FT/SR2*TWOM1
      LSTRTR(2)=-GG*AMZ/CTHW*GUR*SINBE+SR2*FT*MTQ
      LSTLTR(2)=-FT/SR2*AAT
C     First index is for Higgs sector, second is for sfermion
      LSTRTI(1,1)=LSTLTR(1)*COST-LSTRTR(1)*SINT
      LSTRTI(1,2)=LSTLTR(1)*SINT+LSTRTR(1)*COST
      LSTRTI(2,1)=LSTLTR(2)*COST-LSTRTR(2)*SINT
      LSTRTI(2,2)=LSTLTR(2)*SINT+LSTRTR(2)*COST
      LHTRTI(1,1)=COSA*LSTRTI(1,1)-SINA*LSTRTI(2,1)
      LHTRTI(1,2)=COSA*LSTRTI(1,2)-SINA*LSTRTI(2,2)
      LHTRTI(2,1)=SINA*LSTRTI(1,1)+COSA*LSTRTI(2,1)
      LHTRTI(2,2)=SINA*LSTRTI(1,2)+COSA*LSTRTI(2,2)
      LHTRTR(1)=0.
      LHTLTR(1)=FT/SR2*(-TWOM1*COSBE-AAT*SINBE)
      LHTRTR(2)=0.
      LHTLTR(2)=-FT/SR2*(-TWOM1*SINBE+AAT*COSBE)
      LHTRTI(3,1)=-LHTLTR(1)*COST-LHTRTR(1)*SINT
      LHTRTI(3,2)=-LHTLTR(1)*SINT+LHTRTR(1)*COST
      LHTRTI(4,1)=-LHTLTR(2)*COST-LHTRTR(2)*SINT
      LHTRTI(4,2)=-LHTLTR(2)*SINT+LHTRTR(2)*COST
      LHTRBR(1)=0.
      LHTRBL(1)=-FT*(-TWOM1*COSBE-AAT*SINBE)
      LHTRBR(2)=-FT*MBQ*COSBE-FB*MTQ*SINBE
      LHTRBL(2)=FB*(-TWOM1*SINBE+AAT*COSBE)
      LHTRBI(1,1)=LHTRBL(1)*COSB-LHTRBR(1)*SINB
      LHTRBI(1,2)=LHTRBL(1)*SINB+LHTRBR(1)*COSB
      LHTRBI(2,1)=LHTRBL(2)*COSB-LHTRBR(2)*SINB
      LHTRBI(2,2)=LHTRBL(2)*SINB+LHTRBR(2)*COSB
      TEMP=4.*G3**2/3.*(2.*SSG(P,MSS(1),AMTP)
     $+SINT2*SSF(P,MSS(13),0.)+COST2*SSF(P,MSS(12),0.)
     $+SINT2*SSA0(MSS(13))+COST2*SSA0(MSS(12)))
     $+FT**2*(COST2*SSA0(MSS(13))+SINT2*SSA0(MSS(12))
     $+COSB2*SSA0(MSS(10))+SINB2*SSA0(MSS(11)))
     $+4.*GG**2/CHWW2*GUR**2*SSA0(AMZ)
     $+(2./3.)**2*4.*PI/137.036*(SINT2*SSF(P,MSS(13),0.)
     $+COST2*SSF(P,MSS(12),0.))+GG**2/CHWW2*GUR**2
     $*(SINT2*SSF(P,MSS(13),AMZ)+COST2*SSF(P,MSS(12),AMZ))
     $+GGP**2/4.*YUR**2*(SINT2*SSA0(MSS(13))+COST2*SSA0(MSS(12)))
     $+GGP**2/4.*YUR*(3.*(YUL*(SSA0(MSS(2))+SSA0(MSS(8))
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
        TEMP=TEMP+(FT**2*DU(I)-GG**2*GUR/2./CHWW2*C(I))
     $*SSA0(MH0(I))/2.
      ENDDO
      DO I=3,4
        TEMP=TEMP+(FT**2*DD(I)+GG**2*GUR/2./CHWW2*C(I))
     $*SSA0(MHP(I-2))
      ENDDO
      DO I=1,4
          TEMP=TEMP+(LHTRTI(I,1))**2*SSB0(P2,MH0(I),MSS(13))
     $+(LHTRTI(I,2))**2*SSB0(P2,MH0(I),MSS(12))
      ENDDO
      DO I=1,2
        DO II=1,2
          TEMP=TEMP+(LHTRBI(I,II))**2*SSB0(P2,MSS(9+II),MHP(I))
        ENDDO
      ENDDO
      DO I=1,4
        TEMP=TEMP+FTTRR(I)*SSG(P,ABS(MSS(22+I)),AMTP)
     $-2.*GTTRR(I)*ABS(MSS(22+I))*MTQ*SSB0(P2,ABS(MSS(22+I)),AMTP)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FBTRR(I)*SSG(P,ABS(MSS(26+I)),AMBT)
     $-2.*GBTRR(I)*ABS(MSS(26+I))*MBQ*SSB0(P2,ABS(MSS(26+I)),AMBT)
      ENDDO
      PITRTR=REAL(TEMP)/16./PI**2
      RETURN
      END
