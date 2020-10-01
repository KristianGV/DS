#include "PILOT.inc"
        REAL FUNCTION PILLLR(P2,G1,G2,G3,CTHW)
C-----------------------------------------------------------------------
C          PILLLR: off-diagonal tau slpeton self-energy
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
      COMPLEX TEMP,N(4,4),AC0BBL(4),BC0BBL(4),FLLLR(4),GLLLR(4)
     $,AC0BBR(4),BC0BBR(4)
      REAL P2,P,GG,GGP,CTHW,CHWW2,BE,FL,FB,G1,G2,G3
     $,MH0(4),MHP(2),LHLLLI(4,2),LHLRLI(4,2),LHLLN(2),LHLRN(2)
     $,FNLLR(4),VMAT(2,2),UMAT(2,2),AP0BBL(4),BP0BBL(4),AP0BBR(4)
     $,BP0BBR(4),APPTBR(2),BPPTBR(2),ACPTBR(2),BCPTBR(2)
     $,APPTBL(2),BPPTBL(2),ACPTBL(2),BCPTBL(2),LHBRBR(2)
     $,LSBLBL(2),LSBLBR(2),LSBRBR(2),LHBLBL(2),LHBLBR(2),LSBLBI(2,2)
     $,LHBLTL(2),LSBRBI(2,2),LHBRTL(2)
      REAL COST,SINT,COST2,SINT2,COSB,SINB,COSB2,SINB2
     $,COSL,SINL,COSL2,SINL2,THX,THY
      REAL SINA2,COSA2,SINBE2,COSBE2
     $,I3EL,YER,YEL,EEL,EER,SWW2,GEL,GER,SINA,COSA,SINBE,COSBE
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
      FL=MLQ/VDQ
      FB=MBQ/VDQ
      SINA2=SIN(ALFAH)**2
      COSA2=COS(ALFAH)**2
      SINBE2=SIN(BE)**2
      COSBE2=COS(BE)**2
      SINA=SIN(ALFAH)
      COSA=COS(ALFAH)
      SINBE=SIN(BE)
      COSBE=COS(BE)
      MH0(1)=MSS(30)
      MH0(2)=MSS(29)
      MH0(3)=AMZ
      MH0(4)=MSS(31)
      MHP(2)=AMW
      MHP(1)=MSS(32)
      I3EL=-1./2.
      YER=2.
      YEL=-1.
      EEL=-1.
      EER=1.
      SWW2=1-(AMW/AMZ)**2
      GEL=I3EL-EEL*SWW2
      GER=-EER*SWW2
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
      AP0BBL(3)=FL
      DO I=1,4
        DO II=1,4
          AC0BBL(I)=AC0BBL(I)+CONJG(N(I,II))*AP0BBL(II)
          BC0BBL(I)=BC0BBL(I)+N(I,II)*BC0BBL(II)
        ENDDO
      ENDDO
      DO I=1,4
        AP0BBR(I)=0.
        BP0BBR(I)=0.
        AC0BBR(I)=0.
        BC0BBR(I)=0.
      ENDDO
      AP0BBR(1)=GGP/SR2*YER
      BP0BBR(3)=FL
      DO I=1,4
        DO II=1,4
          AC0BBR(I)=AC0BBR(I)+CONJG(N(I,II))*AP0BBR(II)
          BC0BBR(I)=BC0BBR(I)+N(I,II)*BC0BBR(II)
        ENDDO
      ENDDO
      DO I=1,4
        FLLLR(I)=CONJG(AC0BBL(I))*AC0BBR(I)
     $+CONJG(BC0BBL(I))*BC0BBR(I)
        GLLLR(I)=CONJG(BC0BBL(I))*AC0BBR(I)
     $+CONJG(AC0BBL(I))*BC0BBR(I)
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
        APPTBR(I)=0.
        BPPTBR(I)=0.
        ACPTBR(I)=0.
        BCPTBR(I)=0.
      ENDDO
      BPPTBR(1)=-FL
      DO I=1,2
        DO II=1,2
          ACPTBR(I)=ACPTBR(I)+VMAT(I,II)*APPTBR(II)
          BCPTBR(I)=BCPTBR(I)+UMAT(I,II)*BPPTBR(II)
        ENDDO
      ENDDO
      DO I=1,2
        FNLLR(I)=ACPTBL(I)*ACPTBR(I)+BCPTBL(I)*BCPTBR(I)
      ENDDO
      LSBLBL(1)=GG*AMZ/CTHW*GEL*COSBE+SR2*FL*MLQ
      LSBLBR(1)=-FL/SR2*AAL
      LSBLBL(2)=-GG*AMZ/CTHW*GEL*SINBE
      LSBLBR(2)=-FL/SR2*TWOM1
C     First index is for Higgs sector, second is for sfermion
      LSBLBI(1,1)=LSBLBL(1)*COSL-LSBLBR(1)*SINL
      LSBLBI(1,2)=LSBLBL(1)*SINL+LSBLBR(1)*COSL
      LSBLBI(2,1)=LSBLBL(2)*COSL-LSBLBR(2)*SINL
      LSBLBI(2,2)=LSBLBL(2)*SINL+LSBLBR(2)*COSL
      LHLLLI(1,1)=COSA*LSBLBI(1,1)-SINA*LSBLBI(2,1)
      LHLLLI(1,2)=COSA*LSBLBI(1,2)-SINA*LSBLBI(2,2)
      LHLLLI(2,1)=SINA*LSBLBI(1,1)+COSA*LSBLBI(2,1)
      LHLLLI(2,2)=SINA*LSBLBI(1,2)+COSA*LSBLBI(2,2)
      LHBLBL(1)=0.
      LHBLBR(1)=-FL/SR2*(-TWOM1*SINBE-AAL*COSBE)
      LHBLBL(2)=0.
      LHBLBR(2)=-FL/SR2*(-TWOM1*COSBE+AAL*SINBE)
      LHLLLI(3,1)=LHBLBL(1)*COSL-LHBLBR(1)*SINL
      LHLLLI(3,2)=LHBLBL(1)*SINL+LHBLBR(1)*COSL
      LHLLLI(4,1)=LHBLBL(2)*COSL-LHBLBR(2)*SINL
      LHLLLI(4,2)=LHBLBL(2)*SINL+LHBLBR(2)*COSL
      LHBLTL(1)=-GG*AMW/SR2*(COSBE2-SINBE2)+FL*MLQ*COSBE
      LHBLTL(2)=GG*AMW/SR2*2.*COSBE*SINBE-FL*MLQ*SINBE
      LHLLN(1)=LHBLTL(1)
      LHLLN(2)=LHBLTL(2)
      LSBRBR(1)=GG*AMZ/CTHW*GER*COSBE+SR2*FL*MLQ
      LSBLBR(1)=-FL/SR2*AAL
      LSBRBR(2)=-GG*AMZ/CTHW*GER*SINBE
      LSBLBR(2)=-FL/SR2*TWOM1
C     First index is for Higgs sector, second is for sfermion
      LSBRBI(1,1)=LSBLBR(1)*COSL-LSBRBR(1)*SINL
      LSBRBI(1,2)=LSBLBR(1)*SINL+LSBRBR(1)*COSL
      LSBRBI(2,1)=LSBLBR(2)*COSL-LSBRBR(2)*SINL
      LSBRBI(2,2)=LSBLBR(2)*SINL+LSBRBR(2)*COSL
      LHLRLI(1,1)=COSA*LSBRBI(1,1)-SINA*LSBRBI(2,1)
      LHLRLI(1,2)=COSA*LSBRBI(1,2)-SINA*LSBRBI(2,2)
      LHLRLI(2,1)=SINA*LSBRBI(1,1)+COSA*LSBRBI(2,1)
      LHLRLI(2,2)=SINA*LSBRBI(1,2)+COSA*LSBRBI(2,2)
      LHBRBR(1)=0.
      LHBLBR(1)=-FL/SR2*(-TWOM1*SINBE-AAL*COSBE)
      LHBRBR(2)=0.
      LHBLBR(2)=-FL/SR2*(-TWOM1*COSBE+AAL*SINBE)
      LHLRLI(3,1)=-LHBLBR(1)*COSL-LHBRBR(1)*SINL
      LHLRLI(3,2)=-LHBLBR(1)*SINL+LHBRBR(1)*COSL
      LHLRLI(4,1)=-LHBLBR(2)*COSL-LHBRBR(2)*SINL
      LHLRLI(4,2)=-LHBLBR(2)*SINL+LHBRBR(2)*COSL
      LHBRTL(1)=FL*(-TWOM1*SINBE-AAL*COSBE)
      LHBRTL(2)=FL*(-TWOM1*COSBE+AAL*SINBE)
      LHLRN(1)=LHBRTL(1)
      LHLRN(2)=LHBRTL(2)
      TEMP=-FL/2.*(3.*FB*2.*SINB*COSB*(SSA0(MSS(10))-SSA0(MSS(11)))
     $+FL*SINL*COSL*(SSA0(MSS(21))-SSA0(MSS(22))))
     $-GGP**2/4.*YEL*YER*SINL*COSL*(SSA0(MSS(21))-SSA0(MSS(22)))
     $-4.*PI/137.036*SINL*COSL*(SSF(P,MSS(21),0.)
     $-SSF(P,MSS(22),0.))+GG**2/CHWW2*GEL*GER*SINL*COSL
     $*(SSF(P,MSS(21),AMZ)-SSF(P,MSS(22),AMZ))
      DO I=1,4
          TEMP=TEMP+LHLLLI(I,1)*LHLRLI(I,1)
     $*SSB0(P2,MH0(I),MSS(21))+LHLLLI(I,2)*LHLRLI(I,2)
     $*SSB0(P2,MH0(I),MSS(22))
      ENDDO
      DO I=1,2
          TEMP=TEMP+LHLLN(I)*LHLRN(I)*SSB0(P2,MSS(16),MHP(I))
      ENDDO
      DO I=1,4
        TEMP=TEMP+FLLLR(I)*SSG(P,ABS(MSS(22+I)),AMTAU)
     $-2.*GLLLR(I)*ABS(MSS(22+I))*MLQ*SSB0(P2,ABS(MSS(22+I)),AMTAU)
      ENDDO
      DO I=1,2
        TEMP=TEMP+FNLLR(I)*SSG(P,ABS(MSS(26+I)),0.)
      ENDDO
      PILLLR=REAL(TEMP)/16./PI**2
      RETURN
      END
