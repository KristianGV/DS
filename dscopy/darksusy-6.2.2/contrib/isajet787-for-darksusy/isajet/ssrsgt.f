#include "PILOT.inc"
C--------------------------------------------------------------------
      FUNCTION SSRSGT(QS)
C--------------------------------------------------------------------
C
C     Calculate top quark self energy
C     according to Pierce et al. formulae adapted to Isajet
C
C     Modified by Javier 9/2005 /Log threholds already implemented 
C     through RGE decoupling have been substracted 
C     by a redefinition of B1 function
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sspar.inc"
#include "sssm.inc"
#include "sugmg.inc"
#include "sugpas.inc"
      COMPLEX ZAZT1(4),ZAZT2(4),ZBZT1(4),ZBZT2(4),ZI
      COMPLEX ZAWU(2)
      COMPLEX ZAWB1(2),ZAWB2(2),ZBWB1(2),ZBWB2(2)
      COMPLEX*16 SSB0,SSB1F,SSB1,ZZZ,SIGWI,SIGZI
      REAL SR2,PI,E,G,TANB,BETA,SINB,COSB,COSA,SINA,COST,SINT,
     $MG,MT1,MT2,MB1,MB2,GLT,GRT,XM,YM,THX,THY,
     $FT1(4),FT2(4),GT1(4),GT2(4),MZI(4),FT,FB,
     $RSIGT,QS,SSRSGT,ZBW(2),BWP(2),FWB1(2),GWB1(2),FWB2(2),
     $GWB2(2),MWI(2),MW,MZ,MB,MT,SUALFS,COS2W,FAC,
     $COSBE,SINBE,XM3,ST2LP,AT,CF,CA
      REAL*8 LP(11),LPTOT
      REAL*8 REAL8,RSIGWI,RSIGZI

      INTEGER THZ(4),THW(2),I
C
      DATA ZI/(0.,1.)/
C     Recompute weak scale Yukawa couplings including SUSY loops
C     Follow formulae of Pierce et al. NPB491, 3 (1997)
C
      REAL8(ZZZ)=DREAL(ZZZ)
      CA=3.
      CF=4./3.
      SR2=SQRT(2.)
      PI=4*ATAN(1.)
      E=SQRT(4*PI*ALFAEM)
      COS2W=1.-SN2THW
      G=G2
      MW=AMW
      MZ=AMZ
      MB=AMBT
      FAC=16*PI**2
      TANB=VUQ/VDQ
      BETA=ATAN(TANB)
      SINBE=SIN(BETA)
      COSBE=COS(BETA)
      COSA=COS(ALFAH)
      SINA=SIN(ALFAH)
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      MG=ABS(GSS(9))
      XM3=M3Q
      MT=AMTP
      MT1=MSS(12)
      MT2=MSS(13)
      MB1=MSS(10)
      MB2=MSS(11)
      GLT=.5-2*XW/3.
      GRT=2*XW/3.
      FT=MTQ/VUQ
      FB=MBQ/VDQ
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      THW(1)=0
      IF(AMW1SS.LT.0.) THW(1)=1
      THW(2)=0
      IF(AMW2SS.LT.0.) THW(2)=1
 
C     Neutralino-fermion-sfermion couplings
C     a(Pierce)=-beta^* ; b(Pierce)=-alpha^*
C
      SIGZI=(0.,0.)
      DO I=1,4
        THZ(I)=0
        IF (AMZISS(I).LT.0.) THZ(I)=1
        MZI(I)=ABS(AMZISS(I))
        ZAZT1(I)=ZI*(-ZI)**(THZ(I)-1)*(G*ZMIXSS(3,I)+GP*
     $ZMIXSS(4,I)/3.)/SR2*COST+(-ZI)**(THZ(I))*FT*ZMIXSS(1,I)*SINT
        ZAZT2(I)=ZI*(-ZI)**(THZ(I)-1)*(G*ZMIXSS(3,I)+GP*
     $ZMIXSS(4,I)/3.)/SR2*SINT-(-ZI)**(THZ(I))*FT*ZMIXSS(1,I)*COST
        ZBZT1(I)=-ZI*(ZI)**(THZ(I)-1)*4*GP/3./SR2*ZMIXSS(4,I)*SINT
     $-(ZI)**(THZ(I))*FT*ZMIXSS(1,I)*COST
        ZBZT2(I)=ZI*(ZI)**(THZ(I)-1)*4*GP/3./SR2*ZMIXSS(4,I)*COST
     $-(ZI)**(THZ(I))*FT*ZMIXSS(1,I)*SINT
        FT1(I)=ZAZT1(I)*CONJG(ZAZT1(I))+ZBZT1(I)*CONJG(ZBZT1(I))
        GT1(I)=2*REAL(ZAZT1(I)*CONJG(ZBZT1(I)))
        FT2(I)=ZAZT2(I)*CONJG(ZAZT2(I))+ZBZT2(I)*CONJG(ZBZT2(I))
        GT2(I)=2*REAL(ZAZT2(I)*CONJG(ZBZT2(I)))
        SIGZI=SIGZI+FT1(I)*SSB1F(QS,MZI(I),MT1)+GT1(I)*MZI(I)/MTQ*
     $SSB0(QS,MZI(I),MT1)+FT2(I)*SSB1F(QS,MZI(I),MT2)
     $+GT2(I)*MZI(I)/MTQ*SSB0(QS,MZI(I),MT2)
      END DO
      RSIGZI=REAL8(SIGZI)
C
C     Chargino-fermion-sfermion couplings; I labels chargino
C
      ZAWU(1)=ZI*G*SIN(GAMMAL)
      ZAWU(2)=ZI*THX*G*COS(GAMMAL)
      ZBW(1)=-(-1)**THW(1)*FT*COS(GAMMAR)
      ZBW(2)=(-1)**THW(2)*THY*FT*SIN(GAMMAR)
      BWP(1)=-FB*COS(GAMMAL)
      BWP(2)=FB*THX*SIN(GAMMAL)
      SIGWI=(0.,0.)
      MWI(1)=ABS(AMW1SS)
      MWI(2)=ABS(AMW2SS)
      DO I=1,2
        ZAWB1(I)=-ZI*ZAWU(I)*COSB+BWP(I)*SINB
        ZAWB2(I)=-ZI*ZAWU(I)*SINB-BWP(I)*COSB
        ZBWB1(I)=-ZBW(I)*COSB
        ZBWB2(I)=-ZBW(I)*SINB
        FWB1(I)=ZAWB1(I)*CONJG(ZAWB1(I))+ZBWB1(I)*CONJG(ZBWB1(I))
        GWB1(I)=2*REAL(CONJG(ZBWB1(I))*ZAWB1(I))
        FWB2(I)=ZAWB2(I)*CONJG(ZAWB2(I))+ZBWB2(I)*CONJG(ZBWB2(I))
        GWB2(I)=2*REAL(CONJG(ZBWB2(I))*ZAWB2(I))
        SIGWI=SIGWI+FWB1(I)*SSB1F(QS,MWI(I),MB1)+
     $GWB1(I)*MWI(I)/MTQ*SSB0(QS,MWI(I),MB1)+
     $FWB2(I)*SSB1F(QS,MWI(I),MB2)+
     $MWI(I)/MTQ*GWB2(I)*SSB0(QS,MWI(I),MB2)
      END DO
      RSIGWI=REAL8(SIGWI)
      LP(1)=ASMSS/3./PI*(REAL8(SSB1F(QS,MG,MT1))+
     $  REAL8(SSB1F(QS,MG,MT2))-SIN(2*THETAT)*XM3/MTQ*
     $  (REAL8(SSB0(QS,MG,MT1))-REAL8(SSB0(QS,MG,MT2))))
      LP(2)=.5*FT**2*SINA**2*(REAL8(SSB1F(QS,MT,AMHH))+
     $REAL8(SSB0(QS,MT,AMHH)))/FAC
      LP(3)=.5*FT**2*COSA**2*(REAL8(SSB1(QS,MT,AMHL))+
     $REAL8(SSB0(QS,MT,AMHL)))/FAC
      LP(4)=.5*FT**2*COSBE**2*(REAL8(SSB1F(QS,MT,AMHA))-
     $REAL8(SSB0(QS,MT,AMHA)))/FAC
      LP(5)=.5*FT**2*SINBE**2*(REAL8(SSB1(QS,MT,MZ))-
     $REAL8(SSB0(QS,MT,MZ)))/FAC
      LP(6)=.5*((FB**2*SINBE**2+FT**2*COSBE**2)*
     $REAL8(SSB1F(QS,MB,AMHC))+(G2**2+FB**2*COSBE**2+FT**2*SINBE**2)
     $*REAL8(SSB1(QS,MB,MW)))/FAC
      LP(7)=FB**2*COSBE**2*(REAL8(SSB0(QS,MB,AMHC))-
     $REAL8(SSB0(QS,MB,MW)))/FAC
      LP(8)=-(2*E/3.)**2*(5.+3*LOG(QS/MT**2))/FAC
      LP(9)=+G2**2/COS2W*((GLT**2+GRT**2)*REAL8(SSB1(QS,MT,MZ))
     $+4*GLT*GRT*REAL8(SSB0(QS,MT,MZ)))/FAC
      LP(10)=.5*RSIGZI/FAC
      LP(11)=.5*RSIGWI/FAC
      LPTOT=0.D0
      DO I=1,11
        LPTOT=LPTOT+LP(I)
C        WRITE(*,*) 'LP(',I,')=',LP(I)
      END DO
      AT=AAT-MU/TANB
      ST2LP=CF*(ASMSS/4./PI)**2*(47./3.+CF*23./24.+CA*175./72.
     $-4*AT/MSUSY+CF*AT/MSUSY*(7./3.+6*LOG(MT/MSUSY))
     $-8*CA*AT/3./MSUSY)
      SSRSGT=LPTOT
      SSRSGT=SSRSGT-ST2LP
100   RETURN
      END
