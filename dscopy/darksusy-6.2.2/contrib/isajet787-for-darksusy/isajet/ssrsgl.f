#include "PILOT.inc"
C--------------------------------------------------------------------
      FUNCTION SSRSGL(QS)
C--------------------------------------------------------------------
C
C     Calculate tau lepton self energy
C     according to Pierce et al. formulae adapted to Isajet
C
C     Modified by Javier 9/2005 /Log threholds already 
C     implemented through RGE decoupling have been substracted 
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
      COMPLEX ZAZL1(4),ZAZL2(4),ZBZL1(4),ZBZL2(4),ZI
      COMPLEX ZAWL(2)
      COMPLEX ZAWSN(2),ZBWSN(2)
      COMPLEX*16 SSB0,SSB1F,ZZZ,SIGWI,SIGZI
      REAL SR2,PI,E,G,TANB,BETA,SINB,COSB,COSA,SINA,COSL,SINL,
     $MG,ML1,ML2,MSN,GLL,GRL,XM,YM,THX,THY,
     $FL1(4),FL2(4),GL1(4),GL2(4),MZI(4),FL,ML,
     $RSIGL,QS,SSRSGL,BWPP(2),FWSN(2),GWSN(2),
     $MWI(2),MW,MZ,COS2W,FAC,
     $COSBE,SINBE
      REAL*8 LP(11),LPTOT
      REAL*8 REAL8,RSIGWI,RSIGZI

      INTEGER THZ(4),THW(2),I
C
      DATA ZI/(0.,1.)/
C     Recompute weak scale Yukawa couplings including SUSY loops
C     Follow formulae of Pierce et al. NPB491, 3 (1997)
C
      REAL8(ZZZ)=DREAL(ZZZ)
      SR2=SQRT(2.)
      PI=4*ATAN(1.)
      E=SQRT(4*PI*ALFAEM)
      COS2W=1.-SN2THW
      G=G2
      MW=AMW
      MZ=AMZ
      ML=MLQ
      FAC=16*PI**2
      TANB=VUQ/VDQ
      BETA=ATAN(TANB)
      SINBE=SIN(BETA)
      COSBE=COS(BETA)
      COSA=COS(ALFAH)
      SINA=SIN(ALFAH)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      MSN=AMN3SS
      ML1=AML1SS
      ML2=AML2SS
      GLL=-.5+XW
      GRL=-XW
      FL=MLQ/VDQ
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
        ZAZL1(I)=ZI*(-1)*(-ZI)**(THZ(I)-1)*(G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I))/SR2*COSL+(-ZI)**(THZ(I))*FL*ZMIXSS(2,I)*SINL
        ZAZL2(I)=ZI*(-1)*(-ZI)**(THZ(I)-1)*(G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I))/SR2*SINL-(-ZI)**(THZ(I))*FL*ZMIXSS(2,I)*COSL
        ZBZL1(I)=-ZI*(-1)*(ZI)**(THZ(I)-1)*SR2*GP*ZMIXSS(4,I)*SINL
     $  -(ZI)**(THZ(I))*FL*ZMIXSS(2,I)*COSL
        ZBZL2(I)=ZI*(-1)*(ZI)**(THZ(I)-1)*SR2*GP*ZMIXSS(4,I)*COSL
     $  -(ZI)**(THZ(I))*FL*ZMIXSS(2,I)*SINL
        FL1(I)=ZAZL1(I)*CONJG(ZAZL1(I))+ZBZL1(I)*CONJG(ZBZL1(I))
        GL1(I)=2*REAL(ZAZL1(I)*CONJG(ZBZL1(I)))
        FL2(I)=ZAZL2(I)*CONJG(ZAZL2(I))+ZBZL2(I)*CONJG(ZBZL2(I))
        GL2(I)=2*REAL(ZAZL2(I)*CONJG(ZBZL2(I)))
        SIGZI=SIGZI+FL1(I)*SSB1F(QS,MZI(I),ML1)+GL1(I)*MZI(I)/ML*
     $  SSB0(QS,MZI(I),ML1)+FL2(I)*SSB1F(QS,MZI(I),ML2)
     $  +GL2(I)*MZI(I)/ML*SSB0(QS,MZI(I),ML2)
      END DO
      RSIGZI=REAL8(SIGZI)
C
C     Chargino-fermion-sfermion couplings; I labels chargino
C
      ZAWL(1)=ZI*(-1)**THW(1)*G*SIN(GAMMAR)
      ZAWL(2)=ZI*(-1)**THW(2)*THY*G*COS(GAMMAR)
      BWPP(1)=-FL*COS(GAMMAL)
      BWPP(2)=FL*THX*SIN(GAMMAL)
      SIGWI=(0.,0.)
      MWI(1)=ABS(AMW1SS)
      MWI(2)=ABS(AMW2SS)
      DO I=1,2
        ZAWSN(I)=-ZI*ZAWL(I)
        ZBWSN(I)=-BWPP(I)
        FWSN(I)=ZAWSN(I)*CONJG(ZAWSN(I))+ZBWSN(I)*CONJG(ZBWSN(I))
        GWSN(I)=2*REAL(CONJG(ZBWSN(I))*ZAWSN(I))
        SIGWI=SIGWI+FWSN(I)*SSB1F(QS,MWI(I),MSN)+
     $  MWI(I)/MLQ*GWSN(I)*SSB0(QS,MWI(I),MSN)
      END DO
      RSIGWI=REAL8(SIGWI)
      LP(1)=0.
      LP(2)=.5*FL**2*COSA**2*(REAL8(SSB1F(QS,ML,AMHH))+
     $REAL8(SSB0(QS,ML,AMHH)))/FAC
      LP(3)=.5*FL**2*SINA**2*(REAL8(SSB1F(QS,ML,AMHL))+
     $REAL8(SSB0(QS,ML,AMHL)))/FAC
      LP(4)=.5*FL**2*SINBE**2*(REAL8(SSB1F(QS,ML,AMHA))-
     $REAL8(SSB0(QS,ML,AMHA)))/FAC
      LP(5)=.5*FL**2*COSBE**2*(REAL8(SSB1F(QS,ML,MZ))-
     $REAL8(SSB0(QS,ML,MZ)))/FAC
      LP(6)=.5*(FL**2*SINBE**2*REAL8(SSB1F(QS,0.,AMHC))+
     $(G2**2+FL**2*COSBE**2)*REAL8(SSB1F(QS,0.,MW)))/FAC
      LP(7)=0.
      LP(8)=0.
      LP(9)=+G2**2/COS2W*((GLL**2+GRL**2)*REAL8(SSB1F(QS,ML,MZ))
     $+4*GLL*GRL*REAL8(SSB0(QS,ML,MZ)))/FAC
      LP(10)=.5*RSIGZI/FAC
      LP(11)=.5*RSIGWI/FAC
      LPTOT=0.D0
      DO I=1,11
        LPTOT=LPTOT+LP(I)
C        WRITE(6,*) 'LP(',I,')=',LP(I)
      END DO
      SSRSGL=LPTOT
100   RETURN
      END
