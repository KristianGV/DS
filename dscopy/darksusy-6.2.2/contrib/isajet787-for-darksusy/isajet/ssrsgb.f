#include "PILOT.inc"
C--------------------------------------------------------------------
      FUNCTION SSRSGB(QS)
C--------------------------------------------------------------------
C
C     Calculate b quark self energy
C     according to Pierce et al. formulae adapted to Isajet
C
C     Modified by Javier 9/2005 / Log threholds already implemented 
C     through RGE decoupling have been substracted by 
C     a redefinition of B1 function
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sspar.inc"
#include "sssm.inc"
#include "sugmg.inc"
#include "sugpas.inc"
      COMPLEX ZAZB1(4),ZAZB2(4),ZBZB1(4),ZBZB2(4),ZI
      COMPLEX ZAWD(2)
      COMPLEX ZAWT1(2),ZAWT2(2),ZBWT1(2),ZBWT2(2)
      COMPLEX*16 SSB0,SSB1F,ZZZ,SIGWI,SIGZI
      REAL SR2,PI,E,G,TANB,BETA,SINB,COSB,COSA,SINA,COST,SINT,
     $MG,MT1,MT2,MB1,MB2,GLB,GRB,XM,YM,THX,THY,
     $FB1(4),FB2(4),GB1(4),GB2(4),MZI(4),FT,FB,
     $RSIGB,QS,SSRSGB,ZBW(2),BWP(2),FWT1(2),GWT1(2),FWT2(2),
     $GWT2(2),MWI(2),MW,MZ,MB,MT,SUALFS,COS2W,FAC,
     $COSBE,SINBE,XM3
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
      MB=MBQ
C      ASMZ=SUALFS(AMZ**2,.36,AMTP,3)
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
      GLB=-.5+XW/3.
      GRB=-XW/3.
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
        ZAZB1(I)=ZI*(-ZI)**(THZ(I)-1)*(-G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I)/3.)/SR2*COSB+(-ZI)**(THZ(I))*FB*ZMIXSS(2,I)*SINB
        ZAZB2(I)=ZI*(-ZI)**(THZ(I)-1)*(-G*ZMIXSS(3,I)+GP*
     $  ZMIXSS(4,I)/3.)/SR2*SINB-(-ZI)**(THZ(I))*FB*ZMIXSS(2,I)*COSB
        ZBZB1(I)=-ZI*(ZI)**(THZ(I)-1)*(-2)*GP/3./SR2*ZMIXSS(4,I)*SINB
     $  -(ZI)**(THZ(I))*FB*ZMIXSS(2,I)*COSB
        ZBZB2(I)=ZI*(ZI)**(THZ(I)-1)*(-2)*GP/3./SR2*ZMIXSS(4,I)*COSB
     $  -(ZI)**(THZ(I))*FB*ZMIXSS(2,I)*SINB
        FB1(I)=ZAZB1(I)*CONJG(ZAZB1(I))+ZBZB1(I)*CONJG(ZBZB1(I))
        GB1(I)=2*REAL(ZAZB1(I)*CONJG(ZBZB1(I)))
        FB2(I)=ZAZB2(I)*CONJG(ZAZB2(I))+ZBZB2(I)*CONJG(ZBZB2(I))
        GB2(I)=2*REAL(ZAZB2(I)*CONJG(ZBZB2(I)))
        SIGZI=SIGZI+FB1(I)*SSB1F(QS,MZI(I),MB1)+GB1(I)*MZI(I)/MBQ*
     $  SSB0(QS,MZI(I),MB1)+FB2(I)*SSB1F(QS,MZI(I),MB2)
     $  +GB2(I)*MZI(I)/MBQ*SSB0(QS,MZI(I),MB2)
      END DO
      RSIGZI=REAL8(SIGZI)
C
C     Chargino-fermion-sfermion couplings; I labels chargino
C
      ZAWD(1)=ZI*(-1)**THW(1)*G*SIN(GAMMAR)
      ZAWD(2)=ZI*(-1)**THW(2)*THY*G*COS(GAMMAR)
      ZBW(1)=-(-1)**THW(1)*FT*COS(GAMMAR)
      ZBW(2)=(-1)**THW(2)*THY*FT*SIN(GAMMAR)
      BWP(1)=-FB*COS(GAMMAL)
      BWP(2)=FB*THX*SIN(GAMMAL)
      SIGWI=(0.,0.)
      MWI(1)=ABS(AMW1SS)
      MWI(2)=ABS(AMW2SS)
      DO I=1,2
        ZAWT1(I)=-ZI*ZAWD(I)*COST+ZBW(I)*SINT
        ZAWT2(I)=-ZI*ZAWD(I)*SINT-ZBW(I)*COST
        ZBWT1(I)=-BWP(I)*COST
        ZBWT2(I)=-BWP(I)*SINT
        FWT1(I)=ZAWT1(I)*CONJG(ZAWT1(I))+ZBWT1(I)*CONJG(ZBWT1(I))
        GWT1(I)=2*REAL(CONJG(ZBWT1(I))*ZAWT1(I))
        FWT2(I)=ZAWT2(I)*CONJG(ZAWT2(I))+ZBWT2(I)*CONJG(ZBWT2(I))
        GWT2(I)=2*REAL(CONJG(ZBWT2(I))*ZAWT2(I))
        SIGWI=SIGWI+FWT1(I)*SSB1F(QS,MWI(I),MT1)+
     $MWI(I)/MBQ*GWT1(I)*SSB0(QS,MWI(I),MT1)+
     $FWT2(I)*SSB1F(QS,MWI(I),MT2)+
     $MWI(I)/MBQ*GWT2(I)*SSB0(QS,MWI(I),MT2)
      END DO
      RSIGWI=REAL8(SIGWI)
      LP(1)=ASMSS/3./PI*(REAL8(SSB1F(QS,MG,MB1))+
     $  REAL8(SSB1F(QS,MG,MB2))-SIN(2*THETAB)*XM3/MBQ*
     $  (REAL8(SSB0(QS,MG,MB1))-REAL8(SSB0(QS,MG,MB2))))
      LP(2)=.5*FB**2*COSA**2*(REAL8(SSB1F(QS,MB,AMHH))+
     $REAL8(SSB0(QS,MB,AMHH)))/FAC
      LP(3)=.5*FB**2*SINA**2*(REAL8(SSB1F(QS,MB,AMHL))+
     $REAL8(SSB0(QS,MB,AMHL)))/FAC
      LP(4)=.5*FB**2*SINBE**2*(REAL8(SSB1F(QS,MB,AMHA))-
     $REAL8(SSB0(QS,MB,AMHA)))/FAC
      LP(5)=.5*FB**2*COSBE**2*(REAL8(SSB1F(QS,MB,MZ))-
     $REAL8(SSB0(QS,MB,MZ)))/FAC
      LP(6)=.5*((FT**2*COSBE**2+FB**2*SINBE**2)*
     $REAL8(SSB1F(QS,MT,AMHC))+(G2**2+FT**2*SINBE**2+FB**2*COSBE**2)
     $*REAL8(SSB1F(QS,MT,MW)))/FAC
      LP(7)=FT**2*SINBE**2*(REAL8(SSB0(QS,MT,AMHC))-
     $REAL8(SSB0(QS,MT,MW)))/FAC
      LP(8)=0.
      LP(9)=+G2**2/COS2W*((GLB**2+GRB**2)*REAL8(SSB1F(QS,MB,MZ))
     $+4*GLB*GRB*REAL8(SSB0(QS,MB,MZ)))/FAC
      LP(10)=.5*RSIGZI/FAC
      LP(11)=.5*RSIGWI/FAC
      LPTOT=0.D0
      DO I=1,11
        LPTOT=LPTOT+LP(I)
C        WRITE(6,*) 'LP(',I,')=',LP(I)
      END DO
      SSRSGB=LPTOT
100   RETURN
      END
