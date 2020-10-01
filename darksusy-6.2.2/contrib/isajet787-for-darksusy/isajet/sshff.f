#include "PILOT.inc"
      SUBROUTINE SSHFF(IMODEL)
C-----------------------------------------------------------------------
C
C     Calculate all decays higgs -> f fbar, including QCD radiative
C     corrections for quarks.
C
C     Bisset's SETFAC, WDHFFN, QCDRAD
C
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sssm.inc"
#include "sspar.inc"
#include "sstype.inc"
#include "sugpas.inc"
#include "kphggs.inc"
C
      DOUBLE PRECISION PI,SR2,DWID,XG2,MHIH,BETA,BEFAC,ALFAC,MH,MF
     $,MFRUN,FACTOR,ALAM,MF1,MF2,SUM,MF1RUN,MF2RUN,COLOR,TEMP1
     $,QCDFAC
      DOUBLE PRECISION MFIFF(9),MFIF1(6),MFIF2(6)
      DOUBLE PRECISION SSDLAM,SSMQCD,SSHFF1
      REAL WID,COSB,MBMA,TANB,COTB,XI,DELB,DELT,DELL,XMU
      REAL FT,FB,BTA,A1,B1,C1
      INTEGER IH,IDIH,IFF,IDF,ID1,ID2
      INTEGER IDIFF(9),IDIF1(6),IDIF2(6),IMODEL
C
C     Define loop integral for kappas
      XI(A1,B1,C1)=(A1*B1*LOG(A1/B1)+B1*C1*LOG(B1/C1)+C1*A1*LOG(C1/A1))/
     ,((A1-B1)*(B1-C1)*(A1-C1))
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      BETA=ATAN(1./RV2V1)
      BTA=SNGL(BETA)
      COSB=COS(BTA)
      XG2=4.0*PI*ALFAEM/SN2THW
      MBMA=FBMA*COSB*VEV
C
C          Loop over HL, HH, HA and fermions
C
      MFIFF(1)=AME
      IDIFF(1)=IDE
      MFIFF(2)=AMMU
      IDIFF(2)=IDMU
      MFIFF(3)=AMTAU
      IDIFF(3)=IDTAU
      MFIFF(4)=AMDN
      IDIFF(4)=IDDN
      MFIFF(5)=AMST
      IDIFF(5)=IDST
      MFIFF(6)=AMBT
      IDIFF(6)=IDBT
      MFIFF(7)=AMUP
      IDIFF(7)=IDUP
      MFIFF(8)=AMCH
      IDIFF(8)=IDCH
      MFIFF(9)=AMTP
      IDIFF(9)=IDTP
C
      DO 100 IH=1,3
        IF(IH.EQ.1) THEN
          MH=AMHL
          IDIH=ISHL
          BEFAC=COS(BETA)
          ALFAC=SIN(ALFAH)
        ELSEIF(IH.EQ.2) THEN
          MH=AMHH
          IDIH=ISHH
          BEFAC=COS(BETA)
          ALFAC=COS(ALFAH)
        ELSE
          MH=AMHA
          IDIH=ISHA
          BEFAC=1/TAN(BETA)
          ALFAC=1.
        ENDIF
C
C          Down type fermions
C
        DO 110 IFF=1,6
          MF=MFIFF(IFF)
          IDF=IDIFF(IFF)
          FACTOR=1.-4.*MF**2/MH**2
          IF(FACTOR.LE.0) GO TO 110
          FACTOR=SQRT(FACTOR)
          IF(IFF.GE.4) THEN
             COLOR=3.
             MFRUN=SSMQCD(MF,MH)
             IF (IFF.EQ.6.AND.IH.EQ.1) THEN
C              USE ZERWAS MSBAR VALUE MB(100GeV)=2.92
               MFRUN=2.92
             ELSE IF (IFF.EQ.6.AND.IH.GE.2) THEN
                 MFRUN=MBMA
             END IF
             QCDFAC=SSHFF1(MH,MF,IH)
          ELSE
             COLOR=1.
             MFRUN=MF
             QCDFAC=1.
          ENDIF
          DWID=XG2*MFRUN**2*MH*ALFAC**2/(32.*PI*AMW**2*BEFAC**2)
          IF(IH.EQ.1.OR.IH.EQ.2) THEN
            DWID=DWID*FACTOR**3
          ELSEIF(IH.EQ.3) THEN
            DWID=DWID*FACTOR
          ENDIF
          DWID=DWID*COLOR*QCDFAC
          WID=DWID
          CALL SSSAVE(IDIH,WID,IDF,-IDF,0,0,0)
110     CONTINUE
C
C          Up type fermions
C
        IF(IH.EQ.1) THEN
          BEFAC=SIN(BETA)
          ALFAC=COS(ALFAH)
        ELSEIF(IH.EQ.2) THEN
          BEFAC=SIN(BETA)
          ALFAC=SIN(ALFAH)
        ELSE
          BEFAC=TAN(BETA)
          ALFAC=1.
        ENDIF
        DO 120 IFF=7,9
          MF=MFIFF(IFF)
          IDF=IDIFF(IFF)
          FACTOR=1.-4.*MF**2/MH**2
          IF(FACTOR.LE.0) GO TO 120
          FACTOR=SQRT(FACTOR)
          IF (IH.EQ.1.AND.IDF.EQ.8) THEN
C                USE ZERWAS MSBAR VALUE MC(100GeV)=0.62
                 MFRUN=0.62
          END IF
          MFRUN=SSMQCD(MF,MH)
          QCDFAC=SSHFF1(MH,MF,IH)
          DWID=XG2*MFRUN**2*MH*ALFAC**2/(32.*PI*AMW**2*BEFAC**2)
          IF(IH.EQ.1.OR.IH.EQ.2) THEN
            DWID=DWID*FACTOR**3
          ELSEIF(IH.EQ.3) THEN
            DWID=DWID*FACTOR
          ENDIF
          DWID=3.*DWID*QCDFAC
          WID=DWID
          CALL SSSAVE(IDIH,WID,IDF,-IDF,0,0,0)
120     CONTINUE
100   CONTINUE
C
C           HC decays. F1 has Iz=+1/2, F2 has Iz=-1/2
C
      MFIF1(1)=0
      IDIF1(1)=IDNE
      MFIF2(1)=AME
      IDIF2(1)=IDE
      MFIF1(2)=0
      IDIF1(2)=IDNM
      MFIF2(2)=AMMU
      IDIF2(2)=IDMU
      MFIF1(3)=0
      IDIF1(3)=IDNT
      MFIF2(3)=AMTAU
      IDIF2(3)=IDTAU
      MFIF1(4)=AMUP
      IDIF1(4)=IDUP
      MFIF2(4)=AMDN
      IDIF2(4)=IDDN
      MFIF1(5)=AMCH
      IDIF1(5)=IDCH
      MFIF2(5)=AMST
      IDIF2(5)=IDST
      MFIF1(6)=AMTP
      IDIF1(6)=IDTP
      MFIF2(6)=AMBT
      IDIF2(6)=IDBT
      MH=AMHC
C
      DO 200 IFF=1,6
        MF1=MFIF1(IFF)
        MF2=MFIF2(IFF)
        ID1=IDIF1(IFF)
        ID2=IDIF2(IFF)
        SUM=MF1+MF2
        ALAM=SSDLAM(MH**2,MF1**2,MF2**2)
        IF(ALAM.LE.0.OR.SUM.GE.MH) GO TO 200
        IF(IFF.LE.3) THEN
          MF1RUN=MF1
          MF2RUN=MF2
          COLOR=1
        ELSE
          MF1RUN=SSMQCD(MF1,MH)
          IF (MF2.EQ.AMBT.AND.MBMA.NE.0.) THEN
            MF2RUN=MBMA
          ELSE
            MF2RUN=SSMQCD(MF2,MH)
          END IF
          COLOR=3
        ENDIF
        TEMP1=MF1RUN**2*1./TAN(BETA)**2+MF2RUN**2*TAN(BETA)**2
        TEMP1=TEMP1*(MH**2-MF1**2-MF2**2)-4.*MF1**2*MF2**2
        IF (TEMP1.LT.0.0) GO TO 200
        DWID=XG2*COLOR*SQRT(ALAM)*TEMP1/MH**3/(32.0*PI*AMW**2)
        WID=DWID
        CALL SSSAVE(ISHC,WID,ID1,-ID2,0,0,0)
200   CONTINUE
C
C     KAPPA CALCULATION
C
      XMU=-TWOM1
      FT=MTQ/VUQ
      FB=MBQ/VDQ
      COTB=RV2V1
      TANB=1./COTB
      DELB=(2*ALFA3*AMGLSS*XMU/3./PI*XI(AMB1SS**2,AMB2SS**2,AMGLSS**2)+
     ,FT**2*XMU*AAT*XI(AMT1SS**2,AMT2SS**2,XMU**2)/16./PI/PI)*TANB
      DELT=(2*ALFA3*AMGLSS*XMU/3./PI*XI(AMT1SS**2,AMT2SS**2,AMGLSS**2)+
     ,FB**2*XMU*AAB*XI(AMB1SS**2,AMB2SS**2,XMU**2)/16./PI/PI)*COTB
      DELL=-3*ALFA2*ABS(AMW2SS)*XMU*TANB*
     ,XI(AMLLSS**2,AMW2SS**2,XMU**2)/8./PI
      KPB=SIN(BTA+ALFAH)-(TANB-DELB*COTB)*COS(BTA+ALFAH)/(1.+DELB)
      KPL=SIN(BTA+ALFAH)-(TANB-DELL*COTB)*COS(BTA+ALFAH)/(1.+DELL)
      KPT=SIN(BTA+ALFAH)+(COTB-DELT*TANB)*COS(BTA+ALFAH)/(1.+DELT)
      RETURN
      END
