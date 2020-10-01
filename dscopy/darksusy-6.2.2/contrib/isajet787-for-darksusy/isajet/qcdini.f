#include "PILOT.inc"
      SUBROUTINE QCDINI(JIN1,JIN2)
C
C          GENERATE INITIAL-STATE QCD CASCADE USING BACKWARDS
C          EVOLUTION OF GOTTSCHALK AND OF SJOSTRAND.
C
C          IF QCDINI FAILS WHEN ATTEMPTING TO FORCE GL-->QK+QB FOR
C          HEAVY QUARKS, THEN RETURN NJSET=-1.
C
C          VER. 6.40: TRAP W1LIM > 0 TO PREVENT ROUNDING ERRORS.
C
#include "itapes.inc"
#include "idrun.inc"
#include "pinits.inc"
#include "jetpar.inc"
#include "qcdpar.inc"
#include "jetset.inc"
#include "jwork.inc"
#include "jwork2.inc"
#include "const.inc"
#include "primar.inc"
#include "keys.inc"
C
      DIMENSION BOOST1(5),BOOST2(5),B2B1(5),DBL1(5),DBL2(5)
      DIMENSION FXOLD(2),FXNEW(2)
      DIMENSION PJKEEP(5,12),JINS(2),JLIST(16),PFKEEP(5)
#ifdef DOUBLE_X
      DOUBLE PRECISION DBL1,DBL2,DBLM
#endif
C
C          CONVERT IDENT+7 TO JETTYP
      DATA JLIST/13,11,9,7,5,3,0,2,4,6,8,10,12,0,0,1/
      ALAMF(A,B,C)=SQRT((A-B-C)**2-4.*B*C)
C
C          INITIALIZE
C
      JINS(1)=JIN1
      JINS(2)=JIN2
      DO 97 K=1,4
97    PFKEEP(K)=PJSET(K,JIN1)+PJSET(K,JIN2)
C          EXCEPT FOR HIGGS, PFKEEP**2=SHAT
      IF(KEYS(7).OR.KEYS(9)) THEN
        S1KEEP=PFKEEP(4)**2-PFKEEP(1)**2-PFKEEP(2)**2-PFKEEP(3)**2
        PFKEEP(5)=SQRT(S1KEEP)
        PPKEEP=PFKEEP(4)+PFKEEP(3)
        PMKEEP=PFKEEP(4)-PFKEEP(3)
      ELSE
        S1KEEP=SHAT
        PFKEEP(5)=SQRT(S1KEEP)
        IF(PFKEEP(3).GT.0.) THEN
          PPKEEP=PFKEEP(4)+PFKEEP(3)
          PMKEEP=(S1KEEP+PFKEEP(1)**2+PFKEEP(2)**2)/PPKEEP
        ELSE
          PMKEEP=PFKEEP(4)-PFKEEP(3)
          PPKEEP=(S1KEEP+PFKEEP(1)**2+PFKEEP(2)**2)/PMKEEP
        ENDIF
        PFKEEP(4)=.5*(PPKEEP+PMKEEP)
        PFKEEP(3)=.5*(PPKEEP-PMKEEP)
      ENDIF
      DO 98 I=1,NJSET
      DO 98 K=1,5
98    PJKEEP(K,I)=PJSET(K,I)
      NJKEEP=NJSET
      NPASS=0
      NPASS1=0
C
1     CONTINUE
      NPASS1=NPASS1+1
      IF(NPASS1.GT.100) GO TO 9999
      NJSET=NJKEEP
      DO 99 I=1,NJSET
      DO 99 K=1,5
99    PJSET(K,I)=PJKEEP(K,I)
C
      DO 100 K=1,5
100   PFINAL(K)=PFKEEP(K)
      S1=S1KEEP
      PTOTPL=PPKEEP
      PTOTMN=PMKEEP
      TCUT=CUTJET**2
      DO 101 I=1,2
      JI=JINS(I)
      XOLD=(PJSET(4,JI)+ABS(PJSET(3,JI)))/ECM
      JT=JLIST(JTYPE(JI)+7)
      FXOLD(I)=STRUC(XOLD,QSQ,JT,IDIN(I))
101   CONTINUE
C
C          DO FIRST EVOLUTION
      DO 110 I=1,2
      SGN=3-2*I
      JET=10+I
      JI=JINS(I)
      ZMIN=(PJSET(4,JI)+ABS(PJSET(3,JI)))/ECM
      ZMAX=1./(1.+TCUT/S1)
C          DZMAX=1.-ZMAX
      DZMAX=ZMAX*TCUT/S1
      IF(ZMIN.GE.ZMAX) ZMIN=.5*ZMAX
      CALL QCDINT(JI)
      JVIR(I)=JI
110   CONTINUE
C
C          SOLVE INITIAL KINEMATICS
      AM1SQ=PJSET(5,JVIR(1))**2*SIGN(1.,PJSET(5,JVIR(1)))
      AM2SQ=PJSET(5,JVIR(2))**2*SIGN(1.,PJSET(5,JVIR(2)))
      P1PL=(S1+AM1SQ-AM2SQ+ALAMF(S1,AM1SQ,AM2SQ))/(2.*PTOTMN)
      P1MN=AM1SQ/P1PL
      P2MN=(S1+AM2SQ-AM1SQ+ALAMF(S1,AM1SQ,AM2SQ))/(2.*PTOTPL)
      P2PL=AM2SQ/P2MN
      PJSET(3,JVIR(1))=.5*(P1PL-P1MN)
      PJSET(4,JVIR(1))=.5*(P1PL+P1MN)
      PJSET(3,JVIR(2))=.5*(P2PL-P2MN)
      PJSET(4,JVIR(2))=.5*(P2PL+P2MN)
C
C          TEST WHETHER NEW MASS IS PLAUSIBLE
      DO 111 I=1,2
      JI=JINS(I)
      XNEW=(PJSET(4,JI)+ABS(PJSET(3,JI)))/ECM
      IF(XNEW.GE.1.) THEN
        FXNEW(I)=0.
      ELSE
        JT=JLIST(JTYPE(JI)+7)
        FXNEW(I)=STRUC(XNEW,QSQ,JT,IDIN(I))
      ENDIF
111   CONTINUE
      DO 112 I=1,2
      IF(FXNEW(I).LT.FXOLD(I)*RANF()) GO TO 1
112   CONTINUE
C
C          FIND JVIR (SPACE-LIKE PARTON) WITH LARGER (-MASS) FOR NEXT
C          BRANCHING.
10    IF(JDCAY(JVIR(1)).GE.0.AND.JDCAY(JVIR(2)).GE.0) RETURN
      NPASS=NPASS+1
      IF(NPASS.GT.20*NJSET) GO TO 9999
      IF(-PJSET(5,JVIR(1)).GE.-PJSET(5,JVIR(2))) THEN
        IVIR=JVIR(1)
        IVIR2=JVIR(2)
        SGN=+1.
        JET=11
      ELSE
        IVIR=JVIR(2)
        IVIR2=JVIR(1)
        SGN=-1.
        JET=12
      ENDIF
C
      T1=PJSET(5,IVIR)**2
      ZMIN=(PJSET(4,IVIR)+SGN*PJSET(3,IVIR))/ECM
      ZMAX=1./(1.+T1/S1)
      DZMAX=ZMAX*T1/S1
      IF(ZMIN.GE.ZMAX) GO TO 1
C
C          GENERATE Z AND NEW PARTONS.
C          NEWV=SPACELIKE, NEWF=TIMELIKE.
      NEWV=NJSET+1
      NEWF=NJSET+2
      CALL QCDINZ(IVIR)
C
C          IF Z FAILS (BECAUSE OF STRUCTURE FUNCTION) SET NEWV=IVIR,
C          NEWF=NULL AND RE-SOLVE KINEMATICS.
15    IF(.NOT.ZGOOD) THEN
        CALL QCDINT(IVIR)
C
        PP1PL=PJSET(4,IVIR2)+PJSET(3,IVIR2)
        PP1MN=PJSET(4,IVIR2)-PJSET(3,IVIR2)
        AMSQ=PJSET(5,IVIR)**2*SIGN(1.,PJSET(5,IVIR))
        AMPSQ=PJSET(5,IVIR2)**2*SIGN(1.,PJSET(5,IVIR2))
        IF(SGN.GT.0) THEN
          P2PL=(S1-AMSQ-AMPSQ+ALAMF(S1,AMSQ,AMPSQ))/(2.*PP1MN)
          P2MN=AMSQ/P2PL
        ELSE
          P2MN=(S1-AMSQ-AMPSQ+ALAMF(S1,AMSQ,AMPSQ))/(2.*PP1PL)
          P2PL=AMSQ/P2MN
        ENDIF
        PJSET(3,IVIR)=.5*(P2PL-P2MN)
        PJSET(4,IVIR)=.5*(P2PL+P2MN)
C
        NEWV=IVIR
        DO 120 K=1,5
120     PJSET(K,NEWF)=0.
        GO TO 30
      ENDIF
C
C          EVOLVE NEW SPACELIKE PARTON.
      PJSET(5,NEWV)=PJSET(5,IVIR)
      S2=S1/ZZC(IVIR)
      ZMIN=ZMIN/ZZC(IVIR)
      ZMAX=1./(1.+TCUT/S2)
      DZMAX=ZMAX*TCUT/S2
      IF(ZMIN.GE.ZMAX) GO TO 1
      CALL QCDINT(NEWV)
C
C          CALCULATE APPROXIMATE MASS LIMIT AND DO TIMELIKE EVOLUTION.
C          VER. 6.40: TRAP W1LIM < 0 FROM ROUNDING ERRORS.
      W1LIM=T1*(1./(ZZC(IVIR)*(1.+T1/S1))-1.)
      W1LIM=AMIN1(W1LIM,T1)
      PJSET(5,NEWF)=SQRT(ABS(W1LIM))
      JDCAY(NEWF)=-1
20    CALL QCDT(NEWF)
C
C          SOLVE KINEMATICS USING +(PL) AND -(MN) COMPONENTS FOR
C          PJSET(K,NEWV)+PJSET(K,IVIR2)-->PJSET(K,NEWF)+PFINAL
C          STEP 1: SOLVE FOR P2=PJSET(K,NEWV)
      PP1PL=PJSET(4,IVIR2)+PJSET(3,IVIR2)
      PP1MN=PJSET(4,IVIR2)-PJSET(3,IVIR2)
      AMSQ=PJSET(5,NEWV)**2*SIGN(1.,PJSET(5,NEWV))
      AMPSQ=PJSET(5,IVIR2)**2*SIGN(1.,PJSET(5,IVIR2))
      W1=PJSET(5,NEWF)**2
      IF(SGN.GT.0) THEN
        P2PL=(S2-AMSQ-AMPSQ+ALAMF(S2,AMSQ,AMPSQ))/(2.*PP1MN)
        P2MN=AMSQ/P2PL
      ELSE
        P2MN=(S2-AMSQ-AMPSQ+ALAMF(S2,AMSQ,AMPSQ))/(2.*PP1PL)
        P2PL=AMSQ/P2MN
      ENDIF
C
C          STEP 2: SOLVE FOR Q1(K)=PJSET(K,IVIR)
      DEN=P2PL*PP1MN-P2MN*PP1PL
      Q1PL=(+P2PL*(S1+T1-AMPSQ)+PP1PL*(W1+T1-AMSQ))/DEN
      Q1MN=(-P2MN*(S1+T1-AMPSQ)-PP1MN*(W1+T1-AMSQ))/DEN
      WPL=P2PL-Q1PL
      WMN=P2MN-Q1MN
C          CALCULATE TRANSVERSE MOMENTUM AND REJECT IF UNPHYSICAL.
      Q1TR2=T1+Q1PL*Q1MN
      IF(Q1TR2.LT.0.) THEN
        IF(JDCAY(NEWF).EQ.-1) GO TO 20
        ZGOOD=.FALSE.
        GO TO 15
      ENDIF
C
C          DO ONE TIMELIKE BRANCHING TO INSURE CORRECT MASS. MUST FIRST
C          SHIFT NJSET TO PUT DECAY PRODUCTS IN CORRECT PLACE.
      IF(JDCAY(NEWF).EQ.-1) THEN
        NJSET=NJSET+2
        CALL QCDZ(NEWF)
        NJSET=NJSET-2
        Z1=ZZC(NEWF)
        E0=.5*(WPL+WMN)
        P0=SQRT(.25*(WPL-WMN)**2+Q1TR2)
        WM0=PJSET(5,NEWF)
        ZLIM=AMAX1((WM0/(E0+P0))**2,CUTJET/(E0+P0))
        IF(Z1.LE.ZLIM.OR.Z1.GE.1.-ZLIM) GO TO 20
        NEWF1=NEWF+1
        NEWF2=NEWF+2
        JDCAY(NEWF)=NEWF1*JPACK+NEWF2
        CALL QCDT(NEWF1)
        CALL QCDT(NEWF2)
        JORIG(NEWF1)=JPACK*JET+NEWF
        JORIG(NEWF2)=JORIG(NEWF1)
        DO 130 K=1,4
        PJSET(K,NEWF1)=0.
130     PJSET(K,NEWF2)=0.
      ENDIF
C
C          GOOD BRANCHING!
      PHIQ1=2.*PI*RANF()
      Q1TR=SQRT(Q1TR2)

      Q1X=Q1TR*COS(PHIQ1)
      Q1Y=Q1TR*SIN(PHIQ1)
C
      PJSET(1,IVIR)=Q1X
      PJSET(2,IVIR)=Q1Y
      PJSET(3,IVIR)=.5*(Q1PL-Q1MN)
      PJSET(4,IVIR)=.5*(Q1PL+Q1MN)
      JDCAY(IVIR)=JPACK*NEWV+NEWF
C
      PJSET(1,NEWV)=0.
      PJSET(2,NEWV)=0.
      PJSET(3,NEWV)=.5*(P2PL-P2MN)
      PJSET(4,NEWV)=.5*(P2PL+P2MN)
      JORIG(NEWV)=JPACK*JET+IVIR
C
      PJSET(1,NEWF)=-Q1X
      PJSET(2,NEWF)=-Q1Y
      PJSET(3,NEWF)=.5*(WPL-WMN)
      PJSET(4,NEWF)=.5*(WPL+WMN)
      JORIG(NEWF)=JPACK*JET+IVIR
C
C          BOOST ALL FINAL VECTORS (EXCEPT NEW ONES) AND RECALCULATE
C          VIRTUAL MOMENTA.  BOOST IS DETERMINED BY DIFFERENCE OF
C          NEW AND OLD TOTAL FINAL MOMENTA, B2B1=BOOST2-BOOST1.
C
30    CONTINUE
      DO 201 K=1,4
201   BOOST1(K)=PFINAL(K)
      BMASS=PFINAL(5)
      DO 202 K=1,4
202   BOOST2(K)=PJSET(K,NEWV)+PJSET(K,IVIR2)-PJSET(K,NEWF)
C
C          PARAMETERS FOR COMBINED BOOSTS.
#ifdef SINGLE_X
      BDOTB=BOOST1(4)*BOOST2(4)-BOOST1(1)*BOOST2(1)-BOOST1(2)*BOOST2(2)
     $-BOOST1(3)*BOOST2(3)
      DO 203 K=1,4
203   B2B1(K)=BOOST2(K)-BOOST1(K)
#elif defined(DOUBLE_X)
C          DOUBLE PRECISION FOR 32-BIT MACHINES USING 3-VECTORS AND MASS
C          AS EXACT.
      DO 204 K=1,3
      DBL1(K)=BOOST1(K)
204   DBL2(K)=BOOST2(K)
      DBLM=BMASS
      DBL1(4)=DSQRT(DBL1(1)**2+DBL1(2)**2+DBL1(3)**2+DBLM**2)
      DBL2(4)=DSQRT(DBL2(1)**2+DBL2(2)**2+DBL2(3)**2+DBLM**2)
      BDOTB=DBL1(4)*DBL2(4)-DBL1(1)*DBL2(1)-DBL1(2)*DBL2(2)
     $-DBL1(3)*DBL2(3)
      DO 205 K=1,4
205   B2B1(K)=DBL2(K)-DBL1(K)
#endif
      B44=BDOTB/BMASS**2
      BI41=1./BMASS
      BI42=(BDOTB-BMASS**2-B2B1(4)*BMASS)/(BMASS**2*(BOOST2(4)+BMASS))
      B4K1=BI41
      B4K2=(BMASS**2-BDOTB-B2B1(4)*BMASS)/(BMASS**2*(BOOST1(4)+BMASS))
      BIK1=-1./(BMASS*(BOOST1(4)+BMASS))
      BIK2=1./(BMASS*(BOOST2(4)+BMASS))
      BIK3=(BMASS**2-BDOTB)/(BMASS**2*(BOOST1(4)+BMASS)
     $*(BOOST2(4)+BMASS))
C
C          BOOST FINAL JETS
      DO 210 J=1,NJSET
      IF(J.EQ.IVIR.OR.J.EQ.IVIR2) GO TO 210
      IF(PJSET(5,J).LT.0.) GO TO 210
      IF(JDCAY(J).EQ.-1) GO TO 210
      BP1=0.
      BP21=0.
      DO 215 K=1,3
      BP1=BP1+BOOST1(K)*PJSET(K,J)
215   BP21=BP21+B2B1(K)*PJSET(K,J)
      DO 220 K=1,3
220   PJSET(K,J)=PJSET(K,J)
     $+(B2B1(K)*BI41+BOOST2(K)*BI42)*PJSET(4,J)
     $+B2B1(K)*BP1*BIK1+BOOST2(K)*BP21*BIK2+BOOST2(K)*BP1*BIK3
      PJSET(4,J)=B44*PJSET(4,J)+BP21*B4K1+BP1*B4K2
210   CONTINUE
C
C          SET PFINAL TO BOOST2
      DO 230 K=1,4
230   PFINAL(K)=BOOST2(K)
      PFINAL(5)=BMASS
C
C          RESET REMAINING VECTORS
      DO 240 J=NJSET,1,-1
      IF(J.EQ.IVIR.OR.J.EQ.IVIR2) GO TO 240
      IF(PJSET(5,J).GE.0.) GO TO 240
      JX1=JDCAY(J)/JPACK
      JX2=JDCAY(J)-JPACK*JX1
      DO 250 K=1,4
      PJSET(K,J)=PJSET(K,JX1)-PJSET(K,JX2)
250   DBL1(K)=PJSET(K,J)
#ifdef SINGLE_X
      AMJ=SQRT(ABS(DBL1(4)**2-DBL1(1)**2-DBL1(2)**2-DBL1(3)**2))
#elif defined(DOUBLE_X)
      AMJ=DSQRT(ABS(DBL1(4)**2-DBL1(1)**2-DBL1(2)**2-DBL1(3)**2))
#endif
      PJSET(5,J)=-AMJ
240   CONTINUE
C
C          RESET PFINAL, ETC.
#ifdef SINGLE_X
      DO 300 K=1,4
300   PFINAL(K)=PFINAL(K)+PJSET(K,NEWF)
      S1=PFINAL(4)**2-PFINAL(1)**2-PFINAL(2)**2-PFINAL(3)**2
      IF(S1.LT.0.) GO TO 9999
      PFINAL(5)=SQRT(S1)
      PTOTPL=PJSET(4,NEWV)+PJSET(3,NEWV)+PJSET(4,IVIR2)+PJSET(3,IVIR2)
      PTOTMN=PJSET(4,NEWV)-PJSET(3,NEWV)+PJSET(4,IVIR2)-PJSET(3,IVIR2)
#elif defined(DOUBLE_X)
C          NEED DOUBLE PRECISION ON 32-BIT MACHINES
      CALL DBLVEC(PFINAL,DBL1)
      CALL DBLVEC(PJSET(1,NEWF),DBL2)
      DO 300 K=1,4
      DBL1(K)=DBL1(K)+DBL2(K)
300   PFINAL(K)=DBL1(K)
      S1=DBL1(4)**2-DBL1(1)**2-DBL1(2)**2-DBL1(3)**2
      PFINAL(5)=SQRT(S1)
      IF(S1.LT.0.) GO TO 9999
      PFINAL(5)=SQRT(S1)
      PTOTPL=PJSET(4,NEWV)+PJSET(3,NEWV)+PJSET(4,IVIR2)+PJSET(3,IVIR2)
      PTOTMN=PJSET(4,NEWV)-PJSET(3,NEWV)+PJSET(4,IVIR2)-PJSET(3,IVIR2)
#endif
C
C          SET NJSET AND POINTERS IF Z WAS GOOD
      IF(.NOT.ZGOOD) GO TO 10
      NJSET=NJSET+2
      IF(JDCAY(NEWF).GT.0) NJSET=NJSET+2
      JVIR(JET-10)=NEWV
      GO TO 10
C          ERROR -- DISCARD EVENT.
9999  CONTINUE
      WRITE(ITLIS,9998) IEVT
9998  FORMAT(/' ***** ERROR IN QCDINI ... EVENT',I8,' DISCARDED *****')
      NJSET=-1
      RETURN
      END