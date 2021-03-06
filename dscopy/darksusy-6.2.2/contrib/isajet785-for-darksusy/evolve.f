CDECK  ID>, EVOLVE.
      SUBROUTINE EVOLVE
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : 
C-        Call for each process a subroutine to set up
C-        Lorentz frames and perform initial and final QCD jet
C-        evolution in leading-log approximation.
C-
C-   Created  13-AUG-1991   Frank E. Paige,Serban D. Protopopescu
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
      INTEGER MXJETS
      PARAMETER (MXJETS=10)
      COMMON/PJETS/PJETS(5,MXJETS),IDJETS(MXJETS),QWJET(5),IDENTW 
     $,PPAIR(5,4),IDPAIR(4),JPAIR(4),NPAIR,IFRAME(MXJETS)
      SAVE /PJETS/
      INTEGER   IDJETS,IDENTW,IDPAIR,JPAIR,NPAIR,IFRAME
      REAL      PJETS,QWJET,PPAIR
      COMMON/PINITS/PINITS(5,2),IDINIT(2)
      SAVE /PINITS/
      INTEGER   IDINIT
      REAL      PINITS
      INTEGER   MXJSET,JPACK
      PARAMETER (MXJSET=400,JPACK=1000)
      COMMON/JETSET/NJSET,PJSET(5,MXJSET),JORIG(MXJSET),JTYPE(MXJSET),
     $JDCAY(MXJSET)
      SAVE /JETSET/
      INTEGER   NJSET,JORIG,JTYPE,JDCAY
      REAL      PJSET
      COMMON/JWORK/ZZC(MXJSET),JMATCH(MXJSET),TNEW,P1CM(4),
     1J1,J2,J3,J4,J5,E1CM,E2CM,E3CM,E4CM,E5CM
      SAVE /JWORK/
      LOGICAL TNEW
      EQUIVALENCE (J1,JJ(1)),(E1CM,EE(1))
      INTEGER   JMATCH,J1,J2,J3,J4,J5,JJ(5)
      REAL      ZZC,P1CM,E1CM,E2CM,E3CM,E4CM,E5CM,EE(5)
      COMMON/JWORK2/JVIR(2),PFINAL(5),SGN,ZMIN,ZMAX,DZMAX,JET,GLFORC(2),
     $ZGOOD,JIN(400),FXTEST(MXJSET)
      SAVE /JWORK2/
      LOGICAL GLFORC,ZGOOD
      INTEGER   JVIR,JET,JIN
      REAL      PFINAL,SGN,ZMIN,ZMAX,DZMAX,FXTEST
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      COMMON/FRAME/FRAME(5,3),N0JETS,N0W,N0PAIR
      SAVE /FRAME/
      INTEGER   N0JETS,N0W,N0PAIR
      REAL      FRAME
      REAL BP,PINCOM
      INTEGER I,K,J,JJET,IFR
C----------------------------------------------------------------------
C          Initialize
      NJSET=0
      N0JETS=0
      N0W=0
      N0PAIR=0
C
C          Copy momenta from /PINITS/ to /JETSET/
      IF(.NOT.KEYS(2)) THEN
        DO 100 I=1,2
          NJSET=NJSET+1
          JORIG(NJSET)=JPACK*(10+I)
          JTYPE(NJSET)=IDINIT(I)
          JDCAY(NJSET)=JPACK*I+I
          DO 105 K=1,5
105       PJSET(K,NJSET)=PINITS(K,I)
100     CONTINUE
      ENDIF
C
C       Handle each process separately
C
      IF(KEYS(1).OR.KEYS(8)) THEN
        CALL EVOL01
      ELSEIF(KEYS(2)) THEN
        CALL EVOL02
      ELSEIF(KEYS(3)) THEN
        CALL EVOL03
      ELSEIF(KEYS(5)) THEN
        CALL EVOL05
      ELSEIF(KEYS(6).OR.KEYS(10)) THEN
        CALL EVOL06
      ELSEIF(KEYS(7).OR.KEYS(9)) THEN
        CALL EVOL07
      ELSEIF(KEYS(11)) THEN
        CALL EVOL11
      ELSEIF(KEYS(12)) THEN
        CALL EVOL01
      ENDIF
C
      IF(NJSET.LT.0) RETURN
C
C          Boost /JETSET/ partons back to PP COM
C
      DO 500 J=1,NJSET
        JJET=JORIG(J)/JPACK
        IF ( JJET.EQ.0 ) THEN
          IFR=1
        ELSE
          IF(JJET.GT.10) GO TO 500
          IF(IDJETS(JJET).EQ.10.AND.KEYS(6)) GO TO 500
          IFR=IFRAME(JJET)
        ENDIF
        BP=0.
        DO 505 K=1,3
505     BP=BP+FRAME(K,IFR)*PJSET(K,J)
        BP=BP/FRAME(5,IFR)
        DO 510 K=1,3
510     PJSET(K,J)=PJSET(K,J)+FRAME(K,IFR)*PJSET(4,J)/FRAME(5,IFR)
     1  +FRAME(K,IFR)*BP/(FRAME(4,IFR)+FRAME(5,IFR))
        PJSET(4,J)=FRAME(4,IFR)*PJSET(4,J)/FRAME(5,IFR)+BP
500   CONTINUE
C
C          Reset PBEAM
      DO 530 J=1,NJSET
        IF(JDCAY(J).EQ.JPACK*J+J) THEN
          JJET=JORIG(J)/JPACK-10
          PINCOM=.5*(PJSET(4,J)+ABS(PJSET(3,J)))
          PBEAM(JJET)=HALFE-PINCOM
        ENDIF
530   CONTINUE
C
C          Check for zero energy partons
      CALL IRMOV0
C
      RETURN
      END
