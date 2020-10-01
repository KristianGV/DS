CDECK  ID>, ZJJ.
      SUBROUTINE ZJJ
C-----------------------------------------------------------------------
C
C          Use MadGraph/Helas to generate Z + 2 jets after setup by
C          ZJJ0 using cross section routines from MadGraph:
C          ZJJ1:   q1 q1b -> Z q2 q2b, q1 != q2
C          ZJJ2:   g  g   -> Z q2 q2b
C          ZJJ3:   q1 q1b -> Z g  g
C          ZJJ4:   q1 q1b -> Z q1 q1b
C          ZJJ5:   q1 q2  -> Z q1 q2
C          ZJJ6:   q1 q1  -> Z q1 q1
C          ZJJ7:   g  q   -> Z g  q
C
C          Note: The Z is always jet1, but the other two jets are 
C          symmetrized so a symmetry factor of 1/2 is needed for every
C          subprocess. This is included by MadGraph for identical
C          particles!
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
C          Jet limits
      INTEGER MXLIM
      PARAMETER (MXLIM=8)
      INTEGER MXLX12
      PARAMETER (MXLX12=12*MXLIM)
      COMMON/JETLIM/PMIN(MXLIM),PMAX(MXLIM),PTMIN(MXLIM),PTMAX(MXLIM),
     $YJMIN(MXLIM),YJMAX(MXLIM),PHIMIN(MXLIM),PHIMAX(MXLIM),
     $XJMIN(MXLIM),XJMAX(MXLIM),THMIN(MXLIM),THMAX(MXLIM),
     $SETLMJ(12*MXLIM)
      SAVE /JETLIM/
      COMMON/FIXPAR/FIXP(MXLIM),FIXPT(MXLIM),FIXYJ(MXLIM),
     $FIXPHI(MXLIM),FIXXJ(MXLIM),FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      SAVE /FIXPAR/
      COMMON/SGNPAR/CTHS(2,MXLIM),THS(2,MXLIM),YJS(2,MXLIM),XJS(2,MXLIM)  
      SAVE /SGNPAR/
      REAL      PMIN,PMAX,PTMIN,PTMAX,YJMIN,YJMAX,PHIMIN,PHIMAX,XJMIN,
     +          XJMAX,THMIN,THMAX,BLIMS(12*MXLIM),CTHS,THS,YJS,XJS
      LOGICAL SETLMJ
      LOGICAL FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      LOGICAL FIXP,FIXPT,FIXYJ,FIXPHI,FIXXJ
      EQUIVALENCE(BLIMS(1),PMIN(1))
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
      INTEGER   MXJSET,JPACK
      PARAMETER (MXJSET=400,JPACK=1000)
      COMMON/JETSET/NJSET,PJSET(5,MXJSET),JORIG(MXJSET),JTYPE(MXJSET),
     $JDCAY(MXJSET)
      SAVE /JETSET/
      INTEGER   NJSET,JORIG,JTYPE,JDCAY
      REAL      PJSET
      INTEGER   MXPTCL,IPACK
      PARAMETER (MXPTCL=4000,IPACK=10000)
      COMMON/PARTCL/NPTCL,PPTCL(5,MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL)
     1,IDCAY(MXPTCL)
      SAVE /PARTCL/
      INTEGER   NPTCL,IORIG,IDENT,IDCAY
      REAL      PPTCL
      COMMON/PINITS/PINITS(5,2),IDINIT(2)
      SAVE /PINITS/
      INTEGER   IDINIT
      REAL      PINITS
      INTEGER MXJETS
      PARAMETER (MXJETS=10)
      COMMON/PJETS/PJETS(5,MXJETS),IDJETS(MXJETS),QWJET(5),IDENTW 
     $,PPAIR(5,4),IDPAIR(4),JPAIR(4),NPAIR,IFRAME(MXJETS)
      SAVE /PJETS/
      INTEGER   IDJETS,IDENTW,IDPAIR,JPAIR,NPAIR,IFRAME
      REAL      PJETS,QWJET,PPAIR
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
      COMMON/TOTALS/NKINPT,NWGEN,NKEEP,SUMWT,WT
      SAVE /TOTALS/
      INTEGER   NKINPT,NWGEN,NKEEP
      REAL      SUMWT,WT
C          Double precision PJETS; MXJETS defined in /JETLIM/
C          Format matches MadGraph
      COMMON/MGKIN/PJETS8(0:3,MXLIM+2),AMJET8(MXLIM+2)
      REAL*8 PJETS8,AMJET8
      SAVE /MGKIN/
C=====     Begin common blocks used by MadGraph
      REAL*8            GW, GWWA, GWWZ
      COMMON /COUP1/    GW, GWWA, GWWZ
      SAVE /COUP1/
      REAL*8            GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/   GAL,   GAU,   GAD,   GWF
      SAVE /COUP2A/
      REAL*8            GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/   GZN,   GZL,   GZU,   GZD,   G1
      SAVE /COUP2B/
      REAL*8            GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/    GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      SAVE /COUP3/
      COMPLEX*16        GCHF(2,12)
      COMMON /COUP4/    GCHF
      SAVE /COUP4/
      REAL*8            WMASS,WWIDTH,ZMASS,ZWIDTH
      COMMON /VMASS1/   WMASS,WWIDTH,ZMASS,ZWIDTH
      SAVE /VMASS1/
      REAL*8            AMASS,AWIDTH,HMASS,HWIDTH
      COMMON /VMASS2/   AMASS,AWIDTH,HMASS,HWIDTH
      SAVE /VMASS2/
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
      SAVE /FERMIONS/
      REAL*8            GG(2), G
      COMMON /COUPQCD/  GG,    G
      SAVE /COUPQCD/
C=====     End common blocks used by MadGraph
C
C          Running totals for MadGraph cross sections
C          WTTOT8/NWTTOT  = total cross section
C          WTSUM8/NWT8    = channel cross section
C          IFUNC8, IDENT8 = MadGraph function code channel flavors
C
      INTEGER MXSIG8
      PARAMETER (MXSIG8=1000)
      COMMON /MGSIGS/WTTOT8,WTSUM8(MXSIG8),WTMAX8(MXSIG8),NSIG8,
     $NWTTOT,NWT8(MXSIG8),IFUNC8(MXSIG8),IDENT8(MXLIM+2,MXSIG8),
     $ISORT8(MXSIG8)
      REAL*8 WTTOT8,WTSUM8,WTMAX8
      INTEGER NSIG8,NWTTOT,NWT8,IFUNC8,IDENT8,ISORT8
      SAVE /MGSIGS/
C
      INTEGER IMAD(6)
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3)
      EQUIVALENCE (P1(0),PJETS8(0,1))
      EQUIVALENCE (P2(0),PJETS8(0,2))
      EQUIVALENCE (P3(0),PJETS8(0,3))
      EQUIVALENCE (P4(0),PJETS8(0,4))
      EQUIVALENCE (P5(0),PJETS8(0,5))
      REAL QFCN,XX,QQ,RND,RANF,SIG,FJAC,STRUC,ALQCD
      REAL*8 SZJJ1,WT8,TERM,SUM,SZJJ2,SZJJ3,SZJJ4,SIG8,SIGI8
      REAL*8 SZJJ5,SZJJ6,SZJJ7
      INTEGER IQ,IH,ISIG8,IFL1,IFL2,IM1,IM2,IQ1,IQ2,NTRY,I,II,K,IWT8
C
C          Map Jettype/2 to MadGraph
      DATA IMAD/3,4,8,7,12,11/
C
C          Parton distributions
      QFCN(XX,IQ,IH)=STRUC(XX,QQ,IQ,IDIN(IH))/XX
C
C          Begin
C
      NTRY=0
      NJSET=0
      NPTCL=0
C
C          Select process
C
      RND=RANF()
      ISIG8=0
      SIG=0
      DO 10 I=1,NSIG8
        SIG=SIG+WTSUM8(I)/NWT8(I)
10    CONTINUE
      SUM=0
      DO 20 I=1,NSIG8
        II=ISORT8(NSIG8+1-I)
        SUM=SUM+WTSUM8(II)/NWT8(II)
        IF(SUM.GE.RND*SIG) THEN
          ISIG8=II
          GO TO 100
        ENDIF
20    CONTINUE
      WRITE(ITLIS,*) 'ERROR IN ZJJ: NO MODE FOUND'
      STOP99
C
100   CONTINUE
      SIG8=0
      FJAC=UNITS/SCM
      NTRY=NTRY+1
      IF(NTRY.GT.NTRIES) THEN
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NTRY = ',NTRY
        WRITE(ITLIS,*) 'PROCESS WAS ',(IDENT8(K,ISIG8),K=1,5)
        SIGI8=WTSUM8(ISIG8)/NWT8(ISIG8)
        WRITE(ITLIS,*) 'PROCESS SIGMA/MAX = ',SIGI8,WTMAX8(ISIG8)
        WRITE(ITLIS,*) 'CHECK YOUR LIMITS OR INCREASE NTRIES'
        STOP99
      ENDIF
C
C          Cases 1,4: q1 q1b -> z q2 q2b
C
      IF(IFUNC8(ISIG8).EQ.1.OR.IFUNC8(ISIG8).EQ.4) THEN
        AMJET8(3)=ZMASS
        IFL1=IABS(IDENT8(1,ISIG8))
        IM1=IMAD(IFL1)
        IQ1=2*IFL1
        IQ2=IQ1+1
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM1)
        IFL2=IABS(IDENT8(4,ISIG8))
        IM2=IMAD(IFL2)
        AMJET8(4)=FMASS(IM2)
        AMJET8(5)=FMASS(IM2)
        DO 210 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 220
210     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
220     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).GT.0.AND.IDENT8(4,ISIG8).GT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P4,P5,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P4,P5,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).GT.0.AND.IDENT8(4,ISIG8).LT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P5,P4,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P5,P4,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).LT.0.AND.IDENT8(4,ISIG8).GT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P4,P5,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P4,P5,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
        ELSEIF(IDENT8(1,ISIG8).LT.0.AND.IDENT8(4,ISIG8).LT.0) THEN
          IF(IFUNC8(ISIG8).EQ.1) THEN
            TERM=SZJJ1(P1,P2,P3,P5,P4,IM1,IM2)
          ELSE
            TERM=SZJJ4(P1,P2,P3,P5,P4,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ1'
          STOP99
        ENDIF
        SIG8=0.5*TERM
        GO TO 900
      ENDIF
C
C          Case 2: g g -> z q2 q2b
C
      IF(IFUNC8(ISIG8).EQ.2) THEN
        AMJET8(3)=ZMASS
        IFL1=IABS(IDENT8(1,ISIG8))
        AMJET8(1)=0
        AMJET8(2)=0
        IFL2=IABS(IDENT8(4,ISIG8))
        IM2=IMAD(IFL2)
        AMJET8(4)=FMASS(IM2)
        AMJET8(5)=FMASS(IM2)
        DO 310 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 320
310     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
320     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(4,ISIG8).GT.0) THEN
          TERM=SZJJ2(P1,P2,P3,P4,P5,IM2)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
        ELSEIF(IDENT8(4,ISIG8).LT.0) THEN
          TERM=SZJJ2(P1,P2,P3,P5,P4,IM2)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,1,1)*QFCN(X2,1,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ2'
          STOP99
        ENDIF
        SIG8=0.5*TERM
        GO TO 900
      ENDIF
C
C          Case 3: q1 q1b -> z g g
C
      IF(IFUNC8(ISIG8).EQ.3) THEN
        AMJET8(3)=ZMASS
        IFL1=IABS(IDENT8(1,ISIG8))
        IQ1=2*IFL1
        IQ2=IQ1+1
        IM1=IMAD(IFL1)
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM1)
        IFL2=9
        AMJET8(4)=0
        AMJET8(5)=0
        DO 410 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 420
410     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
420     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).GT.0) THEN
          TERM=SZJJ3(P1,P2,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).LT.0) THEN
          TERM=SZJJ3(P2,P1,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ2,1)*QFCN(X2,IQ1,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ3'
          STOP99
        ENDIF
        SIG8=TERM
        GO TO 900
      ENDIF
C
C          Cases 5,6: q1 q2 -> z q1 q2
C
      IF(IFUNC8(ISIG8).EQ.5.OR.IFUNC8(ISIG8).EQ.6) THEN
        IFL1=IABS(IDENT8(1,ISIG8))
        IM1=IMAD(IFL1)
        IFL2=IABS(IDENT8(2,ISIG8))
        IM2=IMAD(IFL2)
        IQ1=2*IFL1
        IQ2=2*IFL2
        IF(IDENT8(1,ISIG8).LT.0) IQ1=IQ1+1
        IF(IDENT8(2,ISIG8).LT.0) IQ2=IQ2+1
        AMJET8(1)=FMASS(IM1)
        AMJET8(2)=FMASS(IM2)
        AMJET8(3)=ZMASS
        AMJET8(4)=FMASS(IM1)
        AMJET8(5)=FMASS(IM2)
        DO 510 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 520
510     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
520     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).EQ.IDENT8(4,ISIG8)) THEN
          IF(IFUNC8(ISIG8).EQ.5) THEN
            TERM=SZJJ5(P1,P2,P3,P4,P5,IM1,IM2)
          ELSE
            TERM=SZJJ6(P1,P2,P3,P4,P5,IM1)
          ENDIF
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).EQ.IDENT8(5,ISIG8)) THEN
          TERM=SZJJ5(P1,P2,P3,P5,P4,IM1,IM2)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ1'
          STOP99
        ENDIF
        SIG8=TERM
        IF(IFL1.NE.IFL2) SIG8=0.5*SIG8
        GO TO 900
      ENDIF
C
C          Case 7: g q -> z g q
C
      IF(IFUNC8(ISIG8).EQ.7) THEN
        IF(IDENT8(1,ISIG8).EQ.9) THEN
          IFL1=IABS(IDENT8(2,ISIG8))
          IM1=IMAD(IFL1)
          AMJET8(1)=0
          AMJET8(2)=FMASS(IM1)
          IQ1=1
          IQ2=2*IFL1
          IF(IDENT8(2,ISIG8).LT.0) IQ2=IQ2+1
        ELSE
          IFL1=IABS(IDENT8(1,ISIG8))
          IM1=IMAD(IFL1)
          AMJET8(1)=FMASS(IM1)
          AMJET8(2)=0
          IQ2=1
          IQ1=2*IFL1
          IF(IDENT8(1,ISIG8).LT.0) IQ1=IQ1+1
        ENDIF
        AMJET8(3)=ZMASS
        IF(IDENT8(4,ISIG8).EQ.9) THEN
          AMJET8(4)=0
          AMJET8(5)=FMASS(IM1)
        ELSE
          AMJET8(4)=FMASS(IM1)
          AMJET8(5)=0
        ENDIF
        DO 610 I=1,NTRIES
          IWT8=I
          CALL MULJET(WT8)
          IF(WT8.GT.0) GO TO 620
610     CONTINUE
        WRITE(ITLIS,*) 'ERROR IN ZJJ: NO PHASE SPACE POINT IN ',
     $  NTRIES,' TRIES'
        STOP99
620     NWTTOT=NWTTOT+IWT8-1
        NWT8(ISIG8)=NWT8(ISIG8)+IWT8-1
        X1=(P1(0)+P1(3))/ECM
        X2=(P2(0)-P2(3))/ECM
        QQ=P3(1)**2+P3(2)**2+P4(1)**2+P4(2)**2+P5(1)**2+
     $  P5(2)**2+AMJET8(3)**2+AMJET8(4)**2+AMJET8(5)**2
C
C          Subcases
C
        IF(IDENT8(1,ISIG8).EQ.9.AND.IDENT8(4,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P1,P2,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(2,ISIG8).EQ.9.AND.IDENT8(4,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P2,P1,P3,P4,P5,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(1,ISIG8).EQ.9.AND.IDENT8(5,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P1,P2,P3,P5,P4,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSEIF(IDENT8(2,ISIG8).EQ.9.AND.IDENT8(5,ISIG8).EQ.9) THEN
          TERM=SZJJ7(P2,P1,P3,P5,P4,IM1)
          TERM=TERM*(4*PI*ALQCD(REAL(QQ)))**2
          TERM=TERM*WT8*FJAC*QFCN(X1,IQ1,1)*QFCN(X2,IQ2,2)
        ELSE
          WRITE(ITLIS,*) 'ERROR IN ZJJ...INVALID FLAVOR FOR ZJJ1'
          STOP99
        ENDIF
        SIG8=0.5*TERM
        GO TO 900
      ENDIF
C
C          Increment totals and test
C
900   WTTOT8=WTTOT8+SIG8
      NWTTOT=NWTTOT+1
      WTSUM8(ISIG8)=WTSUM8(ISIG8)+SIG8
      WTMAX8(ISIG8)=MAX(WTMAX8(ISIG8),SIG8)
      NWT8(ISIG8)=NWT8(ISIG8)+1
      IF(SIG8.LT.RANF()*WTMAX8(ISIG8)) GO TO 100
C
C          Good event
C
      DO 910 I=1,3
        DO 911 K=1,3
          PJETS(K,I)=PJETS8(K,I+2)
911     CONTINUE
        PJETS(4,I)=PJETS8(0,I+2)
        PJETS(5,I)=AMJET8(I+2)
        IDJETS(I)=IDENT8(I+2,ISIG8)
910   CONTINUE
      DO 920 I=1,2
        DO 921 K=1,3
          PINITS(K,I)=PJETS8(K,I)
921     CONTINUE
        PINITS(4,I)=PJETS8(0,I)
        PINITS(5,I)=AMJET8(I)
        IDINIT(I)=IDENT8(I,ISIG8)
920   CONTINUE
C
      QSQ=QQ
      SHAT=(P1(0)+P2(0))**2-(P1(3)+P2(3))**2
      PBEAM(1)=(1.-X1)*HALFE
      PBEAM(2)=(1.-X2)*HALFE
C
C          Set /TOTALS/
C
      NKINPT=NWTTOT
      SUMWT=WTTOT8
C
      RETURN
      END
