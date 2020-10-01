CDECK  ID>, GETTOT.
      SUBROUTINE GETTOT(PRFLAG)
C
C          Calculate total cross section within jet limits.
C          If PRFLAG=.TRUE. print summary.
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/FINAL/NKINF,SIGF,ALUM,ACCEPT,NRECS
      SAVE /FINAL/
      INTEGER   NKINF,NRECS
      REAL      SIGF,ALUM,ACCEPT
      COMMON/TIMES/TIME1,TIME2
      SAVE /TIMES/
      REAL      TIME1,TIME2
      COMMON/TOTALS/NKINPT,NWGEN,NKEEP,SUMWT,WT
      SAVE /TOTALS/
      INTEGER   NKINPT,NWGEN,NKEEP
      REAL      SUMWT,WT
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
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
      COMMON/DYLIM/QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     2  THWMAX,PHWMIN,PHWMAX
     3  ,SETLMQ(12)
      SAVE /DYLIM/
      LOGICAL SETLMQ
      EQUIVALENCE(BLIM1(1),QMIN)
      REAL      QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     +          THWMAX,PHWMIN,PHWMAX,BLIM1(12)
      INTEGER MXTYPE
      PARAMETER (MXTYPE=8)
      COMMON/TYPES/LOC(100),NTYP,NJTTYP(MXTYPE),NWWTYP(2),NWMODE(3)
      COMMON/XTYPES/PARTYP(40),TITLE(10),JETYP(30,MXTYPE),WWTYP(30,2)
     $,WMODES(30,3)
      SAVE /TYPES/,/XTYPES/
      CHARACTER*8 JETYP,WWTYP,TITLE,PARTYP,WMODES
      INTEGER   LOC,NTYP,NJTTYP,NWWTYP,NWMODE
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      COMMON/SEEDg/XSEED
      SAVE /SEEDg/
      CHARACTER*24 XSEED
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/ISLOOP/NEVOLV,NFRGMN,IEVOL,IFRG
      SAVE /ISLOOP/
      INTEGER NEVOLV,NFRGMN,IEVOL,IFRG
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
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C
      REAL DELPHI,SIGF2,FRAC,TMEAN,ALUM2,SIGF3
      LOGICAL PRFLAG
      INTEGER I,II,KK
      REAL TMP
      LOGICAL MGFLAG
      INTEGER L,LINT,LK1,LK2
C
C          Calculate jet cross sections
C
      SIGF=0.
      ALUM=0.
      ACCEPT=0.
      NKINF=NKINPT
C          For 2-body processes we can use the totals.
C          For MadGraph we must sum the partial cross sections.
      MGFLAG=KEYS(12)
      IF(NKINPT.GT.0.AND..NOT.MGFLAG) THEN
        SIGF=SUMWT/NKINPT
        DELPHI=2.*PI
        IF(KEYS(1).OR.KEYS(2).OR.KEYS(5).OR.KEYS(6).OR.KEYS(8)
     $  .OR.KEYS(9)) THEN
          DELPHI=PHIMAX(1)-PHIMIN(1)
        ELSEIF(KEYS(3).AND..NOT.STDDY) THEN
          DELPHI=PHWMAX-PHWMIN
        ENDIF
        SIGF=SIGF*DELPHI/(2.*PI)
      ELSEIF(MGFLAG) THEN
        DO 10 I=1,NSIG8
          SIGF=SIGF+WTSUM8(I)/NWT8(I)
10      CONTINUE
      ENDIF
C
C          Print summary if desired
C
      IF(.NOT.PRFLAG) RETURN
C
C          Print header and title
      WRITE(ITLIS,100)
100   FORMAT('1',30('*')/' *',28X,'*'/
     1' *',5X,'ISAJET RUN SUMMARY',5X,'*'/
     2' *',28X,'*'/1X,30('*')//)
      WRITE(ITLIS,101) TITLE
101   FORMAT(//11X,10A8)
      IF(NKINPT.EQ.0) GO TO 300
C
C          Print cross section
      WRITE(ITLIS,102) NKINPT
102   FORMAT(//' NO. KINEMATIC POINTS GENERATED      =',I13)
      SIGF2=SIGF*NEVOLV*NFRGMN
      WRITE(ITLIS,103) SIGF2
103   FORMAT(//' MONTE CARLO JET CROSS SECTION       =',E13.4,' MB')
      IF(SIGF.EQ.0.) WRITE(ITLIS,111)
111   FORMAT(' CROSS SECTION IS ZERO IF VARIABLES ARE FIXED')
C
C          Print W decay acceptance
      IF(KEYS(3)) THEN
        ACCEPT=FLOAT(NKEEP)/FLOAT(NWGEN)
        WRITE(ITLIS,105) ACCEPT
105     FORMAT(//' ACCEPTANCE FOR W DECAYS             =',E13.4)
      ELSEIF(KEYS(7)) THEN
        ACCEPT=FLOAT(NKEEP)/FLOAT(NWGEN)
        WRITE(ITLIS,106) ACCEPT
106     FORMAT(//' ACCEPTANCE FOR H DECAYS             =',E13.4)
      ENDIF
C
C          Print luminosity
      IF(SIGF.GT.0.) THEN
        ALUM=NEVENT/SIGF
        IF(KEYS(4)) ALUM=NKINPT/SIGF
        WRITE(ITLIS,104) ALUM
104     FORMAT(//' EQUIVALENT INTEGRAL LUMINOSITY      =',E13.4,
     $  ' /MB')
      ENDIF
C
C          Print final multijet cross sections
      IF(KEYS(12)) THEN
        WRITE(ITLIS,9001)
9001    FORMAT(//6X,'FINAL MULTIJET CROSS SECTIONS'/
     $  6X,'PROCESS',18X,'SIGMA',10X,'MAX(SIGMA)')
        DO 992 I=1,NSIG8
          II=ISORT8(I)
          TMP=WTSUM8(II)/NWT8(II)
          WRITE(ITLIS,9002) (IDENT8(KK,II),KK=1,5),TMP,WTMAX8(II)
9002      FORMAT(2X,5I5,2E15.5)
992     CONTINUE
        WRITE(ITLIS,*)
      ENDIF
C
C          Print statistics for multiple evolution and fragmentation
      IF(NEVOLV.GT.1.OR.NFRGMN.GT.1) THEN
        FRAC=FLOAT(IEVGEN)/FLOAT(IEVT)
        WRITE(ITLIS,201) IEVGEN
201     FORMAT(//' NUMBER OF ACCEPTED EVENTS           =',I13)
        WRITE(ITLIS,202) FRAC
202     FORMAT(' FRACTION OF ACCEPTED EVENTS         =',E13.4)
        SIGF3=SIGF2*FRAC
        WRITE(ITLIS,203) SIGF3
203     FORMAT(' CROSS SECTION FOR ACCEPTED EVENTS   =',E13.4)
      ENDIF
C
C          Print mean time per event
300   TMEAN=(TIME2-TIME1)/NEVENT
      WRITE(ITLIS,301) TMEAN
301   FORMAT(//' MEAN TIME PER GENERATED EVENT       =',E13.4,
     $' SEC')
C
C          Print final seed(s)
      CALL RLUXAT(L,LINT,LK1,LK2)
      WRITE(ITLIS,302) LINT,LK1,LK2
302   FORMAT(//' FINAL RANLUX SEEDS                  =',3I12)
C
C          Print run identifier
      WRITE(ITLIS,303) IDG
303   FORMAT(//' END OF ISAJET RUN                   =',2I9.6)
C
C  Print LHE summary: from Azar
      IF (WRTLHE) THEN
        WRITE(LHEOUT,1000)
        WRITE(LHEOUT,1002) NEVENT
        WRITE(LHEOUT,1003) SIGF2*1e9
        WRITE(LHEOUT,1004) ALUM*1e-12
        WRITE(LHEOUT,1009)
      ENDIF
1000  FORMAT('<MCinfo>')
1002  FORMAT(' Number of Events                     :',I13)
1003  FORMAT(' Cross section (pb)                   :',E13.4)
1004  FORMAT(' Equivalent integral luminosity (/fb) :',E13.4)
1009  FORMAT('</MCinfo>')
      RETURN
      END
