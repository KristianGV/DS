#include "PILOT.inc"
      SUBROUTINE ISAWBG
C
C          Write initial record (type 200)
C          Inverse of RDBEG
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "itapes.inc"
#include "dylim.inc"
#include "frgpar.inc"
#include "idrun.inc"
#include "jetlim.inc"
#include "keys.inc"
#include "primar.inc"
#include "qcdpar.inc"
#include "qlmass.inc"
#include "q1q2.inc"
#include "types.inc"
#include "xmssm.inc"
#include "sslun.inc"
#include "jetsig.inc"
C
#include "zvout.inc"
#include "zevel.inc"
#include "final.inc"
#include "sssm.inc"
#include "sstype.inc"
C
      INTEGER NL,IL,ITA
      INTEGER I,K,PDFGUP(2),PDFSUP(2),IDWTUP,IDBM1,IDBM2,IMDL,IMODIN
      INTEGER NOUT
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      REAL SIG
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
C
C       Keep entry point WRBEG for backward compatibility
      ENTRY WRBEG
C
      ITA=IABS(ITEVT)
      NKINF=0
      SIGF=0.
      ALUM=0.
      ACCEPT=0.
      NRECS=0
C     Define some (irrelevant) PDF labels for ISALHE output
      DO I=1,2
      PDFGUP(I)=0
      PDFSUP(I)=0
      END DO
      IDWTUP=3
      LHEOUT=13
#ifdef CDCPACK_X
      CALL ZEROL(ZVOUT,512)
#endif
      CALL ZEROL(ZEVEL,MAXLEN)
      IL=3
      CALL MOVLEI(IDVER,IZEVEL(IL),4)
      IL=IL+4
      CALL MOVLEI(NJET,IZEVEL(IL),7)
      IL=IL+7
      NL=NJET*MXGOQ
      IF(NJET.NE.0) CALL MOVLEL(GOQ(1,1),LZEVEL(IL),NL)
      IL=IL+NL
      CALL MOVLEL(KEYS(1),LZEVEL(IL),10)
      IL=IL+10
      CALL MOVLEV(PMIN(1),ZEVEL(IL),36)
      IL=IL+36
      IF(.NOT.KEYS(3)) GOTO 11
      CALL MOVLEV(QMIN,ZEVEL(IL),12)
      IL=IL+12
   11 CONTINUE
      CALL MOVLEL(GODY(1),LZEVEL(IL),5)
      IL=IL+5
      CALL MOVLEV(PUD,ZEVEL(IL),22)
      IL=IL+22
      CALL MOVLEV(ALAM,ZEVEL(IL),4)
      IL=IL+4
      CALL MOVLEV(AMLEP(6),ZEVEL(IL),3)
      IL=IL+3
      CALL MOVLEI(LOC(1),IZEVEL(IL),100)
      IL=IL+100
      CALL MOVLEL(GOMSSM,LZEVEL(IL),1)
      IL=IL+1
      CALL MOVLEV(XGLSS,ZEVEL(IL),11)
      IL=IL+11
      CALL MOVLEL(GOSUG,LZEVEL(IL),1)
      IL=IL+1
      CALL MOVLEV(XM0SU,ZEVEL(IL),5)
      IL=IL+5
C
      IZEVEL(1)=200
      IZEVEL(2)=1
      CALL BUFOUT(IL)
      IF (WRTLHE) THEN
C     Convert beam type to HEPEVT numbering scheme
      IF (ABS(IDIN(1)).EQ.1120) IDBM1=SIGN(1,IDIN(1))*2212
      IF (ABS(IDIN(1)).EQ.12) IDBM1=SIGN(1,IDIN(1))*11
      IF (ABS(IDIN(2)).EQ.1120) IDBM2=SIGN(1,IDIN(2))*2212
      IF (ABS(IDIN(2)).EQ.12) IDBM2=SIGN(1,IDIN(2))*11
C
      OPEN(UNIT=LHEOUT,FILE='isajet.lhe',STATUS='NEW',FORM='FORMATTED')
C       Write version number
        WRITE(LHEOUT,2001)
        WRITE(LHEOUT,2002)
        WRITE(LHEOUT,2003) IDVER
        WRITE(LHEOUT,2004)
C       Write initialization info
        WRITE(LHEOUT,2005)
c        IF (NSIGS.LT.500) THEN
c       Pythia chokes on more than 500 subprocesses so dump out 500 max
        WRITE(LHEOUT,1001) IDBM1,IDBM2,HALFE,HALFE,PDFGUP,PDFSUP,
     $IDWTUP,NSIGS
        DO I=1,NSIGS
          WRITE(LHEOUT,1002) SIGS(I),0.,1.,INOUT(I)
        END DO
c        ELSE IF (NSIGS.GE.500.AND.GOMSSM) THEN
c        WRITE(LHEOUT,1001) IDBM1,IDBM2,HALFE,HALFE,PDFGUP,PDFSUP,
c     $IDWTUP,1
c        SIG=0.
c        DO I=1,NSIGS
c          SIG=SIG+SIGS(I)
c        END DO
c          WRITE(LHEOUT,1002) SIG,0.,1.,2160
c        END IF
        WRITE(LHEOUT,2006)
C
C     write LHA into LHE file: from Azar
        IF(GOMSSM) THEN
          WRITE(LHEOUT,2007)
            IF(LOC(55).gt.0) THEN
              IMDL=1   ! SUGRA 
              IF((LOC(61)+LOC(62)+LOC(63)+LOC(64)+LOC(65)+LOC(82)).gt.0)
     &          THEN
                IMODIN=3  ! NUSUG
              ELSE
                IMODIN=1  !mSUGRA
              ENDIF
            ELSEIF (LOC(76).gt.0) THEN   ! mAMSB
              IMDL=7
              IMODIN=7
            ELSEIF (LOC(81).gt.0) THEN   ! mAMSB
              IMDL=7
              IMODIN=8
            ELSEIF (LOC(83).gt.0) THEN   ! MMAMSB
              IMDL=9
              IMODIN=9
            ELSEIF (LOC(85).gt.0) THEN   ! HCAMSB
              IMDL=10
              IMODIN=10
            ELSEIF (LOC(88).gt.0) THEN   ! GENMIRAGE
              IMDL=12
              IMODIN=12
            ELSEIF (LOC(89).gt.0) THEN   ! NATAMSB
              IMDL=13
              IMODIN=13
            ELSEIF(LOC(60).gt.0) THEN    ! GMSB
              IMDL=2
              IF(LOC(69).gt.0) THEN
                IMODIN=5  ! non-minimal GMSB
              ELSE
                IMODIN=2
              ENDIF
            ELSEIF(LOC(75).gt.0) THEN    ! SUGRA+RHN
              IMDL=1
              IMODIN=6
            ENDIF
          CALL ISALHA(LHEOUT,IMDL,IMODIN,AMTP)
          CALL WRTQNM(LHEOUT)
          CALL ISALHD(LHEOUT,IDOUT(k),0,NOUT)
          DO K=1,NOUT
            CALL ISALHD(LHEOUT,IDOUT(k),k,NOUT)
          ENDDO
          WRITE(LHEOUT,2008)
        ENDIF
C
      END IF
1001  FORMAT(4X,I5,3X,I5,2X,F12.3,2X,F12.3,5(4X,I5),4X,I5)
1002  FORMAT(2X,E12.6,3X,E12.6,3X,F12.6,5X,I9)
2001  FORMAT('<LesHouchesEvents version="1.0">')
2002  FORMAT('<!--')
2003  FORMAT('File generated with Isajet ',I5)
2004  FORMAT('-->')
2005  FORMAT('<init>')
2006  FORMAT('</init>')
2007  FORMAT('<slha>')
2008  FORMAT('</slha>')
      RETURN
      END
