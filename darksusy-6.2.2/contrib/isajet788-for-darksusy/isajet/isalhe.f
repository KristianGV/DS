#include "PILOT.inc"
      SUBROUTINE ISALHE
C
C     USING NOEVOL AND NOHADR, DECAY SUBPROCESS PARTICLES TO FILL
C     PARTCL COMMON BLOCK. THEN WRITE TO A .lhe FILE,
C     SO EVENT CAN BE PASSED TO OTHER GENERATORS FOR
C     SHOWERING, HADRONIZATION AND UNDERLYING EVENT
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "keys.inc"
#include "partcl.inc"
#include "primar.inc"
#include "jetpar.inc"
#include "pjets.inc"
#include "pinits.inc"
#include "idrun.inc"
#include "sstype.inc"
#include "listss.inc"
#include "xmssm.inc"
#include "sslun.inc"
#include "const.inc"
C
      INTEGER I,IFL1,IFL2,IP1,JET,NFIRST,IP
      INTEGER LISTJ(17),LISTW(4),LISTSM(30),IPAK,ID
      INTEGER ICOLOR(2,100),ISTAT,ITRANS,I1
      INTEGER IF1,IF2,IF3,JSPIN,INDEX,indx1,indx2,indx3,indx4
      INTEGER ND,N1,N2,N3,L1,IMO1,IMO2
      REAL AMASS
C
       INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2
      PARAMETER (MSUPL=-ISUPL)
      PARAMETER (MSDNL=-ISDNL)
      PARAMETER (MSSTL=-ISSTL)
      PARAMETER (MSCHL=-ISCHL)
      PARAMETER (MSBT1=-ISBT1)
      PARAMETER (MSTP1=-ISTP1)
      PARAMETER (MSUPR=-ISUPR)
      PARAMETER (MSDNR=-ISDNR)
      PARAMETER (MSSTR=-ISSTR)
      PARAMETER (MSCHR=-ISCHR)
      PARAMETER (MSBT2=-ISBT2)
      PARAMETER (MSTP2=-ISTP2)
      PARAMETER (MSW1=-ISW1)
      PARAMETER (MSW2=-ISW2)
      PARAMETER (MSNEL=-ISNEL)
      PARAMETER (MSEL=-ISEL)
      PARAMETER (MSNML=-ISNML)
      PARAMETER (MSMUL=-ISMUL)
      PARAMETER (MSNTL=-ISNTL)
      PARAMETER (MSTAU1=-ISTAU1)
      PARAMETER (MSER=-ISER)
      PARAMETER (MSMUR=-ISMUR)
      PARAMETER (MSTAU2=-ISTAU2)
      DATA LISTSM/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,-12,13,-13,
     $14,-14,15,-15,16,-16,10,80,-80,90,81/
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,7,-7,8,-8/
      DATA LISTW/10,80,-80,90/
      DATA IPAK/100/
C     FILL PARTCL FROM JETPAR: FINAL PARTONS
        NPTCL=NJET
        DO 100 I=1,NJET
        PPTCL(1,I)=PT(I)*COS(PHI(I))
        PPTCL(2,I)=PT(I)*SIN(PHI(I))
        PPTCL(3,I)=P(I)*CTH(I)
        IF(KEYS(1)) THEN
          IDENT(I)=LISTJ(JETTYP(I))
        ELSEIF(KEYS(2)) THEN
          IDENT(I)=IDJETS(I)
        ELSEIF(KEYS(5).OR.(KEYS(10).AND.GOMSSM)) THEN
          IDENT(I)=LISTSS(JETTYP(I))
        ELSEIF(KEYS(6)) THEN
          IDENT(I)=LISTW(JETTYP(I))
        ELSEIF(KEYS(8)) THEN
          IF(JETTYP(1).LE.13) THEN
            IFL1=LISTJ(JETTYP(1))
          ELSE
            IFL1=10
          ENDIF
          IF(JETTYP(2).LE.13) THEN
            IFL2=LISTJ(JETTYP(2))
          ELSE
            IFL2=10
          ENDIF
          IDENT(1)=IFL1
          IDENT(2)=IFL2
        ELSEIF(KEYS(10)) THEN
          IDENT(I)=LISTSM(JETTYP(I))
        ENDIF
        PPTCL(5,I)=AMASS(IDENT(I))
        PPTCL(4,I)=SQRT(P(I)**2+PPTCL(5,I)**2)
        IORIG(I)=-(IPACK*I)
        IDCAY(I)=0
100   CONTINUE
C     Implement color connection for 2-> 2 subprocess
      IF (NPTCL.EQ.2) THEN
        CALL COLR22(IDINIT(1),IDINIT(2),IDENT(1),IDENT(2),ICOLOR)
      END IF
C     NOW DECAY FINAL STATE PARTONS
      DO 610 IP=1,NJET
        NFIRST=NPTCL+1
        JET=IP
        CALL DECAY(IP)
c
        ND=NPTCL-(NFIRST-1)
        N1=NFIRST
        N2=NFIRST+1
        L1=IP
        IF (ND.EQ.2) THEN
          CALL COLR12(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,ICOLOR)
        END IF
        IF (ND.EQ.3) THEN
          N3=NFIRST+2
          CALL COLR13(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,
     $                IDENT(N3),N3,ICOLOR)
        END IF
C
        DO 620 IP1=NFIRST,NPTCL
620     IORIG(IP1)=ISIGN(IABS(IORIG(IP1))+IPACK*JET,IORIG(IP1))
610   CONTINUE
C     NOW DECAY THE DECAY PRODUCTS
      IP=NJET+1
700   NFIRST=NPTCL+1
      JET=IABS(IORIG(IP))/IPACK
      CALL DECAY(IP)
c
        ND=NPTCL-(NFIRST-1)
        N1=NFIRST
        N2=NFIRST+1
        L1=IP
        IF (ND.EQ.2) THEN
          CALL COLR12(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,ICOLOR)
        END IF
        IF (ND.EQ.3) THEN
          N3=NFIRST+2
          CALL COLR13(IDENT(IP),L1,IDENT(N1),N1,IDENT(N2),N2,
     $                IDENT(N3),N3,ICOLOR)
        END IF
C
      DO 720 IP1=NFIRST,NPTCL
720     IORIG(IP1)=ISIGN(IABS(IORIG(IP1))+IPACK*JET,IORIG(IP1))
      IP=IP+1
      IF (IP.LE.NPTCL) GO TO 700
C
C     Now output to isajet.lhe
      WRITE(LHEOUT,1001) 
C     Here one needs to invert particle IDs using LISTJ 
C     for idinit or LISTSS for IDENT if running SUSY
C     in order to match up with INOUT reaction code.
C
        ID=IPAK**3*JETTYP(2)+IPAK**2*JETTYP(1)+IPAK*INITYP(2)+INITYP(1)
c     If we have SUSY production, just dump out one type of subprocess, 
C     since Pythia can only handle 500 or less
C        IF (GOMSSM) THEN
C          ID=2160
C        END IF
       WRITE(LHEOUT,1002) NPTCL+2,ID,1.,QSQ,ALFA,ALFQSQ
C     Write out initial state particles
        DO I=1,2
          WRITE(LHEOUT,1003) ITRANS(IDINIT(I),1),-1,0,0,
     $ICOLOR(1,I),ICOLOR(2,I),PINITS(1,I),PINITS(2,I),PINITS(3,I),
     $PINITS(4,I),PINITS(5,I),0.,9.
        END DO
        DO I=1,NPTCL
          IF (IDCAY(I).EQ.0.) THEN 
            ISTAT=1
          ELSE
            ISTAT=2
          END IF
          IF (IORIG(I).EQ.0) ISTAT=-1
          I1=IABS(IORIG(I))
          JET=I1/IPACK
          I1=I1-IPACK*JET
          I1=ISIGN(I1,IORIG(I))
          IF (I.LE.2) THEN
            IMO1=1
            IMO2=2
          ELSE
            IMO1=I1+2
            IMO2=0
          END IF
          WRITE(LHEOUT,1003) ITRANS(IDENT(I),1),ISTAT,IMO1,IMO2,
     $ICOLOR(1,I+2),ICOLOR(2,I+2),PPTCL(1,I),PPTCL(2,I),PPTCL(3,I),
     $PPTCL(4,I),PPTCL(5,I),0.,9.
        END DO
      WRITE(LHEOUT,1004) 
1001  FORMAT('<event>')
1002  FORMAT(4X,I3,4X,I8,3X,F12.5,3X,E12.6,3X,F12.6,3X,F12.6)
1003  FORMAT(6X,I8,3(2X,I4),2(2X,I3),5(2X,E12.6),2(1X,F2.0))
1004  FORMAT('</event>')
      RETURN
      END
