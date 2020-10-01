#include "PILOT.inc"
      SUBROUTINE ELCTRN
C          GENERATE E+ E- ----> QK QB EVENT USING SIGEE CROSS SECTION.
#include "itapes.inc"
#include "jetsig.inc"
#include "eepar.inc"
#include "primar.inc"
#include "pjets.inc"
#include "pinits.inc"
#include "jetpar.inc"
#include "jetlim.inc"
#include "const.inc"
#include "totals.inc"
#include "partcl.inc"
#include "xmssm.inc"
#include "sstype.inc"
      REAL AMQ(2),SSXLAM,RSH,XD,GAM,V,DUMMY
      INTEGER MSUPL,MSDNL,MSSTL,MSCHL,MSBT1,MSTP1,
     $MSUPR,MSDNR,MSSTR,MSCHR,MSBT2,MSTP2,MSW1,MSW2,
     $MSNEL,MSEL,MSNML,MSMUL,MSNTL,MSTAU1,MSER,MSMUR,MSTAU2,IDSS(85)
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
      DIMENSION LISTJ(30)
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,
     111,-11,12,-12,13,-13,14,-14,15,-15,16,-16,10,80,-80,90,81/
      DATA IDSS/0,
     $ISUPL,MSUPL,ISDNL,MSDNL,ISSTL,MSSTL,ISCHL,MSCHL,ISBT1,MSBT1,
     $ISTP1,MSTP1,
     $ISUPR,MSUPR,ISDNR,MSDNR,ISSTR,MSSTR,ISCHR,MSCHR,ISBT2,MSBT2,
     $ISTP2,MSTP2,ISW1,MSW1,ISW2,MSW2,ISZ1,ISZ2,ISZ3,ISZ4,
     $ISNEL,MSNEL,ISEL,MSEL,ISNML,MSNML,ISMUL,MSMUL,
     $ISNTL,MSNTL,ISTAU1,MSTAU1,ISER,MSER,ISMUR,MSMUR,
     $ISTAU2,MSTAU2,
     $9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,-12,13,-13,
     $14,-14,15,-15,16,-16,10,80,-80,90,82,83,84,86,-86/
C          ENTRY
      NPTCL=0
      NREJ=-1
      SIGMA=0.
      NSIGS=0
      DO 10 I=1,MXSIGS
10    SIGS(I)=0.
      WT=1.
C          GENERATE NEXT KINEMATIC POINT
100   CONTINUE
      NREJ=NREJ+1
      IF(NREJ.GT.NTRIES) GO TO 9999
      NKINPT=NKINPT+1
      SUMWT=SUMWT+SIGMA*WT
      IF (IBREM) THEN
        RSH=RSHMIN+(RSHMAX-RSHMIN)*RANF()
        SHAT=RSH**2
        QSQ=SHAT
        XD=(1.-SHAT/SCM)*(-1.+2*RANF())
        X1=(XD+SQRT(XD**2+4*SHAT/SCM))/2.
        X2=X1-XD
      ELSE
        SHAT=SCM
        RSH=SQRT(SHAT)
        X1=1.
        X2=1.
      END IF  
      PHI(1)=PHIMIN(1)+(PHIMAX(1)-PHIMIN(1))*RANF()
      PHI(2)=AMOD(PHI(1)+PI,2.*PI)
      CTH(1)=XJMIN(1)+(XJMAX(1)-XJMIN(1))*RANF()
      CTH(2)=-CTH(1)
      DO 110 I=1,2
      TH(I)=ACOS(CTH(I))
      STH(I)=SIN(TH(I))
      PT(I)=HALFE*STH(I)
      YJ(I)=.5*ALOG((1+CTH(I))/(1-CTH(I)))
      XJ(I)=CTH(I)
      IDINIT(I)=IDIN(I)
      PINITS(1,I)=0.
      PINITS(2,I)=0.
      PINITS(5,I)=AME
110   CONTINUE
C     Set some PINITS parameters
      PINITS(3,1)=X1*HALFE
      PINITS(3,2)=-X2*HALFE
      PINITS(4,1)=SQRT(PINITS(3,1)**2+AME**2)
      PINITS(4,2)=SQRT(PINITS(3,2)**2+AME**2)
C          CALCULATE CROSS SECTION
      IF (GOMSSM) THEN
        CALL SIGSSE
      ELSE
        CALL SIGEE
      END IF
      WT=XJMAX(1)-XJMIN(1)
C          TEST CROSS SECTION
      IF(SIGMA.GT.SGMXEE) SGMXEE=SIGMA
      IF(SIGMA.LT.SGMXEE*RANF()) GO TO 100
      SUMWT=SUMWT+SIGMA*WT
      NKEEP=NKEEP+1
C          SELECT JET TYPES
      SIGINV=1./SIGMA
      TRY=RANF()
      SUM=0.
      DO 200 I=1,NSIGS
      SUM=SUM+SIGS(I)*SIGINV
      IF(SUM.LT.TRY) GO TO 200
C          FIND REACTION
      ISIGS=I
      SIGEVT=SIGS(ISIGS)
      II=INOUT(I)/IOPAK**2
      JETTYP(1)=MOD(II,IOPAK)
      II=II/IOPAK
      JETTYP(2)=MOD(II,IOPAK)
      GO TO 210
200   CONTINUE
      GO TO 9998
C          SET PJETS. RESET P AND PT INCLUDING MASSES.
210   CONTINUE
      IF (GOMSSM) THEN
        AMQ(1)=AMASS(IDSS(JETTYP(1)))
        AMQ(2)=AMASS(IDSS(JETTYP(2)))
      ELSE
        AMQ(1)=AMASS(LISTJ(JETTYP(1)))
        AMQ(2)=AMASS(LISTJ(JETTYP(2)))
      END IF
      PCM=SQRT(SSXLAM(SHAT,AMQ(1)**2,AMQ(2)**2))/2./RSH
      DO 220 I=1,2
      PJETS(1,I)=PCM*STH(I)*COS(PHI(I))
      PJETS(2,I)=PCM*STH(I)*SIN(PHI(I))
      PJETS(3,I)=PCM*CTH(I)
      PJETS(4,I)=SQRT(PCM**2+AMQ(I)**2)
      PJETS(5,I)=AMQ(I)
      IF (GOMSSM) THEN
        IDJETS(I)=IDSS(JETTYP(I))
      ELSE
        IDJETS(I)=LISTJ(JETTYP(I))
      END IF
      P(I)=PCM
      PT(I)=P(I)*STH(I)
220   CONTINUE
C     IF BREMSSTRAHLUNG, THEN BOOST TO LAB FRAME
      IF (IBREM) THEN
        GAM=(X1+X2)*ECM/2./RSH
        V=-SIGN(1.,(X1-X2))*SQRT(ABS(1.-1./GAM)*(1.+1./GAM))
        DO I=1,2
          DUMMY=PJETS(4,I)
          PJETS(4,I)=GAM*(PJETS(4,I)-V*PJETS(3,I))
          PJETS(3,I)=GAM*(PJETS(3,I)-V*DUMMY)
        END DO
      END IF
      RETURN
C          ERROR MESSAGES
9998  CONTINUE
      CALL PRTEVT(0)
      WRITE(ITLIS,1010)
1010  FORMAT(//' ERROR IN ELCTRN...NO GOOD JET TYPES FOUND')
      STOP 99
9999  CONTINUE
      CALL PRTEVT(0)
      WRITE(ITLIS,1020) NTRIES
1020  FORMAT(//' IT IS TAKING MORE THAN',I5,' TRIES TO GENERATE AN',
     $' EVENT. CHECK LIMITS OR INCREASE NTRIES.')
      STOP 99
      END
