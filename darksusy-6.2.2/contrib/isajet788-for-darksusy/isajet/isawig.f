#include "PILOT.inc"
C-----------------------------------------------------------------------
C
C                          HERWIG-ISAJET interface
C
C-----------------------------------------------------------------------
C
C                    VERSION 1.200 of 24th September 2002.
C
C-----------------------------------------------------------------------
C  Subroutines to produce an output file from ISAJET with SUSY masses 
C  and decay modes in a form HERWIG can read in.
C-----------------------------------------------------------------------
C  We also include the code to calculate the R-parity violating decay
C  modes which ISAJET does not include
C-----------------------------------------------------------------------
C  We now include an interface to HDECAY to allow the use of 
C  next-to-leading order Higgs decay rates
C-----------------------------------------------------------------------
C
C  We have made changes to allow the use of the code with either ISAJET
C  7.63/64 or 7.58 this is to enable the Snowmass accord points to be 
C  still be generated with the defined version of ISAJET.
C
C  The default is ISAJET7.64 in order to obtain code compatible with
C  7.63/7.58 the compiler option
C
C  -DISAJET758   for ISAJET7.58
C  -DISAJET763   for ISAJET7.63
C
C  should be used.
C
C-----------------------------------------------------------------------
C
C  We have also made changes to make using HDECAY easier and to allow 
C  the use of the new version of HDECAY.
C
C  The default is to compile a dummy subroutine. If you wish to use
C  HDECAY you should use the following compiler options
C
C -DHDECAY2   for version 2.0 of HDECAY
C -DHDECAY3   for version 3.0 of HDECAY
C
C-----------------------------------------------------------------------
C
C       F. Paige, March 2005:
C  Patchy format; Isajet common blocks included with +CDE statements.
C  OPEN statement removed.
C  Unit number for output passed as first argument, IH, and DATA
C  statement removed.
C  HDECAY version number (0,2,3) passed as second argument, IHDCY.
C  Ignored unless corresponding version selected with Patchy flag.
C
C-----------------------------------------------------------------------
C
C  This block contains the following code
C
C  ISAWIG to output the HERWIG decay table
C
C  Subroutines for R parity
C
C  RPDECY calculates all the decay rates
C  RPINF1 is a function for the amplitude squared terms in 3-body ME's
C  RPINF2 is a function for the interference terms in 3-body ME's
C  RPINT1 is the integrand for RPINF1
C  RPINT2 is the integrand for LR interference terms
C  RPINT3 is the integrand for the remaining terms
C  RPMAIN routine for the user to enter the couplings
C  RPNORM adds R-parity modes to tables, and removes modes < MINBR
C  RPMODA adds an R-parity violating mode to the table
C  RPRATE is a function for the 2-body rates
C  RPRTCH routine to test for negative square roots
C
C  HDECAY interface routines
C
C  HDCYAD is a routine used by the HDECAY interface to add decay modes
C         to the ISAJET decay tables
C  HDCYSY is the interface routine to HDECAY
C
C--13/04/99 Modified to work with ISAJET 7.42  by Peter Richardson 
C--02/07/01 Modified to work with AMSB models in ISAJET 7.48 and 7.51
C                       by Bryan Webber
C--09/04/02 Fixes to various RPV bugs and top decays to charged Higgs
C                       by Peter Richardson
C--19/09/02 Changes to allow use with either ISAJET7.58 or ISAJET7.63
C                       by Peter Richardson
C--19/09/02 Incorperated the HDECAY interface and made changes to allow
C           the use of HDECAY2 or 3 by Peter Richardson
C--24/09/02 Checked the code worked with ISAJET7.64 Peter Richardson
C
C-----------------------------------------------------------------------
C All the RPV rates from JHEP 0004:008,2000 
C                        by H. Dreiner, P. Richardson and M.H. Seymour
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CDECK  ID>, ISAWIG
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author : Bryan Webber, Kosuke Odagiri & Peter Richardson
C-----------------------------------------------------------------------       
      SUBROUTINE ISAWIG(IH,IHDCY,MT,M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,
     &                      MER2,MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR)
C-----------------------------------------------------------------------
C
C  Writes a table of SUSY particle (and top quark) properties
C  and decays in a format suitable for reading by HERWIG.
C  Call at end of SSRUN or SUGRUN.
C  MT = top mass used in ISASUSY or ISASUGRA.
C  Other parameters passed via common /SUGMG/ & /SSMODE/
C
C--18/08/98 modified by BRW to include gravitino (for GMSB option)
C--11/01/99 modified by KO to output mixing matrices in correct format
C--18/01/99 modified by KO to output correct sign of ALFAH
C--02/04/99 modified by PR to include R-parity violation
C--13/04/99 modified by PR to work with ISAJET 7.42
C--23/04/99 modified by PR to inculde KO A and mu terms
C--25/05/99 modified by PR to change R-parity ME code to 300
C--18/06/99 modified by PR to move data statements and rewrite a format
C                             statement caused problems with LINUX 
C--14/08/99 modified by PR to change file format and remove modes < MINBR
C--16/11/99 modified by PR to interface with HDECAY
C--31/03/00 modified by PR to fix a bug in the chargino mixing,
C                             the sfermion mixing angles
C--25/05/00 modified by BRW to include 1st and 2nd generation pseudoscalar
C                              mesons as possible decay products (for AMSB)
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
C
C          ISAJET common blocks and EQUIVALENCE
C
#include "ssmode.inc"
#include "sspar.inc"
      REAL MSS(72)
      EQUIVALENCE (AMGLSS,MSS(1))
C
C  Common Block for the R-parity violating couplings
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
C--Inputs from ISAJET
      REAL MT,M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,
     &                        MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR
C
      REAL GEV2S,WIDTH,LTIM,MASS,THX,THY
      INTEGER IDEC(5),I,J,K,L,M,ID,JD,IH,NSS,NPP
      INTEGER IHDCY
C  SUSY HERWIG names, charge & PDG codes, ISAJET ID & mass codes
      PARAMETER (NSS=65,NPP=NSS+43)
      INTEGER IMISA(NSS),NDEC(NSS),IDHW(NPP),IDISA(NPP)
      DATA (IDISA(I),IMISA(I),I=1,NSS-2)/
     & 22, 4, 21, 2, 23, 6, 24, 8, 25,12, 26,16,-22, 4,-21, 2,
     &-23, 6,-24, 8,-25,12,-26,16, 42, 5, 41, 3, 43, 7, 44, 9,
     & 45,13, 46,17,-42, 5,-41, 3,-43, 7,-44, 9,-45,13,-46,17,
     & 32,18, 31,26, 34,20, 33,27, 36,24, 35,28,-32,18,-31,26,
     &-34,20,-33,27,-36,24,-35,28, 52,19, 51, 0, 54,21, 53, 0,
     & 56,25, 55, 0,-52,19,-51, 0,-54,21,-53, 0,-56,25,-55, 0,
     & 29, 1, 30,31, 40,32, 50,33, 60,34, 39,51, 49,52,-39,51,
     &-49,52, 82,55, 83,56, 84,57, 86,58,-86,58, 91,66/
C  non-SUSY HERWIG names & codes, PDG & ISAJET codes
      DATA (IDHW(I),IDISA(I),I=NSS-1,NPP)/
     &  6,  6, 12, -6,  1,  2,  2,  1,  3,  3,  4,  4,  5,  5,  7, -2,
     &  8, -1,  9, -3, 10, -4, 11, -5,121, 12,122, 11,123, 14,124, 13,
     &125, 16,126, 15,127,-12,128,-11,129,-14,130,-13,131,-16,132,-15,
     & 13,  9, 59, 10,198, 80,199,-80,200, 90,
     & 21, 110, 38, 120, 30,-120, 22, 220, 46, 130, 34,-130,
     & 50, 230, 42,-230, 25, 330,175, 140,140,-140,171, 240,
     &136,-240,179, 340,144,-340,163, 440/
C--constants etc.
      DATA GEV2S/6.582122E-25/
      DATA NDEC/NSS*0/
C Input variable
      CHARACTER *1 HYORN
C  Decide whether or not to use HDECAY
C  Valid version (2 or 3) must be defined with Patchy flag.
#ifdef HDECAY2_X
      IF(IHDCY.EQ.2) THEN
        CALL HDCYSY(M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,
     &                 MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR) 
      ENDIF
#elif defined(HDECAY3_X)
      IF(IHDCY.EQ.3) THEN
        CALL HDCYSY(M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,
     &                 MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR) 
      ENDIF
#endif
C  Calculate the R-parity violating modes
      CALL RPMAIN
C  Add R-parity violating modes if needed and remove all modes less MINBR
      CALL RPNORM
C  Tell HERWIG where to store SUSY particles
      DO I=1,57
        IDHW(I)=400+I
      ENDDO
      DO I=58,62
        IDHW(I)=145+I
      ENDDO
      IDHW(63)=458
C
C  Output SUSY particle + top quark table
C
      WRITE (IH,'(I4)') NSS
      DO I=1,NSS
        IF (I.LT.63) THEN
          MASS=0.
          IF (IMISA(I).NE.0) MASS=MSS(IMISA(I))
          IF (I.GT.49.AND.I.LT.58) MASS=-MASS
        ELSEIF (I.EQ.63) THEN
          MASS=0.
        ELSE
C t and tbar
          MASS=MT
        ENDIF
C  Compute lifetime in ps
        ID=ABS(IDISA(I))
        WIDTH=0.
        DO J=1,NSSMOD
          IF (ISSMOD(J).EQ.ID) THEN
            NDEC(I)=NDEC(I)+1
            WIDTH=WIDTH+GSSMOD(J)
          ENDIF
        ENDDO
        IF (WIDTH.NE.0.) THEN
          LTIM=GEV2S/WIDTH
        ELSE
          LTIM=1E30
        ENDIF
        WRITE (IH,1) IDHW(I),MASS,LTIM
      ENDDO
    1 FORMAT(I5,F12.4,E15.5)
C 
C  Output decay modes
C
      DO I=1,NSS
        WRITE (IH,'(I4)') NDEC(I)
        IF (NDEC(I).NE.0) THEN
          ID=ABS(IDISA(I))
          DO J=1,NSSMOD
            IF (ISSMOD(J).EQ.ID) THEN
C  Translate decay products
              DO K=1,5
                L=JSSMOD(K,J)
                IF (L.EQ.0) THEN
                  IDEC(K)=0
                ELSE
                  IF (IDISA(I).LT.0) L=-L
                  DO M=1,NPP
                    IF (IDISA(M).EQ.L) THEN
                      IDEC(K)=IDHW(M)
                      GO TO 10
                    ENDIF
                  ENDDO
C   Antiparticle=particle
                  L=-L
                  DO M=1,NPP
                    IF (IDISA(M).EQ.L) THEN
                      IDEC(K)=IDHW(M)
                      GO TO 10
                    ENDIF
                  ENDDO
                  PRINT *,' Unknown ISAJET ID =',L
                  IDEC(K)=20
   10             CONTINUE
                ENDIF
              ENDDO
              JD=IDHW(I)
              IF (JD.EQ.6.OR.JD.EQ.12) THEN
C--bug fix 9/04/02 by P.R. for t --> H b
C   Special for t and tbar
                IF (IDEC(1).GT.120.AND.IDEC(1).LT.133.AND.
     &              IDEC(3).NE.0) THEN
C   Leptonic decay: ISAJET order is wrong for M.E. calcn
                  WRITE (IH,11) JD,BSSMOD(J),100,IDEC(2),IDEC(1),
     &                  IDEC(3),0,0
                ELSEIF(IDEC(3).NE.0) THEN
C   Nonleptonic decay
                  WRITE (IH,11) JD,BSSMOD(J),100,(IDEC(K),K=1,5)
C--Higgs decay
                ELSE
                  WRITE (IH,11) JD,BSSMOD(J),0,(IDEC(K),K=1,5)
                ENDIF
C--end of fix
C   RPARITY 3-body matrix elements, all code 200, HERWIG decides what to do
              ELSEIF(JD.GE.449.AND.JD.LE.457.AND.IDEC(1).LE.140.
     &          AND.IDEC(2).LE.140.AND.IDEC(3).LE.140) THEN
                WRITE (IH,11) JD,BSSMOD(J),300,(IDEC(K),K=1,5)
              ELSE
                WRITE (IH,11) JD,BSSMOD(J),0,(IDEC(K),K=1,5)
              ENDIF
   11         FORMAT(I6,E16.8,6I6)
   12         FORMAT(1X,A8,'  -->   ',3(A8,2X),1P,E12.4)
            ENDIF
          ENDDO
        ENDIF
      ENDDO  
C  tan(beta) and alpha, neutralinos, charginos, and sfermion mixings
      WRITE (IH,'(2F16.8)') 1./MSS(30),-MSS(59)
C  neutralino mixing matrix - in ascending mass order (rows)
C             and in the order (bino, w3ino, higgs1, higgs2) (columns)
      WRITE (IH,13) MSS(38),MSS(37),-MSS(36),-MSS(35)
      WRITE (IH,13) MSS(42),MSS(41),-MSS(40),-MSS(39)
      WRITE (IH,13) MSS(46),MSS(45),-MSS(44),-MSS(43)
      WRITE (IH,13) MSS(50),MSS(49),-MSS(48),-MSS(47)
C--Bug fix 22/03/00 by Peter Richardson
      THX=SIGN(1.0D0,1.0D0/TAN(MSS(53)))
      THY=SIGN(1.0D0,1.0D0/TAN(MSS(54)))
c                   WMXVSS(1,1) WMXVSS(1,2)  WMXVSS(2,1) WMXVSS(2,2)
      WRITE (IH,13)     -SIN(MSS(54)),    -COS(MSS(54)),
     &              -THY*COS(MSS(54)),THY*SIN(MSS(54))
c                   WMXUSS(1,1) WMXUSS(1,2)  WMXUSS(2,1) WMXUSS(2,2)
      WRITE (IH,13)      -SIN(MSS(53)),  -COS(MSS(53)),
     &               -THX*COS(MSS(53)),THX*SIN(MSS(53))
C  sfermion mixing A terms and mu for 
      WRITE (IH,'(3F16.8)') -MSS(61),-MSS(63),-MSS(65)
      WRITE (IH,'(3F16.8)')  MSS(60),MSS(62),MSS(64)
      WRITE (IH,'( F16.8)') -MSS(29)
C  R-parity violating couplings
      WRITE(IH,'(L5)') RPARTY
      IF(.NOT.RPARTY) THEN
        WRITE(IH,20) (((LAMDA1(I,J,K),K=1,3),J=1,3),I=1,3)
        WRITE(IH,20) (((LAMDA2(I,J,K),K=1,3),J=1,3),I=1,3)
        WRITE(IH,20) (((LAMDA3(I,J,K),K=1,3),J=1,3),I=1,3)
      ENDIF
      RETURN
 13   FORMAT(4F16.8)
 19   FORMAT(I6)
 20   FORMAT(27E16.8)
      END
CDECK  ID>, RPDECY
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
      SUBROUTINE RPDECY
C----------------------------------------------------------------------- 
C     SUBROUTINE TO CALCULATE ALL THE R-PARITY VIOLATING RATES
C     2-BODY SQUARK AND SLEPTON
C     3-BODY NEUTRALINO, GLUINO AND CHARGINO
C----------------------------------------------------------------------- 
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
C
C          ISAJET common blocks and EQUIVALENCE
C
#include "ssmode.inc"
#include "sspar.inc"
      REAL MSS(72)
      EQUIVALENCE (AMGLSS,MSS(1))
C
C--Common block containing the Rparity violating couplings
C
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
C   
C--Common block containing the Standard Model parameters          
C 
#include "sssm.inc"
C
C--Common Block to contain R-parity violating decay rates
C
      INTEGER NRPMD
      PARAMETER (NRPMD=5000)
      COMMON/RPARRT/NSSMD2,ISSMD2(NRPMD),JSSMD2(5,NRPMD),
     &              GSSMD2(NRPMD),BSSMD2(NRPMD)
      REAL GSSMD2,BSSMD2
      INTEGER ISSMD2,JSSMD2,NSSMD2
      SAVE /RPARRT/
#include "sslun.inc"
      REAL M(4),RESM(6),WIDTH(6),A(6),B(6)
      LOGICAL CRSTRM(6)
C--External isajet gaussian quadrature routine
      REAL RPINF1,RPINF2
      EXTERNAL RPRATE,RPINF1,RPINF2
C--Local Variables
      REAL MQU(6),MLP(3),MSLLT(6),MSLRT(3),MSQLT(6),MSQRT(6),MIXING(6),
     &     RATE,RPRATE,SLLTWD(6),SLRTWD(3),SQLTWD(6),CHARWD,CHARM(2), 
     &     SQRTWD(6),PI,G,ECHAR,CBETA,SBETA,CWEAK,LAMCOL,NEUTWD,  
     &     SWEAK,NPRIME(4,4),TBETA,GLUWD,WMXVSS(2,2),WMXUSS(2,2),
     &     BMIXSS(2,2),TMIXSS(2,2),LMIXSS(2,2),LCHAR(4),RCHAR(4),
     &     MCHAR(2),THX,THY
      INTEGER I,J,K,L,PIN,POUT1,POUT2,POUT3,CHANEL,N,NMSIGN(4),MIX,
     &        CMSIGN(2)
      REAL ZERO,ONE,EPS
      PARAMETER (ZERO=0.0,ONE=1.0,EPS=1E-40)
C--Couplings, etc
      PI = 3.1415926E0
      ECHAR = SQRT(4*PI*ALFAEM) 
      G     = SQRT(4*PI*ALFA2)
      TBETA = 1./RV2V1
      CBETA = COS(ATAN(TBETA))
      SBETA = SIN(ATAN(TBETA))
      CWEAK = COS(ASIN(SQRT(SN2THW)))
      SWEAK = SQRT(SN2THW)
C--Neutralino mixing
      DO I=1,4
        NMSIGN(I) = -INT(AMZISS(I)/ABS(AMZISS(I)))
      ENDDO
      DO I=1,4
        NPRIME(I,1) = ZMIXSS(4,I)*CWEAK+ZMIXSS(3,I)*SWEAK
        NPRIME(I,2) = -ZMIXSS(4,I)*SWEAK+ZMIXSS(3,I)*CWEAK
        NPRIME(I,3) = -ZMIXSS(2,I)
        NPRIME(I,4) = -ZMIXSS(1,I)
      ENDDO
C--Chargino mixing, fixed 01/04/00 PR
      THX=SIGN(1.0D0,1.0D0/TAN(GAMMAL))
      THY=SIGN(1.0D0,1.0D0/TAN(GAMMAR))
      CMSIGN(1) = -INT(AMW1SS/ABS(AMW1SS))
      CMSIGN(2) = -INT(AMW2SS/ABS(AMW2SS))
      WMXVSS(1,1) =    -SIN(GAMMAR)
      WMXVSS(1,2) =    -COS(GAMMAR)
      WMXVSS(2,1) =-THY*COS(GAMMAR)
      WMXVSS(2,2) = THY*SIN(GAMMAR) 
      WMXUSS(1,1) =    -SIN(GAMMAL)
      WMXUSS(1,2) =    -COS(GAMMAL)
      WMXUSS(2,1) =-THX*COS(GAMMAL)
      WMXUSS(2,2) = THX*SIN(GAMMAL)
      CHARM(1) = AMW1SS
      CHARM(2) = AMW2SS
C--Number of R-parity violating modes
      NSSMD2    = 0
C--Set up local mass variables 
      MQU(1)   = AMDN 
      MQU(2)   = AMUP
      MQU(3)   = AMST
      MQU(4)   = AMCH
      MQU(5)   = AMBT
      MQU(6)   = AMTP
      MLP(1)   = AME 
      MLP(2)   = AMMU
      MLP(3)   = AMTAU
      MSLLT(1) = AMELSS
      MSLLT(2) = AMN1SS
      MSLLT(3) = AMMLSS
      MSLLT(4) = AMN2SS
      MSLLT(5) = AML1SS
      MSLLT(6) = AMN3SS
      MSLRT(1) = AMERSS
      MSLRT(2) = AMMRSS
      MSLRT(3) = AML2SS
      MSQLT(1) = AMDLSS
      MSQLT(2) = AMULSS
      MSQLT(3) = AMSLSS
      MSQLT(4) = AMCLSS
      MSQLT(5) = AMB1SS
      MSQLT(6) = AMT1SS
      MSQRT(1) = AMDRSS
      MSQRT(2) = AMURSS
      MSQRT(3) = AMSRSS
      MSQRT(4) = AMCRSS
      MSQRT(5) = AMB2SS
      MSQRT(6) = AMT2SS
C--Scalar top/bottom/tau mixing, bug fix 01/04/00 PR
      TMIXSS(1,1) =  COS(THETAT) 
      TMIXSS(1,2) = -SIN(THETAT)
      TMIXSS(2,1) =  SIN(THETAT)
      TMIXSS(2,2) =  COS(THETAT)
      BMIXSS(1,1) =  COS(THETAB)
      BMIXSS(1,2) = -SIN(THETAB)
      BMIXSS(2,1) =  SIN(THETAB)
      BMIXSS(2,2) =  COS(THETAB)
      LMIXSS(1,1) =  COS(THETAL)
      LMIXSS(1,2) = -SIN(THETAL)
      LMIXSS(2,1) =  SIN(THETAL)
      LMIXSS(2,2) =  COS(THETAL)
C--Now begin the rate calculation
C  Scalar decay rates via Rparity violation
C--CHECKED AGAINST CALCULATIONS PR 1/4/99
C--First the rates of left charged leptons
      DO I=1,3
        DO J=1,3
          DO K=1,3
C--Via LLE
            RATE = ZERO
            IF(I.NE.3) THEN 
              IF(ABS(LAMDA1(J,I,K)).GT.EPS) 
     &          RATE = RPRATE(LAMDA1(J,I,K)**2,
     &                      MSLLT(2*I-1),MLP(K),ZERO)
                CALL RPMODA(RATE,30+2*I,10+2*K,-(9+2*J),0)
            ELSE
C--New for left/right stau mixing
              IF(ABS(LAMDA1(J,I,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA1(J,I,K)**2*LMIXSS(1,1)**2,
     &                      MSLLT(2*I-1),MLP(K),ZERO) 
                CALL RPMODA(RATE,36,10+2*K,-(9+2*J),0)
                RATE = RPRATE(LAMDA1(J,I,K)**2*LMIXSS(1,2)**2,
     &                      MSLRT(I),MLP(K),ZERO)
                CALL RPMODA(RATE,56,10+2*K,-(9+2*J),0)
              ENDIF
            ENDIF
C--Via LUD
            RATE = ZERO
            IF(J.EQ.1) POUT1 = -1
            IF(J.GT.1) POUT1 = -2*J
            IF(K.EQ.1) POUT2 = 2
            IF(K.GT.1) POUT2 = 2*K-1
            IF(I.NE.3) THEN
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) 
     &          RATE = RPRATE(3*LAMDA2(I,J,K)**2,MSLLT(2*I-1),
     &                        MQU(2*J),MQU(2*K-1))
                CALL RPMODA(RATE,30+2*I,POUT1,POUT2,0)
            ELSE
C--New for left/right stau mixing
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                RATE=RPRATE(3*LAMDA2(I,J,K)**2*LMIXSS(1,1)**2,
     &                      MSLLT(2*I-1), MQU(2*J),MQU(2*K-1))
                CALL RPMODA(RATE,30+2*I,POUT1,POUT2,0)
                RATE=RPRATE(3*LAMDA2(I,J,K)**2*LMIXSS(1,2)**2,
     &                      MSLRT(I), MQU(2*J),MQU(2*K-1))
                CALL RPMODA(RATE,56,POUT1,POUT2,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Now right charged leptons via LLE
      DO K=1,3
        DO I=1,3
          DO J=1,3
            RATE = ZERO
            IF(K.NE.3) THEN
              IF(ABS(LAMDA1(I,J,K)).GT.EPS) 
     &          RATE = RPRATE(LAMDA1(I,J,K)**2,MSLRT(K),
     &                     MLP(J),ZERO) 
                CALL RPMODA(RATE,50+2*K,10+2*J,9+2*I,0)
            ELSE
C--New for left/right stau mixing
              IF(ABS(LAMDA1(I,J,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA1(I,J,K)**2*LMIXSS(2,1)**2,
     &                     MSLLT(2*K-1),MLP(J),ZERO) 
                CALL RPMODA(RATE,36,10+2*J,9+2*I,0)
                RATE = RPRATE(LAMDA1(I,J,K)**2*LMIXSS(2,2)**2,
     &                     MSLRT(K),MLP(J),ZERO) 
                CALL RPMODA(RATE,56,10+2*J,9+2*I,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Now sneutrinos 
      DO I=1,3
        DO J=1,3
          DO K=1,3
C--Via LLE
            RATE = ZERO
            IF(ABS(LAMDA1(J,I,K)).GT.EPS) 
     &        RATE = RPRATE(LAMDA1(J,I,K)**2,MSLLT(2*I),
     &                      MLP(J),MLP(K))
              CALL RPMODA(RATE,29+2*I,-(10+2*J),10+2*K,0)
C--Via LUD
            RATE = ZERO
            IF(ABS(LAMDA2(I,J,K)).GT.EPS)
     &         RATE = RPRATE(3*LAMDA2(I,J,K)**2,MSLLT(2*I)
     &                      ,MQU(2*J-1),MQU(2*K-1))
              IF(J.EQ.1) POUT1 = -2
              IF(J.GT.1) POUT1 = 1-2*J
              IF(K.EQ.1) POUT2 = 2
              IF(K.GT.1) POUT2 = 2*K-1
              CALL RPMODA(RATE,29+2*I,POUT1,POUT2,0)
          ENDDO
        ENDDO
      ENDDO  
C--Now left up squarks via LUD
      DO J=1,3
        DO I=1,3
          DO K=1,3
            RATE = ZERO
            IF(J.EQ.1) PIN   = 21
            IF(J.GT.1) PIN   = 20+2*J
            IF(K.EQ.1) POUT1 = 2
            IF(K.GT.1) POUT1 = 2*K-1
            IF(J.NE.3) THEN
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) 
     &           RATE = RPRATE(LAMDA2(I,J,K)**2,MSQLT(2*J),
     &                        MQU(2*K-1),MLP(I))
              CALL RPMODA(RATE,PIN,-(2*I+10),POUT1,0)
            ELSE
C-- New for left/right stop mixing
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA2(I,J,K)**2*TMIXSS(1,1)**2,
     &                        MSQLT(2*J),MQU(2*K-1),MLP(I)) 
                CALL RPMODA(RATE,26,-(2*I+10),POUT1,0)
                RATE = RPRATE(LAMDA2(I,J,K)**2*TMIXSS(1,2)**2,
     &                        MSQRT(2*J),MQU(2*K-1),MLP(I)) 
                CALL RPMODA(RATE,46,-(2*I+10),POUT1,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
C--Now left down squarks via LUD
      DO J=1,3
        DO I=1,3
          DO K=1,3
            RATE = ZERO 
            IF(J.EQ.1) PIN   = 22
            IF(J.GT.1) PIN   = 19+2*J
            IF(K.EQ.1) POUT1 = 2
            IF(K.GT.1) POUT1 = 2*K-1
            IF(J.NE.3) THEN
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) 
     &          RATE = RPRATE(LAMDA2(I,J,K)**2,MSQLT(2*J-1),
     &                        MQU(2*K-1),ZERO)
              CALL RPMODA(RATE,PIN,-9-2*I,POUT1,0)
            ELSE
C--New for left/right sbottom mixing
              IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(1,1)**2,
     &                        MSQLT(2*J-1),MQU(2*K-1),ZERO)
                CALL RPMODA(RATE,25,-9-2*I,POUT1,0)
C--bug fix 9/04/02 by P.R.
                RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(2,1)**2,
     &                        MSQRT(2*J-1),MQU(2*K-1),ZERO)
                CALL RPMODA(RATE,45,-9-2*I,POUT1,0)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
C--Now right up squarks via UDD
      DO I=1,3
        DO J=1,3
          DO K=1,3
            RATE = ZERO 
            IF(I.EQ.1) PIN   = 41
            IF(I.GT.1) PIN   = 40+2*I
            IF(J.EQ.1) POUT1 = -2
            IF(J.GT.1) POUT1 = 1-2*J
            IF(K.EQ.1) POUT2 = -2
            IF(K.GT.1) POUT2 = 1-2*K
            IF(I.NE.3) THEN
              IF(ABS(LAMDA3(I,J,K)).GT.EPS.AND.J.LT.K) 
     &           RATE = RPRATE(2*LAMDA3(I,J,K)**2,MSQRT(2*I),
     &                       MQU(2*K-1),MQU(2*J-1))
              CALL RPMODA(RATE,PIN,POUT2,POUT1,0)
            ELSE
C--New for left right stop mixing
              IF(ABS(LAMDA3(I,J,K)).GT.EPS.AND.J.LT.K) THEN 
                RATE = RPRATE(2*LAMDA3(I,J,K)**2*TMIXSS(2,1)**2,
     &                     MSQLT(2*I),MQU(2*K-1),MQU(2*J-1))
                CALL RPMODA(RATE,26,POUT2,POUT1,0)
                RATE = RPRATE(2*LAMDA3(I,J,K)**2*TMIXSS(2,2)**2,
     &                     MSQRT(2*I),MQU(2*K-1),MQU(2*J-1))
                CALL RPMODA(RATE,46,POUT2,POUT1,0)
              ENDIF       
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Now right down squarks 
      DO I=1,3
        DO J=1,3
          DO K=1,3
C--via LUD 1st rate
            RATE = ZERO 
            IF(I.EQ.1) PIN   = 42
            IF(I.GT.1) PIN   = 39+2*I
            IF(K.EQ.1) POUT1 = 2
            IF(K.GT.1) POUT1 = 2*K-1
            IF(I.NE.3) THEN
              IF(ABS(LAMDA2(J,K,I)).GT.EPS) 
     &          RATE = RPRATE(LAMDA2(J,K,I)**2,MSQRT(2*I-1),
     &                        ZERO,MQU(2*K-1))
              CALL RPMODA(RATE,PIN,9+2*J,POUT1,0)
            ELSE
C--New for left/right sbottom mixing
              IF(ABS(LAMDA2(J,K,I)).GT.EPS) THEN
C--bug fix 09/04/02 by P.R.
                RATE =  RPRATE(LAMDA2(J,K,I)**2*BMIXSS(2,1)**2,
     &                         MSQLT(2*I-1),ZERO,MQU(2*K-1))
                CALL RPMODA(RATE,25,9+2*J,POUT1,0)
                RATE =  RPRATE(LAMDA2(J,K,I)**2*BMIXSS(2,2)**2,
     &                         MSQRT(2*I-1),ZERO,MQU(2*K-1))
                CALL RPMODA(RATE,45,9+2*J,POUT1,0)
              ENDIF
            ENDIF
C--via LUD 2nd rate
              RATE = ZERO
              IF(K.EQ.1) PIN = 42
              IF(K.GT.1) PIN = 39+2*K
              IF(J.EQ.1) POUT1 = 1
              IF(J.GT.1) POUT1 = 2*J 
              IF(K.NE.3) THEN
                IF(ABS(LAMDA2(I,J,K)).GT.EPS)
     &            RATE = RPRATE(LAMDA2(I,J,K)**2,MSQRT(2*K-1),
     &                          MLP(I),MQU(2*J))
                CALL RPMODA(RATE,PIN,10+2*I,POUT1,0) 
              ELSE
C-- New for left/right sbottom mixing
                IF(ABS(LAMDA2(I,J,K)).GT.EPS) THEN
                  RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(2,1)**2,
     &                           MSQLT(2*K-1),MLP(I),MQU(2*J))
                  CALL RPMODA(RATE,25,10+2*I,POUT1,0) 
C--bug fix 09/04/02 by P.R.
                  RATE = RPRATE(LAMDA2(I,J,K)**2*BMIXSS(2,2)**2,
     &                           MSQRT(2*K-1),MLP(I),MQU(2*J))
                  CALL RPMODA(RATE,45,10+2*I,POUT1,0) 
                ENDIF
              ENDIF
C--via UDD
              RATE = ZERO
              IF(I.EQ.1) PIN   = 42
              IF(I.GT.1) PIN   = 39+2*I
              IF(J.EQ.1) POUT1 = -1
              IF(J.GT.1) POUT1 = -2*J
              IF(K.EQ.1) POUT2 = -2
              IF(K.GT.1) POUT2 = 1-2*K
              IF(I.NE.3) THEN
                IF(ABS(LAMDA3(J,K,I)).GT.EPS) 
     &            RATE = RPRATE(2*LAMDA3(J,K,I)**2,MSQRT(2*I-1),
     &                            MQU(2*K-1),MQU(2*J))
                CALL RPMODA(RATE,PIN,POUT1,POUT2,0)
              ELSE
C--New for left right sbottom mixing
                IF(ABS(LAMDA3(J,K,I)).GT.EPS) THEN
                  RATE = RPRATE(2*LAMDA3(J,K,I)**2*BMIXSS(2,1)**2,
     &                            MSQLT(2*I-1),MQU(2*K-1),MQU(2*J))
                  CALL RPMODA(RATE,25,POUT1,POUT2,0)
                  RATE = RPRATE(2*LAMDA3(J,K,I)**2*BMIXSS(2,2)**2,
     &                            MSQRT(2*I-1),MQU(2*K-1),MQU(2*J))
                  CALL RPMODA(RATE,45,POUT1,POUT2,0)
                ENDIF
              ENDIF
          ENDDO
        ENDDO
      ENDDO
C------------END OF SCALAR DECAY RATES
C--Masses and Widths for the 3-body decays
C--First obtain all the widths we need and put in the arrays
      DO I=1,6
        IF(I.LE.3) SLRTWD(I) =  ZERO
        SLLTWD(I) =  ZERO
        SQLTWD(I) =  ZERO
        SQRTWD(I) =  ZERO
      ENDDO
      DO J=1,NSSMOD
        DO I=1,6
          IF(ISSMOD(J).EQ.(20+I)) SQLTWD(I) = SQLTWD(I)+GSSMOD(J)
          IF(ISSMOD(J).EQ.(30+I)) SLLTWD(I) = SLLTWD(I)+GSSMOD(J)
          IF(ISSMOD(J).EQ.(40+I)) SQRTWD(I) = SQRTWD(I)+GSSMOD(J)
          IF(ISSMOD(J).EQ.(50+2*I).AND.I.LE.3) SLRTWD(I) = 
     &                                        SLRTWD(I)+GSSMOD(J)
        ENDDO
      ENDDO
      DO J=1,NSSMD2
        DO I=1,6
          IF(ISSMD2(J).EQ.(20+I)) SQLTWD(I) = SQLTWD(I)+GSSMD2(J)
          IF(ISSMD2(J).EQ.(30+I)) SLLTWD(I) = SLLTWD(I)+GSSMD2(J)
          IF(ISSMD2(J).EQ.(40+I)) SQRTWD(I) = SQRTWD(I)+GSSMD2(J)
          IF(ISSMD2(J).EQ.(50+2*I).AND.I.LE.3) SLRTWD(I) =
     &                                        SLRTWD(I)+GSSMD2(J)
        ENDDO
      ENDDO
C--Change order of down and up to make loops easier
      RATE = SQLTWD(2)
      SQLTWD(2) = SQLTWD(1)
      SQLTWD(1) = RATE
      RATE = SQRTWD(2)
      SQRTWD(2) = SQRTWD(1)
      SQRTWD(1) = RATE
C--Now calculate the rates
C--Neutralino Decay Rates via R-parity violation
      DO CHANEL=1,4
        DO L=1,4
C--First calculate the charges we will need everywhere
          LCHAR(1) = -(ECHAR*NPRIME(L,1)/3.
     &                 +G*(.5-SN2THW/3)*NPRIME(L,2)/CWEAK)
          LCHAR(2) = ECHAR*NPRIME(L,1)*2/3.
     &                 +G*(.5-2*SN2THW/3)*NPRIME(L,2)/CWEAK
          LCHAR(3) = -(ECHAR*NPRIME(L,1)
     &                 +G*(.5-SN2THW)*NPRIME(L,2)/CWEAK)
          LCHAR(4) = G*NPRIME(L,2)/(2*CWEAK)
          RCHAR(1) = (ECHAR*NPRIME(L,1)
     &                -G*SN2THW*NPRIME(L,2)/CWEAK)/3.
          RCHAR(2) = -(ECHAR*NPRIME(L,1)-
     &                  G*SN2THW*NPRIME(L,2)/CWEAK)*2/3.
          RCHAR(3) =  ECHAR*NPRIME(L,1)
     &                -G*SN2THW*NPRIME(L,2)/CWEAK
          RCHAR(4) = ZERO
          MCHAR(1) = G*NPRIME(L,3)/(2*AMW*CBETA)
          MCHAR(2) = G*NPRIME(L,4)/(2*AMW*SBETA)
          DO I=1,3
            DO J=1,3
              DO K=1,3
                DO MIX=1,3
                  MIXING(2*MIX-1) = ONE
                  MIXING(2*MIX)   = ZERO
                  RESM(2*MIX-1)   = ZERO
                  RESM(2*MIX)     = ZERO
                  WIDTH(2*MIX-1)  = ZERO
                  WIDTH(2*MIX)    = ZERO
                ENDDO   
                NEUTWD = ZERO
                M(4) = ABS(AMZISS(L))
C--Charged lepton LQD mode
                IF(CHANEL.EQ.1) THEN
                  WIDTH(1)=SLLTWD(2*I)
                  WIDTH(3)=SQLTWD(2*J)
                  WIDTH(5)=SQRTWD(2*K-1)      
                  M(1) = MLP(I)
                  M(2) = MQU(2*J)
                  M(3) = MQU(2*K-1)
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSQLT(2*J)
                  RESM(5) = MSQRT(2*K-1)
                  IF(I.NE.3) THEN
                    A(1) = MCHAR(1)*M(1)
                    B(1) = LCHAR(3)
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2  
                      MIXING(MIX) = LMIXSS(1,MIX)
                      A(MIX) = MCHAR(1)*M(1)*LMIXSS(1,MIX)+
     &                         RCHAR(3)*LMIXSS(2,MIX)
                      B(MIX) = MCHAR(1)*M(1)*LMIXSS(2,MIX)+
     &                         LCHAR(3)*LMIXSS(1,MIX)
                    ENDDO         
                    RESM(2) = MSLRT(3)
                    WIDTH(2)= SLRTWD(3)
                  ENDIF
                  IF(J.NE.3) THEN
                    A(3) = MCHAR(2)*M(2)
                    B(3) = LCHAR(2)
                  ELSE
C--left/right stop mixing
                    DO MIX=1,2
                      MIXING(2+MIX) = TMIXSS(1,MIX) 
                      A(2+MIX) = MCHAR(2)*M(2)*TMIXSS(1,MIX)+
     &                           RCHAR(2)*TMIXSS(2,MIX)
                      B(2+MIX) = MCHAR(2)*M(2)*TMIXSS(2,MIX)+
     &                           LCHAR(2)*TMIXSS(1,MIX)     
                    ENDDO      
                    RESM(4) = MSQRT(2*J)
                    WIDTH(4) = SQRTWD(2*J)
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(4+MIX) = BMIXSS(2,MIX)
                      A(4+MIX) = MCHAR(1)*M(3)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)     
                    ENDDO  
                    RESM(5) = MSQLT(2*K-1)
                    RESM(6) = MSQRT(2*K-1)    
                    WIDTH(5)=SQLTWD(2*K-1)
                    WIDTH(6)=SQRTWD(2*K-1)      
                  ENDIF
                  LAMCOL = 6*LAMDA2(I,J,K)**2 
                  POUT1 = -10-2*I
                  IF(J.EQ.1) POUT2 = -1
                  IF(J.GT.1) POUT2 = -2*J
                  IF(K.EQ.1) POUT3 = 2
                  IF(K.GT.1) POUT3 = 2*K-1
C--neutrino LQD mode
                ELSEIF(CHANEL.EQ.2) THEN
                  WIDTH(1) = SLLTWD(2*I-1)
                  WIDTH(3) = SQLTWD(2*J-1)
                  WIDTH(5) = SQRTWD(2*K-1) 
                  M(1) = ZERO
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K-1)
                  RESM(1) = MSLLT(2*I)
                  RESM(3) = MSQLT(2*J-1)
C--bug fix 09/04/02 by P.R.
                  RESM(5) = MSQRT(2*K-1)
                  A(1) = ZERO
                  B(1) = LCHAR(4)
                  IF(J.NE.3) THEN
                    A(3) = MCHAR(1)*M(2)
                    B(3) = LCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2 
                      MIXING(2+MIX) = BMIXSS(1,MIX)
                      A(2+MIX) = MCHAR(1)*M(2)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)
                      B(2+MIX) = MCHAR(1)*M(2)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                    ENDDO      
                    RESM(4) = MSQRT(2*J-1)
                    WIDTH(4) = SQRTWD(2*J-1)
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(4+MIX) = BMIXSS(2,MIX) 
                      A(4+MIX) = MCHAR(1)*M(3)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)     
                    ENDDO  
                    RESM(5)= MSQLT(2*K-1)
                    RESM(6)= MSQRT(2*K-1)    
                    WIDTH(5)=SQLTWD(2*K-1)
                    WIDTH(6)=SQRTWD(2*K-1)      
                  ENDIF
                  LAMCOL = 6*LAMDA2(I,J,K)**2
                  POUT1 = -9-2*I 
                  IF(J.EQ.1) POUT2 = -2
                  IF(J.GT.1) POUT2 = 1-2*J
                  IF(K.EQ.1) POUT3 = 2
                  IF(K.GT.1) POUT3 = 2*K-1
C--LLE mode
                ELSEIF(CHANEL.EQ.3) THEN
                  WIDTH(1) = SLLTWD(2*I)
                  WIDTH(3) = SLLTWD(2*J-1)
                  WIDTH(5) = SLRTWD(K)
                  M(1) = MLP(I) 
                  M(2) = ZERO 
                  M(3) = MLP(K)
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSLLT(2*J)
                  RESM(5) = MSLRT(K)
                  A(3) = ZERO
                  B(3) = LCHAR(4)
                  IF(I.NE.3) THEN
                    A(1) = MCHAR(1)*M(1)
                    B(1) = LCHAR(3)
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2  
                      MIXING(MIX) = LMIXSS(1,MIX)
                      A(MIX) =MCHAR(1)*M(1)*LMIXSS(1,MIX)+
     &                        RCHAR(3)*LMIXSS(2,MIX)
                      B(MIX) =MCHAR(1)*M(1)*LMIXSS(2,MIX)+
     &                        LCHAR(3)*LMIXSS(1,MIX)
                    ENDDO         
                    RESM(2) = MSLRT(3)
                    WIDTH(2)= SLRTWD(3)
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(3)
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2
                      MIXING(4+MIX) = LMIXSS(2,MIX)  
                      A(4+MIX) = MCHAR(1)*M(3)*LMIXSS(2,MIX)+
     &                           LCHAR(3)*LMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*LMIXSS(1,MIX)+
     &                           RCHAR(3)*LMIXSS(2,MIX)
                    ENDDO  
                    RESM(5)= MSLLT(2*K-1)
                    RESM(6)= MSLRT(K)    
                    WIDTH(5)=SLLTWD(2*K)
                    WIDTH(6)=SLRTWD(K)      
                  ENDIF
                  LAMCOL = 2*LAMDA1(I,J,K)**2
                  POUT1 = -10-2*I
                  POUT2 =  -9-2*J
                  POUT3 = 10+2*K
C--UDD mode
                ELSEIF(CHANEL.EQ.4) THEN
                  WIDTH(1) = SQRTWD(2*I)
                  WIDTH(3) = SQRTWD(2*J-1)
                  WIDTH(5) = SQRTWD(2*K-1) 
                  M(1) = MQU(2*I)
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K-1)
                  RESM(1) = MSQRT(2*I)
                  RESM(3) = MSQRT(2*J-1)
                  RESM(5) = MSQRT(2*K-1)
                  IF(I.NE.3) THEN
                    A(1) = MCHAR(2)*M(1)
                    B(1) = RCHAR(2)
                  ELSE
C--left/right stop mixing
                    DO MIX=1,2  
                      MIXING(MIX) = TMIXSS(2,MIX)
                      A(MIX) = MCHAR(2)*M(1)*TMIXSS(2,MIX)+
     &                         LCHAR(2)*TMIXSS(1,MIX)
                      B(MIX) = MCHAR(2)*M(1)*TMIXSS(1,MIX)+
     &                         RCHAR(2)*TMIXSS(2,MIX)
                    ENDDO  
                    RESM(1) = MSQLT(2*I)
                    RESM(2) = MSQRT(2*I)
                    WIDTH(1) = SQLTWD(2*I)
                    WIDTH(2) = SQRTWD(2*I)
                  ENDIF
                  IF(J.NE.3) THEN
                    A(3) = MCHAR(1)*M(2)
                    B(3) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(2+MIX) = BMIXSS(2,MIX)
                      A(2+MIX) = MCHAR(1)*M(2)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(2+MIX) = MCHAR(1)*M(2)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)  
                    ENDDO     
                    RESM(3)= MSQLT(2*J-1)
                    RESM(4)= MSQRT(2*J-1)    
                    WIDTH(3)=SQLTWD(2*J-1)
                    WIDTH(4)=SQRTWD(2*J-1)  
                  ENDIF
                  IF(K.NE.3) THEN
                    A(5) = MCHAR(1)*M(3)
                    B(5) = RCHAR(1)
                  ELSE
C--left/right sbottom mixing
                    DO MIX=1,2  
                      MIXING(4+MIX) = BMIXSS(2,MIX)
                      A(4+MIX) = MCHAR(1)*M(3)*BMIXSS(2,MIX)+
     &                           LCHAR(1)*BMIXSS(1,MIX)
                      B(4+MIX) = MCHAR(1)*M(3)*BMIXSS(1,MIX)+
     &                           RCHAR(1)*BMIXSS(2,MIX)     
                    ENDDO  
                    RESM(5)= MSQLT(2*K-1)
                    RESM(6)= MSQRT(2*K-1)    
                    WIDTH(5)=SQLTWD(2*K-1)
                    WIDTH(6)=SQRTWD(2*K-1)      
                  ENDIF
                  IF(I.EQ.1) POUT1 = 1
                  IF(I.GT.1) POUT1 = 2*I
                  IF(J.EQ.1) POUT2 = 2
                  IF(J.GT.1) POUT2 = 2*J-1
                  IF(K.EQ.1) POUT3 = 2
                  IF(K.GT.1) POUT3 = 2*K-1
                  IF(J.LT.K) THEN
                    LAMCOL = 12*LAMDA3(I,J,K)**2
                  ELSE
                    LAMCOL = ZERO
                  ENDIF
                ENDIF
                IF(NMSIGN(L).LT.0) THEN
                  IF(CHANEL.EQ.4) THEN
                    DO N=1,6
                      B(N) = -B(N)
                    ENDDO
                  ELSE
                    DO N=1,2
                      A(N)   = -A(N)
                      A(N+2) = -A(N+2)
                      B(N+4) = -B(N+4)
                    ENDDO
                  ENDIF
                ENDIF
C--Decide whether to remove diagrams
                DO N=1,6
                  CRSTRM(N) = .FALSE.
                ENDDO
                IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 10
                CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                       AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
                CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                       AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                       .OR.(ABS(RESM(2)).LT.EPS)
                CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                       AND.ABS(RESM(3)).GT.(M(1)+M(3)))
                CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                       AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                       .OR.(ABS(RESM(4)).LT.EPS)
                CRSTRM(5) = (M(4).GT.(M(3)+ABS(RESM(5))).
     &                       AND.ABS(RESM(5)).GT.(M(1)+M(2)))
                CRSTRM(6) = (M(4).GT.(M(3)+ABS(RESM(6))).
     &                       AND.ABS(RESM(6)).GT.(M(1)+M(2)))
     &                       .OR.(ABS(RESM(6)).LT.EPS)
C--Calculation of the rate
                NEUTWD = ZERO
                IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &          CRSTRM(4).AND.CRSTRM(5).AND.CRSTRM(6)).
     &          OR.LAMCOL.LT.EPS) GOTO 10
C--first the diagram squared pieces
                IF(.NOT.CRSTRM(1)) NEUTWD = NEUTWD+ MIXING(1)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
                IF(.NOT.CRSTRM(2)) NEUTWD = NEUTWD+MIXING(2)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
                IF(.NOT.CRSTRM(3)) NEUTWD = NEUTWD+MIXING(3)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
                IF(.NOT.CRSTRM(4)) NEUTWD = NEUTWD+ MIXING(4)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
                IF(.NOT.CRSTRM(5)) NEUTWD = NEUTWD+MIXING(5)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(5),RESM(5),A(5),B(5))
                IF(.NOT.CRSTRM(6)) NEUTWD = NEUTWD+MIXING(6)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(6),RESM(6),A(6),B(6))
C--now for the light/heavy interference due left/right mixing
                IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) NEUTWD=NEUTWD+
     &          MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &           WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
                IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) NEUTWD=NEUTWD+
     &          MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &          WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
                IF(.NOT.CRSTRM(5).AND..NOT.CRSTRM(6)) NEUTWD=NEUTWD+
     &          MIXING(5)*MIXING(6)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),
     &          WIDTH(6),RESM(5),RESM(6),A(5),A(6),B(5),B(6),1)
C--now for the true interference terms
                IF(.NOT.CRSTRM(1)) THEN
                  IF(.NOT.CRSTRM(3)) NEUTWD=NEUTWD-MIXING(1)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &            RESM(1),A(1),A(3),B(1),B(3),2)
                  IF(.NOT.CRSTRM(4)) NEUTWD=NEUTWD-MIXING(1)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &            RESM(1),A(1),A(4),B(1),B(4),2)
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(1)*MIXING(5)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(1),RESM(5),
     &            RESM(1),A(1),A(5),B(1),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(1)*MIXING(6)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(1),RESM(6),
     &            RESM(1),A(1),A(6),B(1),B(6),2)
                ENDIF
                IF(.NOT.CRSTRM(2)) THEN
                  IF(.NOT.CRSTRM(3)) NEUTWD=NEUTWD-MIXING(2)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &            RESM(2),A(2),A(3),B(2),B(3),2)
                  IF(.NOT.CRSTRM(4)) NEUTWD=NEUTWD-MIXING(2)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &            RESM(2),A(2),A(4),B(2),B(4),2)
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(2)*MIXING(5)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(2),RESM(5),
     &            RESM(2),A(2),A(5),B(2),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(2)*MIXING(6)*
     &            RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(2),RESM(6),
     &            RESM(2),A(2),A(6),B(2),B(6),2)
                ENDIF
                IF(.NOT.CRSTRM(3)) THEN
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(3)*MIXING(5)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(3),RESM(5),
     &            RESM(3),A(3),A(5),B(3),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(3)*MIXING(6)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(3),RESM(6),
     &            RESM(3),A(3),A(6),B(3),B(6),2)
                ENDIF
                IF(.NOT.CRSTRM(4)) THEN
                  IF(.NOT.CRSTRM(5)) NEUTWD=NEUTWD-MIXING(4)*MIXING(5)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(4),RESM(5),
     &            RESM(4),A(4),A(5),B(4),B(5),2)
                  IF(.NOT.CRSTRM(6)) NEUTWD=NEUTWD-MIXING(4)*MIXING(6)*
     &            RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(4),RESM(6),
     &            RESM(4),A(4),A(6),B(4),B(6),2)
                ENDIF
                NEUTWD = LAMCOL*NEUTWD/((2*PI)**3*32*M(4)**3)
C--Output the rate and particles to decay tables
 10             IF(NEUTWD.GT.EPS) THEN
                  CALL RPMODA(NEUTWD,20+10*L,POUT1,POUT2,POUT3)
                  CALL RPMODA(NEUTWD,20+10*L,-POUT1,-POUT2,-POUT3)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Gluino decay rates all the info we need is already in the arrays
C--Easier as only three gluino channels. Just calculate the total rate
C--and handle the colour structure inside HERWIG
C--Now calculate the rates
C--First the LQD rates
      DO CHANEL=1,2
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO MIX=1,2
                MIXING(2*MIX-1) = ONE
                MIXING(2*MIX)   = ZERO
                RESM(2*MIX-1)   = ZERO
                RESM(2*MIX)     = ZERO
                WIDTH(2*MIX-1)  = ZERO
                WIDTH(2*MIX)    = ZERO 
              ENDDO
C--charged lepton, ubar, down
              IF(CHANEL.EQ.1) THEN
                WIDTH(1)=SQLTWD(2*J)
                WIDTH(3)=SQRTWD(2*K-1)      
                M(1) = MLP(I)
                M(2) = MQU(2*J)
                M(3) = MQU(2*K-1)
                M(4) = AMGLSS
                RESM(1) = MSQLT(2*J)
                RESM(3) = MSQRT(2*K-1)
                IF(J.NE.3) THEN
                  A(1) =ZERO 
                  B(1) = -ONE 
                ELSE
C--left/right stop mixing
                  DO MIX=1,2
                    MIXING(MIX) = TMIXSS(1,MIX)  
                    A(MIX) = TMIXSS(2,MIX)
                    B(MIX) =-TMIXSS(1,MIX)
                  ENDDO    
C--bug fix 09/04/02 by P.R.  
                  RESM(2) = MSQRT(2*J)
                  WIDTH(2) = SQRTWD(2*J)
                ENDIF
                IF(K.NE.3) THEN
                  A(3) = ZERO
                  B(3) = ONE
                ELSE
C--left/right sbottom mixing
                  MIXING(3) = BMIXSS(2,1) 
                  MIXING(4) = BMIXSS(2,2)
                  DO MIX=1,2 
                    MIXING(2+MIX) = BMIXSS(2,MIX) 
                    A(2+MIX) = -BMIXSS(1,MIX) 
                    B(2+MIX) =  BMIXSS(2,MIX)
                  ENDDO 
                  RESM(3)= MSQLT(2*K-1)
                  RESM(4)= MSQRT(2*K-1)    
                  WIDTH(3)=SQLTWD(2*K-1)
                  WIDTH(4)=SQRTWD(2*K-1)      
                ENDIF
                LAMCOL = ALFA3*4*PI*LAMDA2(I,J,K)**2 
                POUT1 = -10-2*I
                IF(J.EQ.1) POUT2 = -1
                IF(J.GT.1) POUT2 = -2*J
                IF(K.EQ.1) POUT3 = 2
                IF(K.GT.1) POUT3 = 2*K-1
C--neutrino, dbar, down.
              ELSEIF(CHANEL.EQ.2) THEN
                WIDTH(1) = SQLTWD(2*J-1)
                WIDTH(3) = SQRTWD(2*K-1) 
                M(1) = ZERO
                M(2) = MQU(2*J-1)
                M(3) = MQU(2*K-1)
                M(4) = AMGLSS
                RESM(1) = MSQLT(2*J-1)
C--bug fix 09/04/02 by P.R. 
                RESM(3) = MSQRT(2*K-1)
                IF(J.NE.3) THEN
                  A(1) =  ZERO
                  B(1) =-ONE
                ELSE
C--left/right sbottom mixing
                  DO MIX=1,2
                    MIXING(MIX) = BMIXSS(1,MIX) 
                    A(MIX) = BMIXSS(2,MIX)
                    B(MIX) = -BMIXSS(1,MIX) 
                  ENDDO      
                  RESM(2) = MSQRT(2*J-1)
                  WIDTH(2) = SQRTWD(2*J-1)
                ENDIF
                IF(K.NE.3) THEN
                  A(3) = ZERO
                  B(3) = ONE
                ELSE
C--left/right sbottom mixing
                  DO MIX=1,2 
                    MIXING(2+MIX) = BMIXSS(2,MIX)  
                    A(2+MIX) = -BMIXSS(1,MIX) 
                    B(2+MIX) =  BMIXSS(2,MIX)
                  ENDDO  
                  RESM(3)= MSQLT(2*K-1)
                  RESM(4)= MSQRT(2*K-1)    
                  WIDTH(3)=SQLTWD(2*K-1)
                  WIDTH(4)=SQRTWD(2*K-1)      
                ENDIF
                LAMCOL = ALFA3*4*PI*LAMDA2(I,J,K)**2
                POUT1 = -9-2*I 
                IF(J.EQ.1) POUT2 = -2
                IF(J.GT.1) POUT2 = 1-2*J
                IF(K.EQ.1) POUT3 = 2
                IF(K.GT.1) POUT3 = 2*K-1
              ENDIF
C--Decide whether to remove diagrams
              DO N=1,4
                CRSTRM(N) = .FALSE.
              ENDDO
C--bug fix 09/04/02 by P.R.
              CRSTRM(1) = (M(4).GT.(M(2)+ABS(RESM(1))).
     &                     AND.ABS(RESM(1)).GT.(M(1)+M(3)))
              CRSTRM(2) = (M(4).GT.(M(2)+ABS(RESM(2))).
     &                       AND.ABS(RESM(2)).GT.(M(1)+M(3)))
     &                       .OR.(ABS(RESM(2)).LT.EPS)
              CRSTRM(3) = (M(4).GT.(M(3)+ABS(RESM(3))).
     &                       AND.ABS(RESM(3)).GT.(M(1)+M(2)))
              CRSTRM(4) = (M(4).GT.(M(3)+ABS(RESM(4))).
     &                       AND.ABS(RESM(4)).GT.(M(1)+M(2)))
     &                       .OR.(ABS(RESM(4)).LT.EPS)
C--Calculation of the rate
              GLUWD = 0
              IF((CRSTRM(1).AND.CRSTRM(2).AND.
     &            CRSTRM(3).AND.CRSTRM(4)).OR.
     &           LAMCOL.LT.EPS.OR.(M(1)+M(2)+M(3)).GT.M(4)) GOTO 20
C--First the amplitude square pieces
              IF(.NOT.CRSTRM(1)) GLUWD = GLUWD+MIXING(1)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(1),RESM(1),A(1),B(1))
              IF(.NOT.CRSTRM(2)) GLUWD = GLUWD+ MIXING(2)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(2),RESM(2),A(2),B(2))
              IF(.NOT.CRSTRM(3)) GLUWD = GLUWD+MIXING(3)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(3),RESM(3),A(3),B(3))
              IF(.NOT.CRSTRM(4)) GLUWD = GLUWD+MIXING(4)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(4),RESM(4),A(4),B(4))
C--now for the light/heavy interference due left/right mixing
              IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) GLUWD=GLUWD+
     &        MIXING(1)*MIXING(2)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(1),
     &        WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
              IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) GLUWD=GLUWD+
     &        MIXING(3)*MIXING(4)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(3),
     &        WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
C--now for the true interference terms
              IF(.NOT.CRSTRM(3)) THEN
                IF(.NOT.CRSTRM(3)) GLUWD=GLUWD-MIXING(1)*MIXING(3)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &          RESM(1),A(1),A(3),B(1),B(3),2)
                IF(.NOT.CRSTRM(4)) GLUWD=GLUWD-MIXING(1)*MIXING(4)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &          RESM(1),A(1),A(4),B(1),B(4),2)
              ENDIF
              IF(.NOT.CRSTRM(2)) THEN
                IF(.NOT.CRSTRM(3)) GLUWD=GLUWD-MIXING(2)*MIXING(3)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &          RESM(2),A(2),A(3),B(2),B(3),2)
                IF(.NOT.CRSTRM(4)) GLUWD=GLUWD-MIXING(2)*MIXING(4)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &          RESM(2),A(2),A(4),B(2),B(4),2)
              ENDIF
              GLUWD = LAMCOL*GLUWD
              GLUWD = GLUWD/((2*PI)**3*32*M(4)**3)
C--Output the rate and particles to decay tables
 20           IF(GLUWD.GT.EPS) THEN
                CALL RPMODA(GLUWD,29,POUT1,POUT2,POUT3)
                CALL RPMODA(GLUWD,29,-POUT1,-POUT2,-POUT3)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Now the gluino decay via UDD, there are three diagrams before LR mixing
C--i.e. gluino -> u d d
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO MIX=1,3
              MIXING(2*MIX-1) = ONE
              MIXING(2*MIX)   = ZERO
              RESM(2*MIX-1)   = ZERO
              RESM(2*MIX)     = ZERO
              WIDTH(2*MIX-1)  = ZERO
              WIDTH(2*MIX)    = ZERO 
            ENDDO
            WIDTH(1) = SQRTWD(2*I)
            WIDTH(3) = SQRTWD(2*J-1)
            WIDTH(5) = SQRTWD(2*K-1) 
            M(1) = MQU(2*I)
            M(2) = MQU(2*J-1)
            M(3) = MQU(2*K-1)
            M(4) = AMGLSS
            RESM(1) = MSQRT(2*I)
            RESM(3) = MSQRT(2*J-1)
            RESM(5) = MSQRT(2*K-1)
            IF(I.NE.3) THEN
              A(1) = ZERO
              B(1) = ONE
            ELSE
C--left/right stop mixing
              DO MIX=1,2
                MIXING(MIX) = TMIXSS(2,MIX)  
                A(MIX) = -TMIXSS(1,MIX)
                B(MIX) = TMIXSS(2,MIX)
              ENDDO  
              RESM(1) = MSQLT(2*I)
              RESM(2) = MSQRT(2*I)
              WIDTH(1) = SQLTWD(2*I)
              WIDTH(2) = SQRTWD(2*I)
            ENDIF
            IF(J.NE.3) THEN
              A(3) = ZERO
              B(3) = ONE
            ELSE
C--left/right sbottom mixing
              DO MIX=1,2  
                MIXING(2+MIX) = BMIXSS(2,MIX)
                A(2+MIX) =-BMIXSS(1,MIX)
                B(2+MIX) = BMIXSS(2,MIX)
              ENDDO     
              RESM(3)= MSQLT(2*J-1)
              RESM(4)= MSQRT(2*J-1)    
              WIDTH(3)=SQLTWD(2*J-1)
              WIDTH(4)=SQRTWD(2*J-1)  
            ENDIF
            IF(K.NE.3) THEN
              A(5) = ZERO
              B(5) = ONE 
            ELSE
C--left/right sbottom mixing
            DO MIX=1,2  
              MIXING(4+MIX)=BMIXSS(2,MIX)
              A(4+MIX) =-BMIXSS(1,MIX)
              B(4+MIX) = BMIXSS(2,MIX)
            ENDDO  
              RESM(5)= MSQLT(2*K-1)
              RESM(6)= MSQRT(2*K-1)    
              WIDTH(5)=SQLTWD(2*K-1)
              WIDTH(6)=SQRTWD(2*K-1)      
            ENDIF
            IF(I.EQ.1) POUT1 = 1
            IF(I.GT.1) POUT1 = 2*I
            IF(J.EQ.1) POUT2 = 2
            IF(J.GT.1) POUT2 = 2*J-1
            IF(K.EQ.1) POUT3 = 2
            IF(K.GT.1) POUT3 = 2*K-1
            IF(J.LT.K) THEN
              LAMCOL = 2*ALFA3*4*PI*LAMDA3(I,J,K)**2
            ELSE
              LAMCOL = 0
            ENDIF
C--Decide whether to remove diagrams
            DO N=1,6
              CRSTRM(N) = .FALSE.
            ENDDO
            GLUWD = ZERO
            IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 30
            CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                  AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
            CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                  AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                  .OR.(ABS(RESM(2)).LT.EPS)
            CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                  AND.ABS(RESM(3)).GT.(M(1)+M(3)))
            CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                  AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                  .OR.(ABS(RESM(4)).LT.EPS)
            CRSTRM(5) = (M(4).GT.(M(3)+ABS(RESM(5))).
     &                  AND.ABS(RESM(5)).GT.(M(1)+M(2)))
            CRSTRM(6) = (M(4).GT.(M(3)+ABS(RESM(6))).
     &                  AND.ABS(RESM(6)).GT.(M(1)+M(2)))
     &                  .OR.(ABS(RESM(6)).LT.EPS)
C--Calculation of the rate
            IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &      CRSTRM(4).AND.CRSTRM(5).AND.CRSTRM(6)).
     &      OR.LAMCOL.LT.EPS) GOTO 30
C--first the diagram squared pieces
            IF(.NOT.CRSTRM(1)) GLUWD = GLUWD+ MIXING(1)**2*
     &      RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
            IF(.NOT.CRSTRM(2)) GLUWD = GLUWD+MIXING(2)**2*
     &      RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
            IF(.NOT.CRSTRM(3)) GLUWD = GLUWD+MIXING(3)**2*
     &      RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
            IF(.NOT.CRSTRM(4)) GLUWD = GLUWD+ MIXING(4)**2*
     &      RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
            IF(.NOT.CRSTRM(5)) GLUWD = GLUWD+MIXING(5)**2*
     &      RPINF1(M(1),M(2),M(3),M(4),WIDTH(5),RESM(5),A(5),B(5))
            IF(.NOT.CRSTRM(6)) GLUWD = GLUWD+MIXING(6)**2*
     &      RPINF1(M(1),M(2),M(3),M(4),WIDTH(6),RESM(6),A(6),B(6))
C--now for the light/heavy interference due left/right mixing
            IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) GLUWD=GLUWD+
     &      MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &      WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
            IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) GLUWD=GLUWD+
     &      MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &      WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
            IF(.NOT.CRSTRM(5).AND..NOT.CRSTRM(6)) GLUWD=GLUWD+
     &      MIXING(5)*MIXING(6)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),
     &      WIDTH(6),RESM(5),RESM(6),A(5),A(6),B(5),B(6),1)
C--now for the true interference terms
            IF(.NOT.CRSTRM(1)) THEN
              IF(.NOT.CRSTRM(3)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(3)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &        RESM(1),A(1),A(3),B(1),B(3),2)
              IF(.NOT.CRSTRM(4)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(4)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &        RESM(1),A(1),A(4),B(1),B(4),2)
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(5)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(1),RESM(5),
     &        RESM(1),A(1),A(5),B(1),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(1)*MIXING(6)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(1),RESM(6),
     &        RESM(1),A(1),A(6),B(1),B(6),2)
            ENDIF
            IF(.NOT.CRSTRM(2)) THEN
              IF(.NOT.CRSTRM(3)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(3)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &        RESM(2),A(2),A(3),B(2),B(3),2)
              IF(.NOT.CRSTRM(4)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(4)*
     &        RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &        RESM(2),A(2),A(4),B(2),B(4),2)
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(5)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(2),RESM(5),
     &        RESM(2),A(2),A(5),B(2),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(2)*MIXING(6)*
     &        RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(2),RESM(6),
     &        RESM(2),A(2),A(6),B(2),B(6),2)
            ENDIF
            IF(.NOT.CRSTRM(3)) THEN
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(3)*MIXING(5)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(3),RESM(5),
     &        RESM(3),A(3),A(5),B(3),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(3)*MIXING(6)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(3),RESM(6),
     &        RESM(3),A(3),A(6),B(3),B(6),2)
            ENDIF
            IF(.NOT.CRSTRM(4)) THEN
              IF(.NOT.CRSTRM(5)) GLUWD=GLUWD+0.5*MIXING(4)*MIXING(5)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(4),RESM(5),
     &        RESM(4),A(4),A(5),B(4),B(5),2)
              IF(.NOT.CRSTRM(6)) GLUWD=GLUWD+0.5*MIXING(4)*MIXING(6)*
     &        RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(4),RESM(6),
     &        RESM(4),A(4),A(6),B(4),B(6),2)
            ENDIF
            GLUWD = LAMCOL*GLUWD
            GLUWD = GLUWD/((2*PI)**3*32*M(4)**3)
C--Output the rate and particles to decay tables
 30         IF(GLUWD.GT.EPS) THEN
              CALL RPMODA(GLUWD,29,POUT1,POUT2,POUT3)
              CALL RPMODA(GLUWD,29,-POUT1,-POUT2,-POUT3)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C--Decay modes of the chargino 
C--we will do the one diagram via LLE and LQD modes first
      MCHAR(1) = 1./(SQRT(2.)*AMW*CBETA)
      MCHAR(2) = 1./(SQRT(2.)*AMW*SBETA)
      DO L=1,2
        DO CHANEL=1,3
          DO I=1,3
            DO J=1,3
              DO K=1,3
                MIXING(1) = ONE
                MIXING(2) = ZERO
                DO MIX=1,2
                  RESM(MIX) = ZERO
                  WIDTH(MIX) = ZERO
                ENDDO
                IF(CHANEL.EQ.1) THEN
C--LLE decay to charged lepton, neutrino, antineutrino
                  M(1) = ZERO
                  M(2) = MLP(J)
                  M(3) = ZERO
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLRT(K)
                  WIDTH(1) = SLRTWD(K)
                  IF(K.NE.3) THEN
                    A(1) = -WMXUSS(L,2)*MLP(K)*MCHAR(1)
                    B(1) = ZERO
                  ELSE
C--left/right stau mixing
                    DO MIX=1,2  
                      MIXING(MIX) = LMIXSS(2,MIX)
                      A(MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &             WMXUSS(L,2)*MLP(K)*MCHAR(1)*LMIXSS(2,MIX)
                      B(MIX) = ZERO
                    ENDDO     
                    RESM(1)= MSLLT(2*K-1)
                    RESM(2)= MSLRT(K)    
                    WIDTH(1)=SLLTWD(2*K)
                    WIDTH(2)=SLRTWD(K) 
                  ENDIF
                  LAMCOL = G**2*LAMDA1(I,J,K)**2
                  POUT1 = -9-2*I 
                  POUT2 = -10-2*J  
                  POUT3 = 9+2*K
                ELSEIF(CHANEL.EQ.2) THEN
C--LQD decay to antineutrino, dbar, up
                  M(1) = ZERO
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSQRT(2*K-1)
                  WIDTH(1)= SQRTWD(2*K-1)
                  IF(K.NE.3) THEN
                    A(1) = -WMXUSS(L,2)*MQU(2*K-1)*MCHAR(1)
                    B(1) = 0
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = BMIXSS(2,MIX)
                      A(MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &            WMXUSS(L,2)*BMIXSS(2,MIX)*MQU(2*K-1)*MCHAR(1)
                      B(MIX) = -MQU(2*K)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                         MCHAR(2) 
                      IF(CMSIGN(L).LT.0) B(MIX)=-B(MIX)
                    ENDDO
                    RESM(1)= MSQLT(2*K-1)
                    RESM(2)= MSQRT(2*K-1)    
                    WIDTH(1)=SQLTWD(2*K-1)
                    WIDTH(2)=SQRTWD(2*K-1)
                  ENDIF
                  LAMCOL = 3.*G**2*LAMDA2(I,J,K)**2
                  POUT1 = -9-2*I 
                  POUT2 = 1-2*J
                  POUT3 = 2*K
                  IF(J.EQ.1) POUT2 = -2
                  IF(K.EQ.1) POUT3 = 1
                ELSEIF(CHANEL.EQ.3) THEN
C--LQD decay to charged lepton, ubar, up
                  M(1) = MLP(I)
                  M(2) = MQU(2*J)
                  M(3) = MQU(2*K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSQRT(2*K-1)
                  WIDTH(1)= SQRTWD(2*K-1)
                  IF(K.NE.3) THEN
                    A(1) = -WMXUSS(L,2)*MQU(2*K-1)*MCHAR(1)
                    B(1) = ZERO
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = BMIXSS(2,MIX)
                      A(MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &            WMXUSS(L,2)*BMIXSS(2,MIX)*MQU(2*K-1)*MCHAR(1)
                      B(MIX) = -MQU(2*K)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                         MCHAR(2)
                      IF(CMSIGN(L).LT.0) B(MIX)=-B(MIX)
                    ENDDO
                    RESM(1)= MSQLT(2*K-1)
                    RESM(2)= MSQRT(2*K-1)    
                    WIDTH(1)=SQLTWD(2*K-1)
                    WIDTH(2)=SQRTWD(2*K-1) 
                  ENDIF
                  LAMCOL = 3.*G**2*LAMDA2(I,J,K)**2
                  POUT1 = -10-2*I 
                  POUT2 = -2*J
                  POUT3 = 2*K
                  IF(J.EQ.1) POUT2 = -1
                  IF(K.EQ.1) POUT3 = 1
                ENDIF  
                CHARWD = ZERO
                IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 40
C--Decide whether to remove diagrams
                DO N=1,2
                  CRSTRM(N) = .FALSE.
                ENDDO
                CRSTRM(1) = (M(4).GT.(M(3)+ABS(RESM(1))).
     &                      AND.ABS(RESM(1)).GT.(M(1)+M(2))) 
                CRSTRM(2) = (M(4).GT.(M(3)+ABS(RESM(2))).
     &                      AND.ABS(RESM(2)).GT.(M(1)+M(2))) 
     &                     .OR.(ABS(RESM(2)).LT.EPS)
C--Calculation of the rate
               IF((CRSTRM(1).AND.CRSTRM(2)).OR.LAMCOL.LT.EPS) GOTO 40
C--first the diagram squared pieces
                IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(1),RESM(1),A(1),B(1))
                IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &          RPINF1(M(1),M(2),M(3),M(4),WIDTH(2),RESM(2),A(2),B(2))
C--now for the light/heavy interference due left/right mixing
                IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &          MIXING(1)*MIXING(2)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(1),
     &          WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
                CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
 40             IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,
     &                                          POUT1,POUT2,POUT3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Now for the two diagram LLE and LQD modes, these tend to have
C--higher branching ratios
      DO L=1,2
        DO CHANEL=1,4
          DO I=1,3
            DO J=1,3
              DO K=1,3
                DO N=1,6
                  A(N) = ZERO
                  B(N) = ZERO
                  RESM(N) = ZERO
                  WIDTH(N) = ZERO
                ENDDO
                DO MIX=1,2
                  MIXING(2*MIX-1) = 1.
                  MIXING(2*MIX) = ZERO
                ENDDO
C--LLE to neutrino, neutrino, charged lepton
                IF(CHANEL.EQ.1) THEN
                  M(1) = ZERO
                  M(2) = ZERO
                  M(3) = MLP(K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSLLT(2*J-1)
                  WIDTH(1)= SLLTWD(2*I)
                  WIDTH(3)= SLLTWD(2*J)
                  IF(I.NE.3) THEN
                    B(1) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = LMIXSS(1,MIX)
                      B(MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &                WMXUSS(L,2)*LMIXSS(2,MIX)*MLP(I)*MCHAR(1)
                    ENDDO
                    RESM(2)= MSLRT(I)    
                    WIDTH(2)=SLRTWD(I) 
                  ENDIF
                  IF(J.NE.3) THEN
                    B(3) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(2+MIX) = LMIXSS(1,MIX)
                      B(2+MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &                WMXUSS(L,2)*LMIXSS(2,MIX)*MLP(J)*MCHAR(1)
                    ENDDO
                    RESM(4)= MSLRT(J)    
                    WIDTH(4)=SLRTWD(J) 
                  ENDIF
                  IF(I.GT.J) THEN
                    LAMCOL = G**2*LAMDA1(I,J,K)**2
                  ELSE
                    LAMCOL = ZERO
                  ENDIF
                  POUT1 = 9+2*I 
                  POUT2 = 9+2*J
                  POUT3 = -10-2*K
C--LLE +ve lepton, +ve lepton, -ve lepton 
                ELSEIF(CHANEL.EQ.2) THEN
                  M(1) = MLP(I)
                  M(2) = MLP(J)
                  M(3) = MLP(K)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I)
                  RESM(3) = MSLLT(2*J)
                  WIDTH(1)= SLLTWD(2*I-1)
                  WIDTH(3)= SLLTWD(2*J-1)
                  A(1) = -MLP(I)*WMXUSS(L,2)*MCHAR(1)
                  B(1) = WMXVSS(L,1)
                  A(3) = -MLP(J)*WMXUSS(L,2)*MCHAR(1)
                  B(3) = WMXVSS(L,1)
                  IF(CMSIGN(L).LT.0) THEN
                    B(1) = -B(1)
                    B(3) = -B(3)
                  ENDIF
                  IF(I.GT.J) THEN
                    LAMCOL = G**2*LAMDA1(I,J,K)**2
                  ELSE
                    LAMCOL = ZERO
                  ENDIF
                  POUT1 = -10-2*I 
                  POUT2 = -10-2*J
                  POUT3 =  10+2*K
C--LQD to charged lepton, dbar, down
                ELSEIF(CHANEL.EQ.3) THEN
                  M(1) = MLP(I)
                  M(2) = MQU(2*J-1)
                  M(3) = MQU(2*K-1)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I)
                  RESM(3) = MSQLT(2*J)
                  WIDTH(1)= SLLTWD(2*I-1)
                  WIDTH(3)= SQLTWD(2*J)
                  A(1) = -MLP(I)*WMXUSS(L,2)*MCHAR(1)
                  B(1) =  WMXVSS(L,1)
                  IF(J.NE.3) THEN
                    A(3) = -MQU(2*J-1)*WMXUSS(L,2)*MCHAR(1)
                    B(3) = WMXVSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(2+MIX) = TMIXSS(1,MIX)
                      A(2+MIX) = -MQU(2*J-1)*WMXUSS(L,2)*MCHAR(1)
     &                           *TMIXSS(1,MIX)
                      B(2+MIX) = WMXVSS(L,1)*TMIXSS(1,MIX)-
     &                WMXVSS(L,2)*TMIXSS(2,MIX)*MQU(2*J)*MCHAR(2)
                    ENDDO
                    RESM(4)= MSQRT(2*J)    
                    WIDTH(4)=SQRTWD(2*J) 
                  ENDIF
                  IF(CMSIGN(L).LT.0) THEN
                    B(1) = -B(1)
                    B(3) = -B(3)
                    B(4) = -B(4)
                  ENDIF
                  LAMCOL = 3*G**2*LAMDA2(I,J,K)**2
                  POUT1 = -10-2*I 
                  POUT2 = 1-2*J
                  POUT3 = 2*K-1
                  IF(J.EQ.1) POUT2 = -2 
                  IF(K.EQ.1) POUT3 = 2
C--LQD to neutrino, up, dbar
                ELSEIF(CHANEL.EQ.4) THEN
                  M(1) = ZERO
                  M(2) = MQU(2*J)
                  M(3) = MQU(2*K-1)
                  M(4) = ABS(CHARM(L))
                  RESM(1) = MSLLT(2*I-1)
                  RESM(3) = MSQLT(2*J-1)
                  WIDTH(1)= SLLTWD(2*I)
                  WIDTH(3)= SQLTWD(2*J-1)
                  IF(I.NE.3) THEN
                    B(1) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(MIX) = LMIXSS(1,MIX)
                      B(MIX) = WMXUSS(L,1)*LMIXSS(1,MIX)-
     &                WMXUSS(L,2)*LMIXSS(2,MIX)*MLP(I)*MCHAR(1)
                    ENDDO
                    RESM(2)= MSLRT(I)    
                    WIDTH(2)=SLRTWD(I) 
                  ENDIF
                  IF(J.NE.3) THEN
                    A(3) = -MQU(2*J)*WMXVSS(L,2)*MCHAR(2)
                    B(3) = WMXUSS(L,1)
                  ELSE
                    DO MIX=1,2
                      MIXING(2+MIX) = BMIXSS(1,MIX)
                      A(2+MIX) = -MQU(2*J)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                         MCHAR(2)
                      B(2+MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &            WMXUSS(L,2)*BMIXSS(2,MIX)*MQU(2*J-1)*MCHAR(1)
                    ENDDO
                    RESM(4)= MSQRT(2*J-1)    
                    WIDTH(4)=SQRTWD(2*J-1) 
                  ENDIF
                  IF(CMSIGN(L).LT.0) THEN
                    DO N=1,4
                      A(N) = -A(N)
                    ENDDO
                  ENDIF
                  LAMCOL = 3*G**2*LAMDA2(I,J,K)**2
                  POUT1 = 9+2*I 
                  POUT2 = 2*J
                  POUT3 = 1-2*K
                  IF(J.EQ.1) POUT2 = 1 
                  IF(K.EQ.1) POUT3 = -2
                ENDIF
C--Decide whether to remove diagrams
                CHARWD =ZERO
                IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 50
                DO N=1,4
                  CRSTRM(N) = .FALSE.
                ENDDO
                CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                       AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
                CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                       AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                       .OR.(ABS(RESM(2)).LT.EPS)
                CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                       AND.ABS(RESM(3)).GT.(M(1)+M(3)))
                CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                       AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                       .OR.(ABS(RESM(4)).LT.EPS)
C--Calculation of the rate
                CHARWD = ZERO
                IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &          CRSTRM(4)).OR.LAMCOL.LT.EPS) GOTO 50
C--first the diagram squared pieces
                IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
                IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &          RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
                IF(.NOT.CRSTRM(3)) CHARWD = CHARWD+MIXING(3)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
                IF(.NOT.CRSTRM(4)) CHARWD = CHARWD+ MIXING(4)**2*
     &          RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
C--now for the light/heavy interference due left/right mixing
                IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &          MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &           WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
                IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) CHARWD=CHARWD+
     &          MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &          WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
C--now for the true interference terms
                IF(.NOT.CRSTRM(1)) THEN
                  IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(1)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &            RESM(1),A(1),A(3),B(1),B(3),2)
                  IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(1)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &            RESM(1),A(1),A(4),B(1),B(4),2)
                ENDIF
                IF(.NOT.CRSTRM(2)) THEN
                  IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(2)*MIXING(3)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &            RESM(2),A(2),A(3),B(2),B(3),2)
                  IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(2)*MIXING(4)*
     &            RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &            RESM(2),A(2),A(4),B(2),B(4),2)
                ENDIF
C--final factors
                CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
 50             IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,
     &                                          POUT1,POUT2,POUT3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Finally the BV chargino
C--first the decay to two up type and one down type quark
      DO L=1,2
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO MIX=1,2
                MIXING(MIX) = LAMDA3(J,I,K)
                MIXING(2+MIX)= LAMDA3(I,J,K)
              ENDDO
              M(1) = MQU(2*I)
              M(2) = MQU(2*J)
              M(3) = MQU(2*K-1)
              M(4) = ABS(CHARM(L))
              RESM(1) = MSQRT(2*I-1)
              RESM(3) = MSQRT(2*J-1)
              WIDTH(1)= SQRTWD(2*I-1)
              WIDTH(3)= SQRTWD(2*J-1)
              RESM(2) =ZERO
              WIDTH(2)=ZERO
              RESM(4) =ZERO
              WIDTH(4)=ZERO
              IF(I.NE.3) THEN
                A(1)=-WMXUSS(L,2)*MCHAR(1)*MQU(2*I-1)
                B(1)= ZERO
                MIXING(2) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX) = BMIXSS(2,MIX)*MIXING(MIX)
                  A(MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &         WMXUSS(L,2)*MCHAR(1)*MQU(2*I-1)*BMIXSS(2,MIX)
                  B(MIX) = -MQU(2*I)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                     MCHAR(2)
                ENDDO
                RESM(1) = MSQLT(2*I-1)
                RESM(2) = MSQRT(2*I-1)
                WIDTH(1)= SQLTWD(2*I-1)
                WIDTH(2)= SQRTWD(2*I-1)
              ENDIF
              IF(J.NE.3) THEN
                A(3) = -MQU(2*J-1)*WMXUSS(L,2)*MCHAR(1)
                B(3) = ZERO
                MIXING(4) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(2+MIX) = BMIXSS(2,MIX)*MIXING(MIX+2)
                  A(2+MIX) = WMXUSS(L,1)*BMIXSS(1,MIX)-
     &         WMXUSS(L,2)*MCHAR(1)*MQU(2*J-1)*BMIXSS(2,MIX)
                  B(2+MIX) = -MQU(2*J)*WMXVSS(L,2)*BMIXSS(1,MIX)*
     &                     MCHAR(2)
                ENDDO 
                RESM(3) = MSQLT(2*J-1)
                RESM(4) = MSQRT(2*J-1)
                WIDTH(3)= SQLTWD(2*J-1)
                WIDTH(4)= SQRTWD(2*J-1)
              ENDIF
              IF(CMSIGN(L).LT.0) THEN
                DO N=1,4
                  B(N) = -B(N)
                ENDDO
              ENDIF
              IF(I.LE.J) THEN
                LAMCOL = 6*G**2
              ELSE
                LAMCOL = ZERO
              ENDIF
              POUT1 = 2*I
              POUT2 = 2*J
              POUT3 = 2*K-1
              IF(I.EQ.1) POUT1 = 1
              IF(J.EQ.1) POUT2 = 1 
              IF(K.EQ.1) POUT3 = 2
C--Decide whether to remove diagrams
              DO N=1,4
                CRSTRM(N) = .FALSE.
              ENDDO
              CHARWD = ZERO
              IF(M(4).LT.(M(1)+M(2)+M(3)).OR.LAMCOL.LT.EPS) GOTO 60
              CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                     AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
              CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                     AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                     .OR.(ABS(RESM(2)).LT.EPS)
              CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                     AND.ABS(RESM(3)).GT.(M(1)+M(3)))
              CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                     AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                     .OR.(ABS(RESM(4)).LT.EPS)
C--Calculation of the rate
              CHARWD = ZERO
              IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &        CRSTRM(4))) GOTO 60
C--first the diagram squared pieces
              IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
              IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
              IF(.NOT.CRSTRM(3)) CHARWD = CHARWD+MIXING(3)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
              IF(.NOT.CRSTRM(4)) CHARWD = CHARWD+ MIXING(4)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
C--now for the light/heavy interference due left/right mixing
              IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &        MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &        WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
              IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) CHARWD=CHARWD+
     &        MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &        WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
C--now for the true interference terms
              IF(.NOT.CRSTRM(1)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(1)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &          RESM(1),A(1),A(3),B(1),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(1)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &          RESM(1),A(1),A(4),B(1),B(4),2)
              ENDIF
              IF(.NOT.CRSTRM(2)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD+MIXING(2)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &          RESM(2),A(2),A(3),B(2),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD+MIXING(2)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &          RESM(2),A(2),A(4),B(2),B(4),2)
              ENDIF
              CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
C--Identical particle symmetry factor
              IF(I.EQ.J) CHARWD = CHARWD/2.0
C--Output the rate and particles to decay tables
 60           IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,POUT1,
     &                                                    POUT2,POUT3)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C--Now the decay to three down quarks
      DO L=1,2
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO MIX=1,2
                MIXING(MIX)   = LAMDA3(I,J,K)
                MIXING(MIX+2) = LAMDA3(J,K,I)
                MIXING(MIX+4) = LAMDA3(K,I,J)
              ENDDO   
              DO MIX=1,3
                RESM(2*MIX-1)   = ZERO
                RESM(2*MIX)     = ZERO
                WIDTH(2*MIX-1)  = ZERO
                WIDTH(2*MIX)    = ZERO 
              ENDDO
              CHARWD = ZERO
              M(1) = MQU(2*I-1)
              M(2) = MQU(2*J-1)
              M(3) = MQU(2*K-1)
              M(4) = ABS(CHARM(L))
              RESM(1) = MSQRT(2*I)
              RESM(3) = MSQRT(2*J)
              RESM(5) = MSQRT(2*K)
              WIDTH(1)= SQRTWD(2*I)
              WIDTH(3)= SQRTWD(2*J)
              WIDTH(5)= SQRTWD(2*K)
              IF(I.NE.3) THEN
                A(1) = -WMXVSS(L,2)*MQU(2*I)*MCHAR(2)
                B(1) = ZERO
                MIXING(2) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX) = TMIXSS(2,MIX)*MIXING(MIX)
                  A(MIX) =WMXVSS(L,1)*TMIXSS(1,MIX)-
     &     WMXVSS(L,2)*MQU(2*I)*MCHAR(2)*TMIXSS(2,MIX)
                  B(MIX) =-MQU(2*I-1)*TMIXSS(1,MIX)*WMXUSS(L,2)*
     &                     MCHAR(1)
                ENDDO
C--bug fix 09/04/02 by P.R.
                RESM(1) = MSQLT(2*I)
                RESM(2) = MSQRT(2*I)
                WIDTH(1)= SQLTWD(2*I)
                WIDTH(2)= SQRTWD(2*I)
              ENDIF
              IF(J.NE.3) THEN
                A(3) = -WMXVSS(L,2)*MQU(2*J)*MCHAR(2)
                B(3) = ZERO 
                MIXING(4) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX+2) = TMIXSS(2,MIX)*MIXING(MIX+2)
                  A(MIX+2) =WMXVSS(L,1)*TMIXSS(1,MIX)-
     &     WMXVSS(L,2)*MQU(2*J)*MCHAR(2)*TMIXSS(2,MIX)
                  B(MIX+2) =-MQU(2*I-1)*TMIXSS(1,MIX)*WMXUSS(L,2)*
     &                     MCHAR(1)
                ENDDO
C--bug fix 09/04/02 by P.R.
                RESM(3) = MSQLT(2*J)
                RESM(4) = MSQRT(2*J)
                WIDTH(3)= SQLTWD(2*J)
                WIDTH(4)= SQRTWD(2*J)
              ENDIF
              IF(K.NE.3) THEN
                A(5) = -WMXVSS(L,2)*MQU(2*K)*MCHAR(2)
                B(5) = ZERO
C--bug fix 09/04/02 by P.R.
                MIXING(6) = ZERO
              ELSE
                DO MIX=1,2
                  MIXING(MIX+4) = TMIXSS(2,MIX)*MIXING(MIX+4)
                  A(MIX+4) =WMXVSS(L,1)*TMIXSS(1,MIX)-
     &     WMXVSS(L,2)*MQU(2*K)*MCHAR(2)*TMIXSS(2,MIX)
                  B(MIX+4) =-MQU(2*K-1)*TMIXSS(1,MIX)*WMXUSS(L,2)*
     &                     MCHAR(1)
                ENDDO
C--bug fix 09/04/02 by P.R.
                RESM(5) = MSQLT(2*K)
                RESM(6) = MSQRT(2*K)
                WIDTH(5)= SQLTWD(2*K)
                WIDTH(6)= SQRTWD(2*K)
              ENDIF
              IF(CMSIGN(L).LT.0) THEN
                DO N=1,6
                  A(N) = -A(N)
                ENDDO
              ENDIF
              IF(K.LE.J.AND.J.LE.I) THEN
                LAMCOL = 6*G**2
              ELSE
                LAMCOL = ZERO
              ENDIF
              POUT1 = 1-2*I
              POUT2 = 1-2*J
              POUT3 = 1-2*K
              IF(I.EQ.1) POUT1 = -2
              IF(J.EQ.1) POUT2 = -2
              IF(K.EQ.1) POUT3 = -2
C--Decide whether to remove diagrams
              DO N=1,6
                CRSTRM(N) = .FALSE.
              ENDDO
              IF(M(4).LT.(M(1)+M(2)+M(3))) GOTO 70
              CRSTRM(1) = (M(4).GT.(M(1)+ABS(RESM(1))).
     &                     AND.ABS(RESM(1)).GT.(M(2)+M(3))) 
              CRSTRM(2) = (M(4).GT.(M(1)+ABS(RESM(2))).
     &                     AND.ABS(RESM(2)).GT.(M(2)+M(3))) 
     &                     .OR.(ABS(RESM(2)).LT.EPS)
              CRSTRM(3) = (M(4).GT.(M(2)+ABS(RESM(3))).
     &                     AND.ABS(RESM(3)).GT.(M(1)+M(3)))
              CRSTRM(4) = (M(4).GT.(M(2)+ABS(RESM(4))).
     &                     AND.ABS(RESM(4)).GT.(M(1)+M(3)))
     &                     .OR.(ABS(RESM(4)).LT.EPS)
              CRSTRM(5) = (M(4).GT.(M(3)+ABS(RESM(5))).
     &                     AND.ABS(RESM(5)).GT.(M(1)+M(2)))
              CRSTRM(6) = (M(4).GT.(M(3)+ABS(RESM(6))).
     &                     AND.ABS(RESM(6)).GT.(M(1)+M(2)))
     &                     .OR.(ABS(RESM(6)).LT.EPS)
C--Calculation of the rate
              CHARWD = ZERO
              IF((CRSTRM(1).AND.CRSTRM(2).AND.CRSTRM(3).AND.
     &        CRSTRM(4).AND.CRSTRM(5).AND.CRSTRM(6)).
     &        OR.LAMCOL.LT.EPS) GOTO 70
C--first the diagram squared pieces
              IF(.NOT.CRSTRM(1)) CHARWD = CHARWD+ MIXING(1)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(1),RESM(1),A(1),B(1))
              IF(.NOT.CRSTRM(2)) CHARWD = CHARWD+MIXING(2)**2*
     &        RPINF1(M(2),M(3),M(1),M(4),WIDTH(2),RESM(2),A(2),B(2))
              IF(.NOT.CRSTRM(3)) CHARWD = CHARWD+MIXING(3)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(3),RESM(3),A(3),B(3))
              IF(.NOT.CRSTRM(4)) CHARWD = CHARWD+ MIXING(4)**2*
     &        RPINF1(M(1),M(3),M(2),M(4),WIDTH(4),RESM(4),A(4),B(4))
              IF(.NOT.CRSTRM(5)) CHARWD = CHARWD+MIXING(5)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(5),RESM(5),A(5),B(5))
              IF(.NOT.CRSTRM(6)) CHARWD = CHARWD+MIXING(6)**2*
     &        RPINF1(M(1),M(2),M(3),M(4),WIDTH(6),RESM(6),A(6),B(6))
C--now for the light/heavy interference due left/right mixing
              IF(.NOT.CRSTRM(1).AND..NOT.CRSTRM(2)) CHARWD=CHARWD+
     &        MIXING(1)*MIXING(2)*RPINF2(M(2),M(3),M(1),M(4),WIDTH(1),
     &        WIDTH(2),RESM(1),RESM(2),A(1),A(2),B(1),B(2),1)
              IF(.NOT.CRSTRM(3).AND..NOT.CRSTRM(4)) CHARWD=CHARWD+
     &        MIXING(3)*MIXING(4)*RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),
     &        WIDTH(4),RESM(3),RESM(4),A(3),A(4),B(3),B(4),1)
              IF(.NOT.CRSTRM(5).AND..NOT.CRSTRM(6)) CHARWD=CHARWD+
     &        MIXING(5)*MIXING(6)*RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),
     &        WIDTH(6),RESM(5),RESM(6),A(5),A(6),B(5),B(6),1)
C--now for the true interference terms
              IF(.NOT.CRSTRM(1)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD-MIXING(1)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(1),RESM(3),
     &          RESM(1),A(1),A(3),B(1),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD-MIXING(1)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(1),RESM(4),
     &          RESM(1),A(1),A(4),B(1),B(4),2)
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(1)*MIXING(5)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(1),RESM(5),
     &          RESM(1),A(1),A(5),B(1),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(1)*MIXING(6)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(1),RESM(6),
     &          RESM(1),A(1),A(6),B(1),B(6),2)
              ENDIF
              IF(.NOT.CRSTRM(2)) THEN
                IF(.NOT.CRSTRM(3)) CHARWD=CHARWD-MIXING(2)*MIXING(3)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(3),WIDTH(2),RESM(3),
     &          RESM(2),A(2),A(3),B(2),B(3),2)
                IF(.NOT.CRSTRM(4)) CHARWD=CHARWD-MIXING(2)*MIXING(4)*
     &          RPINF2(M(1),M(3),M(2),M(4),WIDTH(4),WIDTH(2),RESM(4),
     &          RESM(2),A(2),A(4),B(2),B(4),2)
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(2)*MIXING(5)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(5),WIDTH(2),RESM(5),
     &          RESM(2),A(2),A(5),B(2),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(2)*MIXING(6)*
     &          RPINF2(M(1),M(2),M(3),M(4),WIDTH(6),WIDTH(2),RESM(6),
     &          RESM(2),A(2),A(6),B(2),B(6),2)
              ENDIF
              IF(.NOT.CRSTRM(3)) THEN
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(3)*MIXING(5)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(3),RESM(5),
     &          RESM(3),A(3),A(5),B(3),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(3)*MIXING(6)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(3),RESM(6),
     &          RESM(3),A(3),A(6),B(3),B(6),2)
              ENDIF
              IF(.NOT.CRSTRM(4)) THEN
                IF(.NOT.CRSTRM(5)) CHARWD=CHARWD-MIXING(4)*MIXING(5)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(5),WIDTH(4),RESM(5),
     &          RESM(4),A(4),A(5),B(4),B(5),2)
                IF(.NOT.CRSTRM(6)) CHARWD=CHARWD-MIXING(4)*MIXING(6)*
     &          RPINF2(M(2),M(1),M(3),M(4),WIDTH(6),WIDTH(4),RESM(6),
     &          RESM(4),A(4),A(6),B(4),B(6),2)
              ENDIF
              CHARWD = LAMCOL*CHARWD/((2*PI)**3*32*M(4)**3)
C--Identical particle symmetry factor
              IF(I.EQ.J.OR.I.EQ.K.OR.J.EQ.K) CHARWD = CHARWD/2.0
C--Output the rate and particles to decay tables
 70           IF(CHARWD.GT.EPS) CALL RPMODA(CHARWD,29+10*L,POUT1,
     &                                              POUT2,POUT3)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      END
CDECK  ID>, RPINF1
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINF1(MI1,MI2,MI3,MI4,WIDTH,RESM,AIN,BIN)
C-----------------------------------------------------------------------
C      FUNCTION TO RETURN THE INTEGRATED AMPLITUDE SQUARE PIECE OF A
C      THREE BODY DECAY MATRIX ELEMENT
C-----------------------------------------------------------------------
       DOUBLE PRECISION LOW,UPP
       REAL MI1,MI2,MI3,MI4,WIDTH,RESM,AIN,BIN,RPINF1
       DOUBLE PRECISION SSDINT,RPINT1
       EXTERNAL SSDINT,RPINT1
C--common block to pass the masses,etc
       COMMON /SQUTRM/ M(4),GAM,MR,A,B
       DOUBLE PRECISION M,GAM,MR,A,B
C--Set the masses and couplings in the integration routine
       M(1) = DBLE(MI1)
       M(2) = DBLE(MI2)
       M(3) = DBLE(MI3)
       M(4) = DBLE(MI4)
       A    = DBLE(AIN)
       B    = DBLE(BIN)
       GAM  = DBLE(WIDTH)
       MR   = DBLE(ABS(RESM))
C--Perform the smoothing
       LOW = ATAN(((M(1)+M(2))**2-MR**2)/(GAM*MR))/(GAM*MR)
       UPP = ATAN(((M(4)-M(3))**2-MR**2)/(GAM*MR))/(GAM*MR)
C--Do the outer integral
       RPINF1 = 0.5E0*REAL(SSDINT(LOW,RPINT1,UPP))
       END
CDECK  ID>, RPINF2
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINF2(MI1,MI2,MI3,MI4,WIDTH1,WIDTH2,RESM1,RESM2,
     &                 AIN1,AIN2,BIN1,BIN2,TYPE)
C-----------------------------------------------------------------------
C      FUNCTION TO CALCULATE THE INTERFERENCE TERMS FOR 3-BODY
C      RPARITY VIOLATING DECAY RATES
C-----------------------------------------------------------------------
       IMPLICIT NONE
       DOUBLE PRECISION LOW,UPP
       REAL MI1,MI2,MI3,MI4,WIDTH1,WIDTH2,RESM1,RESM2,AIN1,AIN2
     &      ,BIN1,BIN2,RPINF2
       INTEGER I,TYPE
       DOUBLE PRECISION SSDINT,RPINT2,RPINT3
       EXTERNAL SSDINT,RPINT2,RPINT3
C--common block to pass the masses,etc
       COMMON /INFTRM/ M(4),MSQ(4),GAM(2),MR(2),MRSQ(2),A(2),B(2)
       DOUBLE PRECISION M,GAM,MR,MRSQ,A,B,MSQ
C--Set the masses and couplings in the integration routine
       M(1)   = DBLE(MI1)
       M(2)   = DBLE(MI2)
       M(3)   = DBLE(MI3)
       M(4)   = DBLE(MI4)
       A(1)   = DBLE(AIN1)
       A(2)   = DBLE(AIN2)
       B(1)   = DBLE(BIN1)
       B(2)   = DBLE(BIN2)
       GAM(1) = DBLE(WIDTH1)
       GAM(2) = DBLE(WIDTH2)
       MR(1)  = DBLE(RESM1)
       MR(2)  = DBLE(RESM2)
       DO I=1,4
         MSQ(I)=M(I)**2
       ENDDO
       DO I=1,2
         MRSQ(I)=MR(I)**2
       ENDDO
C--Perform the smoothing
       LOW = ATAN(((M(1)+M(2))**2-MR(1)**2)/(GAM(1)*ABS(MR(1))))
     &       /(GAM(1)*ABS(MR(1)))
       UPP = ATAN(((M(4)-M(3))**2-MR(1)**2)/(GAM(1)*ABS(MR(1))))
     &       /(GAM(1)*ABS(MR(1)))
C--Do the outer integral
       IF(TYPE.EQ.1) RPINF2 = REAL(SSDINT(LOW,RPINT2,UPP))
       IF(TYPE.EQ.2) RPINF2 = REAL(SSDINT(LOW,RPINT3,UPP))
       END
CDECK  ID>, RPINT1
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINT1(RHO)
C-----------------------------------------------------------------------
C      INTEGRAND FOR THE AMPLITUDE SQUARE PIECE OF THE THREE BODY DECAY
C-----------------------------------------------------------------------
       DOUBLE PRECISION RHO,X,RPINT1,M23MIN,M23MAX,E2STAR,E3STAR,
     &                  Y,LIMIT(2),E2MIMA,E3MIMA
       INTEGER          I
C--common block to pass the masses,etc
       COMMON /SQUTRM/ M(4),GAM,MR,A,B
       DOUBLE PRECISION M,GAM,MR,A,B,RPRTCH
       EXTERNAL RPRTCH
C--Calculate the change of variables
       X = MR**2+GAM*MR*TAN(GAM*MR*RHO)
C--Evaulate limits on the inner integral
       E2STAR = (X-M(1)**2+M(2)**2)/(2*SQRT(X))
       E3STAR = (M(4)**2-X-M(3)**2)/(2*SQRT(X))
       E2MIMA = RPRTCH(M(2),E2STAR)
       E3MIMA = RPRTCH(M(3),E3STAR)
       M23MAX = (E2STAR+E3STAR)**2-(E2MIMA-E3MIMA)**2
       M23MIN = (E2STAR+E3STAR)**2-(E2MIMA+E3MIMA)**2
C--Do the inner integral
       Y = M23MIN
       DO I=1,2
          LIMIT(I) =         Y*(X - M(1)**2 - M(2)**2)*
     -  (4*A*B*M(3)*M(4) + (A**2 + B**2)*(-X + M(3)**2 + M(4)**2))
          Y=M23MAX
       ENDDO
       RPINT1 = LIMIT(2)-LIMIT(1)
       END
CDECK  ID>, RPINT2
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINT2(RHO)
C-----------------------------------------------------------------------
C      FUNCTION FOR THE INTEGRAND FOR THE LIGHT/HEAVY INTERFERENCE
C      TERM IN A 3-BODY R-PARITY VIOLATING DECAY
C-----------------------------------------------------------------------
       DOUBLE PRECISION RHO,X,RPINT2,M23MIN,M23MAX,E2STAR,E3STAR,
     &                  Y,LIMIT(2),E2MIMA,E3MIMA
       INTEGER          I
C--common block to pass the masses,etc
       COMMON /INFTRM/ M(4),MSQ(4),GAM(2),MR(2),MRSQ(2),A(2),B(2)
       DOUBLE PRECISION M,MSQ,GAM,MR,MRSQ,A,B,RPRTCH
       EXTERNAL RPRTCH
C--Calculate the change of variables
       X = MR(1)**2+GAM(1)*ABS(MR(1))*TAN(GAM(1)*ABS(MR(1))*RHO)
C--Evaulate limits on the inner integral
       E2STAR = (X-M(1)**2+M(2)**2)/(2*SQRT(X))
       E3STAR = (M(4)**2-X-M(3)**2)/(2*SQRT(X))
       E2MIMA = RPRTCH(M(2),E2STAR)
       E3MIMA = RPRTCH(M(3),E3STAR)
       M23MAX = (E2STAR+E3STAR)**2-(E2MIMA-E3MIMA)**2
       M23MIN = (E2STAR+E3STAR)**2-(E2MIMA+E3MIMA)**2
C--Do the inner intergral
       Y = M23MIN
       DO I=1,2
          LIMIT(I) =         (Y*(X - MSQ(1) - MSQ(2))*
     -    (2*(A(2)*B(1) + A(1)*B(2))*M(3)*M(4) + 
     -      (A(1)*A(2) + B(1)*B(2))*(-X + MSQ(3) + MSQ(4)))*
     -    (GAM(1)*GAM(2)*MR(1)*MR(2) + (X - MR(1)**2)*(X - MR(2)**2)))
     -   /(GAM(2)**2*MR(2)**2 + (X - MR(2)**2)**2)
          Y=M23MAX
       ENDDO
       RPINT2 = LIMIT(2)-LIMIT(1)
       END
CDECK  ID>, RPINT3
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINT3(RHO)
C-----------------------------------------------------------------------
C      FUNCTION FOR INTEGRAND OF THE NORMAL INTERFERENCE TERM IN A 
C      3-BODY R-PARITY VIOLATING DECAY 
C-----------------------------------------------------------------------
       IMPLICIT NONE
       DOUBLE PRECISION RHO,X,RPINT3,M23MIN,M23MAX,E2STAR,E3STAR,
     &                  Y,LIMIT(2),E2MIMA,E3MIMA
       INTEGER          I
C--common block to pass the masses,etc
       COMMON /INFTRM/ M(4),MSQ(4),GAM(2),MR(2),MRSQ(2),A(2),B(2)
       DOUBLE PRECISION M,MSQ,GAM,MR,MRSQ,A,B,RPRTCH
       EXTERNAL RPRTCH
C--Calculate the change of variables
       X = MR(1)**2+GAM(1)*ABS(MR(1))*TAN(GAM(1)*ABS(MR(1))*RHO)
C--Evaulate limits on the inner integral
       E2STAR = (X-M(1)**2+M(2)**2)/(2*SQRT(X))
       E3STAR = (M(4)**2-X-M(3)**2)/(2*SQRT(X))
       E2MIMA = RPRTCH(M(2),E2STAR)
       E3MIMA = RPRTCH(M(3),E3STAR)
       M23MAX = (E2STAR+E3STAR)**2-(E2MIMA-E3MIMA)**2
       M23MIN = (E2STAR+E3STAR)**2-(E2MIMA+E3MIMA)**2
C--Do the inner integral
       Y = M23MIN
       DO I=1,2
          LIMIT(I) =Y*(X*B(1)*B(2)+A(1)*A(2)*M(1)*M(3)
     -      +A(1)*B(2)*M(1)*M(4))*
     -   (X-MRSQ(1))+(ATAN((Y-MRSQ(2))/(GAM(2)*SQRT(MRSQ(2))))*
     -     (X*A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2) + 
     -       X*A(2)*B(1)*GAM(1)*M(3)*M(4)*MR(1)*MR(2) - 
     -       X**2*B(1)*B(2)*GAM(2)*MRSQ(2) - 
     -       X*A(1)*A(2)*GAM(2)*M(1)*M(3)*MRSQ(2) - 
     -       X*A(1)*B(2)*GAM(2)*M(1)*M(4)*MRSQ(2) + 
     -       X*B(1)*B(2)*GAM(1)*MR(1)*MR(2)*MRSQ(2) + 
     -       A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2)*MRSQ(2) + 
     -       A(1)*B(2)*GAM(1)*M(1)*M(4)*MR(1)*MR(2)*MRSQ(2) + 
     -       X*B(1)*B(2)*GAM(2)*MRSQ(1)*MRSQ(2) + 
     -       A(1)*A(2)*GAM(2)*M(1)*M(3)*MRSQ(1)*MRSQ(2) + 
     -       A(1)*B(2)*GAM(2)*M(1)*M(4)*MRSQ(1)*MRSQ(2) - 
     -       A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2)*MSQ(1) - 
     -       A(2)*B(1)*GAM(1)*M(3)*M(4)*MR(1)*MR(2)*MSQ(1) - 
     -       A(1)*B(2)*GAM(1)*M(1)*M(4)*MR(1)*MR(2)*MSQ(2) - 
     -       A(2)*B(1)*GAM(1)*M(3)*M(4)*MR(1)*MR(2)*MSQ(2) - 
     -       A(1)*A(2)*GAM(1)*M(1)*M(3)*MR(1)*MR(2)*MSQ(3) - 
     -       A(1)*B(2)*GAM(1)*M(1)*M(4)*MR(1)*MR(2)*MSQ(3) - 
     -       B(1)*B(2)*GAM(1)*MR(1)*MR(2)*MSQ(1)*MSQ(3) - 
     - B(1)*B(2)*GAM(1)*MR(1)*MR(2)*MSQ(2)*MSQ(4)))/SQRT(MRSQ(2))+ 
     -  (LOG(Y**2 - 2*Y*MRSQ(2) + GAM(2)**2*MRSQ(2) + MRSQ(2)**2)*
     -     (X**2*A(1)*A(2)*M(1)*M(3) + X**2*A(2)*B(1)*M(3)*M(4) + 
     -       X*B(1)*B(2)*GAM(1)*GAM(2)*MR(1)*MR(2) + 
     -       A(1)*A(2)*GAM(1)*GAM(2)*M(1)*M(3)*MR(1)*MR(2) + 
     -       A(1)*B(2)*GAM(1)*GAM(2)*M(1)*M(4)*MR(1)*MR(2) - 
     - X*A(1)*A(2)*M(1)*M(3)*MRSQ(1)-X*A(2)*B(1)*M(3)*M(4)*MRSQ(1)+ 
     -       X**2*B(1)*B(2)*MRSQ(2) + X*A(1)*A(2)*M(1)*M(3)*MRSQ(2)+ 
     -    X*A(1)*B(2)*M(1)*M(4)*MRSQ(2)-X*B(1)*B(2)*MRSQ(1)*MRSQ(2)- 
     -       A(1)*A(2)*M(1)*M(3)*MRSQ(1)*MRSQ(2) - 
     -       A(1)*B(2)*M(1)*M(4)*MRSQ(1)*MRSQ(2) - 
     - X*A(1)*A(2)*M(1)*M(3)*MSQ(1) - X*A(2)*B(1)*M(3)*M(4)*MSQ(1) + 
     -       A(1)*A(2)*M(1)*M(3)*MRSQ(1)*MSQ(1) + 
     -       A(2)*B(1)*M(3)*M(4)*MRSQ(1)*MSQ(1) - 
     - X*A(1)*B(2)*M(1)*M(4)*MSQ(2) - X*A(2)*B(1)*M(3)*M(4)*MSQ(2) + 
     -       A(1)*B(2)*M(1)*M(4)*MRSQ(1)*MSQ(2) + 
     -       A(2)*B(1)*M(3)*M(4)*MRSQ(1)*MSQ(2) - 
     - X*A(1)*A(2)*M(1)*M(3)*MSQ(3) - X*A(1)*B(2)*M(1)*M(4)*MSQ(3) + 
     -     A(1)*A(2)*M(1)*M(3)*MRSQ(1)*MSQ(3) + 
     - A(1)*B(2)*M(1)*M(4)*MRSQ(1)*MSQ(3)-X*B(1)*B(2)*MSQ(1)*MSQ(3) + 
     -  B(1)*B(2)*MRSQ(1)*MSQ(1)*MSQ(3) - X*B(1)*B(2)*MSQ(2)*MSQ(4) + 
     - B(1)*B(2)*MRSQ(1)*MSQ(2)*MSQ(4)))/2.0D0
          Y=M23MAX
       ENDDO
       RPINT3 = LIMIT(2)-LIMIT(1)
       END
CDECK  ID>, RPMAIN
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C-----------------------------------------------------------------------
      SUBROUTINE RPMAIN 
C-----------------------------------------------------------------------
C     MAIN R-PARITY SUBROUTINE
C-----------------------------------------------------------------------
C     RPARTY   = logical, .TRUE. is consereved/.FALSE. is violated
C     LAMDA1   = LLE couplings
C     LAMDA2   = LUD couplings
C     LAMDA3   = UDD couplings
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
      CHARACTER*2 RYORN
C--Decide if want Rparity violating couplings
C-FP  Assume no Rparity violation, just return
      RPARTY=.TRUE.
      RETURN
      PRINT 100
 100  FORMAT(' R PARITY VIOLATION (Y/N)')
      READ(*,'(A1)')RYORN
      RPARTY = RYORN.NE.'Y'.AND.RYORN.NE.'y'
      IF(RPARTY) RETURN
C--If decide want rparity violation input the couplings
      PRINT 150
 150  FORMAT('Do you want lambda couplings(y/n)?')
      READ (*,'(A1)') RYORN
      IF(RYORN.EQ.'Y'.OR.RYORN.EQ.'y') THEN
        PRINT 200
 200    FORMAT(
     &      'ENTER LAMDBA IN ORDER 121,122,123,131,132,133,231,232,233')
        READ*,LAMDA1(1,2,1),LAMDA1(1,2,2),LAMDA1(1,2,3),LAMDA1(1,3,1),
     &        LAMDA1(1,3,2),LAMDA1(1,3,3),LAMDA1(2,3,1),LAMDA1(2,3,2),
     &        LAMDA1(2,3,3)
C--use the antisymmetry to find the rest
        DO K=1,3
          DO I=1,3
            DO J=1,3
              IF(I.EQ.J) LAMDA1(I,J,K) = 0.0E0
              IF(I.GT.J) LAMDA1(I,J,K) = -LAMDA1(J,I,K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO I=1,3
          DO J=1,3
            DO K=1,3
              LAMDA1(I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      PRINT 250
 250  FORMAT('Do you want lambda` couplings(y/n)?')
      READ(*,'(A1)') RYORN
      IF(RYORN.EQ.'Y'.OR.RYORN.EQ.'y') THEN
        PRINT 300
        PRINT 310
        PRINT 320
        PRINT 330
 300    FORMAT('ENTER LAMBDA` IN THE ORDER')
 310    FORMAT('111,112,113,121,122,123,131,132,133,')
 320    FORMAT('211,212,213,221,222,223,231,232,233,')
 330    FORMAT('311,312,313,321,322,323,331,332,333,')
        READ*,(((LAMDA2(I,J,K) ,K=1,3),J=1,3),I=1,3)
      ELSE
        DO I=1,3
          DO J=1,3
            DO K=1,3
              LAMDA2(I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      PRINT 350
 350  FORMAT('Do you want lambda`` couplings(y/n)?')
      READ (*,'(A1)') RYORN
      IF(RYORN.EQ.'Y'.OR.RYORN.EQ.'y') THEN
        PRINT 400
 400    FORMAT(
     &      'ENTER LAMDBA`` ORDER 112,113,123,212,213,223,312,313,323')
        READ*,LAMDA3(1,1,2),LAMDA3(1,1,3),LAMDA3(1,2,3),LAMDA3(2,1,2),
     &        LAMDA3(2,1,3),LAMDA3(2,2,3),LAMDA3(3,1,2),LAMDA3(3,1,3),
     &        LAMDA3(3,2,3)
C--use the antisymmetry to find the rest
        DO I=1,3
          DO J=1,3
            DO K=1,3
              IF(J.EQ.K) LAMDA3(I,J,K) = 0.0E0
              IF(J.GT.K) LAMDA3(I,J,K) = -LAMDA3(I,K,J)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO I=1,3
          DO J=1,3
            DO K=1,3
              LAMDA3(I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      CALL RPDECY
      END
CDECK  ID>, RPMODA
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
      SUBROUTINE RPMODA(RATE,PIN,POUT1,POUT2,POUT3)
C-----------------------------------------------------------------------
C     SUBROUTINE TO ADDED MODE TO DECAY TABLE
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NRPMD
      PARAMETER (NRPMD=5000)
      COMMON/RPARRT/NSSMD2,ISSMD2(NRPMD),JSSMD2(5,NRPMD),
     &              GSSMD2(NRPMD),BSSMD2(NRPMD)
      INTEGER NSSMD2,ISSMD2,JSSMD2,PIN,POUT1,POUT2,POUT3
      REAL    GSSMD2,BSSMD2,RATE
      SAVE /RPARRT/
      IF(RATE.GT.1E-40.AND.NSSMD2+1.LE.NRPMD) THEN
        NSSMD2=NSSMD2+1
        ISSMD2(NSSMD2)   = PIN
        JSSMD2(1,NSSMD2) = POUT1
        JSSMD2(2,NSSMD2) = POUT2
        JSSMD2(3,NSSMD2) = POUT3 
        JSSMD2(4,NSSMD2) = 0 
        JSSMD2(5,NSSMD2) = 0
        GSSMD2(NSSMD2) = RATE
        BSSMD2(NSSMD2) = 0
      ELSEIF(NSSMD2+1.GT.NRPMD) THEN
        print *,'TOO MANY R-PARITY VIOLATING MODES'
        print *,'INCREASE NRPMD AND RERUN'
        STOP
      ENDIF
      END
CDECK  ID>, RPNORM
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson 
C----------------------------------------------------------------------- 
      SUBROUTINE RPNORM  
C----------------------------------------------------------------------- 
C     SUBROUTINE TO REMOVE ALL MODES WITH BRANCHING RATIO LESS THAN 1E-5
C     AND TO ADD THE RPARITY VIOLATING MODES IF NEEDED
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C--Common block containing the SUSY parameters
C
#include "ssmode.inc"
C
C--Common Block to contain R-parity violating decay rates
C
      INTEGER NRPMD
      PARAMETER (NRPMD=5000)
      COMMON/RPARRT/NSSMD2,ISSMD2(NRPMD),JSSMD2(5,NRPMD),
     &              GSSMD2(NRPMD),BSSMD2(NRPMD)
      REAL GSSMD2,BSSMD2
      INTEGER ISSMD2,JSSMD2,NSSMD2
      SAVE /RPARRT/
#include "sslun.inc"
C  Common Block for the R-parity violating couplings
      COMMON/RSLASH/LAMDA1(3,3,3),LAMDA2(3,3,3),LAMDA3(3,3,3),RPARTY
      LOGICAL RPARTY
      REAL LAMDA1,LAMDA2,LAMDA3
      SAVE /RSLASH/
C--Local varaibles
      INTEGER SUSYLP,SUSYMN(9),SUSYMX(9),SUSYSP(9),I,J,NRMMDS,NRMMDR,
     &        ANUMRM,L,K
      REAL RATE,ZERO,RRATE,RPBR,EPS,RMBRRT,MINBR
      PARAMETER(ZERO=0.,MINBR=1E-5,EPS=1E-40)  
      DATA SUSYMN /31,51,21,41,29,82,86,30,39/
      DATA SUSYMX /36,56,26,46,29,84,86,60,49/
      DATA SUSYSP / 1, 1, 1, 1, 1, 1, 1,10,10/
C--Now we need to recalculate rates,etc.
      DO SUSYLP = 1,9
        DO 50 I=SUSYMN(SUSYLP),SUSYMX(SUSYLP),SUSYSP(SUSYLP) 
C--If Rparity conserved just remove modes less then MINBR
          IF(RPARTY) THEN 
            RMBRRT = ZERO
            NRMMDS = 0
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I.AND.BSSMOD(J).LT.MINBR) THEN
                RMBRRT = RMBRRT + BSSMOD(J)
                ISSMOD(J) = 100000
                NRMMDS = NRMMDS+1 
              ENDIF
            ENDDO
C--Now remove the modes we have removed, and renormalise BR's to 1
            IF(NRMMDS.GT.0.OR.NRMMDR.GT.0) THEN
              ANUMRM = 0
              RMBRRT = 1/(1-RMBRRT)
              DO J=1,NSSMOD
 10             IF(ISSMOD(J+ANUMRM).EQ.100000) ANUMRM=ANUMRM+1
                ISSMOD(J) = ISSMOD(J+ANUMRM)
                DO L=1,5
                  JSSMOD(L,J)=JSSMOD(L,J+ANUMRM)
                ENDDO
                GSSMOD(J)=GSSMOD(J+ANUMRM)
                BSSMOD(J)=BSSMOD(J+ANUMRM)
                IF(ISSMOD(J).EQ.100000) GOTO 10
                IF(ISSMOD(J).EQ.I) THEN
                  BSSMOD(J) = BSSMOD(J)*RMBRRT
                  GSSMOD(J) = GSSMOD(J)*RMBRRT
                ENDIF
              ENDDO
              ANUMRM = 0
              NSSMOD = NSSMOD-NRMMDS
            ENDIF
          ELSE
C--If R-paroty violated need to calculate rate and renormalise
C--Obtain the MSSM rate
            RATE  = ZERO
            RRATE = ZERO
            RPBR  = ZERO
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I) THEN
                RATE = RATE + GSSMOD(J)
              ENDIF
            ENDDO
C--Calculate the Rparity violating rate      
            DO J=1,NSSMD2
              IF(ISSMD2(J).EQ.I) THEN
                RRATE = RRATE + GSSMD2(J)
              ENDIF
            ENDDO        
            IF(RRATE.LT.EPS) GOTO 50
            RPBR = RATE/(RATE+RRATE)
            RATE = RATE+RRATE
C--Reset MSSM rates and branching ratios
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I) BSSMOD(J) = BSSMOD(J)*RPBR
            ENDDO
C--Calculate Rparity violating branching ratios
            DO J=1,NSSMD2
              IF(ISSMD2(J).EQ.I) BSSMD2(J) = GSSMD2(J)/RATE         
            ENDDO
C--Now remove any modes of the particle with branching ratio
C--less than MINBR
            RMBRRT = ZERO
            NRMMDS = 0
            NRMMDR = 0
            DO J=1,NSSMOD
              IF(ISSMOD(J).EQ.I.AND.BSSMOD(J).LT.MINBR) THEN
                RMBRRT = RMBRRT + BSSMOD(J)
                ISSMOD(J) = 100000
                NRMMDS = NRMMDS+1 
              ENDIF
            ENDDO
            IF(RRATE.LT.EPS) GOTO 50
            DO J=1,NSSMD2
              IF(ISSMD2(J).EQ.I.AND.BSSMD2(J).LT.MINBR) THEN
                RMBRRT = RMBRRT + BSSMD2(J)
                ISSMD2(J) = 100000
                NRMMDR = NRMMDR +1
              ENDIF
            ENDDO  
C--Now remove the modes we have removed, and renormalise BR's to 1
            IF(NRMMDS.GT.0.OR.NRMMDR.GT.0) THEN
              ANUMRM = 0
              RMBRRT = 1/(1-RMBRRT)
              DO J=1,NSSMOD
 20             IF(ISSMOD(J+ANUMRM).EQ.100000) ANUMRM=ANUMRM+1
                ISSMOD(J) = ISSMOD(J+ANUMRM)
                DO L=1,5
                  JSSMOD(L,J)=JSSMOD(L,J+ANUMRM)
                ENDDO
                GSSMOD(J)=GSSMOD(J+ANUMRM)
                BSSMOD(J)=BSSMOD(J+ANUMRM)
                IF(ISSMOD(J).EQ.100000) GOTO 20
                IF(ISSMOD(J).EQ.I) THEN
                  BSSMOD(J) = BSSMOD(J)*RMBRRT
                  GSSMOD(J) = GSSMOD(J)*RMBRRT
                ENDIF
              ENDDO
              ANUMRM = 0
              DO J=1,NSSMD2
 30             IF(ISSMD2(J+ANUMRM).EQ.100000) ANUMRM=ANUMRM+1
                  ISSMD2(J) = ISSMD2(J+ANUMRM)
                  DO L=1,5
                    JSSMD2(L,J)=JSSMD2(L,J+ANUMRM)
                  ENDDO
                  GSSMD2(J)=GSSMD2(J+ANUMRM)
                  BSSMD2(J)=BSSMD2(J+ANUMRM)
                  IF(ISSMD2(J).EQ.100000) GOTO 30
                  IF(ISSMD2(J).EQ.I) THEN
                  BSSMD2(J) = BSSMD2(J)*RMBRRT
                  GSSMD2(J) = GSSMD2(J)*RMBRRT
               ENDIF
              ENDDO
              NSSMOD = NSSMOD-NRMMDS
              NSSMD2 = NSSMD2-NRMMDR
            ENDIF
          ENDIF
 50     CONTINUE
      ENDDO
      IF(RPARTY) RETURN
C--NOW WE NEED TO COMBINE THE RPARITY AND MSSM MODES
      IF((NSSMOD+NSSMD2).GT.MXSS) THEN      
        WRITE(LOUT,*) 'WARNING REMOVING',NSSMOD+NSSMD2-MXSS,'MODES'
        print *,'WARNING EXCEEDS ISAJET NO OF MODES'
        print *,'REMOVING THE',NSSMOD+NSSMD2-MXSS,
     &          'MODES WITH LOWEST BRANCHING RATIO'
        print *,'RECOMMEND YOU RERUN WITH HIGHER MXSS'
      ENDIF
 110   IF((NSSMOD+NSSMD2).GT.MXSS) THEN
        RMBRRT = 1
        NRMMDS = 0
        DO J=1,NSSMOD
          IF(BSSMOD(J).LT.RMBRRT) THEN
            RMBRRT = BSSMOD(J)
            NRMMDS = J
          ENDIF
        ENDDO
        DO J=1,NSSMD2
          IF(BSSMD2(J).LT.RMBRRT) THEN
            RMBRRT = BSSMD2(J)
            NRMMDS = J+MXSS
          ENDIF
        ENDDO
C--remove mode with lowest branching ratio and rescale BR's
        RMBRRT = 1/(1-RMBRRT)
        IF(NRMMDS.LE.MXSS) NRMMDR = ISSMOD(NRMMDS)
        IF(NRMMDS.GT.MXSS) NRMMDR = ISSMD2(NRMMDS-MXSS)
        DO J=1,NSSMOD
          IF(ISSMOD(J).EQ.NRMMDR) THEN
            GSSMOD(J) = RMBRRT*GSSMOD(J)
            BSSMOD(J) = RMBRRT*BSSMOD(J)
          ENDIF
        ENDDO
        DO J=1,NSSMD2
          IF(ISSMD2(J).EQ.NRMMDR) THEN
            GSSMD2(J) = RMBRRT*GSSMD2(J)
            BSSMD2(J) = RMBRRT*BSSMD2(J)
          ENDIF
        ENDDO
        IF(NRMMDS.LE.MXSS) THEN
          DO J=NRMMDS,NSSMOD-1
            ISSMOD(J) = ISSMOD(J+1)
            DO L=1,5
              JSSMOD(L,J) = JSSMOD(L,J+1)
            ENDDO
            GSSMOD(J)=GSSMOD(J+1)
            BSSMOD(J)=BSSMOD(J+1)
          ENDDO
          NSSMOD = NSSMOD-1
        ELSE
          DO J=(NRMMDS-MXSS),NSSMD2-1
            ISSMD2(J) = ISSMD2(J+1)
            DO L=1,5
              JSSMD2(L,J) = JSSMD2(L,J+1)
            ENDDO
            GSSMD2(J) = GSSMD2(J+1)
            BSSMD2(J) = BSSMD2(J+1)
          ENDDO
          NSSMD2 = NSSMD2-1
        ENDIF   
      ENDIF   
      IF((NSSMD2+NSSMOD).GT.MXSS) GOTO 110
C--Now less than maximum so add the R parity violating modes
      DO J=1,NSSMD2
        NSSMOD=NSSMOD+1
        ISSMOD(NSSMOD) = ISSMD2(J)
        DO K=1,5
          JSSMOD(K,NSSMOD) = JSSMD2(K,J)
        ENDDO
        GSSMOD(NSSMOD) = GSSMD2(J)
        BSSMOD(NSSMOD) = BSSMD2(J)
        MSSMOD(NSSMOD) = 0
      ENDDO
      END
CDECK  ID>, RPRATE
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
      FUNCTION RPRATE(LAMCOL,M1,M2,M3)
C-----------------------------------------------------------------------
C     FUNCTION TO CALCULATE A 2 BODY R-PARITY VIOLATING DECAY RATE
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL LAMCOL,M1,M2,M3,RPRATE,PCM,PI
      PI = 3.1415926E0
      RPRATE = 0
      IF(M1.LT.(M2+M3)) RETURN
      RPRATE = LAMCOL*(M1**2-M2**2-M3**2)
      PCM = SQRT((M1**2-(M2+M3)**2)*(M1**2-(M2-M3)**2))/(2*M1)
      RPRATE = RPRATE*PCM/(8*PI*M1**2)
      END
CDECK  ID>, RPRTCH
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPRTCH(MASS,ENERGY)
C-----------------------------------------------------------------------
C      AVOID ERRORS DUE TO TAKING SQRT OF SMALL NEGATIVE NUMBERS
C      OCCURS DUE ROUNDING ERRORS
C-----------------------------------------------------------------------
       DOUBLE PRECISION RPRTCH,ENERGY,MASS
#include "sslun.inc"
       RPRTCH =ENERGY**2-MASS**2
       IF(RPRTCH.LT.0) THEN
         IF(RPRTCH/(ENERGY**2+MASS**2).LT.1D-20) THEN
           RPRTCH = 0.0D0
         ELSE
           WRITE(LOUT,*) 'WARNING SQRT OF NEGATIVE NUMBER',RPRTCH
           RPRTCH = 0.0D0
         ENDIF
       ENDIF
       RPRTCH = SQRT(RPRTCH)
       END
C-----------------------------------------------------------------------
C
C                         ISAWIG-HDECAY interface
C
C-----------------------------------------------------------------------
C  This block contains the code to interface HDECAY with ISAWIG so that
C  HDECAY can be used to give NLO Higgs widths and branching ratios
C-----------------------------------------------------------------------
C  This code is designed to work with HDECAY2.0/3.0
C  if you are using a later version of HDECAY there may be problems
CDECK  ID>, HDCYAD
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :  Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HDCYAD(HIGGS,WIDTH,BRFRAC,DECAY)
C-----------------------------------------------------------------------
C--Subroutine to add Higgs modes to the ISAJET decay tables
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Common block containing ISAJET decay modes
#include "ssmode.inc"
C
      DOUBLE PRECISION WIDTH,BRFRAC
      INTEGER HIGGS,DECAY(3),J
      NSSMOD = NSSMOD+1
      IF(NSSMOD.GT.MXSS) THEN
        WRITE(*,*) 'TOO MANY MODES'
        WRITE(*,*) 'RERUN WITH INCREASED MXSS'
        WRITE(*,*) 'STOPPING'
        STOP
      ENDIF
      ISSMOD(NSSMOD) = HIGGS
      DO J=1,3
        JSSMOD(J,NSSMOD) = DECAY(J)
      ENDDO
      DO J=4,5
        JSSMOD(J,NSSMOD) = 0
      ENDDO
      GSSMOD(NSSMOD) = WIDTH*BRFRAC
      BSSMOD(NSSMOD) = BRFRAC
      MSSMOD(NSSMOD) = 0
      END 
C
C   N.B. THERE ARE THREE VERSIONS OF THIS ROUTINE
C   ONE FOR HDECAY 2.0, ONE FOR HDECAY3.0 AND ONE DUMMY
C
#ifdef HDECAY3_X
CDECK  ID>, HDCYSY
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :  Peter Richardson based on HDECAY main program by M. Spira
C-----------------------------------------------------------------------
      SUBROUTINE HDCYSY(M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,
     &                     MER2,MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR) 
C-----------------------------------------------------------------------
C  Subroutine to interface with ISAWIG for MSSM Higgs Decays
C  for use with version 3 of HDECAY
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--varaibles passed from ISAWIG
      REAL M1,M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,MQL2,MUR2,MDR2,
     &     MEL,MER,MSQ,MUR,MDR
C--declarations for HDECAY variables
      DOUBLE PRECISION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     &                 GMSU(2),GMSD(2),GMSE(2),GMSN(2),
     &     AXMPL,AXMGD,HLBRGD,HABRGD,HHBRGD,HCBRGD,BHLSUSY,
     &     BHLST,BHLSB,BHLSTAU,BHHSUSY,BHHST,BHHSB,BHHSTAU,
     &     BHASUSY,BHCSUSY,BHCSTB,HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),
     .          HHBRSN(4,4),HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2),
     &     HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .     HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .     HABRSB,HHBRSQ,HHBRSQT,HCBRSQ,
     .     HCBRSQT,HLBRSQ,HLBRSQT,HLBRS,HLBRH,HLWDTH,
     &     HHBRS,HHBRH,HHWDTH,
     &     HCBRS,HCWDTH,SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .     SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH,HABRS,HAWDTH,AMGLU,
     &     AMQL1,AMUR1,AMDR1,AMEL1,AMER1,AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .     XLBB,XHBB,XABB,GAMC0,GAMT0,GAMT1,GAMW,GAMZ,GAT,
     &     GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .     GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .     GLPM,GHPM,B,A,XLAMBDA,AMC0,AMB0,AMT0,
     .     XMST1,XMST2,AMS,AMC,AMB,AMT,AMSB,GF,ALPH,AMTAU,AMMUON,AMZ,
     &     AMW,VUS,VCB,VUB,AMA,AML,AMH,AMCH,AMAR,AMSM,
     &     AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2,
     &     PI,ALSMZ,RVUB,TGBET,AMABEG,AMAEND,ACC,XITLA_HDEC
      INTEGER IGOLD,NFGG,N0,IHIGGS,NNLO,IPOLE,IMODEL,IONSH,IONWZ,
     &     IOFSUSY,NLOOP,NBER
      EXTERNAL XITLA_HDEC
C--HDECAY Common block
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2 
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/SFER1ST_HDEC/AMQL1,AMUR1,AMDR1,AMEL1,AMER1

      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ


      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/HABRS(10),HAWDTH
      COMMON/WIDTHHL_HDEC/HLBRS(11),HLBRH(3),HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRS(11),HHBRH(4),HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRS(9),HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD
      COMMON/WISFER_HDEC/BHLSUSY(7),BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .                   BHHSUSY(7),BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .                   BHASUSY(3),BHCSUSY(4),BHCSTB(2,2)
      COMMON/SMASS_HDEC/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
C--Common blocks from ISAJET which contains the masses and rates
C--and ISAJET SM parameters
C--Common block containing ISAJET masses and mixings
#include "sspar.inc"
C--Common block containing ISAJET SM paramaters
#include "sssm.inc"
C--Common block containing ISAJET decay modes
#include "ssmode.inc"
C--ISAJET SUGRA common block
#include "sugmg.inc"
C--local variables 
      INTEGER NOHGMD,I,SMDECS(3,11),SMDECA(3,10),TPDEC(3,2),HGDECH(3,5),
     &        CHRDEC(3,2,2),NETDEC(3,4,4),SMDECC(3,9),CNDEC(3,2,4),J,
     &        JMAX,SFRDC(3,15),SFRDCA(3,3),TAUDEC(3,2,2),BOTDEC(3,2,2),
     &        TOPDEC(3,2,2),SFRDCC(3,6),BTDEC(3,2,2),MODE,TOTMOD
      DOUBLE PRECISION EPS,TEMPWD
      REAL SSPOLE
      PARAMETER(EPS=1.0D-20)
      REAL EPS2
      PARAMETER (EPS2=1E-3)
      EXTERNAL SSPOLE
C--particle content of the various decay modes
      DATA SMDECS /  5, -5,  0, 16,-16,  0, 14,-14,  0,  3, -3,  0,
     &               4, -4,  0,  6, -6,  0,  9,  9,  0, 10, 10,  0,
     &              10, 90,  0, 80,-80,  0, 90, 90,  0/  
      DATA SMDECA /  5, -5,  0, 16,-16,  0, 14,-14,  0,  3, -3,  0,
     &               4, -4,  0,  6, -6,  0,  9,  9,  0, 10, 10,  0,
     &              10, 90,  0, 90, 82,  0/
      DATA TPDEC  /-80,  6, -5, 80, -6,  5/
      DATA HGDECH / 82, 82,  0, 84, 84,  0, 90, 84,  0, 80,-86,  0,
     &             -80, 86,  0/           
      DATA CHRDEC / 39,-39,  0, 39,-49,  0, 49,-39,  0, 49,-49,  0/
      DATA NETDEC / 30, 30,  0, 30, 40,  0, 30, 50,  0, 30, 60,  0,
     &              40, 30,  0, 40, 40,  0, 40, 50,  0, 40, 60,  0,
     &              50, 30,  0, 50, 40,  0, 50, 50,  0, 50, 60,  0,
     &              60, 30,  0, 60, 40,  0, 60, 50,  0, 60, 60,  0/
      DATA SMDECC / -5,  4,  0,-16, 15,  0,-14, 13,  0, -5,  1,  0,
     &              -3,  1,  0, -3,  4,  0, -5,  6,  0, 82, 80,  0,
     &              84, 80,  0/           
      DATA CNDEC /  39, 30,  0, 49, 30,  0, 39, 40,  0, 49, 40,  0,
     &              39, 50,  0, 49, 50,  0, 39, 60,  0, 49, 60,  0/
      DATA SFRDC /  31,-31,  0, 33,-33,  0, 35,-35,  0, 32,-32,  0,
     &              34,-34,  0, 52,-52,  0, 54,-54,  0, 21,-21,  0,
     &              24,-24,  0, 41,-41,  0, 44,-44,  0, 22,-22,  0,
     &              23,-23,  0, 42,-42,  0, 43,-43,  0/   
      DATA SFRDCA/  36,-56,  0, 25,-45,  0, 26,-46,  0/
      DATA TAUDEC/  36,-36,  0, 36,-56,  0, 56,-36,  0, 56,-56,  0/
      DATA BOTDEC/  25,-25,  0, 25,-45,  0, 45,-25,  0, 45,-45,  0/
      DATA TOPDEC/  26,-26,  0, 26,-46,  0, 46,-26,  0, 46,-46,  0/
      DATA SFRDCC/  31,-32,  0, 33,-34,  0, 35,-36,  0, 35,-56,  0,
     &              21,-22,  0, 24,-23,  0/
      DATA BTDEC /  26,-25,  0, 26,-45,  0, 46,-25,  0, 46,-45,  0/ 
C--define PI
      PI = 4*DATAN(1D0)
C--Output the authorship of HDECAY and warnings about the spectrum
      WRITE(*,*) 'Deleting ISAJET masses and Higgs decays'
      WRITE(*,*) 'Using HDECAY masses and Higgs decay modes'
      WRITE(*,*) 'HDECAY by A. Djouadi, J. Kalinowski, M. Spira'
      WRITE(*,*) 'Please refer to  Comput.Phys.Commun.108:56-74,1998'
C--Set the options we need
C--calculate all the SUSY Higgs decay modes
      IHIGGS = 5
C   IMODEL: USE SPECIFIC SUBROUTINE FOR MSSM HIGSS MASSES AND COUPLINGS
C           =1: CARENA ET AL., NUCL. PHYS. B461 (1996) 407 (SUBHPOLE)
C           =2: CARENA ET AL., PHYS. LETT. B355 (1995) 209 (SUBH)
C           =3: HABER ET AL.
C           =4: HEINEMEYER ET AL., HEP-PH/0002213 (FEYNHIGGSFAST1.2.2)
C
 500  WRITE(*,*) 
     &     'USE SPECIFIC SUBROUTINE FOR MSSM HIGSS MASSES AND COUPLINGS'
      WRITE(*,*) 
     &     '1: CARENA ET AL., NUCL. PHYS. B461 (1996) 407 (SUBHPOLE)'
      WRITE(*,*) 
     &     '2: CARENA ET AL., PHYS. LETT. B355 (1995) 209 (SUBH)'
      WRITE(*,*) 
     &     '3: HABER ET AL'
      WRITE(*,*) 
     &     '4: HEINEMEYER ET AL., HEP-PH/0002213 (FEYNHIGGSFAST1.2.2)'
      READ(*,*) IMODEL
      IF(IMODEL.LT.1.OR.IMODEL.GT.4) GOTO 500
      IF(IMODEL.EQ.3)THEN
       WRITE(6,*)'MU (UP TO THE SIGN) WILL BE IDENTIFIED WITH M_SQ...'
      ENDIF
C--if using SSRUN
      IF(M2.GE.1.E19) M2 =ALFA2I*SSPOLE(AMGLSS,AMZI**2,-ALFA3I)/ALFA3I
C--check the soft terms are the same for the 1st and 2nd generations
      IF(MEL2.LE.1.E19) THEN
        IF(ABS(MEL1-MEL2).GT.EPS2) WRITE(*,1000) 'left slepton'
        IF(ABS(MER1-MER2).GT.EPS2) WRITE(*,1000) 'right slepton'
        IF(ABS(MQL1-MQL2).GT.EPS2) WRITE(*,1000) 'left squark'
        IF(ABS(MUR1-MUR2).GT.EPS2) WRITE(*,1000) 'right up squark'
        IF(ABS(MDR1-MDR2).GT.EPS2) WRITE(*,1000) 'left up squark'
      ENDIF

C--Now set the parameters for HDECAY equal to the relevant ones from ISAJET
C--first the Standard Model Masses and Widths
      AMS    = DBLE(AMSTI)
      AMC    = DBLE(AMCHI)
      AMB    = DBLE(AMBTI)
      AMT    = DBLE(AMTPI)
      AMTAU  = DBLE(AMTAUI)
      AMMUON = DBLE(AMMUI)
      AMZ    = DBLE(AMZI)
      AMW    = DBLE(AMWI)
      GAMW   = DBLE(GAMWI)
      GAMZ   = DBLE(GAMZI)
C--Standard Model couplings
      ALSMZ  = DBLE(ALFA3I)
      ALPH   = 137.0359895D0
      GF     = SQRT(2.0D0)*PI/(2.0D0*ALPH*AMW**2*DBLE(SN2THWI))
C--Hard wire in the various CKM Matrix Elements(perhaps we should read these in)
      VUS    = 0.2205D0
      VCB    = 0.04D0
      RVUB   = 0.08D0
C--Now the SUSY parameters
      TGBET  = 1.0D0/DBLE(RV2V1)
      AMABEG = DBLE(AMHA)
      AMAEND = AMABEG
      AMU    = -DBLE(TWOM1)
      AM2    = DBLE(M2) 
      AMGLUINO = DBLE(AMGLSS)
      AMGLU = AMGLUINO
      AMEL1  = DBLE(MEL1)
      AMER1  = DBLE(MER1)
      AMQL1  = DBLE(MQL1)
      AMUR1  = DBLE(MUR1)
      AMDR1  = DBLE(MDR1)
      AMEL   = DBLE(MEL)
      AMER   = DBLE(MER)
      AMSQ   = DBLE(MSQ)
      AMUR   = DBLE(MUR)
      AMDR   = DBLE(MDR)
      AL     = DBLE(AAL)
      AU     = DBLE(AAT)
      AD     = DBLE(AAB)   
C--Options for HDECAY
 501  WRITE(*,*) 'OPTION FOR NNLO POLE MASS'
      WRITE(*,*) 
     &     '1: USE O(ALPHA_S) FORMULA FOR POLE MASS --> MSBAR MASS'
      WRITE(*,*) 
     &     '2: USE O(ALPHA_S**2) FORMULA FOR POLE MASS --> MSBAR MASS'
      READ(*,*) NNLO
      IF(NNLO.LT.1.OR.NNLO.GT.2) GOTO 501
 502  WRITE(*,*) 'OPTION FOR OFF_SHELL DECAY'
      WRITE(*,*) '0: INCLUDE OFF_SHELL DECAYS H,A --> T*T*, A --> Z*H,'
      WRITE(*,*) '   H --> W*H+,Z*A, H+ --> W*A, W*H, T*B'
      WRITE(*,*) '1: EXCLUDE THE OFF-SHELL DECAYS ABOVE'
      READ(*,*) IONSH
      IF(IONSH.GT.1.OR.IONSH.LT.0) GOTO 502
 503  WRITE(*,*) 'OPTION FOR GAUGE BOSON DECAYS'
      WRITE(*,*)
     &     '0: INCLUDE DOUBLE OFF-SHELL PAIR DECAYS PHI --> W*W*,Z*Z*'
      WRITE(*,*) 
     &     '1: INCLUDE ONLY SINGLE OFF-SHELL DECAYS PHI --> W*W,Z*Z'
      READ(*,*) IONWZ
      IF(IONWZ.GT.1.OR.IONWZ.LT.0) GOTO 503
 504  WRITE(*,*) 'OPTION FOR HIGGS MASSES'
      WRITE(*,*) '0: COMPUTES RUNNING HIGGS MASSES (FASTER)'
      WRITE(*,*) '1: COMPUTES POLE HIGGS MASSES'
      READ(*,*) IPOLE
      IF(IPOLE.GT.1.OR.IPOLE.LT.0) GOTO 504
 505  WRITE(*,*) 'OPTION FOR HIGGS DECAYS TO SUSY PARTICLES'
      WRITE(*,*) '0: INCLUDE DECAYS (AND LOOPS) INTO SUSY PARTICLES'
      WRITE(*,*) '1: EXCLUDE DECAYS (AND LOOPS) INTO SUSY PARTICLES'
      READ(*,*) IOFSUSY
      IF(IOFSUSY.GT.1.OR.IOFSUSY.LT.0) GOTO 505
 506  WRITE(*,*) 
     &     'NUMBER OF LIGHT FLAVORS INCLUDED IN THE GLUONIC DECAYS'
      WRITE(*,*) 'PHI --> GG* --> GQQ (3,4 OR 5)'
      READ(*,*)  NFGG
      IF(NFGG.LT.3.OR.NFGG.GT.5) GOTO 506
      VUB=RVUB*VCB
      ALPH=1.D0/ALPH
      AMSB = AMS
      IGOLD = 0
      AXMPL = 2.4E18
      AXMGD = AMGVSS
C--Now code for the HDECAY package
C--initialisation and checks
C--CHECK NFGG
      IF(NFGG.GT.5.OR.NFGG.LT.3)THEN
       WRITE(6,*)'NF-GG NOT VALID. TAKING THE DEFAULT NF-GG = 3....'
       NFGG = 3
      ENDIF
      AMC0=AMC
      AMB0=AMB
      AMT0=AMT
      ACC=1.D-8
      NLOOP=2
      XLAMBDA=XITLA_HDEC(NLOOP,ALSMZ,ACC)
      N0=5
      CALL ALSINI_HDEC(ACC)
C--INITIALIZE COEFFICIENTS FOR POLYLOGARITHMS
      NBER = 18
      CALL BERNINI_HDEC(NBER)
C
C--WRITE THE INPUT PARAMTERS 

      WRITE(6,8)'HIGGS    = ',IHIGGS
      WRITE(6,8)'MODEL    = ',IMODEL
      WRITE(6,9)'TGBET    = ',TGBET
      WRITE(6,9)'MABEG    = ',AMABEG
      WRITE(6,9)'MAEND    = ',AMAEND
      WRITE(6,9)'ALS(MZ)  = ',ALSMZ
      WRITE(6,9)'MSBAR(1) = ',AMS
      WRITE(6,9)'MC       = ',AMC
      WRITE(6,9)'MB       = ',AMB
      WRITE(6,9)'MT       = ',AMT
      WRITE(6,9)'MTAU     = ',AMTAU
      WRITE(6,9)'MMUON    = ',AMMUON
      WRITE(6,9)'ALPH     = ',1.D0/ALPH
      WRITE(6,9)'GF       = ',GF
      WRITE(6,9)'GAMW     = ',GAMW
      WRITE(6,9)'GAMZ     = ',GAMZ
      WRITE(6,9)'MZ       = ',AMZ
      WRITE(6,9)'MW       = ',AMW
      WRITE(6,9)'VUS      = ',VUS
      WRITE(6,9)'VCB      = ',VCB
      WRITE(6,9)'VUB/VCB  = ',RVUB
      WRITE(6,9)'MU       = ',AMU
      WRITE(6,9)'M2       = ',AM2
      WRITE(6,9)'MEL1      = ',AMEL1
      WRITE(6,9)'MER1      = ',AMER1
      WRITE(6,9)'MQL1      = ',AMQL1
      WRITE(6,9)'MUR1      = ',AMUR1
      WRITE(6,9)'MDR1      = ',AMDR1
      WRITE(6,9)'MEL      = ',AMEL
      WRITE(6,9)'MER      = ',AMER
      WRITE(6,9)'MSQ      = ',AMSQ
      WRITE(6,9)'MUR      = ',AMUR
      WRITE(6,9)'MDR      = ',AMDR
      WRITE(6,9)'AL       = ',AL
      WRITE(6,9)'AU       = ',AU
      WRITE(6,9)'AD       = ',AD
      WRITE(6,8)'NNLO (M) = ',NNLO
      WRITE(6,8)'ON-SHELL = ',IONSH
      WRITE(6,8)'ON-SH-WZ = ',IONWZ
      WRITE(6,8)'OFF-SUSY = ',IOFSUSY
      WRITE(6,8)'IPOLE    = ',IPOLE 
      WRITE(6,8)'NF-GG    = ',NFGG
      WRITE(6,9)'LAMBDA_5 = ',XLAMBDA
7     FORMAT(A11,I7)
8     FORMAT(A11,I4)
9     FORMAT(A11,G15.6)
C--set the pseudoscalar mass
      AMAR = AMABEG
      AMSM = AMAR
      AMA  = AMAR
C--calclaute the decay modes
      CALL HDEC(TGBET)
C--Now we need to enter the HDECAY modes in the ISAJET event record and
C--overwrite the masses of the HIGGS
C--Remove the Higgs decay modes from the ISAJET decay tables
      NOHGMD = 0
      DO I=1,NSSMOD
        IF(NOHGMD.NE.0) THEN
          ISSMOD(I-NOHGMD) = ISSMOD(I)
          DO J=1,5
            JSSMOD(J,I-NOHGMD) = JSSMOD(J,I)
          ENDDO
          GSSMOD(I-NOHGMD) = GSSMOD(I)
          BSSMOD(I-NOHGMD) = BSSMOD(I)
          MSSMOD(I-NOHGMD) = MSSMOD(I)
        ENDIF
        IF(ABS(ISSMOD(I)).EQ.86.OR.
     & (ISSMOD(I).GE.82.AND.ISSMOD(I).LE.84)) NOHGMD = NOHGMD+1
      ENDDO
      NSSMOD = NSSMOD-NOHGMD
C--Now write out the HDECAY modes into the ISAJET tables
C--this has to be done line by line
C--Decay to Standard Model particles of the lightest Higgs
      DO 10 I=1,11
        IF(HLBRS(I).LT.EPS) GOTO 10
        IF(I.NE.6.OR.AML.GT.2.0D0*AMT) THEN
          CALL HDCYAD(82,HLWDTH,HLBRS(I),SMDECS(1,I))
        ELSE
C--Decay to off mass shell top quark
          DO J=1,2
            CALL HDCYAD(82,HLWDTH,HLBRS(I)/2.0D0,TPDEC(1,J))
          ENDDO
        ENDIF
 10   CONTINUE
C--Decay to charginos of the lightest Higgs
      DO I=1,2
        DO 11 J=1,2
          IF(HLBRSC(I,J).LT.EPS) GOTO 11
          CALL HDCYAD(82,HLWDTH,HLBRSC(I,J),CHRDEC(1,I,J))
 11     CONTINUE
      ENDDO
C--Decay to neutralinos of the lightest scalar Higgs
      DO I=1,4
        DO 12 J=1,I
          TEMPWD = HLBRSN(I,J)
          IF(I.NE.J) TEMPWD = TEMPWD+ HLBRSN(J,I)
          IF(TEMPWD.LT.EPS) GOTO 12
          CALL HDCYAD(82,HLWDTH,TEMPWD,NETDEC(1,I,J))
 12     CONTINUE
      ENDDO
C--Decay to sfermions of the lightest scalar Higgs
C--first the decays to the first 2 generations(degenerate)
      TOTMOD = 0
      MODE = 0
      DO 13 I=1,7
        TOTMOD = TOTMOD+MODE
        MODE = 2
        IF(I.EQ.1) MODE = 3
        IF(BHLSUSY(I).LT.EPS) GOTO 13
        TEMPWD = BHLSUSY(I)/MODE   
        DO J=1,MODE
          CALL HDCYAD(82,HLWDTH,TEMPWD,SFRDC(1,TOTMOD+J))
        ENDDO     
 13   CONTINUE
C--now to the 3rd generation
      DO I=1,2
        DO 14 J=1,2
          IF(BHLSTAU(I,J).GE.EPS) 
     &      CALL HDCYAD(82,HLWDTH,BHLSTAU(I,J),TAUDEC(1,I,J))
          IF(BHLSB(I,J).GE.EPS) 
     &      CALL HDCYAD(82,HLWDTH,BHLSB(I,J),BOTDEC(1,I,J))
          IF(BHLST(I,J).GE.EPS) 
     &      CALL HDCYAD(82,HLWDTH,BHLST(I,J),TOPDEC(1,I,J))
 14     CONTINUE
      ENDDO
C--Decay to Standard Model particles of the heavy scalar Higgs
      DO 20 I=1,11
        IF(HHBRS(I).LT.EPS) GOTO 20
        IF(I.NE.6.OR.AMH.GT.2.0D0*AMT) THEN
          CALL HDCYAD(83,HHWDTH,HHBRS(I),SMDECS(1,I))
        ELSE
C--Decay to off mass shell top quark
          DO J=1,2
            CALL HDCYAD(83,HHWDTH,HHBRS(I)/2.0D0,TPDEC(1,J))
          ENDDO
        ENDIF
 20   CONTINUE 
C--Decay to charginos of the heaviest scalar Higgs
      DO I=1,2
        DO 21 J=1,2
          IF(HHBRSC(I,J).LT.EPS) GOTO 21
          CALL HDCYAD(83,HHWDTH,HHBRSC(I,J),CHRDEC(1,I,J))
 21     CONTINUE
      ENDDO  
C--Decay to neutralinos of the heaviest scalar Higgs
      DO I=1,4
        DO 22 J=1,I
          TEMPWD = HHBRSN(I,J)
          IF(I.NE.J) TEMPWD = TEMPWD+ HHBRSN(J,I)
          IF(TEMPWD.LT.EPS) GOTO 22
          CALL HDCYAD(83,HHWDTH,TEMPWD,NETDEC(1,I,J))
 22     CONTINUE
      ENDDO   
C--Decay to sfermions of the heaviest scalar Higgs
C--1st two generations
      TOTMOD = 0
      MODE   = 0
      DO 23 I=1,7
        TOTMOD = TOTMOD+MODE
        MODE = 2
        IF(I.EQ.1) MODE = 3
        IF(BHHSUSY(I).LT.EPS) GOTO 23
        TEMPWD = BHHSUSY(I)/MODE     
        DO J=1,MODE
          CALL HDCYAD(83,HHWDTH,TEMPWD,SFRDC(1,TOTMOD+J))
        ENDDO   
 23   CONTINUE
C--now to the 3rd generation
      DO I=1,2
        DO 24 J=1,2
          IF(BHHSTAU(I,J).GE.EPS) 
     &      CALL HDCYAD(83,HHWDTH,BHHSTAU(I,J),TAUDEC(1,I,J))
          IF(BHHSB(I,J).GE.EPS) 
     &      CALL HDCYAD(83,HHWDTH,BHHSB(I,J),BOTDEC(1,I,J))
          IF(BHHST(I,J).GE.EPS) 
     &      CALL HDCYAD(83,HHWDTH,BHHST(I,J),TOPDEC(1,I,J))
 24     CONTINUE
      ENDDO
C--Decay to Higgs of heavy scalar Higgs
      DO 25 I=1,4
        IF(HHBRH(I).LT.EPS) GOTO 25
        IF(I.NE.4) THEN
          JMAX = 1
          TEMPWD = HHBRH(I)
        ELSE
          JMAX = 2
          TEMPWD = HHBRH(I)/2.0D0
        ENDIF
        DO J=1,JMAX
          CALL HDCYAD(83,HHWDTH,TEMPWD,HGDECH(1,I+J-1))
        ENDDO
 25   CONTINUE
C--Decay to Standard Model particles of the pseudoscalar Higgs
      DO 30 I=1,10
        IF(HABRS(I).LT.EPS) GOTO 30
        IF(I.NE.6.OR.AMA.GT.2.0D0*AMT) THEN
          CALL HDCYAD(84,HAWDTH,HABRS(I),SMDECA(1,I))
        ELSE
C--Decay to off mass shell top quark
          DO J=1,2
            CALL HDCYAD(84,HAWDTH,HABRS(I)/2.0D0,TPDEC(1,J))
          ENDDO
        ENDIF
 30   CONTINUE  
C--Decay to charginos of the pseudoscalar Higgs
      DO I=1,2
        DO 31 J=1,2
          IF(HABRSC(I,J).LT.EPS) GOTO 31
          CALL HDCYAD(84,HAWDTH,HABRSC(I,J),CHRDEC(1,I,J))
 31     CONTINUE
      ENDDO    
C--Decay to neutralinos of the pseudoscalar Higgs
      DO I=1,4
        DO 32 J=1,I
          TEMPWD = HABRSN(I,J)
          IF(I.NE.J) TEMPWD = TEMPWD+ HABRSN(J,I)
          IF(TEMPWD.LT.EPS) GOTO 32
          CALL HDCYAD(84,HAWDTH,TEMPWD,NETDEC(1,I,J))
 32     CONTINUE
      ENDDO   
C--Decays to sfermions of the pseudoscalar Higgs
      DO 33 I=1,3
        IF(BHASUSY(I).GE.EPS) 
     &    CALL HDCYAD(84,HAWDTH,BHASUSY(I),SFRDCA(1,I))
 33   CONTINUE
C--Decay to Standard Model particles of the charged Higgs
      DO 40 I=1,9
        IF(HCBRS(I).LT.EPS) GOTO 40
        CALL HDCYAD(86,HCWDTH,HCBRS(I),SMDECC(1,I))
 40   CONTINUE
C--Decay to charginos and neutralino of the charged Higgs
      DO I=1,2
        DO 41 J=1,4
        IF(HCBRSU(I,J).LT.EPS) GOTO 41 
        CALL HDCYAD(86,HCWDTH,HCBRSU(I,J),CNDEC(1,I,J))
 41     CONTINUE
      ENDDO
C--Decay to first 2 genration sfermions and stau
      TOTMOD = 0
      MODE = 0
      DO 42 I=1,4
        TOTMOD = TOTMOD+MODE
        MODE = 1
        IF(I.EQ.2.OR.I.EQ.4) MODE = 2
        IF(BHCSUSY(I).LT.EPS) GOTO 42
        TEMPWD = BHCSUSY(I)/MODE
        DO J=1,MODE
          CALL HDCYAD(86,HCWDTH,TEMPWD,SFRDCC(1,TOTMOD+J))
        ENDDO
 42   CONTINUE
C--Now the decay to sbottom and stop
      DO I = 1,2
        DO 43 J=1,2
          IF(BHCSTB(I,J).LT.EPS) GOTO 43
          CALL HDCYAD(86,HCWDTH,BHCSTB(I,J),BTDEC(1,I,J))
 43     CONTINUE
      ENDDO
C--Now reset the ISAJET HIGGS MASSES with the HDECAY ones
      AMHL = REAL(AML)
      AMHH = REAL(AMH)
      AMHA = REAL(AMA)
      AMHC = REAL(AMCH)
      RETURN
 1000 FORMAT('WARNING ',A15,
     &        ' SOFT TERM IS NOT THE SAME FOR 1st AND 2nd GENERATION')
      END
#endif
#ifdef HDECAY2_X
CDECK  ID>, HDCYSY
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :  Peter Richardson based on HDECAY main program by M. Spira
C-----------------------------------------------------------------------
      SUBROUTINE HDCYSY(M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,
     &                     MER2,MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR) 
C-----------------------------------------------------------------------
C  Subroutine to interface with ISAWIG for MSSM Higgs Decays
C  for use with version 2.0 of HDECAY
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--varaibles passed from ISAWIG
      REAL M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,MQL2,MUR2,MDR2,
     &     MEL,MER,MSQ,MUR,MDR
C--declarations for HDECAY variables
      DOUBLE PRECISION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .                 GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DOUBLE PRECISION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .                 HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .                 HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2)
      DOUBLE PRECISION AMT,AMSB,AMMUON,AMZ,
     &                 AMW,VUS,VCB,VUB,AMSM,AMA,AML,AMH,AMCH,AMAR,AMEL,
     &                 AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2,AMQL1,AMUR1,
     &                 AMDR1,AMEL1,AMER1,GAMC0,GAMT0,GAMT1,GAMW,GAMZ,
     &                 GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,GHHH,GLLL,GHLL,
     &                 GLHH,GHAA,GLAA,GLVV,GHVV,GLPM,GHPM,B,A,XLAMBDA,
     &                 AMC0,AMB0,AMT0,SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,
     &                 SMBRT,SMBRG,SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH,
     &                 HABRS,HAWDTH,HHBRS,HHBRH,HHWDTH,HCBRS,HCWDTH,
     &                 HLBRS,HLBRH,HLWDTH,HLBRCHT,HHBRCHT,HABRCHT,AMB,
     &                 HLBRNET,HHBRNET,HABRNET,HCBRCNT,HLBRSL,HHBRSL,
     &                 HCBRSL,HABRSL,HABRST,HABRSB,HHBRSQ,XITLA,GF,ALPH,
     &                 HHBRSQT,HCBRSQ,HCBRSQT,HLBRSQ,HLBRSQT,AMS,AMC,
     &                 BHLSUSY,BHLST,BHLSB,BHLSTAU,BHHSUSY,BHHST,
     &                 BHHSB,BHHSTAU,BHASUSY,BHCSUSY,AMTAU,BHCSTB,
     &                 PI,ALSMZ,RVUB,TGBET,AMABEG,AMAEND,ACC                 
      INTEGER N0,IHIGGS,NNLO,IPOLE,IONSH,IONWZ,IOFSUSY,NFGG,NLOOP,NBER
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/STRANGE/AMSB
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR/VUS,VCB,VUB
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/SFER1ST/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG/IHIGGS,NNLO,IPOLE
      COMMON/ONSHELL/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH/NFGG
      COMMON/WIDTHSM/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,SMBRGA,
     .               SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA/HABRS(10),HAWDTH
      COMMON/WIDTHHL/HLBRS(11),HLBRH(3),HLWDTH
      COMMON/WIDTHHH/HHBRS(11),HHBRH(4),HHWDTH
      COMMON/WIDTHHC/HCBRS(9),HCWDTH
      COMMON/WISUSY/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,HCBRSU,
     .              HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,HABRNET,
     .              HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,HABRSB,
     .              HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER/BHLSUSY(7),BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSUSY(7),BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASUSY(3),BHCSUSY(4),BHCSTB(2,2)
      COMMON/SMASS/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
C--Common blocks from ISAJET which contains the masses and rates
C--and ISAJET SM parameters
C--Common block containing ISAJET masses and mixings
#include "sspar.inc"
C--Common block containing ISAJET SM paramaters
#include "sssm.inc"
C--Common block containing ISAJET decay modes
#include "ssmode.inc"
C--ISAJET SUGRA common block
#include "sugmg.inc"
C--local variables 
      INTEGER NOHGMD,I,SMDECS(3,11),SMDECA(3,10),TPDEC(3,2),HGDECH(3,5),
     &        CHRDEC(3,2,2),NETDEC(3,4,4),SMDECC(3,9),CNDEC(3,2,4),J,
     &        JMAX,SFRDC(3,15),SFRDCA(3,3),TAUDEC(3,2,2),BOTDEC(3,2,2),
     &        TOPDEC(3,2,2),SFRDCC(3,6),BTDEC(3,2,2),MODE,TOTMOD
      DOUBLE PRECISION EPS,TEMPWD
      REAL SSPOLE
      PARAMETER(EPS=1.0D-20)
      REAL EPS2
      PARAMETER (EPS2=1E-3)
      EXTERNAL SSPOLE
C--particle content of the various decay modes
      DATA SMDECS /  5, -5,  0, 16,-16,  0, 14,-14,  0,  3, -3,  0,
     &               4, -4,  0,  6, -6,  0,  9,  9,  0, 10, 10,  0,
     &              10, 90,  0, 80,-80,  0, 90, 90,  0/  
      DATA SMDECA /  5, -5,  0, 16,-16,  0, 14,-14,  0,  3, -3,  0,
     &               4, -4,  0,  6, -6,  0,  9,  9,  0, 10, 10,  0,
     &              10, 90,  0, 90, 82,  0/
      DATA TPDEC  /-80,  6, -5, 80, -6,  5/
      DATA HGDECH / 82, 82,  0, 84, 84,  0, 90, 84,  0, 80,-86,  0,
     &             -80, 86,  0/           
      DATA CHRDEC / 39,-39,  0, 39,-49,  0, 49,-39,  0, 49,-49,  0/
      DATA NETDEC / 30, 30,  0, 30, 40,  0, 30, 50,  0, 30, 60,  0,
     &              40, 30,  0, 40, 40,  0, 40, 50,  0, 40, 60,  0,
     &              50, 30,  0, 50, 40,  0, 50, 50,  0, 50, 60,  0,
     &              60, 30,  0, 60, 40,  0, 60, 50,  0, 60, 60,  0/
      DATA SMDECC / -5,  4,  0,-16, 15,  0,-14, 13,  0, -5,  1,  0,
     &              -3,  1,  0, -3,  4,  0, -5,  6,  0, 82, 80,  0,
     &              84, 80,  0/           
      DATA CNDEC /  39, 30,  0, 49, 30,  0, 39, 40,  0, 49, 40,  0,
     &              39, 50,  0, 49, 50,  0, 39, 60,  0, 49, 60,  0/
      DATA SFRDC /  31,-31,  0, 33,-33,  0, 35,-35,  0, 32,-32,  0,
     &              34,-34,  0, 52,-52,  0, 54,-54,  0, 21,-21,  0,
     &              24,-24,  0, 41,-41,  0, 44,-44,  0, 22,-22,  0,
     &              23,-23,  0, 42,-42,  0, 43,-43,  0/   
      DATA SFRDCA/  36,-56,  0, 25,-45,  0, 26,-46,  0/
      DATA TAUDEC/  36,-36,  0, 36,-56,  0, 56,-36,  0, 56,-56,  0/
      DATA BOTDEC/  25,-25,  0, 25,-45,  0, 45,-25,  0, 45,-45,  0/
      DATA TOPDEC/  26,-26,  0, 26,-46,  0, 46,-26,  0, 46,-46,  0/
      DATA SFRDCC/  31,-32,  0, 33,-34,  0, 35,-36,  0, 35,-56,  0,
     &              21,-22,  0, 24,-23,  0/
      DATA BTDEC /  26,-25,  0, 26,-45,  0, 46,-25,  0, 46,-45,  0/ 
C--define PI
      PI = 4*DATAN(1D0)
C--Output the authorship of HDECAY and warnings about the spectrum
      WRITE(*,*) 'Deleting ISAJET masses and Higgs decays'
      WRITE(*,*) 'Using HDECAY masses and Higgs decay modes'
      WRITE(*,*) 'HDECAY by A. Djouadi, J. Kalinowski, M. Spira'
      WRITE(*,*) 'Please refer to  Comput.Phys.Commun.108:56-74,1998'
C--Set the options we need
C--calculate all the SUSY Higgs decay modes
      IHIGGS = 5
C--if using SSRUN
      IF(M2.GE.1.E19) M2 =ALFA2I*SSPOLE(AMGLSS,AMZI**2,-ALFA3I)/ALFA3I
C--check the soft terms are the same for the 1st and 2nd generations
      IF(MEL2.LE.1.E19) THEN
        IF(ABS(MEL1-MEL2).GT.EPS2) WRITE(*,1000) 'left slepton'
        IF(ABS(MER1-MER2).GT.EPS2) WRITE(*,1000) 'right slepton'
        IF(ABS(MQL1-MQL2).GT.EPS2) WRITE(*,1000) 'left squark'
        IF(ABS(MUR1-MUR2).GT.EPS2) WRITE(*,1000) 'right up squark'
        IF(ABS(MDR1-MDR2).GT.EPS2) WRITE(*,1000) 'left up squark'
      ENDIF
C--Now set the parameters for HDECAY equal to the relevant ones from ISAJET
C--first the Standard Model Masses and Widths
      AMS    = DBLE(AMSTI)
      AMC    = DBLE(AMCHI)
      AMB    = DBLE(AMBTI)
      AMT    = DBLE(AMTPI)
      AMTAU  = DBLE(AMTAUI)
      AMMUON = DBLE(AMMUI)
      AMZ    = DBLE(AMZI)
      AMW    = DBLE(AMWI)
      GAMW   = DBLE(GAMWI)
      GAMZ   = DBLE(GAMZI)
C--Standard Model couplings
      ALSMZ  = DBLE(ALFA3I)
      ALPH   = 137.0359895D0
      GF     = SQRT(2.0D0)*PI/(2.0D0*ALPH*AMW**2*DBLE(SN2THWI))
C--Hard wire in the various CKM Matrix Elements(perhaps we should read these in)
      VUS    = 0.2205D0
      VCB    = 0.04D0
      RVUB   = 0.08D0
C--Now the SUSY parameters
      TGBET  = 1.0D0/DBLE(RV2V1)
      AMABEG = DBLE(AMHA)
      AMAEND = AMABEG
      AMU    = -DBLE(TWOM1)
      AM2    = DBLE(M2) 
      AMEL1  = DBLE(MEL1)
      AMER1  = DBLE(MER1)
      AMQL1  = DBLE(MQL1)
      AMUR1  = DBLE(MUR1)
      AMDR1  = DBLE(MDR1)
      AMEL   = DBLE(MEL)
      AMER   = DBLE(MER)
      AMSQ   = DBLE(MSQ)
      AMUR   = DBLE(MUR)
      AMDR   = DBLE(MDR)
      AL     = DBLE(AAL)
      AU     = DBLE(AAT)
      AD     = DBLE(AAB)   
C--Options for HDECAY
      NNLO   = 1
      IONSH  = 0
      IONWZ  = 0
      IPOLE  = 1
      IOFSUSY= 0
      NFGG   = 5
      VUB=RVUB*VCB
      ALPH=1.D0/ALPH
      AMSB = AMS
9     FORMAT(A11,G15.6)
C--Now code for the HDECAY package
C--initialisation and checks
C--CHECK NFGG
      IF(NFGG.GT.5.OR.NFGG.LT.3)THEN
       WRITE(6,*)'NF-GG NOT VALID. TAKING THE DEFAULT NF-GG = 3....'
       NFGG = 3
      ENDIF
      AMC0=AMC
      AMB0=AMB
      AMT0=AMT
      ACC=1.D-8
      NLOOP=2
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      N0=5
      CALL ALSINI(ACC)
C--INITIALIZE COEFFICIENTS FOR POLYLOGARITHMS
      NBER = 18
      CALL BERNINI(NBER)
C--set the pseudoscalar mass
      AMAR = AMABEG
      AMSM = AMAR
      AMA  = AMAR
C--SUSY MASS and coupling calculations
      CALL SUSYCP(TGBET)
C--calclaute the decay modes
      CALL HDEC(TGBET)
C--Now we need to enter the HDECAY modes in the ISAJET event record and
C--overwrite the masses of the HIGGS
C--Remove the Higgs decay modes from the ISAJET decay tables
      NOHGMD = 0
      DO I=1,NSSMOD
        IF(NOHGMD.NE.0) THEN
          ISSMOD(I-NOHGMD) = ISSMOD(I)
          DO J=1,5
            JSSMOD(J,I-NOHGMD) = JSSMOD(J,I)
          ENDDO
          GSSMOD(I-NOHGMD) = GSSMOD(I)
          BSSMOD(I-NOHGMD) = BSSMOD(I)
          MSSMOD(I-NOHGMD) = MSSMOD(I)
        ENDIF
        IF(ABS(ISSMOD(I)).EQ.86.OR.
     & (ISSMOD(I).GE.82.AND.ISSMOD(I).LE.84)) NOHGMD = NOHGMD+1
      ENDDO
      NSSMOD = NSSMOD-NOHGMD
C--Now write out the HDECAY modes into the ISAJET tables
C--this has to be done line by line
C--Decay to Standard Model particles of the lightest Higgs
      DO 10 I=1,11
        IF(HLBRS(I).LT.EPS) GOTO 10
        IF(I.NE.6.OR.AML.GT.2.0D0*AMT) THEN
          CALL HDCYAD(82,HLWDTH,HLBRS(I),SMDECS(1,I))
        ELSE
C--Decay to off mass shell top quark
          DO J=1,2
            CALL HDCYAD(82,HLWDTH,HLBRS(I)/2.0D0,TPDEC(1,J))
          ENDDO
        ENDIF
 10   CONTINUE
C--Decay to charginos of the lightest Higgs
      DO I=1,2
        DO 11 J=1,2
          IF(HLBRSC(I,J).LT.EPS) GOTO 11
          CALL HDCYAD(82,HLWDTH,HLBRSC(I,J),CHRDEC(1,I,J))
 11     CONTINUE
      ENDDO
C--Decay to neutralinos of the lightest scalar Higgs
      DO I=1,4
        DO 12 J=1,I
          TEMPWD = HLBRSN(I,J)
          IF(I.NE.J) TEMPWD = TEMPWD+ HLBRSN(J,I)
          IF(TEMPWD.LT.EPS) GOTO 12
          CALL HDCYAD(82,HLWDTH,TEMPWD,NETDEC(1,I,J))
 12     CONTINUE
      ENDDO
C--Decay to sfermions of the lightest scalar Higgs
C--first the decays to the first 2 generations(degenerate)
      TOTMOD = 0
      MODE = 0
      DO 13 I=1,7
        TOTMOD = TOTMOD+MODE
        MODE = 2
        IF(I.EQ.1) MODE = 3
        IF(BHLSUSY(I).LT.EPS) GOTO 13
        TEMPWD = BHLSUSY(I)/MODE   
        DO J=1,MODE
          CALL HDCYAD(82,HLWDTH,TEMPWD,SFRDC(1,TOTMOD+J))
        ENDDO     
 13   CONTINUE
C--now to the 3rd generation
      DO I=1,2
        DO 14 J=1,2
          IF(BHLSTAU(I,J).GE.EPS) 
     &      CALL HDCYAD(82,HLWDTH,BHLSTAU(I,J),TAUDEC(1,I,J))
          IF(BHLSB(I,J).GE.EPS) 
     &      CALL HDCYAD(82,HLWDTH,BHLSB(I,J),BOTDEC(1,I,J))
          IF(BHLST(I,J).GE.EPS) 
     &      CALL HDCYAD(82,HLWDTH,BHLST(I,J),TOPDEC(1,I,J))
 14     CONTINUE
      ENDDO
C--Decay to Standard Model particles of the heavy scalar Higgs
      DO 20 I=1,11
        IF(HHBRS(I).LT.EPS) GOTO 20
        IF(I.NE.6.OR.AMH.GT.2.0D0*AMT) THEN
          CALL HDCYAD(83,HHWDTH,HHBRS(I),SMDECS(1,I))
        ELSE
C--Decay to off mass shell top quark
          DO J=1,2
            CALL HDCYAD(83,HHWDTH,HHBRS(I)/2.0D0,TPDEC(1,J))
          ENDDO
        ENDIF
 20   CONTINUE 
C--Decay to charginos of the heaviest scalar Higgs
      DO I=1,2
        DO 21 J=1,2
          IF(HHBRSC(I,J).LT.EPS) GOTO 21
          CALL HDCYAD(83,HHWDTH,HHBRSC(I,J),CHRDEC(1,I,J))
 21     CONTINUE
      ENDDO  
C--Decay to neutralinos of the heaviest scalar Higgs
      DO I=1,4
        DO 22 J=1,I
          TEMPWD = HHBRSN(I,J)
          IF(I.NE.J) TEMPWD = TEMPWD+ HHBRSN(J,I)
          IF(TEMPWD.LT.EPS) GOTO 22
          CALL HDCYAD(83,HHWDTH,TEMPWD,NETDEC(1,I,J))
 22     CONTINUE
      ENDDO   
C--Decay to sfermions of the heaviest scalar Higgs
C--1st two generations
      TOTMOD = 0
      MODE   = 0
      DO 23 I=1,7
        TOTMOD = TOTMOD+MODE
        MODE = 2
        IF(I.EQ.1) MODE = 3
        IF(BHHSUSY(I).LT.EPS) GOTO 23
        TEMPWD = BHHSUSY(I)/MODE     
        DO J=1,MODE
          CALL HDCYAD(83,HHWDTH,TEMPWD,SFRDC(1,TOTMOD+J))
        ENDDO   
 23   CONTINUE
C--now to the 3rd generation
      DO I=1,2
        DO 24 J=1,2
          IF(BHHSTAU(I,J).GE.EPS) 
     &      CALL HDCYAD(83,HHWDTH,BHHSTAU(I,J),TAUDEC(1,I,J))
          IF(BHHSB(I,J).GE.EPS) 
     &      CALL HDCYAD(83,HHWDTH,BHHSB(I,J),BOTDEC(1,I,J))
          IF(BHHST(I,J).GE.EPS) 
     &      CALL HDCYAD(83,HHWDTH,BHHST(I,J),TOPDEC(1,I,J))
 24     CONTINUE
      ENDDO
C--Decay to Higgs of heavy scalar Higgs
      DO 25 I=1,4
        IF(HHBRH(I).LT.EPS) GOTO 25
        IF(I.NE.4) THEN
          JMAX = 1
          TEMPWD = HHBRH(I)
        ELSE
          JMAX = 2
          TEMPWD = HHBRH(I)/2.0D0
        ENDIF
        DO J=1,JMAX
          CALL HDCYAD(83,HHWDTH,TEMPWD,HGDECH(1,I+J-1))
        ENDDO
 25   CONTINUE
C--Decay to Standard Model particles of the pseudoscalar Higgs
      DO 30 I=1,10
        IF(HABRS(I).LT.EPS) GOTO 30
        IF(I.NE.6.OR.AMA.GT.2.0D0*AMT) THEN
          CALL HDCYAD(84,HAWDTH,HABRS(I),SMDECA(1,I))
        ELSE
C--Decay to off mass shell top quark
          DO J=1,2
            CALL HDCYAD(84,HAWDTH,HABRS(I)/2.0D0,TPDEC(1,J))
          ENDDO
        ENDIF
 30   CONTINUE  
C--Decay to charginos of the pseudoscalar Higgs
      DO I=1,2
        DO 31 J=1,2
          IF(HABRSC(I,J).LT.EPS) GOTO 31
          CALL HDCYAD(84,HAWDTH,HABRSC(I,J),CHRDEC(1,I,J))
 31     CONTINUE
      ENDDO    
C--Decay to neutralinos of the pseudoscalar Higgs
      DO I=1,4
        DO 32 J=1,I
          TEMPWD = HABRSN(I,J)
          IF(I.NE.J) TEMPWD = TEMPWD+ HABRSN(J,I)
          IF(TEMPWD.LT.EPS) GOTO 32
          CALL HDCYAD(84,HAWDTH,TEMPWD,NETDEC(1,I,J))
 32     CONTINUE
      ENDDO   
C--Decays to sfermions of the pseudoscalar Higgs
      DO 33 I=1,3
        IF(BHASUSY(I).GE.EPS) 
     &    CALL HDCYAD(84,HAWDTH,BHASUSY(I),SFRDCA(1,I))
 33   CONTINUE
C--Decay to Standard Model particles of the charged Higgs
      DO 40 I=1,9
        IF(HCBRS(I).LT.EPS) GOTO 40
        CALL HDCYAD(86,HCWDTH,HCBRS(I),SMDECC(1,I))
 40   CONTINUE
C--Decay to charginos and neutralino of the charged Higgs
      DO I=1,2
        DO 41 J=1,4
        IF(HCBRSU(I,J).LT.EPS) GOTO 41 
        CALL HDCYAD(86,HCWDTH,HCBRSU(I,J),CNDEC(1,I,J))
 41     CONTINUE
      ENDDO
C--Decay to first 2 genration sfermions and stau
      TOTMOD = 0
      MODE = 0
      DO 42 I=1,4
        TOTMOD = TOTMOD+MODE
        MODE = 1
        IF(I.EQ.2.OR.I.EQ.4) MODE = 2
        IF(BHCSUSY(I).LT.EPS) GOTO 42
        TEMPWD = BHCSUSY(I)/MODE
        DO J=1,MODE
          CALL HDCYAD(86,HCWDTH,TEMPWD,SFRDCC(1,TOTMOD+J))
        ENDDO
 42   CONTINUE
C--Now the decay to sbottom and stop
      DO I = 1,2
        DO 43 J=1,2
          IF(BHCSTB(I,J).LT.EPS) GOTO 43
          CALL HDCYAD(86,HCWDTH,BHCSTB(I,J),BTDEC(1,I,J))
 43     CONTINUE
      ENDDO
C--Now reset the ISAJET HIGGS MASSES with the HDECAY ones
      AMHL = REAL(AML)
      AMHH = REAL(AMH)
      AMHA = REAL(AMA)
      AMHC = REAL(AMCH)
      RETURN
 1000 FORMAT('WARNING ',A15,
     &        ' SOFT TERM IS NOT THE SAME FOR 1st AND 2nd GENERATION')
      END
#endif
#ifndef HDECAY2
#ifndef HDECAY3
CDECK  ID>, HDCYSY
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :  Peter Richardson 
C-----------------------------------------------------------------------
      SUBROUTINE HDCYSY(M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,
     &                     MER2,MQL2,MUR2,MDR2,MEL,MER,MSQ,MUR,MDR) 
C-----------------------------------------------------------------------
C  Dummy version of subroutine to interface HDECAY with ISAWIG for MSSM
C  Higgs Decays
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--variables passed from ISAWIG
      REAL M2,MEL1,MER1,MQL1,MUR1,MDR1,MEL2,MER2,MQL2,MUR2,MDR2,
     &     MEL,MER,MSQ,MUR,MDR
      WRITE (6,10)
   10 FORMAT(/10X,'HDECAY CALLED BUT NOT LINKED')
      STOP
      END
#endif
#endif
