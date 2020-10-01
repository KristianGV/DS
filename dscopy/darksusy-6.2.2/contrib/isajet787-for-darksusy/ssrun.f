#include "PILOT.inc"
      PROGRAM SSRUN
C-----------------------------------------------------------------------
C
C     Main program for ISASUSY
C
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sssm.inc"
#include "sspar.inc"
#include "ssmode.inc"
#include "sstype.inc"
#include "isapw.inc"
#include "xmssm.inc"
C
      REAL XMG,XMU,XMHA,XTANB,XMQ1,XMDR,XMUR,XML1,XMER
     $,XMQ2,XMSR,XMCR,XML2,XMMR,XMQ3,XMBR,XMTR,XML3,XMLR
     $,XAT,XAB,XAL,XM1,XM2,XMT
      REAL QSUSY,ASMB,MBMB,ASMT,MTMT,SUALFS,PI,GG
      REAL AMASS
      REAL RSEE,PLEM,PLEP
      DOUBLE PRECISION SSMQCD
      INTEGER J,K,NOGOOD,IMHL,IMHC,IMSQ
      CHARACTER*60 FNAME,FNWIG,FNLHA
      CHARACTER*5 LBLIN,LBLOUT(5),SSID
      CHARACTER*40 V,VISAJE
      INTEGER NOUT,IALLOW,IITEST
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      CHARACTER*30 ISAPW2
      SAVE ISAPW2
      LOGICAL GOWIG,GOLHA
      INTEGER LWIG,ILHA
      REAL MU2,MU1,AMGLMZ,SSPOLE
      INTEGER ISATLS,IMODEL,IMODIN
C
C          Isatools common blocks and variables
C
      COMMON/SUGRED/OMGH2,SIGMA,XFREEZ,NSTEPS,FFF_V
      REAL OMGH2,SIGMA,XFREEZ,FFF_V
      INTEGER NSTEPS
      REAL ALEMIGM2,BFBSG,ALEMI
      COMMON/SUGRES/SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
      REAL*8 SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
      SAVE/SUGRES/
C-FP  INTEGER INUHM
      REAL*8 DAMU,DBFBSG
      REAL BRBS,BRBD
      INTEGER IRED,IRES,IAMU,IBSG,IBLL
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
C          ISAPW2 is used to check whether ALDATA is loaded
      DATA ISAPW2/'ALDATA REQUIRED BY FORTRAN G,H'/
      DATA ILHA/11/
C
C          Initialize
C
      IF(ISAPW1.NE.ISAPW2) THEN
        PRINT*, ' ERROR: BLOCK DATA ALDATA HAS NOT BEEN LOADED.'
        PRINT*, ' ISAJET CANNOT RUN WITHOUT IT.'
        PRINT*, ' PLEASE READ THE FINE MANUAL FOR ISAJET.'
        STOP99
      ENDIF
C
      LOUT=1
      NOGOOD=0
      PRINT 1000
1000  FORMAT(' ENTER output file name (in single quotes)')
      READ*, FNAME
      OPEN(1,FILE=FNAME,STATUS='NEW',FORM='FORMATTED')
C          Print version
      V=VISAJE()
      WRITE(LOUT,1001) V
1001  FORMAT(' ',44('*')/' *',42X,'*'/
     C  ' * ',A40,' *'/
     C  ' *',42X,'*'/' ',44('*')/)
C
C          Open LHA file
      FNLHA=''
      GOLHA=.FALSE.
      PRINT*,'ENTER SUSY Les Houches Accord filename [/ for none]:'
      READ*,FNLHA
      IF(FNLHA.NE.'') THEN
        GOLHA=.TRUE.
        IMODEL=0
        IMODIN=11
        OPEN(ILHA,FILE=FNLHA,STATUS='NEW',FORM='FORMATTED')
      ENDIF

C          Open Isawig file
      FNWIG=''
      GOWIG=.FALSE.
      LWIG=2
      PRINT 1010
1010  FORMAT(' ENTER Isawig (Herwig interface) filename [/ for none]:')
      READ*,FNWIG
      IF(FNWIG.NE.'') THEN
        GOWIG=.TRUE.
        OPEN(LWIG,FILE=FNWIG,STATUS='NEW',FORM='FORMATTED')
      ENDIF
      GORGE=.FALSE.
C
      PRINT 1100
1100  FORMAT(' ENTER M(TP)')
      READ*, XMT
      PRINT 1200
1200  FORMAT(' ENTER M(GLSS), MU, M(A), TAN(BETA)')
      READ*, XMG,XMU,XMHA,XTANB
      PRINT 1250
1250  FORMAT(' ENTER M(Q1), M(D1), M(U1), M(L1), M(E1)')
      READ*, XMQ1,XMDR,XMUR,XML1,XMER
      PRINT 1300
1300  FORMAT(' ENTER M(Q3), M(D3), M(U3), M(L3), M(E3), A_T, A_B, A_L')
      READ*, XMQ3,XMBR,XMTR,XML3,XMLR,XAT,XAB,XAL
      XMQ2=1.E20
      XMSR=1.E20
      XMCR=1.E20
      XML2=1.E20
      XMMR=1.E20
      PRINT 1400
1400  FORMAT(' ENTER OPTIONAL 2ND GEN MASSES (/ FOR DEFAULT):')
      PRINT 1401 
1401  FORMAT(' ENTER M(Q2), M(D2), M(U2), M(L2), M(E2)')
      READ*, XMQ2,XMSR,XMCR,XML2,XMMR
      XM1=1.E20
      XM2=1.E20
      PRINT 1500
1500  FORMAT(' ENTER OPTIONAL GAUGINO MASSES M1, M2 (/ FOR DEFAULT):')
      READ*, XM1,XM2
      AMGVSS=1.E20
      PRINT 1501
1501  FORMAT(' ENTER OPTIONAL GRAVITINO MASS (/ FOR DEFAULT):')
      READ*, AMGVSS
C
C          Calculate...
C
C          First calculate fermion masses at QSUSY
      PI=4.*ATAN(1.)
      QSUSY=SQRT(XMQ3*XMTR)
      ALQCD4=0.177
      AMBT=AMASS(5)
      AMTP=XMT
      ASMB=SUALFS(AMBT**2,.36,AMTP,3)
      MBMB=AMBT*(1.-4*ASMB/3./PI)
      MBQ=SNGL(SSMQCD(DBLE(MBMB),DBLE(QSUSY)))
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SNGL(SSMQCD(DBLE(MTMT),DBLE(QSUSY)))
      MLQ=1.7463
C     For MSSM solution TANBQ=TANB; for SUGRA, TANBQ=/ TANB
      ALFAEM=1./128.
      SN2THW=.232
      AMW=80.423
      GG=SQRT(4*PI*ALFAEM/SN2THW)
      VUQ=SQRT(2*AMW**2/GG**2/(1.+1./XTANB**2))
      VDQ=VUQ/XTANB
C
      CALL SSMSSM(XMG,XMU,XMHA,XTANB,XMQ1,XMDR,XMUR,
     $XML1,XMER,XMQ2,XMSR,XMCR,XML2,XMMR,XMQ3,XMBR,XMTR,XML3,
     $XMLR,XAT,XAB,XAL,XM1,XM2,XMT,IALLOW,0,IMHL,IMHC,IMSQ)
C
C      IF (MHPNEG.EQ.1) THEN
C        PRINT*, 'BAD POINT: M(H_P)^2<0!!'
C        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!!'
C        NOGOOD=3
C      END IF
      IF(IMHL.EQ.1) THEN
        NOGOOD=8
        PRINT*, 'BAD POINT: MHL^2<0!'
        WRITE(LOUT,*) 'BAD POINT: MHL^2<0!'
      ENDIF
      IF(IMHC.EQ.1) THEN
        NOGOOD=8
        PRINT*, 'BAD POINT: M(H_C)^2<0!'
        WRITE(LOUT,*) 'BAD POINT: M(H_C)^2<0!'
      ENDIF
      IF(NOGOOD.NE.0) STOP99
C
C          Test parameters
C
      IF(IALLOW.NE.0) THEN
        WRITE(LOUT,2001)
2001    FORMAT(//' MSSM WARNING: Z1SS IS NOT LSP')
      ENDIF
C
      CALL SSTEST(IALLOW)
      IITEST=IALLOW/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2002)
2002    FORMAT(' MSSM WARNING: Z -> Z1SS Z1SS EXCEEDS BOUND')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2004)
2004    FORMAT(' MSSM WARNING: Z -> CHARGINOS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2008)
2008    FORMAT(' MSSM WARNING: Z -> Z1SS Z2SS TO BIG')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2016)
2016    FORMAT(' MSSM WARNING: Z -> SQUARKS, SLEPTONS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2032)
2032    FORMAT(' MSSM WARNING: Z -> Z* HL0 EXCEEDS BOUND')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
        WRITE(LOUT,2064)
2064    FORMAT(' MSSM WARNING: Z -> HL0 HA0 ALLOWED')
      ENDIF
C
C          Print results. Note decay modes need not be together, so
C          need to select each parent particle separately.
C
C
C          Print ISAJET MSSMi equivalent input
C
      WRITE(LOUT,3000)
3000  FORMAT(/' ISAJET equivalent input:')
      WRITE(LOUT,3001) XMG,XMU,XMHA,XTANB
3001  FORMAT(' MSSMA: ',4F8.2)
      WRITE(LOUT,3002) XMQ1,XMDR,XMUR,XML1,XMER
3002  FORMAT(' MSSMB: ',5F8.2)
      WRITE(LOUT,3003) XMQ3,XMBR,XMTR,XML3,XMLR,XAT,XAB,XAL
3003  FORMAT(' MSSMC: ',8F8.2)
      IF (XML2.LT.1.E19) THEN
      WRITE(LOUT,3004) XMQ2,XMSR,XMCR,XML2,XMMR
3004  FORMAT(' MSSMD: ',5F8.2)
      END IF
      IF (ABS(XM1).LT.1.E19) THEN
      WRITE(LOUT,3005) XM1,XM2
3005  FORMAT(' MSSME: ',2F8.2)
      END IF

      WRITE(LOUT,3006) XMT,ALFAEM,SN2THW,ALFA3
3006  FORMAT(/
     $' M(TP)     =',F10.3/
     $' ALPHAEM   =',F10.5,'   SIN2(THW) =',F10.5,'   ALPHA3 =',F10.5)
C
C
C          Print mass spectrum from SSPAR
C
      WRITE(LOUT,3007) AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,
     $AMB1SS,AMB2SS,AMT1SS,AMT2SS,
     $AMN1SS,AMELSS,AMERSS,
     $AMN3SS,AML1SS,AML2SS,
     $AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,
     $AMW1SS,AMW2SS,
     $AMHL,AMHH,AMHA,AMHC
3007  FORMAT(/' ISAJET masses (with signs):'/
     $' M(GL)  =',F9.2/
     $' M(UL)  =',F9.2,'   M(UR)  =',F9.2,'   M(DL)  =',F9.2,
     $'   M(DR) =',F9.2/
     $' M(B1)  =',F9.2,'   M(B2)  =',F9.2,'   M(T1)  =',F9.2,
     $'   M(T2) =',F9.2/
     $' M(SN)  =',F9.2,'   M(EL)  =',F9.2,'   M(ER)  =',F9.2/
     $' M(NTAU)=',F9.2,'   M(TAU1)=',F9.2,'   M(TAU2)=',F9.2/
     $' M(Z1)  =',F9.2,'   M(Z2)  =',F9.2,'   M(Z3)  =',F9.2,
     $'   M(Z4) =',F9.2/
     $' M(W1)  =',F9.2,'   M(W2)  =',F9.2/
     $' M(HL)  =',F9.2,'   M(HH)  =',F9.2,'   M(HA)  =',F9.2,
     $'   M(H+) =',F9.2)
      WRITE(LOUT,3008) THETAT,THETAB,THETAL,ALFAH
3008  FORMAT(/,' theta_t=',F9.4,'   theta_b=',F9.4,
     $'   theta_l=',F9.4,'   alpha_h=',F9.4)
C
C     Write out chargino /neutralino masses/eigenvectors
C
      WRITE(LOUT,3009) AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS
3009  FORMAT(/' NEUTRALINO MASSES (SIGNED) =',4F10.3)
      DO 100 J=1,4
        WRITE(LOUT,3010) J,(ZMIXSS(K,J),K=1,4)
3010    FORMAT(' EIGENVECTOR ',I1,'       =',4F10.5)
100   CONTINUE
      WRITE(LOUT,3011) AMW1SS,AMW2SS
3011  FORMAT(/' CHARGINO MASSES (SIGNED)  =',2F10.3)
      WRITE(LOUT,3012) GAMMAL,GAMMAR
3012  FORMAT(' GAMMAL, GAMMAR             =',2F10.5/)

C
      WRITE(LOUT,3600)
3600  FORMAT(' PARENT --> DAUGHTERS',18X,'WIDTH',10X,
     $'BRANCHING RATIO'/)
C          Write all modes
      DO 200 J=1,NOUT
        CALL SSPRT(IDOUT(J))
200   CONTINUE
C
C          Write optional LHA file
C
      IF(GOLHA) THEN
        XGLSS=XMG
        XMUSS=XMU
        XHASS=XMHA
        XTBSS=XTANB
        XQ1SS=XMQ1
        XDRSS=XMDR
        XURSS=XMUR
        XL1SS=XML1
        XERSS=XMER
        XQ3SS=XMQ3
        XBRSS=XMBR
        XTRSS=XMTR
        XL3SS=XML3
        XTARSS=XMLR
        XATSS=XAT
        XABSS=XAB
        XATASS=XAL
        XQ2SS=XMQ2
        XSRSS=XMSR
        XCRSS=XMCR
        XL2SS=XML2
        XMRSS=XMMR
        IF(ABS(XM2).GT.1.E19.AND.ABS(XM1).GT.1.E19) THEN
          AMGLMZ=SSPOLE(SIGN(1.,XMG)*AMGLSS,AMGLSS**2,-ALFA3)
           XM2=-ALFA2*AMGLMZ/ALFA3
           XM1=5*SN2THW/3./(1.-SN2THW)*XM2
        ENDIF
        XM1SS=XM1
        XM2SS=XM2
        CALL ISALHA(ILHA,IMODEL,IMODIN,XMT)
C       write LHD header
        CALL ISALHD(ILHA,IDOUT(J),0,NOUT)        
        DO 210 J=1,NOUT
          CALL ISALHD(ILHA,IDOUT(J),J,NOUT)
210     CONTINUE
      ENDIF

C
C          Write optional Isawig file
C          C.f. SUGRUN. Arguments are not actually used??
C
      IF(GOWIG) THEN
        IF(ABS(XM2).LT.1.E19.AND.ABS(XM1).LT.1.E19) THEN
           MU2=-XM2
           MU1=-XM1
        ELSE
          AMGLMZ=SSPOLE(SIGN(1.,XMG)*AMGLSS,AMGLSS**2,-ALFA3)
           MU2=-ALFA2*AMGLMZ/ALFA3
           MU1=5*SN2THW/3./(1.-SN2THW)*MU2
        ENDIF
        IF(XMQ2.LT.1.0E19) XMQ2=XMQ1
        IF(XMSR.LT.1.0E19) XMSR=XMDR
        IF(XMCR.LT.1.0E19) XMCR=XMUR
        IF(XML2.LT.1.0E19) XML2=XML1
        IF(XMMR.LT.1.0E19) XMMR=XMER
        CALL ISAWIG(LWIG,0,XMT,MU2,XML3,XMLR,
     $  XMQ1,XMUR,XMDR,
     $  XML1,XMER,XMQ2,
     $  XMCR,XMSR,XML3,
     $  XMLR,XMQ3,XMTR,
     $  XMBR)
      ENDIF
C
C          Optionally, compute e+e- cross sections
C
#ifdef PRTEESIG_X
      IF(GOLHA) THEN
        RSEE=0.
        PLEM=0.
        PLEP=0.
        PRINT 3700
3700    FORMAT(' Print e+e- xsecs? Enter Ecm,Pol-,Pol+ [/ for none]:')
        PRINT 3701 
3701    FORMAT(' ENTER RSEE,PLEM,PLEP')
        READ*, RSEE,PLEM,PLEP
        IF (RSEE.GT.0.) THEN
          CALL ISASEE(RSEE,PLEM,PLEP,ILHA)
        END IF
      ENDIF
#endif
C
      STOP
      END
