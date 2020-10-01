#include "PILOT.inc"
      PROGRAM SUGRUN
C
C     Main program to calculate MSSM input parameters for ISAJET
C     from renormalization group equations and supergravity.
C     All external names are of the form SUxxxx.
C     Must link with block data ALDATA.
C
C     Includes optional link to ISATOOLS, which requires libisared.a.
C     Make this from isared.tar; see the Makefile for instructions.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sspar.inc"
#include "sstype.inc"
#include "sugmg.inc"
#include "sugxin.inc"
#include "sugpas.inc"
#include "sugnu.inc"
#include "isapw.inc"
      CHARACTER*80 FNAME,FNLHA,FNWIG
      LOGICAL GOLHA,GOWIG
      INTEGER ILHA,IWIG,IMHL,IMHC,IMSQ
      REAL M0,MHF,A0,TANB,SGNMU,MT,XLAMGM,XMESGM,XN5GM,AMPL,XCMGV
      INTEGER NSTEP,IMODEL,INUSUG,IMODIN,ISATLS,IRGEFL
      INTEGER K,NOUT,IALLOW,IITEST,J,II
      CHARACTER*40 VERSN,VISAJE
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      CHARACTER*30 ISAPW2
      SAVE ISAPW2
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
      COMMON/RGEFNM/FNRGE
      CHARACTER*128 FNRGE
      REAL*8 DAMU,DBFBSG
      REAL BRBS,BRBD
      REAL R,BFBTN
      INTEGER IRED,IRES,IAMU,IBSG,IBLL,IBTN
      REAL RSEE,PLEM,PLEP
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
      DATA AMPL/2.4E18/
C          ISAPW2 is used to check whether ALDATA is loaded
      DATA ISAPW2/'ALDATA REQUIRED BY FORTRAN G,H'/
      DATA ILHA/11/,IWIG/12/
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
      NSTEP=1000
      XNRIN(2)=1.E20
      XSUGIN(7)=0.
      INUHM=0
      IDTERM=0
C
C          Open files
C
      PRINT*,'ENTER output filename in single quotes:'
      READ*,FNAME
      OPEN(1,FILE=FNAME,STATUS='NEW',FORM='FORMATTED')
      FNLHA=''
      GOLHA=.FALSE.
      PRINT*,'ENTER SUSY Les Houches Accord filename [/ for none]:'
      READ*,FNLHA
      IF(FNLHA.NE.'') THEN
        GOLHA=.TRUE.
        OPEN(ILHA,FILE=FNLHA,STATUS='NEW',FORM='FORMATTED')
      ENDIF
      FNWIG=''
      GOWIG=.FALSE.
      PRINT*,'ENTER Isawig (Herwig interface) filename [/ for none]:'
      READ*,FNWIG
      IF(FNWIG.NE.'') THEN
        GOWIG=.TRUE.
        OPEN(IWIG,FILE=FNWIG,STATUS='NEW',FORM='FORMATTED')
      ENDIF
C
      PRINT*,'ENTER 1 for mSUGRA:'
      PRINT*,'ENTER 2 for mGMSB:'
      PRINT*,'ENTER 3 for non-universal SUGRA:'
      PRINT*,'ENTER 4 for SUGRA with truly unified gauge couplings:'
      PRINT*,'ENTER 5 for non-minimal GMSB:'
      PRINT*,'ENTER 6 for SUGRA+right-handed neutrino:'
      PRINT*,'ENTER 7 for minimal anomaly-mediated SUSY breaking:'
      PRINT*,'ENTER 8 for non-minimal AMSB:'
      PRINT*,'ENTER 9 for mixed moduli-AMSB:'
      PRINT*,'ENTER 10 for Hypercharged-AMSB:'
      PRINT*,'ENTER 11 for NUHM from D-term:'
      PRINT*,'ENTER 12 for general mirage mediation (GMM):'
      READ*,IMODIN
      IMODEL=IMODIN
      IF (IMODEL.EQ.4) THEN
        IAL3UN=1
        IMODEL=1
      END IF
      IF (IMODEL.EQ.1.OR.IMODEL.EQ.3.OR.IMODEL.EQ.6) THEN
        PRINT*,'ENTER M_0, M_(1/2), A_0, tan(beta), sgn(mu), M_t:'
        READ*,M0,MHF,A0,TANB,SGNMU,MT
        IF (IMODEL.EQ.6) THEN
          IMODEL=1
          PRINT*,' ENTER M(nu_3)[=0], M_Majorana, A_N, M(NRSS)'
          PRINT*,'If M(nu_3)=0, do top-down with fn-ft unification'
          READ*,XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4)
          GO TO 15
        END IF
        IF (IMODEL.EQ.3) THEN
          IMODEL=1
10        PRINT*,' ENTER 1,...,5 for NUSUGx keyword; 0 to continue:'
          PRINT*,' NUSUG1 = GUT scale gaugino masses'
          PRINT*,' NUSUG2 = GUT scale A terms'
          PRINT*,' NUSUG3 = GUT scale Higgs masses'
          PRINT*,' NUSUG4 = GUT scale 1st/2nd generation masses'
          PRINT*,' NUSUG5 = GUT scale 3rd generation masses'
          PRINT*,' ENTER 6 to activate right-hand neutrino'
          PRINT*,' ENTER 7 to enter alternate high scale Q_max.ne.M_GUT'
          PRINT*,' ENTER 8 for mu, mA input at weak scale (NUHM model)'
          READ*,INUSUG
          IF (INUSUG.EQ.0) THEN
            GO TO 15
          ELSE IF (INUSUG.EQ.1) THEN
            PRINT*,'Enter GUT scale M_1, M_2, M_3:'
            READ*,XNUSUG(1),XNUSUG(2),XNUSUG(3)
C            IF (XNUSUG(3).LE.0.) THEN
C            PRINT*, ' NEGATIVE M_3 IS NOT ALLOWED'
C            STOP 99
C            END IF
          ELSE IF (INUSUG.EQ.2) THEN
            PRINT*,'Enter GUT scale A_t, A_b, A_tau:'
            READ*,XNUSUG(6),XNUSUG(5),XNUSUG(4)
          ELSE IF (INUSUG.EQ.3) THEN
            PRINT*,'Enter GUT scale m_Hd, m_Hu:'
            READ*,XNUSUG(7),XNUSUG(8)
          ELSE IF (INUSUG.EQ.4) THEN
            PRINT*,'Enter GUT scale M(ul), M(dr), M(ur), M(el), M(er):'
            READ*,XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10),XNUSUG(9)
          ELSE IF (INUSUG.EQ.5) THEN
            PRINT*,'Enter GUT scale M(tl), M(br), M(tr), M(Ll), M(Lr):'
            READ*,XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15),XNUSUG(14)
          ELSE IF (INUSUG.EQ.6) THEN
            PRINT*,' ENTER M(nu_3), M_Majorana, A_N, M(NRSS)'
            READ*,XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4)
          ELSE IF (INUSUG.EQ.7) THEN
            PRINT*,' ENTER Q_max high scale for SUSY BCs'
            READ*,XSUGIN(7)
          ELSE IF (INUSUG.EQ.8) THEN
            PRINT*,' ENTER mu(Q), mA(Q)'
            READ*,XNUSUG(19),XNUSUG(20)
            MU=XNUSUG(19)
            AMHA=XNUSUG(20)
            TWOM1=-MU
            INUHM=1
          END IF
          GO TO 10
        END IF
      ELSE IF (IMODEL.EQ.2.OR.IMODEL.EQ.5) THEN
          PRINT*,'ENTER Lambda, M_mes, N_5, tan(beta), sgn(mu), ',
     $    'M_t, C_gv:'
          READ*,M0,MHF,A0,TANB,SGNMU,MT,XCMGV
          XGMIN(7)=XCMGV
          XGMIN(8)=1.
          AMGVSS=M0*MHF*XCMGV/SQRT(3.)/AMPL
          IF (IMODEL.EQ.5) THEN
            IMODEL=2
            PRINT*,'Rsl = factor multiplying gaugino masses at M_mes'
            PRINT*,'dmH_d^2, dmH_u^2 = Higgs mass**2 shifts at M_mes'
            PRINT*,'d_Y = mass**2 shifts proportional to Y at M_mes'
            PRINT*,'n5_1,n5_2,n5_3 = n5 values for U(1),SU(2),SU(3)'
            PRINT*,'ENTER Rsl, dmH_d^2, dmH_u^2, d_Y, n5_1, n5_2, n5_3'
            READ*,XGMIN(8),XGMIN(9),XGMIN(10),XGMIN(11),XGMIN(12),
     $      XGMIN(13),XGMIN(14)
            END IF
      ELSE IF (IMODEL.EQ.7) THEN
        PRINT*,'ENTER M_0, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ*,M0,MHF,TANB,SGNMU,MT
        A0=0.
        DO 101 II=1,7
101     XAMIN(II)=1.
      ELSE IF (IMODEL.EQ.8) THEN
        PRINT*,'ENTER M_0, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ*,M0,MHF,TANB,SGNMU,MT
        A0=0.
        PRINT*,'ENTER cQ, cD, cU, cL, cE, cHd, cHu:'
        READ*,
     $XAMIN(1),XAMIN(2),XAMIN(3),XAMIN(4),XAMIN(5),XAMIN(6),XAMIN(7)
        IMODEL=7
      ELSE IF (IMODEL.EQ.9) THEN
        PRINT*,'ENTER alpha, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ*,M0,MHF,TANB,SGNMU,MT
        A0=0.
C          Set defaults
        DO 102 II=1,7
          XAMIN(II)=0
102     CONTINUE
        DO 103 II=8,10
          XAMIN(II)=1
103     CONTINUE
        PRINT*,'ENTER moduli weights nQ, nD, nU, nL, nE, nHd, nHu ',
     $  '[/ for all 0]:'
        READ*,
     $XAMIN(1),XAMIN(2),XAMIN(3),XAMIN(4),XAMIN(5),XAMIN(6),XAMIN(7)
        PRINT*,'ENTER moduli parameters L1, L2, L3 [/ for all 1]:'
        READ*,XAMIN(8),XAMIN(9),XAMIN(10)
      ELSE IF (IMODEL.EQ.10) THEN
        PRINT*,'ENTER alpha, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ*,XAMIN(11),MHF,TANB,SGNMU,MT
        M0=0.
        A0=0.
      ELSE IF (IMODEL.EQ.11) THEN
        PRINT*,'ENTER M_16, M_(1/2), A_0, tan(beta), mu, mA, M_t:'
        READ*,M0,MHF,A0,TANB,XNUSUG(19),XNUSUG(20),MT
        SGNMU=SIGN(1.,XNUSUG(19))
        MU=XNUSUG(19)
        AMHA=XNUSUG(20)
        TWOM1=-MU
        INUHM=1
        IDTERM=1
        IMODEL=1
      ELSE IF (IMODEL.EQ.12) THEN
        PRINT*,'ENTER al,m32,c_m,c_m3,a_3,tan(beta),sgn(mu),m_t,INUHM:'
        PRINT*,'INUHM=0 for cHu,cHd inputs or 1 for mu, mA inputs:'
        READ*,XAMIN(11),MHF,XAMIN(1),XAMIN(5),XAMIN(4),
     $TANB,SGNMU,MT,INUHM
        IF (INUHM.EQ.0) THEN
          PRINT*,'ENTER c_Hu, c_Hd:'
          READ*,XAMIN(2),XAMIN(3)
        ELSE IF (INUHM.EQ.1) THEN
          PRINT*,'ENTER mu, mA:'
          READ*,XNUSUG(19),XNUSUG(20)
          MU=XNUSUG(19)
          AMHA=XNUSUG(20)
          TWOM1=-MU
        END IF
        IMODEL=12
      ELSE
        PRINT*,'Invalid model choice.'
        STOP99
      END IF
C
C          Solve RG equations
C
15    CALL SUGRA(M0,MHF,A0,TANB,SGNMU,MT,IMODEL)
C
C          Print results
C
      
      VERSN=VISAJE()
      WRITE(LOUT,20) VERSN
20    FORMAT(' ',44('*')/' *',42X,'*'/
     $  ' * ',A40,' *'/
     $  ' *',42X,'*'/' ',44('*')/)
      IF (NOGOOD.EQ.1) THEN
        PRINT*, 'BAD POINT: TACHYONIC PARTICLES!'
        WRITE(LOUT,*) 'BAD POINT: TACHYONIC PARTICLES!'
      ELSE IF (NOGOOD.EQ.2) THEN
        PRINT*, 'BAD POINT: NO EW SYMMETRY BREAKING!'
        WRITE(LOUT,*) 'BAD POINT: NO EW SYMMETRY BREAKING!'
      ELSE IF (NOGOOD.EQ.3) THEN
        PRINT*, 'BAD POINT: M(H_P)^2<0!'
        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!'
      ELSE IF (NOGOOD.EQ.4) THEN
        PRINT*, 'BAD POINT: YUKAWA>10!'
        WRITE(LOUT,*) 'BAD POINT: YUKAWA>10!'
      ELSE IF (NOGOOD.EQ.5.AND.IMODEL.EQ.1) THEN
        PRINT*, 'SUGRA BAD POINT: Z1SS NOT LSP!'
        WRITE(LOUT,*) 'SUGRA BAD POINT: Z1SS NOT LSP!'
      ELSE IF (NOGOOD.EQ.7) THEN
        PRINT*, 'BAD POINT: XT EWSB BAD!'
        WRITE(LOUT,*) 'BAD POINT: XT EWSB BAD!'
      ELSE IF (NOGOOD.EQ.8) THEN
        PRINT*, 'BAD POINT: MHL^2<0!'
        WRITE(LOUT,*) 'BAD POINT: MHL^2<0!'
      ELSE IF (NOGOOD.EQ.-1) THEN
        PRINT*, 'BAD POINT: NO RGE SOLUTION FOUND'
        WRITE(LOUT,*) 'BAD POINT: NO RGE SOLUTION FOUND'
      END IF
      IF (MHPNEG.EQ.1) THEN
        PRINT*, 'BAD POINT: M(H_P)^2<0!!'
        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!!'
        NOGOOD=3
      END IF
      IF (MHLNEG.EQ.1) THEN
        PRINT*, 'BAD POINT: M(H_L)^2<0!!'
        WRITE(LOUT,*) 'BAD POINT: M(H_L)^2<0!!'
        NOGOOD=8
      END IF
      IF (MHCNEG.EQ.1) THEN
        PRINT*, 'BAD POINT: M(H^+)^2<0!!'
        WRITE(LOUT,*) 'BAD POINT: M(H^+)^2<0!!'
        NOGOOD=8
      END IF
      IF (MSQNEG.EQ.1) THEN
        PRINT*, 'BAD POINT: M(3RD)^2<0!!'
        WRITE(LOUT,*) 'BAD POINT: M(3RD)^2<0!!'
        NOGOOD=9
      END IF
      IF(NOGOOD.NE.0) STOP99
      IF(ITACHY.NE.0) THEN
        WRITE(LOUT,*) 'WARNING: TACHYONIC SLEPTONS AT GUT SCALE'
        WRITE(LOUT,*) '         POINT MAY BE INVALID'
      ENDIF
CCC      IF (IGUTST.EQ.1) THEN
CCC        PRINT*, 'WARNING: GUT INSTABILITY IN NUHM MODEL'
CCC      END IF
C
C          Print selected model and results
C
      IF(IMODIN.EQ.1) WRITE(LOUT,1001)
1001  FORMAT(//' Minimal supergravity (mSUGRA) model:'/)
      IF(IMODIN.EQ.2) WRITE(LOUT,1002)
1002  FORMAT(//' Minimal gauge mediated (GMSB) model:'/)
      IF(IMODIN.EQ.3) WRITE(LOUT,1003)
1003  FORMAT(//' Non-universal supergravity model:'/)
      IF(IMODIN.EQ.4) WRITE(LOUT,1004)
1004  FORMAT(//' Supergravity model with truly unified couplings:'/)
      IF(IMODIN.EQ.5) WRITE(LOUT,1005)
1005  FORMAT(//' Non-minimal gauge mediated (GMSB) model:'/)
      IF(IMODIN.EQ.6) WRITE(LOUT,1006)
1006  FORMAT(//' Supergravity model with right-handed neutrinos:'/)
      IF(IMODIN.EQ.7) WRITE(LOUT,1007)
1007  FORMAT(//' Anomaly-mediated SUSY breaking model:'/)
      IF(IMODIN.EQ.8) WRITE(LOUT,1008)
1008  FORMAT(//' Non-minimal anomaly-mediated SUSY breaking model:'/)
      IF(IMODIN.EQ.9) WRITE(LOUT,1009)
1009  FORMAT(//' Mixed modulus-AMSB SUSY breaking model:'/)
      IF(IMODIN.EQ.10) WRITE(LOUT,1010)
1010  FORMAT(//' Hypercharged-AMSB SUSY breaking model:'/)
      IF(IMODIN.EQ.11) WRITE(LOUT,1011)
1011  FORMAT(//' NUHM model with D-term splitting (NUHMDT):'/)
      IF(IMODIN.EQ.12) WRITE(LOUT,1012)
1012  FORMAT(//' general mirage mediation model (GMM or GMMprime):'/)
C
C          Calculate all masses and decay modes
C
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ MT,IALLOW,IMODEL,IMHL,IMHC,IMSQ)
C
C          Execute Isatools
C
#ifdef ISATOOLS_X
      ISATLS=0
      IRED=0
      IRES=0
      IAMU=0
      IBSG=0
      IBLL=0
      IBTN=0
      PRINT*,'Run Isatools? Choose 2=all, 1=some, 0=none:'
      READ*,ISATLS
      IF(ISATLS.EQ.2) THEN
        IRED=1
        IRES=1
        IAMU=1
        IBSG=1
        IBLL=1
        IBTN=1
      ELSE IF (ISATLS.EQ.1) THEN
        PRINT*,'Select desired ISATools packages'
        PRINT*,'Neutralino Relic Density [1/0]:'
        READ*,IRED
        PRINT*,'Neutralino DD rates [1/0]:'
        READ*,IRES
        PRINT*,'Muon (g-2)/2 [1/0]:'
        READ*,IAMU
        PRINT*,'b->s gamma branching fraction [1/0]:'
        READ*,IBSG
        PRINT*,'B_s->ll branching fractions [1/0]:'
        READ*,IBLL
        PRINT*,'B_u->tau+nu_tau branching fractions [1/0]:'
        READ*,IBTN
      ENDIF  
      DAMU=0.
      BFBSG=0.
      OMGH2=0.
      BRBS=0.
      BRBD=0.
      BFBTN=0.
      SIGMA0PROT=0.
      SIGMA0NEUT=0.
      SIGMASPROT=0.
      SIGMASNEUT=0.
C          g_mu - 2
      IF (IAMU.EQ.1) THEN
        ALEMI=128.
      CALL ISAAMU(RV2V1,ALEMI,GAMMAL,GAMMAR,TWOM1,AAL,SQRT(GSS(16)),
     $SQRT(GSS(15)),AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,AMW1SS,AMW2SS,AMN2SS,
     $              ZMIXSS,0,DAMU)
      ENDIF
C          B -> s gamma
      IF (IBSG.EQ.1) THEN
        CALL ISABSG(IMODEL,M0,MHF,A0,DBFBSG,0)
        BFBSG=SNGL(DBFBSG)
      ENDIF
C          Bu -> tau+nu_tau
      IF (IBTN.EQ.1) THEN
        CALL ISABTN(TANB,MSS(32),MU,MSS(1),MSS(10),MSS(11),R,BFBTN)
      ENDIF
C          Dark matter cross sections
      IF (IRES.EQ.1) THEN
        CALL ISARES(0)
      ENDIF
C          Relic density calculation -- requires isared.tar
      IF (IRED.EQ.1) THEN
        CALL ISARED(0)
      ENDIF
C          B_s -> mu mu (BRBS) and B -> tau tau (BRBD)
      IF (IBLL.EQ.1) THEN
        CALL ISABMM(MT,TANB,MSS(29),MSS(30),MSS(31),MSS(1),
     $              MSS(10),MSS(11),GSS(11),THETAB,
     $              MSS(12),MSS(13),GSS(12),THETAT,
     $              MU,GSS(8),ALFAH,MSS(6),MSS(8),BRBS,BRBD)
      ENDIF
#endif
C
      CALL SUGPRT(IMODEL,IMODIN)
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
2008    FORMAT(' MSSM WARNING: Z -> Z1SS Z2SS TOO BIG')
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
C          Isatools output
C
#ifdef ISATOOLS_X
      BFBSG=1.0E-4*BFBSG
      WRITE(LOUT,3500) DAMU,BFBSG,OMGH2,FFF_V,BRBS,BRBD,BFBTN,
     $SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
3500  FORMAT(//' Output from ISATOOLS:'/
     $' Delta a_mu       = ',E12.5,'    BF(b->s gamma)    = ',E12.5/
     $' Omega h^2        = ',E12.5,'    <s.v>(v=0)[cm^3/s]= ',E12.5/
     $' BF(Bs -> mu mu)  = ',E12.5,'    BF(B -> tau tau)  = ',E12.5/
     $' BF(Bu -> tau nu) = ',E12.5/
     $' LSP-nucleon spin independent (SI) and dependent (SD) sigmas:'/
     $' sigma(p,SI)[pb]  = ',E12.5,'    sigma(n,SI)[pb]   = ',E12.5/
     $' sigma(p,SD)[pb]  = ',E12.5,'    sigma(n,SD)[pb]   = ',E12.5)
#endif
      WRITE(LOUT,3600)
3600  FORMAT(//' ISASUSY decay modes:'/
     $' Parent --> daughters',18X,'Width',10X,'Branching ratio'/)
C          Write all modes
      DO 200 J=1,NOUT
        CALL SSPRT(IDOUT(J))
200   CONTINUE
C
C          Make optional output files
C
      IF(GOLHA) THEN
        CALL ISALHA(ILHA,IMODEL,IMODIN,MT)
C       write LHD header
        CALL ISALHD(ILHA,IDOUT(J),0,NOUT)        
        DO 210 J=1,NOUT
          CALL ISALHD(ILHA,IDOUT(J),J,NOUT)
210     CONTINUE
      ENDIF
      IF(GOWIG) THEN
        CALL ISAWIG(IWIG,0,MT,GSS(8),SQRT(GSS(16)),SQRT(GSS(15)),
     &              SQRT(GSS(19)),SQRT(GSS(18)),SQRT(GSS(17)),
     &              SQRT(GSS(16)),SQRT(GSS(15)),SQRT(GSS(19)),
     &              SQRT(GSS(18)),SQRT(GSS(17)),SQRT(GSS(21)),
     &              SQRT(GSS(20)),SQRT(GSS(24)),SQRT(GSS(23)),
     &              SQRT(GSS(22)))
      ENDIF
C
#ifdef ISAFLAVR_X
      FNRGE=''
      IRGEFL=0
      PRINT*,'To run RGEFLAV, enter filename Prefix [/ for none]:'
      PRINT*,'RGEFLAV will open file Prefix.rgein, and print to files'
      PRINT*,'Prefix.weakout, Prefix.gutout, Prefix.sqm2u, Prefix.sqm2d'
      READ*,FNRGE
      IF(FNRGE.NE.'') THEN
        IRGEFL=1
      ENDIF
      IF(IRGEFL.EQ.1) THEN
        CALL RGEFLAV(M0,MHF,A0,TANB,SGNMU,MT)
      END IF
#endif
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
