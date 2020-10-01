#include "PILOT.inc"
#ifdef ISATOOLS_X
C--------------------------------------------------------------------
      subroutine ISABSG(IMODEL,M0,MHF,A0,BFBSG,IWR)
C--------------------------------------------------------------------
C    Calculates branching ratio of b->s\gamma decay.
C
C    Note: This version uses just running MBQ.
C   
C    Ref: H.Anlauf hep-ph/9406286;
C         C.Greub, T.Hurth and D.Wyler  PRD54 (1996);
C         A.Masiero et al  NPB353 (1991) 591-649;
C         M.Ciuchini et al  hep-ph/9304257
C
C    Created by H.Baer and M.Brhlik
C    Modified:
c      03/13/07 by Azar Mustafayev  - converted to double precision
c               and improved interface with ISASUGRA
c      05/29/07 by Azar Mustafayev  - neutrino sector added
c      12/13/07 by Azar Mustafayev  - fixes to speedup
c
c  NB: This version is design to work with ISAJET 7.75
c
C--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IMODEL,IWR
      REAL M0,MHF,A0
      REAL*8 BFBSG
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RSIGT,RSIGB,RSIGL
c  GISA(i) - values of RGE parameters at MZ in DRbar:
C     GISA( 1) = g_1        GISA( 2) = g_2        GISA( 3) = g_3
C     GISA( 4) = y_tau      GISA( 5) = y_b        GISA( 6) = y_t
C     GISA( 7) = M_1        GISA( 8) = M_2        GISA( 9) = M_3
C     GISA(10) = A_tau      GISA(11) = A_b        GISA(12) = A_t
C     GISA(13) = M_hd^2     GISA(14) = M_hu^2     GISA(15) = M_er^2
C     GISA(16) = M_el^2     GISA(17) = M_dnr^2    GISA(18) = M_upr^2
C     GISA(19) = M_upl^2    GISA(20) = M_taur^2   GISA(21) = M_taul^2
C     GISA(22) = M_btr^2    GISA(23) = M_tpr^2    GISA(24) = M_tpl^2
C     GISA(25) = mu         GISA(26) = B          GISA(27) = Y_N
C     GISA(28) = M_nr       GISA(29) = A_n        GISA(30) = vdq
C     GISA(31) = vuq
c
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
c     RSIGT,RSIGB,RSIGL - radiative corrections to top, bottom and tau
c                         Yukawas at MSUSY
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RSIGT,RSIGB,RSIGL
      SAVE /BSG/
c     
      REAL*8 G(157)
c   G(i) - values of RGE parameters at MZ in DRbar:
C     G(1) = g_1
C     G(2) = g_2
C     G(3) = g_3
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
C     G(31) = M_1
C     G(32) = M_2
C     G(33) = M_3
C     G(34)-G(42) = h_u
C     G(43)-G(51) = h_d
C     G(52)-G(60) = h_e
C     G(61)= mu
C     G(62)= B*mu
C     G(63)= m^2_Hu
C     G(64)= m^2_Hd
C     G(65)-G(73) = m^2_Q
C     G(74)-G(82) = m^2_L
C     G(83)-G(91) = m^2_u
C     G(92)-G(100)= m^2_d
C     G(101)-G(109)=m^2_e
C     G(110)= v_u
C     G(111)= v_d
C     G(112)-G(120) = Y_nu  
C     G(121)-G(129) = M_RHN 
C     G(130)-G(138) = h_nu
C     G(139)-G(147) = m^2_nu
C     G(148)-G(156) = \kappa
C     G(157) = \lambda   - SM quartic higgs coupling
      REAL*8 CKM(3,3),YU(3,3),YD(3,3),YE(3,3),HU(3,3),HD(3,3),HE(3,3),
     &       xYU(3,3),xYD(3,3),xHU(3,3),xHD(3,3),GF(157),xYUMT(3,3),
     &       VULMT(3,3),VURMT(3,3),
     &       YUGUT(3,3),YDGUT(3,3),YEGUT(3,3),YNGUT(3,3)
      REAL SQ,SQOLD,SFTMT
      REAL*8 HUGUT(9),HDGUT(9),HEGUT(9)
      REAL*8 CSM(9),CH1(9),CH2(19),
     $     CC111(12),CC112(12),CC113(12),
     $     CC121(12),CC122(12),CC123(12),
     $     CC211(17),CC212(17),CC213(17),
     $     CC221(17),CC222(17),CC223(17),
     $     DSSM(9),DSH1(9),DSH2(9),
     $     DSC111(9),DSC112(9),DSC113(9),
     $     DSC121(9),DSC122(9),DSC123(9),
     $     DSC211(9),DSC212(9),DSC213(9),
     $     DSC221(9),DSC222(9),DSC223(9),
     $     WSM(27),WHP1(27),WHP2(57),
     $     WSC111(36),WSC112(36),WSC113(36),
     $     WSC121(36),WSC122(36),WSC123(36),
     $     WSC211(51),WSC212(51),WSC213(51),
     $     WSC221(51),WSC222(51),WSC223(51),
     $     CI1(8),WCI1(24),CI2(8),WCI2(24),
     $     GK(2,3),HK(2,3)
      REAL*8 GAM(6,6),W2(471)

      REAL*8 AGE(4,6,6),AHE(4,6,6)
      INTEGER IGF(157)
      EXTERNAL RGE157,GAMMASM,GAMMAHP,GAMMAC1,GAMMAC2,GAMMAWB1,GAMMAWB2
      REAL SUALFE,SUALFS
      REAL*8 DDILOG,EI,BI,CI,GES,FES
C
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),IRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER IRHN
      SAVE /BSGDEC/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
c       MCHA - chargino masses
c       MSQU - up-squark masses
c       MSQD - down-squark masses
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX
      SAVE /CHRGN/        
      COMMON /GLNN/ MCH0(4),MDSB(6)
c       MCH0 - neutralino masses
c       MDSB - 
      REAL*8 MCH0,MDSB
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
c        M1,M2,M3  - gaugino masses at MZ
c        ABOT,ATOP,ATAU  - soft trilinear scalar couplings at MZ 
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON /G3/ G3
c        G3 = g_3 - strong coupling at current stage of RGE evolution
      REAL*8 G3
      SAVE /G3/
      COMMON /SCALEQ/ Q
      REAL*8 Q
      SAVE /SCALEQ/
C
      REAL*8 PI,SR2,C12,C13,C23,
     &       XLAMGM,XMESGM,XN5GM,XLM,THRF,THRG,BETA,SINB,COSB,VEV,
     &       MELMZ,MMUMZ,MUPMZ,MDOMZ,MCHMZ,MSTMZ,FELMZ,FMUMZ,FUPMZ,
     &       FDOMZ,FCHMZ,FSTMZ,TZ,TGUT,T,DT,QOLD,XC,DY,BLHAT,BBHAT,
     &       BTHAT,QNEW,COS2B
      REAL*8 VULQ(3,3),VURQ(3,3),VDLQ(3,3),VDRQ(3,3)
      REAL*8 FTAGUT,FBGUT,FTGUT,FELGUT,FMUGUT,FDOGUT,FSTGUT,FUPGUT,
     &       FCHGUT
      REAL*8 X,Y,FN1,FN2,FN3,FN4,YPS,SCALE,XWA,XCB,ALES,ALEL,G3MHP,MSQ,
     &       MCH,G3SS,QMAX,SWI1,SWI2,SWI3,CER,CEG,ALFA2,RALPH,FACG7,
     &       FACG8,FACN7,FACN8,SINTW2I,XTW,ASMB,AEMB,G2,GEF,BQLOG,ZI,
     &       LOGZI,RER2,RER7,RER8,REDE,IMR2,IMR7,IMR8,IMDE,EF,GVIRT,
     &       GAMMASG,gammasg0,GAMMASL,GBREMS,XWB
      REAL*8 c3smnai,c2smnai,c3hpnai,c3c11nai,c3c12nai,c3c13nai,
     &       c3c23nai,C2SMEFF,C3SMEFF,C2HPEFF,C3HPEFF,C2C11EFF,C3C11EFF,
     &       C2C12EFF,C3C12EFF,C2C13EFF,C3C13EFF,C2C21EFF,C3C21EFF,
     &       C2C22EFF,C3C22EFF,C2C23EFF,C3C23EFF,C3GS1EFF,C2GS1EFF,
     &       C3GS2EFF,C2GS2EFF,C3GS3EFF,C2GS3EFF,C3GS4EFF,C2GS4EFF,
     &       C3GS5EFF,C2GS5EFF,C3GS6EFF,C2GS6EFF,C3NS1EFF,C2NS1EFF,
     &       C3NS2EFF,C2NS2EFF,C3NS3EFF,C2NS3EFF,C3NS4EFF,C2NS4EFF,
     &       C3NS5EFF,C2NS5EFF,C3NS6EFF,C2NS6EFF,C3CHEFF,C3GLEFF,
     &       C3NEEFF,C2CHEFF,C2GLEFF,C2NEEFF,C3C1EFF,C3C2EFF,
     &       C3C1NAI,C3C2NAI,C3CHNAI,C3GLNAI,C3NENAI,c3c21nai,c3c22nai
      REAL*8 BRSL,COBA,BRSG,BRSG0,brup,brdo,BFBSG0
      REAL XQ2,XMT
      INTEGER NSTEP,I,J,K,L,II
      DATA NSTEP/30000/
      
      real*8 ESQ(3,3),ESD(3,3),EAD(3,3),EMD(3,3)
      
      PI=4.d0*DATAN(1.d0)
      SR2=SQRT(2.d0)
      
c...Set experimentaly measured parameters
      BRSL=0.104    ! semileptonic branching ratio
      COBA=0.95     ! CKM factor

c...Reset the branching fraction             
      BFBSG=-1.d0

c...Print logo
      IF (IWR.EQ.1) THEN
        PRINT*,'%%%%%%%%%%%%%%%%%%% ISABSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      ENDIF

c...Transfer input parameters
      CALL SUG2BSG

C...DEFINE CKM MATRIX
      C12=SQRT(1.d0-S12**2)
      C13=SQRT(1.d0-S13**2)
      C23=SQRT(1.d0-S23**2)
      CKM(1,1)= C12*C13
      CKM(1,2)= S12*C13
      CKM(1,3)= S13
      CKM(2,1)=-S12*C23-C12*S23*S13
      CKM(2,2)= C12*C23-S12*S23*S13
      CKM(2,3)= S23*C13
      CKM(3,1)= S12*S23-C12*C23*S13
      CKM(3,2)=-C12*S23-S12*C23*S13
      CKM(3,3)= C23*C13
      
C---------------------------------------------
      CALL CHARGINO(GK,HK)

C...Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLAMGM=M0
        XMESGM=MHF
        XN5GM=A0
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
      END IF

C...Assemble (s)quark matrices and transform into correct basis if necessary
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      COS2B=COS(2.d0*BETA)
      VEV=SQRT(V**2+VP**2)
      DO I=1,3
        DO J=1,3
	  xYU(I,J)=0.d0
	  xYD(I,J)=0.d0
	  xHU(I,J)=0.d0
	  xHD(I,J)=0.d0
          YU(I,J)=0.d0
          YD(I,J)=0.d0
          YE(I,J)=0.d0
          HU(I,J)=0.d0
          HD(I,J)=0.d0
          HE(I,J)=0.d0
        ENDDO
      ENDDO
      
      xYU(1,1)= 0.0027 /SINB/VEV 
      xYU(2,2)= 0.620  /SINB/VEV
C      xYU(3,3)= GISA(6)
      xYU(3,3)= 0.D0
      DO I=1,3
        DO J=1,3
          DO K=1,3
            YU(I,J)=YU(I,J)+CKM(K,I)*xYU(K,J)
          ENDDO
        ENDDO
      ENDDO
      HU(3,3)= GISA(12)*GISA(6)
      YD(1,1)= 0.0058 /COSB/VEV 
      YD(2,2)= 0.113  /COSB/VEV 
      YD(3,3)= GISA(5)
      HD(3,3)= GISA(11)*GISA(5)
C...(s)lepton sector
      YE(1,1)= 0.00051/COSB/VEV 
      YE(2,2)= 0.106  /COSB/VEV
      YE(3,3)= GISA(4)
      HE(3,3)= GISA(10)*GISA(4)

C...Set RGE parameters at M_Z in DRbar      
c   NOTE difference in Yukawa and trilinear definitions!      
      G(1)=GISA(1)
      G(2)=GISA(2)
      G(3)=GISA(3)

      K=0
      DO I=1,3
        DO J=1,3
          G(4+K) = YU(J,I)
          G(13+K)= YD(J,I)
          G(22+K)= YE(J,I)
          G(34+K)= HU(J,I)
          G(43+K)= HD(J,I)
          G(52+K)= HE(J,I)
          K=K+1
        ENDDO
      ENDDO
      
      G(31)=GISA(7)
      G(32)=GISA(8)
      G(33)=GISA(9)
      
      DO I=65,157
        G(I)=0.d0
      ENDDO
      G(61)=GISA(25)
      G(62)=GISA(26)*GISA(25)
      G(63)=GISA(14)
      G(64)=GISA(13)
      
      G(65)=GISA(19)
      G(69)=GISA(19)
      G(73)=GISA(24)
      
      G(74)=GISA(16)
      G(78)=GISA(16)
      G(82)=GISA(21)
      
      G(83)=GISA(18)
      G(87)=GISA(18)
      G(91)=GISA(23)
      
      G(92)=GISA(17)
      G(96)=GISA(17)
      G(100)=GISA(22)
      
      G(101)=GISA(15)
      G(105)=GISA(15)
      G(109)=GISA(20)
      
      G(110)=GISA(31)
      G(111)=GISA(30)
      
      G(120)=GISA(27)
      G(129)=GISA(28)
      G(138)=GISA(29)*GISA(27)
      G(147)=0.d0
      G(156)=0.d0
      G(157)=0.5
C-------------------------------------------------------------------
C----- RUN RGEs FROM MZ TO MGUT -------------------------------------
C-------------------------------------------------------------------
      TZ=DLOG(MZ/MGUT)
      TGUT=0.d0
      DT=(TGUT-TZ)/FLOAT(NSTEP)
      DO II=1,NSTEP
        T=TZ+DT*FLOAT(II-1)
        QOLD=Q
        Q=MGUT*EXP(T)
c...Correct up Yukawa matrix at M_top
        IF (QOLD.LE.MT.AND.Q.GT.MT) THEN
          SFTMT=SNGL(FTMT)
	  CALL TACTIV(4,SFTMT,G)
	ENDIF
C...Implement sparticle threshold corrections at Q=HIGFRZ
        IF (QOLD.LE.MSUSY.AND.Q.GT.MSUSY) THEN
	  CALL YUKDIAG(G,4,xYU,VULQ,VURQ)
	  xYU(3,3)=xYU(3,3)/(1.d0-DBLE(RSIGT))
	  CALL AYUKDIAG(G,4,xYU,VULQ,VURQ)
	  CALL YUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  xYD(3,3)=xYD(3,3)/(1.d0-DBLE(RSIGB))
	  CALL AYUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  G(30)=G(30)/(1.d0-DBLE(RSIGL))
      CALL GLUNENO(G,GAM,AGE,AHE)
        END IF
c...neutrino sector
	IF (IRHN.GT.0) THEN
	  IF (QOLD.LE.DABS(MRNDEC(3)).AND.Q.GT.DABS(MRNDEC(3))) THEN
            G(120)=GISA(27)
            G(129)=GISA(28)
            G(138)=GISA(29)*GISA(27)
	  ENDIF
        ENDIF
        CALL DRKSTP(157,DT,T,G,RGE157,W2)   
      END DO
C----------------------------------------------------------------------
C----- SET GUT-scale boundary conditions ------------------------------
C----------------------------------------------------------------------

c...Yukawa couplings      
      CALL VEC2MAT(G,4,YUGUT,-1)
      CALL VEC2MAT(G,13,YDGUT,-1)
      CALL VEC2MAT(G,22,YEGUT,-1)
      IF (IRHN.EQ.1)  CALL VEC2MAT(G,112,YNGUT,-1)
      
      CALL BSGGUT(YUGUT,YDGUT,YEGUT,YNGUT,G,IMODEL,M0,MHF,A0)
c...reset freezeout flags
      DO I=1,157
        IF(I.LT.148) THEN 
	  IGF(I)=0
	ELSE
	  IGF(I)=1
	ENDIF
	GF(I)=0.d0
      ENDDO
C----------------------------------------------------------------------
C-------- SET VALUES OF WILSON COEFFICIENTS ---------------------------
C----------------------------------------------------------------------
      CALL WILSON(GK,HK,
     &  	  CSM,CH1,CH2,CC111,CC112,CC113,CC121,CC122,CC123,
     &  	  CC211,CC212,CC213,CC221,CC222,CC223,
     &  	  DSSM,DSH1,DSH2,DSC111,DSC112,DSC113,DSC121,
     &  	  DSC122,DSC123,DSC211,DSC212,DSC213,DSC221,
     &  	  DSC222,DSC223)
C
      c3smnai=  -1./2.*(CSM(1)+DSSM(1))
     $          +(CSM(3)+DSSM(3))
     $          -1./2.*(CSM(4)+DSSM(4)+2./9.d0)
     $          -1./4.*(CSM(5)+DSSM(5)-7./9.d0)
     $          +1./4.*(CSM(7)+DSSM(7)+1.d0)
     $          -1./4.*(CSM(9)+DSSM(9)+9.d0)          
      c2smnai=  -1./2.*(CSM(1)+DSSM(1))
     $          +(CSM(2)+DSSM(2))
     $          -1./2.*(CSM(4)+DSSM(4)+2./9.d0)
     $          -1./4.*(CSM(5)+DSSM(5)-7./9.d0)
     $          +1./4.*(CSM(7)+DSSM(7)+1.d0)         
      c3hpnai=  -1./2.d0*(CH1(1)+DSH1(1))
     $          +(CH1(3)+DSH1(3))
     $          -1./2.d0*(CH1(4)+DSH1(4))
     $          -1./4.d0*(CH1(5)+DSH1(5))
     $          +1./4.d0*(CH1(7)+DSH1(7))
     $          -1./4.d0*(CH1(9)+DSH1(9))        
      c3c11nai= -1./2.d0*(CC211(1)+DSC211(1))                                
     $          +(CC211(3)+DSC211(3))
     $          -1./2.d0*(CC211(4)+DSC211(4))
     $          -1./4.d0*(CC211(5)+DSC211(5))
     $          +1./4.d0*(CC211(7)+DSC211(7))
     $          -1./4.d0*(CC211(9)+DSC211(9))                  
      c3c12nai= -1./2.d0*(CC212(1)+DSC212(1))
     $          +(CC212(3)+DSC212(3))
     $          -1./2.d0*(CC212(4)+DSC212(4))
     $          -1./4.d0*(CC212(5)+DSC212(5))
     $          +1./4.d0*(CC212(7)+DSC212(7))
     $          -1./4.d0*(CC212(9)+DSC212(9))                  
      c3c13nai= -1./2.d0*(CC213(1)+DSC213(1))
     $          +(CC213(3)+DSC213(3))
     $          -1./2.d0*(CC213(4)+DSC213(4))
     $          -1./4.d0*(CC213(5)+DSC213(5))
     $          +1./4.d0*(CC213(7)+DSC213(7))
     $          -1./4.d0*(CC213(9)+DSC213(9))                    
      c3c21nai= -1./2.d0*(CC221(1)+DSC221(1))
     $          +(CC221(3)+DSC221(3))
     $          -1./2.d0*(CC221(4)+DSC221(4))
     $          -1./4.d0*(CC221(5)+DSC221(5))
     $          +1./4.d0*(CC221(7)+DSC221(7))
     $          -1./4.d0*(CC221(9)+DSC221(9))                         
      c3c22nai= -1./2.d0*(CC222(1)+DSC222(1))
     $          +(CC222(3)+DSC222(3))
     $          -1./2.d0*(CC222(4)+DSC222(4))
     $          -1./4.d0*(CC222(5)+DSC222(5))
     $          +1./4.d0*(CC222(7)+DSC222(7))
     $          -1./4.d0*(CC222(9)+DSC222(9))                               
      c3c23nai= -1./2.d0*(CC223(1)+DSC223(1))
     $          +(CC223(3)+DSC223(3))
     $          -1./2.d0*(CC223(4)+DSC223(4))
     $          -1./4.d0*(CC223(5)+DSC223(5))
     $          +1./4.d0*(CC223(7)+DSC223(7))
     $          -1./4.d0*(CC223(9)+DSC223(9))                           
C----------------------------------------------------------------------
C------- NOW WE RUN from MGUT past MZ to Q=1 GeV ----------------------
C----------------------------------------------------------------------
      TZ=LOG(1.d0/MGUT)
      DT=TZ/FLOAT(NSTEP)
C
      DO 400 II=1,NSTEP+2   
        T=    DT*FLOAT(II-1)
        Q   = MGUT*EXP(T)
        QOLD= MGUT*EXP(T-dt)
        QNEW= MGUT*EXP(T+dt)
C
        RMIX(1)=0.d0
        RMIX(2)=0.d0 
        SMIX(1)=0.d0
        SMIX(2)=0.d0  
        SMIX(3)=0.d0
        SMIX(4)=0.d0  
        SMIX(5)=0.d0
        SMIX(6)=0.d0  
        ABMIX(1)=0.d0
        ABMIX(2)=0.d0      
C
C-------- Evolution of MSSM parameters --------------------------------
C
        CALL DRKSTP(157,DT,T,G,RGE157,W2)
C...TEST YUKAWA DIVERGENCE
        DO I=4,30
	  IF (G(i).gt.5.d0) THEN
            PRINT*,'ISABSG: NON-PERTURBATIVE YUKAWA'
	    PRINT*,'G(',i,')>= 5:',G(i),' Q=',Q
          END IF
	ENDDO
c...Decouple top quark in  up Yukawa elements at M_top
        IF (Q.LT.MT.AND.QOLD.GE.MT) THEN
	  CALL YUKDIAG(G,4,xYUMT,VULMT,VURMT)
	  xYUMT(3,3)=0.d0
	  CALL AYUKDIAG(G,4,xYUMT,VULMT,VURMT)
	ENDIF
C...Implement sparticle threshold corrections at Q=HIGFRZ
        IF (QOLD.GT.MSUSY.AND.Q.LE.MSUSY) THEN
      CALL GLUNENO(G,GAM,AGE,AHE)
	  CALL YUKDIAG(G,4,xYU,VULQ,VURQ)
	  xYU(3,3)=xYU(3,3)*(1.d0-DBLE(RSIGT))
	  CALL AYUKDIAG(G,4,xYU,VULQ,VURQ)
	  CALL YUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  xYD(3,3)=xYD(3,3)*(1.d0-DBLE(RSIGB))
	  CALL AYUKDIAG(G,13,xYD,VDLQ,VDRQ)
	  G(30)=G(30)*(1.d0-DBLE(RSIGL))
        END IF
c...decoupling in neutrino sector
        IF (QNEW.LT.DABS(MRNDEC(3))) THEN
          G(120)=0.d0
          G(129)=0.d0
          G(138)=0.d0
	ENDIF
C
C------- Evolution of Wilson coefficients -----------------------------
C
        G3=G(3)
        ALES=G(3)**2/4.d0/PI
        ALEL=DBLE(SUALFE(SNGL(MW**2)))
C
c ...  SM contribution
C
       IF(Q.LT.MT.AND.Q.GT.MW) THEN
        T=DT*FLOAT(II-1)
        CALL DRKSTP(9,DT,T,CSM,GAMMASM,WSM)
       ENDIF

       IF(Q.GT.MW.AND.QNEW.LT.MW) THEN  
         CSM(1)=CSM(1)+DSSM(1)
         CSM(2)=CSM(2)+DSSM(2)
         CSM(3)=CSM(3)+DSSM(3)       
         CSM(4)=CSM(4)+2./9.d0+DSSM(4)
         CSM(5)=CSM(5)-7./9.d0+DSSM(5)
         CSM(6)=CSM(6)+2./9.d0+DSSM(6)
         CSM(7)=CSM(7)+1.d0+DSSM(7)
         CSM(8)=CSM(8)-3./2.d0+DSSM(8)
         CSM(9)=CSM(9)+9.d0+DSSM(9)             
       ENDIF
        C2SMEFF=-1./2.d0*CSM(1)+CSM(2)-1./2.d0*CSM(4)-1./4.d0*CSM(5)
     $          +1./4.d0*CSM(7)            
        C3SMEFF=-1./2.d0*CSM(1)+CSM(3)-1./2.d0*CSM(4)-1./4.d0*CSM(5)
     $          +1./4.d0*CSM(7)-1./4.d0*CSM(9) 
C
C ...  CH (charged higgs) contribution  
C
       IF(Q.GT.MHPLUS.AND.QNEW.LT.MHPLUS) THEN
        G3MHP=G(3)
        X=(MT/MHPLUS)**2
        CH2(11)=-1./2.d0*X/TANB**2
        CH2(17)= 16.*PI**2*X/G3MHP**2
       ENDIF
C       
       IF(MHPLUS.LT.MT) THEN
        IF(Q.LT.MT.AND.Q.GT.MW) THEN
         T=DT*FLOAT(II-1)
         CALL DRKSTP(9,DT,T,CH1,GAMMASM,WHP1)     
         IF(Q.GT.MHPLUS.AND.QNEW.LE.MHPLUS) THEN
          CH1(1)=CH1(1)+DSH1(1)   
          CH1(2)=CH1(2)+DSH1(2)
          CH1(3)=CH1(3)+DSH1(3)
          CH1(4)=CH1(4)+DSH1(4)
          CH1(5)=CH1(5)+DSH1(5)
          CH1(6)=CH1(6)+DSH1(6)
          CH1(7)=CH1(7)+DSH1(7)
          CH1(8)=CH1(8)+DSH1(8)
          CH1(9)=CH1(9)+DSH1(9)    
         ENDIF
         C2HPEFF=-1./2.d0*CH1(1)+CH1(2)-1./2.d0*CH1(4)-1./4.d0*CH1(5)
     $           +1./4.d0*CH1(7)      
         C3HPEFF=-1./2.d0*CH1(1)+CH1(3)-1./2.d0*CH1(4)-1./4.d0*CH1(5)
     $           +1./4.d0*CH1(7)-1./4.d0*CH1(9)       
        ENDIF
       ELSE
        IF(Q.LT.MHPLUS.AND.Q.GT.MT) THEN
         T=DT*FLOAT(II-1)
         CALL DRKSTP(19,DT,T,CH2,GAMMAHP,WHP2)     
          SMIX(1)=CH2(10)
          SMIX(2)=CH2(11)
          SMIX(3)=CH2(12)
          SMIX(4)=CH2(13)
          SMIX(5)=CH2(14)
          SMIX(6)=CH2(15)
C
         IF(Q.GT.MT.AND.QNEW.LE.MT) THEN
          CH2(1)=CH2(1)+DSH2(1) 
          CH2(2)=CH2(2)+DSH2(2)
          CH2(3)=CH2(3)+DSH2(3)
          CH2(4)=CH2(4)+DSH2(4)
          CH2(5)=CH2(5)+DSH2(5)
          CH2(6)=CH2(6)+DSH2(6)
          CH2(7)=CH2(7)+DSH2(7)
          CH2(8)=CH2(8)+DSH2(8)
          CH2(9)=CH2(9)+DSH2(9)    
         ENDIF
        ENDIF                           
C       
       IF(Q.LT.MT.AND.Q.GT.MW) THEN
         T=DT*FLOAT(II-1)
         CALL DRKSTP(9,DT,T,CH2,GAMMASM,WHP2)    
C
         C2HPEFF=-1./2.d0*CH2(1)+CH2(2)-1./2.d0*CH2(4)-1./4.d0*CH2(5)
     $           +1./4.*CH2(7)       
         C3HPEFF=-1./2.d0*CH2(1)+CH2(3)-1./2.d0*CH2(4)-1./4.d0*CH2(5)
     $           +1./4.d0*CH2(7)-1./4.d0*CH2(9) 
        
	ENDIF
       ENDIF
C
C ...  chargino(1)-squark(1) contribution  
C
       MSQ=MSQU(1)
       MCH=MCHA(1)    
C
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(1))**2
        Y=MCHA(1)/MBQ
        CC111(10)= GK(1,1)*HK(1,1)*16.*PI**2*X*Y/G3SS**2
        CC111(11)= GK(1,1)*HK(1,1)*16.*PI**2*X*Y/G3SS**2
        CC111(12)=-GK(1,1)*GK(1,1)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(1).GT.MCHA(1)) THEN
         IF(MCHA(1).GT.(MW+0.5)) THEN
           QMAX=MCH
         ELSE
           QMAX=MW
         ENDIF
         IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN  
           T=DT*FLOAT(II-1)
           CALL DRKSTP(12,DT,T,CC111,GAMMAC1,WSC111)        
         ENDIF
         IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
	   CC111(1)=CC111(1)+DSC111(1)        
           CC111(2)=CC111(2)+DSC111(2)
           CC111(3)=CC111(3)+DSC111(3)
           CC111(4)=CC111(4)+DSC111(4)
           CC111(5)=CC111(5)+DSC111(5)
           CC111(6)=CC111(6)+DSC111(6)
           CC111(7)=CC111(7)+DSC111(7)
           CC111(8)=CC111(8)+DSC111(8)
           CC111(9)=CC111(9)+DSC111(9)
 	 ENDIF
         IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
           T=DT*FLOAT(II-1)
           CALL DRKSTP(9,DT,T,CC111,GAMMASM,WSC111)       
         ENDIF
c
         C2C11EFF=-1./2.d0*CC111(1)+CC111(2)
     $            -1./2.d0*CC111(4)-1./4.d0*CC111(5)
     $            +1./4.d0*CC111(7)       
         C3C11EFF=-1./2.d0*CC111(1)+CC111(3)
     $            -1./2.d0*CC111(4)-1./4.d0*CC111(5)
     $            +1./4.d0*CC111(7)-1./4.d0*CC111(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC211(13)
        RMIX(2)=RMIX(2)+SWI1*CC211(14)       
        IF(MSQU(1).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(2)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                
        ABMIX(1)=SWI2*CC212(13)+SWI3*CC213(13)
        ABMIX(2)=SWI2*CC212(14)+SWI3*CC213(14)
C 
        IF(MSQU(1).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC211,GAMMAC2,WSC211)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
         CC211(1)=CC211(1)+DSC211(1)        
         CC211(2)=CC211(2)+DSC211(2)
         CC211(3)=CC211(3)+DSC211(3)
         CC211(4)=CC211(4)+DSC211(4)
         CC211(5)=CC211(5)+DSC211(5)
         CC211(6)=CC211(6)+DSC211(6)
         CC211(7)=CC211(7)+DSC211(7)
         CC211(8)=CC211(8)+DSC211(8)
         CC211(9)=CC211(9)+DSC211(9)        
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC211,GAMMASM,WSC211)      
        ENDIF                  
         C2C11EFF=-1./2.d0*CC211(1)+CC211(2)
     $            -1./2.d0*CC211(4)-1./4.d0*CC211(5)
     $            +1./4.d0*CC211(7)      
         C3C11EFF=-1./2.d0*CC211(1)+CC211(3)
     $            -1./2.d0*CC211(4)-1./4.d0*CC211(5)
     $            +1./4.d0*CC211(7)-1./4.d0*CC211(9)       
       ENDIF
C
C ...  chargino(1)-squark(2) contribution  
C
       MSQ=MSQU(2)
       MCH=MCHA(1)    
c
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(2))**2
        Y=MCHA(1)/MBQ
        CC112(10)= GK(1,2)*HK(1,2)*16.*PI**2*X*Y/G3SS**2
        CC112(11)= GK(1,2)*HK(1,2)*16.*PI**2*X*Y/G3SS**2
        CC112(12)=-GK(1,2)*GK(1,2)*16.*PI**2*X/G3SS**2
       ENDIF  
C
       IF(MSQU(2).GT.MCHA(1)) THEN
        IF(MCHA(1).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC112,GAMMAC1,WSC112)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
         CC112(1)=CC112(1)+DSC112(1)        
         CC112(2)=CC112(2)+DSC112(2)
         CC112(3)=CC112(3)+DSC112(3)
         CC112(4)=CC112(4)+DSC112(4)
         CC112(5)=CC112(5)+DSC112(5)
         CC112(6)=CC112(6)+DSC112(6)
         CC112(7)=CC112(7)+DSC112(7)
         CC112(8)=CC112(8)+DSC112(8)
         CC112(9)=CC112(9)+DSC112(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC112,GAMMASM,WSC112)      
        ENDIF
         C2C12EFF=-1./2.d0*CC112(1)+CC112(2)
     $            -1./2.d0*CC112(4)-1./4.d0*CC112(5)
     $            +1./4.d0*CC112(7)       
         C3C12EFF=-1./2.d0*CC112(1)+CC112(3)
     $            -1./2.d0*CC112(4)-1./4.d0*CC112(5)
     $            +1./4.d0*CC112(7)-1./4.d0*CC112(9)
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC212(13)
        RMIX(2)=RMIX(2)+SWI1*CC212(14)
        IF(MSQU(2).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                
        ABMIX(1)=SWI2*CC211(13)+SWI3*CC213(13)
        ABMIX(2)=SWI2*CC211(14)+SWI3*CC213(14)
C 
        IF(MSQU(2).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC212,GAMMAC2,WSC212)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC212(1)=CC212(1)+DSC212(1)        
         CC212(2)=CC212(2)+DSC212(2)
         CC212(3)=CC212(3)+DSC212(3)
         CC212(4)=CC212(4)+DSC212(4)
         CC212(5)=CC212(5)+DSC212(5)
         CC212(6)=CC212(6)+DSC212(6)
         CC212(7)=CC212(7)+DSC212(7)
         CC212(8)=CC212(8)+DSC212(8)
         CC212(9)=CC212(9)+DSC212(9)        
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC212,GAMMASM,WSC212)      
        ENDIF                  
         C2C12EFF=-1./2.d0*CC212(1)+CC212(2)
     $            -1./2.d0*CC212(4)-1./4.d0*CC212(5)
     $            +1./4.d0*CC212(7)       
         C3C12EFF=-1./2.d0*CC212(1)+CC212(3)
     $            -1./2.d0*CC212(4)-1./4.d0*CC212(5)
     $            +1./4.d0*CC212(7)-1./4.d0*CC212(9) 
       ENDIF
C  
C ...  chargino(1)-squark(3) contribution  
C
       MSQ=MSQU(3)
       MCH=MCHA(1)    
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(3))**2
        Y=MCHA(1)/MBQ
        CC113(10)= HK(1,3)*16.*PI**2*X*Y/G3SS**2
        CC113(11)= HK(1,3)*16.*PI**2*X*Y/G3SS**2
        CC113(12)=-GK(1,3)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(3).GT.MCHA(1)) THEN
        IF(MCHA(1).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC113,GAMMAC1,WSC113)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC113(1)=CC113(1)+DSC113(1)        
         CC113(2)=CC113(2)+DSC113(2)
         CC113(3)=CC113(3)+DSC113(3)
         CC113(4)=CC113(4)+DSC113(4)
         CC113(5)=CC113(5)+DSC113(5)
         CC113(6)=CC113(6)+DSC113(6)
         CC113(7)=CC113(7)+DSC113(7)
         CC113(8)=CC113(8)+DSC113(8)
         CC113(9)=CC113(9)+DSC113(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC113,GAMMASM,WSC113)      
        ENDIF
         C2C13EFF=-1./2.d0*CC113(1)+CC113(2)
     $            -1./2.d0*CC113(4)-1./4.d0*CC113(5)
     $            +1./4.d0*CC113(7)       
         C3C13EFF=-1./2.d0*CC113(1)+CC113(3)
     $            -1./2.d0*CC113(4)-1./4.d0*CC113(5)
     $            +1./4.d0*CC113(7)-1./4.d0*CC113(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC213(13)
        RMIX(2)=RMIX(2)+SWI1*CC213(14)
        IF(MSQU(3).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(1).AND.Q.GT.MSQU(2)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                
        ABMIX(1)=SWI2*CC211(13)+SWI3*CC212(13)
        ABMIX(2)=SWI2*CC211(14)+SWI3*CC212(14)
C 
        IF(MSQU(3).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC213,GAMMAC2,WSC213)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LT.QMAX) THEN
         CC213(1)=CC213(1)+DSC213(1)        
         CC213(2)=CC213(2)+DSC213(2)
         CC213(3)=CC213(3)+DSC213(3)
         CC213(4)=CC213(4)+DSC213(4)
         CC213(5)=CC213(5)+DSC213(5)
         CC213(6)=CC213(6)+DSC213(6)
         CC213(7)=CC213(7)+DSC213(7)
         CC213(8)=CC213(8)+DSC213(8)
         CC213(9)=CC213(9)+DSC213(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC213,GAMMASM,WSC213)      
        ENDIF                  
         C2C13EFF=-1./2.d0*CC213(1)+CC213(2)
     $            -1./2.d0*CC213(4)-1./4.d0*CC213(5)
     $            +1./4.d0*CC213(7)       
         C3C13EFF=-1./2.d0*CC213(1)+CC213(3)
     $            -1./2.d0*CC213(4)-1./4.d0*CC213(5)
     $            +1./4.d0*CC213(7)-1./4.d0*CC213(9)       
       ENDIF
C
C ...  chargino(2)-squark(1) contribution  
C
       MSQ=MSQU(1)
       MCH=MCHA(2)    
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(1))**2
        Y=MCHA(2)/MBQ
        CC121(10)= GK(2,1)*HK(2,1)*16.*PI**2*X*Y/G3SS**2
        CC121(11)= GK(2,1)*HK(2,1)*16.*PI**2*X*Y/G3SS**2
        CC121(12)=-GK(2,1)*GK(2,1)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(1).GT.MCHA(2)) THEN
        IF(MCHA(2).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC121,GAMMAC1,WSC121)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC121(1)=CC121(1)+DSC121(1)        
         CC121(2)=CC121(2)+DSC121(2)
         CC121(3)=CC121(3)+DSC121(3)
         CC121(4)=CC121(4)+DSC121(4)
         CC121(5)=CC121(5)+DSC121(5)
         CC121(6)=CC121(6)+DSC121(6)
         CC121(7)=CC121(7)+DSC121(7)
         CC121(8)=CC121(8)+DSC121(8)
         CC121(9)=CC121(9)+DSC121(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC121,GAMMASM,WSC121)      
        ENDIF
         C2C21EFF=-1./2.d0*CC121(1)+CC121(2)
     $            -1./2.d0*CC121(4)-1./4.d0*CC121(5)
     $            +1./4.d0*CC121(7)       
         C3C21EFF=-1./2.d0*CC121(1)+CC121(3)
     $            -1./2.d0*CC121(4)-1./4.d0*CC121(5)
     $            +1./4.d0*CC121(7)-1./4.d0*CC121(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC221(13)
        RMIX(2)=RMIX(2)+SWI1*CC221(14)
        IF(MSQU(1).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(2)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                    
        ABMIX(1)=SWI2*CC222(13)+SWI3*CC223(13)
        ABMIX(2)=SWI2*CC222(14)+SWI3*CC223(14)
C 
        IF(MSQU(1).GT.MW) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC221,GAMMAC2,WSC221)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC221(1)=CC221(1)+DSC221(1)        
         CC221(2)=CC221(2)+DSC221(2)
         CC221(3)=CC221(3)+DSC221(3)
         CC221(4)=CC221(4)+DSC221(4)
         CC221(5)=CC221(5)+DSC221(5)
         CC221(6)=CC221(6)+DSC221(6)
         CC221(7)=CC221(7)+DSC221(7)
         CC221(8)=CC221(8)+DSC221(8)
         CC221(9)=CC221(9)+DSC221(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC221,GAMMASM,WSC221)      
        ENDIF                  
         C2C21EFF=-1./2.d0*CC221(1)+CC221(2)
     $            -1./2.d0*CC221(4)-1./4.d0*CC221(5)
     $            +1./4.d0*CC221(7)      
         C3C21EFF=-1./2.d0*CC221(1)+CC221(3)
     $            -1./2.d0*CC221(4)-1./4.d0*CC221(5)
     $            +1./4.d0*CC221(7)-1./4.d0*CC221(9)       
       ENDIF
C
C ...  chargino(2)-squark(2) contribution  
C
       MSQ=MSQU(2)
       MCH=MCHA(2)    
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(2))**2
        Y=MCHA(2)/MBQ
        CC122(10)= GK(2,2)*HK(2,2)*16.*PI**2*X*Y/G3SS**2
        CC122(11)= GK(2,2)*HK(2,2)*16.*PI**2*X*Y/G3SS**2
        CC122(12)=-GK(2,2)*GK(2,2)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(2).GT.MCHA(2)) THEN
        IF(MCHA(2).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC122,GAMMAC1,WSC122)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC122(1)=CC122(1)+DSC122(1)        
         CC122(2)=CC122(2)+DSC122(2)
         CC122(3)=CC122(3)+DSC122(3)
         CC122(4)=CC122(4)+DSC122(4)
         CC122(5)=CC122(5)+DSC122(5)
         CC122(6)=CC122(6)+DSC122(6)
         CC122(7)=CC122(7)+DSC122(7)
         CC122(8)=CC122(8)+DSC122(8)
         CC122(9)=CC122(9)+DSC122(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC122,GAMMASM,WSC122)      
C
        ENDIF
         C2C22EFF=-1./2.d0*CC122(1)+CC122(2)
     $            -1./2.d0*CC122(4)-1./4.d0*CC122(5)
     $            +1./4.d0*CC122(7)           
         C3C22EFF=-1./2.d0*CC122(1)+CC122(3)
     $            -1./2.d0*CC122(4)-1./4.d0*CC122(5)
     $            +1./4.d0*CC122(7)-1./4.d0*CC122(9)           
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC222(13)
        RMIX(2)=RMIX(2)+SWI1*CC222(14)
        IF(MSQU(2).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(3)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                    
        ABMIX(1)=SWI2*CC221(13)+SWI3*CC223(13)
        ABMIX(2)=SWI2*CC221(14)+SWI3*CC223(14)
C 
        IF(MSQU(2).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC222,GAMMAC2,WSC222)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC222(1)=CC222(1)+DSC222(1)        
         CC222(2)=CC222(2)+DSC222(2)
         CC222(3)=CC222(3)+DSC222(3)
         CC222(4)=CC222(4)+DSC222(4)
         CC222(5)=CC222(5)+DSC222(5)
         CC222(6)=CC222(6)+DSC222(6)
         CC222(7)=CC222(7)+DSC222(7)
         CC222(8)=CC222(8)+DSC222(8)
         CC222(9)=CC222(9)+DSC222(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC222,GAMMASM,WSC222)      
        ENDIF                  
         C2C22EFF=-1./2.d0*CC222(1)+CC222(2)
     $            -1./2.d0*CC222(4)-1./4.d0*CC222(5)
     $            +1./4.d0*CC222(7)       
         C3C22EFF=-1./2.d0*CC222(1)+CC222(3)
     $            -1./2.d0*CC222(4)-1./4.d0*CC222(5)
     $            +1./4.d0*CC222(7)-1./4.d0*CC222(9)       
       ENDIF
C
C ...  chargino(2)-squark(3) contribution  
C
       MSQ=MSQU(3)
       MCH=MCHA(2)   
       IF(Q.GT.MSQ.AND.QNEW.LE.MSQ) THEN
        G3SS=G(3)
        X=(MW/MSQU(3))**2
        Y=MCHA(2)/MBQ
        CC123(10)= HK(2,3)*16.*PI**2*X*Y/G3SS**2
        CC123(11)= HK(2,3)*16.*PI**2*X*Y/G3SS**2
        CC123(12)=-GK(2,3)*16.*PI**2*X/G3SS**2
       ENDIF    
C
       IF(MSQU(3).GT.MCHA(2)) THEN
        IF(MCHA(2).GT.(MW+0.5)) THEN
         QMAX=MCH
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MSQ.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(12,DT,T,CC123,GAMMAC1,WSC123)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC123(1)=CC123(1)+DSC123(1)        
         CC123(2)=CC123(2)+DSC123(2)
         CC123(3)=CC123(3)+DSC123(3)
         CC123(4)=CC123(4)+DSC123(4)
         CC123(5)=CC123(5)+DSC123(5)
         CC123(6)=CC123(6)+DSC123(6)
         CC123(7)=CC123(7)+DSC123(7)
         CC123(8)=CC123(8)+DSC123(8)
         CC123(9)=CC123(9)+DSC123(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC123,GAMMASM,WSC123)      
        ENDIF
         C2C23EFF=-1./2.d0*CC123(1)+CC123(2)
     $            -1./2.d0*CC123(4)-1./4.d0*CC123(5)
     $            +1./4.d0*CC123(7)       
         C3C23EFF=-1./2.d0*CC123(1)+CC123(3)
     $            -1./2.d0*CC123(4)-1./4.d0*CC123(5)
     $            +1./4.d0*CC123(7)-1./4.d0*CC123(9)       
       ELSE 
        IF(Q.LT.MCH.AND.Q.GT.MSQ) THEN
         SWI1=1.
        ELSE
         SWI1=0.
        ENDIF              
        RMIX(1)=RMIX(1)+SWI1*CC223(13)
        RMIX(2)=RMIX(2)+SWI1*CC223(14)
        IF(MSQU(3).GT.MT) THEN
         SMIX(1)=0.
         SMIX(2)=0.
         SMIX(3)=0.
         SMIX(4)=0.
         SMIX(5)=0.
         SMIX(6)=0.
        ENDIF
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(1)) THEN
         SWI2=1.
        ELSE
         SWI2=0.
        ENDIF                
        IF(Q.LT.MCHA(2).AND.Q.GT.MSQU(2)) THEN
         SWI3=1.
        ELSE
         SWI3=0.
        ENDIF                    
        ABMIX(1)=SWI2*CC221(13)+SWI3*CC222(13)
        ABMIX(2)=SWI2*CC221(14)+SWI3*CC222(14)
C 
        IF(MSQU(3).GT.(MW+0.5)) THEN
         QMAX=MSQ
        ELSE
         QMAX=MW
        ENDIF
        IF(Q.LT.MCH.AND.Q.GT.QMAX) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(17,DT,T,CC223,GAMMAC2,WSC223)      
        ENDIF
        IF(Q.GT.QMAX.AND.QNEW.LE.QMAX) THEN
         CC223(1)=CC223(1)+DSC223(1)        
         CC223(2)=CC223(2)+DSC223(2)
         CC223(3)=CC223(3)+DSC223(3)
         CC223(4)=CC223(4)+DSC223(4)
         CC223(5)=CC223(5)+DSC223(5)
         CC223(6)=CC223(6)+DSC223(6)
         CC223(7)=CC223(7)+DSC223(7)
         CC223(8)=CC223(8)+DSC223(8)
         CC223(9)=CC223(9)+DSC223(9)
        ENDIF
        IF(QMAX.GT.MW.AND.Q.LT.QMAX.AND.Q.GT.MW) THEN
          T=DT*FLOAT(II-1)
          CALL DRKSTP(9,DT,T,CC223,GAMMASM,WSC223)      
        ENDIF                  
         C2C23EFF=-1./2.d0*CC223(1)+CC223(2)
     $            -1./2.d0*CC223(4)-1./4.d0*CC223(5)
     $            +1./4.d0*CC223(7)       
         C3C23EFF=-1./2.d0*CC223(1)+CC223(3)
     $            -1./2.d0*CC223(4)-1./4.d0*CC223(5)
     $            +1./4.d0*CC223(7)-1./4.d0*CC223(9)       
       ENDIF
c
c ... gluino and neutralinos contributions ...
c      GLUNENO must be recalled at MSUSY to set squark mass matrices
c
       IF(Q.GT.MW.AND.QNEW.LE.MW) THEN  

       CER=4./3.d0
       CEG=3.d0
       ALFA2=ALFAEM/SN2THW
       RALPH=ALES/ALFA2
       FACG7= 6.d0*(-1./3.d0)*RALPH*CER/(CKM(3,2)*CKM(3,3))
       FACG8=-RALPH/(CKM(3,2)*CKM(3,3))
       FACN7= 3.d0*(-1./3.d0)/(CKM(3,2)*CKM(3,3))
       FACN8=-1.d0/(CKM(3,2)*CKM(3,3))
C
C ...  gluino-dsquark1 ...
C              
       X=MGLU**2/MDSB(1)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS1EFF=FACG7*(MW/MDSB(1))**2
     &          *(GAM(3,1)*GAM(2,1)*FN2
     $           -GAM(6,1)*GAM(2,1)*MGLU/MBQ*FN4)
       C2GS1EFF=FACG8*(MW/MDSB(1))**2
     &          *(GAM(3,1)*GAM(2,1)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,1)*GAM(2,1)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark2 ...
C             
       X=MGLU**2/MDSB(2)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS2EFF=FACG7*(MW/MDSB(2))**2
     &          *(GAM(3,2)*GAM(2,2)*FN2
     $           -GAM(6,2)*GAM(2,2)*MGLU/MBQ*FN4)
       C2GS2EFF=FACG8*(MW/MDSB(2))**2
     &          *(GAM(3,2)*GAM(2,2)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,2)*GAM(2,2)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark3 ...
C       
       X=MGLU**2/MDSB(3)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS3EFF=FACG7*(MW/MDSB(3))**2
     &          *(GAM(3,3)*GAM(2,3)*FN2
     $           -GAM(6,3)*GAM(2,3)*MGLU/MBQ*FN4)
       C2GS3EFF=FACG8*(MW/MDSB(3))**2
     &          *(GAM(3,3)*GAM(2,3)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,3)*GAM(2,3)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark4 ...
C       
       X=MGLU**2/MDSB(4)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS4EFF=FACG7*(MW/MDSB(4))**2
     &          *(GAM(3,4)*GAM(2,4)*FN2
     $           -GAM(6,4)*GAM(2,4)*MGLU/MBQ*FN4)
       C2GS4EFF=FACG8*(MW/MDSB(4))**2
     &          *(GAM(3,4)*GAM(2,4)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,4)*GAM(2,4)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark5 ...
C       
       X=MGLU**2/MDSB(5)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS5EFF=FACG7*(MW/MDSB(5))**2
     &          *(GAM(3,5)*GAM(2,5)*FN2
     $           -GAM(6,5)*GAM(2,5)*MGLU/MBQ*FN4)
       C2GS5EFF=FACG8*(MW/MDSB(5))**2
     &          *(GAM(3,5)*GAM(2,5)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,5)*GAM(2,5)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  gluino-dsquark6 ...
C       
       X=MGLU**2/MDSB(6)**2
       CALL FUNS(X,FN1,FN2,FN3,FN4) 
       C3GS6EFF=FACG7*(MW/MDSB(6))**2
     &          *(GAM(3,6)*GAM(2,6)*FN2
     $           -GAM(6,6)*GAM(2,6)*MGLU/MBQ*FN4)
       C2GS6EFF=FACG8*(MW/MDSB(6))**2
     &          *(GAM(3,6)*GAM(2,6)*(-CEG*FN1+(2*CER-CEG)*FN2)
     $           -GAM(6,6)*GAM(2,6)*MGLU/MBQ*(-CEG*FN3+(2*CER-CEG)*FN4))
C
C ...  neutralinos-dsquark1 ...
C              
       C3NS1EFF=0.
       C2NS1EFF=0.
       DO K=1,4                    ! loop through neutralinos
         X=MCH0(K)**2/MDSB(1)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS1EFF=C3NS1EFF
     $            +FACN7*(MW/MDSB(1))**2
     $             *(2.*AGE(K,3,1)*AGE(K,2,1)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,1)-AHE(K,3,1))
     $               *AGE(K,2,1)*MCH0(K)/MBQ*FN4)  
         C2NS1EFF=C2NS1EFF
     $            +FACN8*(MW/MDSB(1))**2
     $             *(2.*AGE(K,3,1)*AGE(K,2,1)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,1)-AHE(K,3,1))
     $               *AGE(K,2,1)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark2 ...
C              
       C3NS2EFF=0.
       C2NS2EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(2)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS2EFF=C3NS2EFF
     $            +FACN7*(MW/MDSB(2))**2
     $             *(2.*AGE(K,3,2)*AGE(K,2,2)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,2)-AHE(K,3,2))
     $               *AGE(K,2,2)*MCH0(K)/MBQ*FN4)  
         C2NS2EFF=C2NS2EFF
     $            +FACN8*(MW/MDSB(2))**2
     $             *(2.*AGE(K,3,2)*AGE(K,2,2)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,2)-AHE(K,3,2))
     $               *AGE(K,2,2)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark3 ...
C              
       C3NS3EFF=0.
       C2NS3EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(3)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS3EFF=C3NS3EFF
     $            +FACN7*(MW/MDSB(3))**2
     $             *(2.*AGE(K,3,3)*AGE(K,2,3)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,3)-AHE(K,3,3))
     $               *AGE(K,2,3)*MCH0(K)/MBQ*FN4)  
         C2NS3EFF=C2NS3EFF
     $            +FACN8*(MW/MDSB(3))**2
     $             *(2.*AGE(K,3,3)*AGE(K,2,3)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,3)-AHE(K,3,3))
     $               *AGE(K,2,3)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark4 ...
C              
       C3NS4EFF=0.
       C2NS4EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(4)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS4EFF=C3NS4EFF
     $            +FACN7*(MW/MDSB(4))**2
     $             *(2.*AGE(K,3,4)*AGE(K,2,4)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,4)-AHE(K,3,4))
     $               *AGE(K,2,4)*MCH0(K)/MBQ*FN4)  
         C2NS4EFF=C2NS4EFF
     $            +FACN8*(MW/MDSB(4))**2
     $             *(2.*AGE(K,3,4)*AGE(K,2,4)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,4)-AHE(K,3,4))
     $               *AGE(K,2,4)*MCH0(K)/MBQ*FN4)  
       ENDDO
c
C ...  neutralinos-dsquark5 ...
C              
       C3NS5EFF=0.
       C2NS5EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(5)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS5EFF=C3NS5EFF
     $            +FACN7*(MW/MDSB(5))**2
     $             *(2.*AGE(K,3,5)*AGE(K,2,5)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,5)-AHE(K,3,5))
     $               *AGE(K,2,5)*MCH0(K)/MBQ*FN4)  
         C2NS5EFF=C2NS5EFF
     $            +FACN8*(MW/MDSB(5))**2
     $             *(2.*AGE(K,3,5)*AGE(K,2,5)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,5)-AHE(K,3,5))
     $               *AGE(K,2,5)*MCH0(K)/MBQ*FN4)  
       ENDDO
C
C ...  neutralinos-dsquark6 ...
C              
       C3NS6EFF=0.
       C2NS6EFF=0.
       DO K=1,4
         X=MCH0(K)**2/MDSB(6)**2
         CALL FUNS(X,FN1,FN2,FN3,FN4) 
         C3NS6EFF=C3NS6EFF
     $            +FACN7*(MW/MDSB(6))**2
     $             *(2.*AGE(K,3,6)*AGE(K,2,6)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,6)-AHE(K,3,6))
     $               *AGE(K,2,6)*MCH0(K)/MBQ*FN4)  
         C2NS6EFF=C2NS6EFF
     $            +FACN8*(MW/MDSB(6))**2
     $             *(2.*AGE(K,3,6)*AGE(K,2,6)*FN2
     $              -SQRT(2.)*(SQRT(2.)*AGE(K,6,6)-AHE(K,3,6))
     $               *AGE(K,2,6)*MCH0(K)/MBQ*FN4)  
       ENDDO
      ENDIF
C
      SINTW2I=1./SN2THW
      
      XTW=MT**2/MW**2
C
c ... here we switch from Anlauf's and Ciuchini's 
c     normalization to that of Greub:
C     Ci(Greub)=Ci(Ciuchini,Buras)                for i=1,...,6
c     C7(Greub)=-Qb*C7(Anlauf)=1/3*C7(Anlauf)
c     C8(Greub)=-C8(Anlauf)
c ...
c
      IF(Q.GT.MW.AND.QNEW.LE.MW) THEN 
       CI1(1)= ALES/4./PI*11./2.d0
       CI1(2)= 1.d0-ALES/4./PI*11./6.d0
       CI1(3)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
     $        +ALEL/6.d0/PI*SINTW2I*(2.*BI(XTW)+CI(XTW))
       CI1(4)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI1(5)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
       CI1(6)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI1(7)= 1./3.d0*(C3SMEFF+C3HPEFF
     $                  +C3C11EFF+C3C12EFF+C3C13EFF
     $                  +C3C21EFF+C3C22EFF+C3C23EFF
     $                  +C3GS1EFF+C3GS2EFF+C3GS3EFF
     $                  +C3GS4EFF+C3GS5EFF+C3GS6EFF
     $                  +C3NS1EFF+C3NS2EFF+C3NS3EFF
     $                  +C3NS4EFF+C3NS5EFF+C3NS6EFF)
       CI1(8)=-1.d0*(C2SMEFF+C2HPEFF
     $               +C2C11EFF+C2C12EFF+C2C13EFF
     $               +C2C21EFF+C2C22EFF+C2C23EFF
     $               +C2GS1EFF+C2GS2EFF+C2GS3EFF
     $               +C2GS4EFF+C2GS5EFF+C2GS6EFF
     $               +C2NS1EFF+C2NS2EFF+C2NS3EFF
     $               +C2NS4EFF+C2NS5EFF+C2NS6EFF)
c
       CI2(1)= ALES/4.d0/PI*11./2.d0
       CI2(2)= 1.d0-ALES/4.d0/PI*11./6.d0
       CI2(3)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
     $        +ALEL/6.d0/PI*SINTW2I*(2.*BI(XTW)+CI(XTW))
       CI2(4)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI2(5)=-ALES/24.d0/PI*(EI(XTW)-2./3.d0)
       CI2(6)= ALES/8.d0/PI*(EI(XTW)-2./3.d0)
       CI2(7)= 1./3.d0*(C3SMEFF+C3HPEFF
     $                  +C3C11EFF+C3C12EFF+C3C13EFF
     $                  +C3C21EFF+C3C22EFF+C3C23EFF
     $                  +C3GS1EFF+C3GS2EFF+C3GS3EFF
     $                  +C3GS4EFF+C3GS5EFF+C3GS6EFF
     $                  +C3NS1EFF+C3NS2EFF+C3NS3EFF
     $                  +C3NS4EFF+C3NS5EFF+C3NS6EFF)
       CI2(8)=-1.d0*(C2SMEFF+C2HPEFF
     $               +C2C11EFF+C2C12EFF+C2C13EFF
     $               +C2C21EFF+C2C22EFF+C2C23EFF
     $               +C2GS1EFF+C2GS2EFF+C2GS3EFF
     $               +C2GS4EFF+C2GS5EFF+C2GS6EFF
     $               +C2NS1EFF+C2NS2EFF+C2NS3EFF
     $               +C2NS4EFF+C2NS5EFF+C2NS6EFF)
c
        C3CHEFF=C3C11EFF+C3C12EFF+C3C13EFF
     $         +C3C21EFF+C3C22EFF+C3C23EFF
        C3GLEFF=C3GS1EFF+C3GS2EFF+C3GS3EFF
     $         +C3GS4EFF+C3GS5EFF+C3GS6EFF
        C3NEEFF=C3NS1EFF+C3NS2EFF+C3NS3EFF
     $         +C3NS4EFF+C3NS5EFF+C3NS6EFF
c
        C2CHEFF=C2C11EFF+C2C12EFF+C2C13EFF
     $         +C2C21EFF+C2C22EFF+C2C23EFF
        C2GLEFF=C2GS1EFF+C2GS2EFF+C2GS3EFF
     $         +C2GS4EFF+C2GS5EFF+C2GS6EFF
        C2NEEFF=C2NS1EFF+C2NS2EFF+C2NS3EFF
     $         +C2NS4EFF+C2NS5EFF+C2NS6EFF
      ENDIF
C
C           
      IF(Q.LT.MW) THEN
        T=DT*FLOAT(II-1)
        CALL DRKSTP(8,DT,T,CI1,GAMMAWB1,WCI1)
        T=DT*FLOAT(II-1)
        CALL DRKSTP(8,DT,T,CI2,GAMMAWB2,WCI2)
C
C------- Computation of b->s\gamma branching ratio --------------------
C
       IF(Q.GT.(2.d0*MB)) THEN
         XQ2 = SNGL(Q**2)
	 XMT = SNGL(MT)
         ASMB= DBLE(SUALFS(XQ2,.36,XMT,3))
         AEMB= DBLE(SUALFE(XQ2))
         ALFA2=AEMB/SN2THW
         G2=SQRT(4.d0*PI*ALFA2)
         GEF=SQRT(2.d0)/8.d0*G2**2/MW**2      
       ENDIF
       BQLOG=LOG(MB/Q)
       ZI=(MC/MB)**2
       LOGZI=LOG(ZI)
       RER2=2./243.d0
     &      *(-833.d0+144.*PI**2*ZI**(3./2.)
     $        +(1728.d0-180.*PI**2-1296.d0*1.20206
     $          +(1296.d0-324.*PI**2)*LOGZI+108.d0*LOGZI**2
     $          +36.d0*LOGZI**3)*ZI
     $        +(648.d0+72.d0*PI**2+(432.d0-216.d0*PI**2)*LOGZI
     $          +36.d0*LOGZI**3)*ZI**2
     $        +(-54.d0-84.d0*PI**2+1092.d0*LOGZI
     $          -756.d0*LOGZI**2)*ZI**3)
       IMR2=16.*PI/81.d0
     &      *(-5.d0+(45.d0-3.d0*PI**2+9.*LOGZI+9.*LOGZI**2)*ZI
     $        +(-3.d0*PI**2+9.d0*LOGZI**2)*ZI**2
     $        +(28.d0-12.d0*LOGZI)*ZI**3)
       RER7= 8./9.d0*(4.d0-PI**2)
       IMR7= 0.d0
       RER8=-4./27.d0*(-33.d0+2.*PI**2)
       IMR8=-4./27.d0*(-6.d0*PI)
       REDE=ASMB/4.d0/PI
     $      *((416./81.d0*BQLOG+RER2)*CI1(2)
     $        +(32./3.d0*BQLOG+RER7)*CI1(7)
     $        -(32./9.d0*BQLOG+RER8)*CI1(8))
       IMDE=ASMB/4.d0/PI
     $      *(IMR2*CI1(2)+IMR7*CI1(7)+IMR8*CI1(8))
       EF=(1.d0-8./3.d0*ASMB/PI)
c
C...Compute bremsstrahlung corrections
c
       CALL BREMS(AEMB,ASMB,MS,MC,MB,Q,CI1,GBREMS)
c
c...Implement virtual QCD corrections according to Eq(5.6) of Greub.
c
       GVIRT=AEMB/32.d0/PI**4*EF*(CI2(7)**2+2.d0*CI2(7)*REDE)
c     
       GAMMASG=GVIRT+GBREMS
       gammasg0=AEMB/32.d0/PI**4*EF*CI1(7)**2
c
c...Evaluate semileptonic decay width as in Eq(5.9) of Greub.
c
       GAMMASL=1./192.d0/PI**3*GES(MC/MB)
     $         *(1.d0-2./3.d0/PI*ASMB*FES(MC/MB))
c
c...Compute b->s\gamma branching ratio using Eq(5.8) of Greub.
c
       BRSG = COBA*GAMMASG/GAMMASL*BRSL
       BRSG0= COBA*GAMMASG0/GAMMASL*BRSL
c
C...Evaluate theoretical uncertainties from scale variation
c
       if(q.Gt.(2.*MB).AND.QNEW.LT.(2.*MB)) then
         brup=BRSG*1.d+4
       ENDIF
       if(q.Gt.(1./2.*MB).AND.QNEW.LT.(1./2.*MB)) then
         brdo=BRSG*1.d+4
       ENDIF
c
C...Computation of the final result
c
       IF(q.Gt.MB.AND.QNEW.LT.MB) THEN
         BFBSG =BRSG*1.d+4
         BFBSG0=BRSG0*1.d+4
c       write(6,*) 'q=',q,ci1(7),ci2(7),rede,gvirt,gbrems
c       write(6,*) 'ci=',ci1(2),ci1(7),ci1(8),rer2,rer7,rer8,asmb
c       write(6,*) 'brsg,brsg0=',brsg,brsg0
        ENDIF
      ENDIF
c---------------------------------------------------------------
400   CONTINUE 
c
        C3C1NAI=C3C11NAI+C3C12NAI+C3C13NAI
        C3C2NAI=C3C21NAI+C3C22NAI+C3C23NAI
        C3CHNAI=C3C11NAI+C3C12NAI+C3C13NAI
     $         +C3C21NAI+C3C22NAI+C3C23NAI
        C3GLNAI=C3GS1EFF+C3GS2EFF+C3GS3EFF
     $         +C3GS4EFF+C3GS5EFF+C3GS6EFF
        C3NENAI=C3NS1EFF+C3NS2EFF+C3NS3EFF
     $         +C3NS4EFF+C3NS5EFF+C3NS6EFF
C
        C3C1EFF=C3C11EFF+C3C12EFF+C3C13EFF
        C3C2EFF=C3C21EFF+C3C22EFF+C3C23EFF
        C3CHEFF=C3C11EFF+C3C12EFF+C3C13EFF
     $         +C3C21EFF+C3C22EFF+C3C23EFF
        C3GLEFF=C3GS1EFF+C3GS2EFF+C3GS3EFF
     $         +C3GS4EFF+C3GS5EFF+C3GS6EFF
        C3NEEFF=C3NS1EFF+C3NS2EFF+C3NS3EFF
     $         +C3NS4EFF+C3NS5EFF+C3NS6EFF
c
      IF (IWR.EQ.1) THEN
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'SM:          '
        WRITE(6,*)'QCD corrected=',c3smeff,' naive=',C3SMNAI                
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'higgs:       '
        WRITE(6,*)'QCD corrected=',c3hpeff,' naive=',C3HPNAI      
        WRITE(6,*)'charged Higgs mass=',mhplus
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(1)-stop(1): '
        WRITE(6,*)'QCD corrected=',c3c11eff,' naive=',C3C11NAI       
        WRITE(6,*)'chargino mass=',mcha(1),' squark mass=',msqu(1)  
        WRITE(6,*)'G=',gk(1,1),' H=',hk(1,1)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(1)-stop(2): '
        WRITE(6,*)'QCD corrected=',c3c12eff,' naive=',C3C12NAI        
        WRITE(6,*)'chargino mass=',mcha(1),' squark mass=',msqu(2)  
        WRITE(6,*)'G=',gk(1,2),' H=',hk(1,2)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(1)-up+charm squarks: '
        WRITE(6,*)'QCD corrected=',c3c13eff,' naive=',C3C13NAI       
        WRITE(6,*)'chargino mass=',mcha(1),' squark mass=',msqu(3)  
        WRITE(6,*)'G=',gk(1,3),' H=',hk(1,3)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(2)-stop(1): '
        WRITE(6,*)'QCD corrected=',c3c21eff,' naive=',C3C21NAI       
        WRITE(6,*)'chargino mass=',mcha(2),' squark mass=',msqu(1)  
        WRITE(6,*)'G=',gk(2,1),' H=',hk(2,1)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'chargino(2)-stop(2): '
        WRITE(6,*)'QCD corrected=',c3c22eff,' naive=',C3C22NAI        
        WRITE(6,*)'chargino mass=',mcha(2),' squark mass=',msqu(2)  
        WRITE(6,*)'G=',gk(2,2),' H=',hk(2,2)
	WRITE(6,*)'____________________________________________'
	WRITE(6,*)'chargino(2)-up+charm squarks: '
        WRITE(6,*)'QCD corrected=',c3c23eff,' naive=',C3C23NAI       
        WRITE(6,*)'chargino mass=',mcha(2),' squark mass=',msqu(3)  
        WRITE(6,*)'G=',gk(2,3),' H=',hk(2,3)
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'sum of chargino-squarks contributions: '
        WRITE(6,*)'chargino(1)'
        WRITE(6,*)'QCD corrected=',c3c1eff,' naive=',C3C1NAI   
        WRITE(6,*)'chargino(2)'
        WRITE(6,*)'QCD corrected=',c3c2eff,' naive=',C3C2NAI    
        WRITE(6,*)'total'
        WRITE(6,*)'QCD corrected=',c3cheff,' naive=',C3ChNAI    
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'sum of gluino-squarks contributions: '
        WRITE(6,*)'QCD corrected=',1./3.*c3gleff,' naive=',
     $              1./3.*C3glNAI    
        WRITE(6,*)'____________________________________________'
        WRITE(6,*)'sum of neutralino-squarks contributions: '
        WRITE(6,*)'QCD corrected=',c3neeff,' naive=',C3neNAI    
        WRITE(6,*)'____________________________________________'
        write(6,*) 'BFBSG,BFBSG0=',BFBSG,BFBSG0,'  x 10^-4'
      END IF
C--------------------------------------------------------------
1000  continue              
      
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE RGE157(T,G,F)
C----------------------------------------------------------------------
C
C  Right hand side of the full set of 157 renormalization 
c  group equations
C              dG_i/dT = F_i(G), i=1..157
c  Sparticles decouple at muliple scales according to Castano et al.
C  This is SURG157 subrotine of ISAJET-M with some minor modifications.
c
C  NOTE: This sibroutine follows convention in which fermion masses 
c  	 are proportional to Yukawas, i.e. m ~ Y, while the rest 
c  	 uses convention with fermion masses proportional to transposed
c  	 Yukawas in gauge eigenbasis, i.e.  m ~ Y^T,
C  	 One has to transpose Yukasas and SSB trilinears when passing
C  	 into or from this subroutine.
C
c  Ref.: Martin & Vaughn PRD 50, 2282 (1994);
c        Castano, Piard, Ramond PRD 49, 4882 (1994)
c        with errors in Eqs (B15) and B(16) fixed;
c        S.Antusch et al hep-ph/0501272;
C        Casas & Ibarra hep-ph/0103065;
c        Casas et al PRD 63,097302 (2001).
c
c  Author: Azar Mustafayev
c  Created: 05/29/07
c  Modified: 13/12/07 by Azar Mustafayev 
c             - to speedup matrix multiplications moved to subroutines
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C  Local variables:
C     YU  = Y_u - up Yukawa coupling matrix
c     YD  = Y_d - down Yukawa coupling matrix
C     YE  = Y_e - lepton Yukawa coupling matrix 
C     YN  = Y_nu- neutrino Yukawa coupling matrix
C     MRHN= M_RHN - RH neutrino mass matrix 
C     KA  = \kappa- effective dim-5 operator from integration out RH neutrinos 
C     TU  = h_u - up soft-breaking trilinear coupling matrix
C     TD  = h_d - down soft-breaking trilinear coupling matrix
c     TE  = h_e - lepton soft-breaking trilinear coupling matrix 
C     TN  = h_nu- neutrino soft-breaking trilinear coupling matrix 
C     M2Q  = m^2_Q  - squark doublet mass^2 matrix
C     M2L  = m^2_L  - slepton doublet mass^2 matrix
C     M2U  = m^2_u  - up-squark mass^2 matrix
C     M2D  = m^2_d  - down-squark mass^2 matrix
C     M2E  = m^2_e  - slepton singlet mass^2 matrix
C     M2N  = m^2_nu - sneutrino singlet mass^2 matrix
C     YUD = (Y_u)^dagger - Hermitian conjugated Y_u
C     YDD = (Y_d)^dagger 
C     YED = (Y_e)^dagger 
C     YND = (Y_nu)^dagger 
C     TUD = (h_u)^dagger 
C     TDD = (h_d)^dagger 
C     TED = (h_e)^dagger 
C     TND = (h_nu)^dagger 
C
C     G(1) = g_1
C     G(2) = g_2
C     G(3) = g_3
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
C     G(31) = M_1
C     G(32) = M_2
C     G(33) = M_3
C     G(34)-G(42) = h_u
C     G(43)-G(51) = h_d
C     G(52)-G(60) = h_e
C     G(61)= mu
C     G(62)= B*mu
C     G(63)= m^2_Hu
C     G(64)= m^2_Hd
C     G(65)-G(73) = m^2_Q
C     G(74)-G(82) = m^2_L
C     G(83)-G(91) = m^2_u
C     G(92)-G(100)= m^2_d
C     G(101)-G(109)=m^2_e
C     G(110)= v_u
C     G(111)= v_d
C     G(112)-G(120) = Y_nu  
C     G(121)-G(129) = M_RHN 
C     G(130)-G(138) = h_nu
C     G(139)-G(147) = m^2_nu
C     G(148)-G(156) = \kappa
C     G(157) = \lambda   - SM quartic higgs coupling
C
C     X - auxiliary 250x3x3 matrix
C     X(200)-(250) - used for neutrino sector
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      IMPLICIT NONE
      REAL*8 T,G(157),F(157)
      REAL*8 TRACE,TR3X3
c
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),IRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER IRHN
      SAVE /BSGDEC/
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
c        M1,M2,M3  - gaugino masses at MZ
c        ABOT,ATOP,ATAU  - soft trilinear scalar couplings at MZ 
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
c
      REAL*8 PI,FAC,S,X(250,3,3),F2(157),B1,B2,B3,
     &       BETA,SINB,COSB,Q,SP,SIG1,SIG2,SIG3,
     &       T1(3,3),T2(3,3),T3(3,3),T4(3,3),T5(3,3),T6(3,3),I3(3,3),
     &       YU(3,3),YD(3,3),YE(3,3),YN(3,3),TU(3,3),TD(3,3),TE(3,3),
     &       TN(3,3),M2Q(3,3),M2L(3,3),M2U(3,3),M2D(3,3),M2E(3,3),
     &       M2N(3,3),YUD(3,3),YDD(3,3),YED(3,3),YND(3,3),TUD(3,3),
     &       TDD(3,3),TED(3,3),TND(3,3),MRHN(3,3),KA(3,3)
      INTEGER STPSHU,STPHD,STPSHL,STPSHH,
     &        STPHL,STPHH,STPSB,STPSW,STPSG,
     &        STPSQ(3),STPSU(3),STPSD(3),STPSL(3),STPSE(3),ITWOLP
C   Step functions (=1 above and =0 below the scale):
C     STPSHU  = theta_sHu  - up higgsino
C     STPSHD  = theta_sHd  - down higgsino
C     STPSHL  = theta_sh   - light higgsino
C     STPSHH  = theta_sH   - heavy higgsino
C     STPHL   = theta_h    - light Higgs
C     STPHH   = theta_H    - heavy Higgs
C     STPSB   = theta_sB   - bino
C     STPSW   = theta_sW   - wino
C     STPSG   = theta_sg   - gluino
C     STPSQ(i)= theta_sQi  - i-th squark doublet
C     STPSU(i)= theta_sUi  - i-th up squark singlet
C     STPSD(i)= theta_sDi  - i-th down squark singlet
C     STPSL(i)= theta_sLi  - i-th slepton doublet
C     STPSE(i)= theta_sEi  - i-th slepton singlet
C     
C     ITWOLP  - switch for the second loop
c              =0 - 1-loop
c              =1 - only gauge and Yukawas at 2-loop
c              =2 - full 2-loop
      INTEGER NSQ,NSU,NSD,NSL,NSE,NSH,NH,NU,NE,ND,NN
      INTEGER I,J,K,L,M,N
      DATA I3/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix

      
      PI=4.d0*DATAN(1.d0)
      FAC=16.d0*PI**2

      ITWOLP=2    ! full 2-loop
      
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=SQRT(1.d0-SINB**2)
      Q=MGUT*DEXP(T)

c-----assemble matrices------------------------------------------------      
      K=0
      DO I=1,3
        DO J=1,3
	  YU(I,J)=G(4+K)
	  YD(I,J)=G(13+K)
	  YE(I,J)=G(22+K)
	  YN(I,J)=G(112+K)
	  MRHN(I,J)=G(121+K)
	  KA(I,J)=G(148+K)
	  TU(I,J)=G(34+K)
	  TD(I,J)=G(43+K)
	  TE(I,J)=G(52+K)
	  TN(I,J)=G(130+K)
	  M2Q(I,J)=G(65+K)
	  M2L(I,J)=G(74+K)
	  M2U(I,J)=G(83+K)
	  M2D(I,J)=G(92+K)
	  M2E(I,J)=G(101+K)
	  M2N(I,J)=G(139+K)
	  YUD(J,I)=YU(I,J)
	  YDD(J,I)=YD(I,J)
	  YED(J,I)=YE(I,J)
	  YND(J,I)=YN(I,J)
	  TUD(J,I)=TU(I,J)
	  TDD(J,I)=TD(I,J)
	  TED(J,I)=TE(I,J)
	  TND(J,I)=TN(I,J)
	  K=K+1
        ENDDO
      ENDDO
C----Reset auxilliary matrices
      DO K=1,250
        DO I=1,3
          DO J=1,3
            X(K,I,J)=0.d0
	  ENDDO
	ENDDO
      ENDDO
c---- CALCULATE 1-LOOP THRESHOLD EFFECTS ----------
      ND=3
      NE=3
      NN=3
      IF (Q.GE.MT) THEN
        NU=3
      ELSE
        NU=2
      END IF

      DO i=1,3
	IF (Q.GE.MSQDEC(i)) THEN
	  STPSQ(i)=1
        ELSE
	  STPSQ(i)=0
        ENDIF
	IF (Q.GE.MSLDEC(i)) THEN
	  STPSL(i)=1
        ELSE
	  STPSL(i)=0
        ENDIF
	IF (Q.GE.MSUDEC(i)) THEN
	  STPSU(i)=1
        ELSE
	  STPSU(i)=0
        ENDIF
	IF (Q.GE.MSDDEC(i)) THEN
	  STPSD(i)=1
        ELSE
	  STPSD(i)=0
        ENDIF
	IF (Q.GE.MSEDEC(i)) THEN
	  STPSE(i)=1
        ELSE
	  STPSE(i)=0
        ENDIF
      ENDDO
      NSU=STPSU(1)+STPSU(2)+STPSU(3)
      NSD=STPSD(1)+STPSD(2)+STPSD(3)
      NSQ=STPSQ(1)+STPSQ(2)+STPSQ(3)
      NSE=STPSE(1)+STPSE(2)+STPSE(3)
      NSL=STPSL(1)+STPSL(2)+STPSL(3)
      
      IF (Q.GE.DBLE(ABS(MU))) THEN
	STPSHL=1
	STPSHH=1
      ELSE 
	STPSHL=0
	STPSHH=0
      END IF
      NSH=STPSHL+STPSHH

      STPHL=1
      IF (Q.GT.MHPLUS) THEN
        STPHH=1
      ELSE
        STPHH=0
      END IF
      NH=STPHL+STPHH

      IF (Q.GE.MGLU) THEN
        STPSG=1
      ELSE
        STPSG=0
      END IF
	
      IF (Q.GE.M1) THEN
        STPSB=1
      ELSE
        STPSB=0
      END IF
      
      IF (Q.GE.M2) THEN
        STPSW=1
      ELSE
        STPSW=0
      END IF
	    
c-----evaluate auxiliary matrices-------------------------------------      
      CALL MPROD2X(X,1,YUD,YU)
      CALL MPROD2X(X,2,YDD,YD)
      CALL MPROD2X(X,3,YED,YE)
      CALL MPROD2X(X,11,YUD,TU)
      CALL MPROD2X(X,16,YDD,TD)
      CALL MPROD2X(X,17,YED,TE)
      CALL MPROD2X(X,24,TUD,TU)
      CALL MPROD2X(X,29,TDD,TD)
      CALL MPROD2X(X,30,TED,TE)
      CALL MPROD2X(X,37,TU,TUD)
      CALL MPROD2X(X,41,TD,TDD)
      CALL MPROD2X(X,45,TE,TED)
      CALL MPROD2X(X,46,YU,YUD)
      CALL MPROD2X(X,47,YD,YDD)
      CALL MPROD2X(X,48,YE,YED)
      CALL MPROD3X(X,4,YU,YUD,YU)
      CALL MPROD3X(X,5,YU,YDD,YD)
      CALL MPROD3X(X,6,YD,YDD,YD)
      CALL MPROD3X(X,7,YD,YUD,YU)
      CALL MPROD3X(X,8,YE,YED,YE)
      CALL MPROD3X(X,9,TU,YUD,YU)
      CALL MPROD3X(X,10,TU,YDD,YD)
      CALL MPROD3X(X,12,YU,YUD,TU)
      CALL MPROD3X(X,13,YU,YDD,TD)
      CALL MPROD3X(X,14,TD,YDD,YD)
      CALL MPROD3X(X,15,TD,YUD,YU)
      CALL MPROD3X(X,18,YD,YDD,TD)
      CALL MPROD3X(X,19,YD,YUD,TU)
      CALL MPROD3X(X,20,TE,YED,YE)
      CALL MPROD3X(X,21,YE,YED,TE)
      CALL MPROD3X(X,22,M2Q,YUD,YU)
      CALL MPROD3X(X,23,YUD,M2U,YU)
      CALL MPROD3X(X,25,M2Q,YDD,YD)
      CALL MPROD3X(X,26,YDD,M2D,YD)
      CALL MPROD3X(X,27,M2L,YED,YE)
      CALL MPROD3X(X,28,YED,M2E,YE)
      CALL MPROD3X(X,31,YUD,YU,M2Q)
      CALL MPROD3X(X,32,YDD,YD,M2Q)
      CALL MPROD3X(X,33,YED,YE,M2L)
      CALL MPROD3X(X,34,M2U,YU,YUD)
      CALL MPROD3X(X,35,YU,M2Q,YUD)
      CALL MPROD3X(X,36,YU,YUD,M2U)
      CALL MPROD3X(X,38,M2D,YD,YDD)
      CALL MPROD3X(X,39,YD,M2Q,YDD)
      CALL MPROD3X(X,40,YD,YDD,M2D)
      CALL MPROD3X(X,42,M2E,YE,YED)
      CALL MPROD3X(X,43,YE,M2L,YED)
      CALL MPROD3X(X,44,YE,YED,M2E)
C
C  1-loop part
C

c...gauge couplings
      B1=2.d0*(17.*NU/12.d0+5.*ND/12.d0+5.*NE/4.d0+1.*NN/4.d0)/5.d0
     &   +1.*NSQ/30.d0+4.*NSU/15.d0+1.*NSD/15.d0+1.*NSL/10.d0
     &   +1.*NSE/5.d0+1.*NSH/5.d0+1.*NH/10.d0
      B2=-22./3.d0+0.5d0*(NU+ND)+1.d0*(NE+NN)/6.d0
     $   +1.*NSQ/2.d0+1.*NSL/6.d0+1.*NSH/3.d0+1.*NH/6.d0+4.*STPSW/3.d0
      B3=2.*(NU+ND)/3.d0+1.*NSQ/3.d0+1.*NSU/6.d0+1.*NSD/6.d0
     $   +2.d0*STPSG-11.d0

      F(1)=G(1)**3*B1/FAC
      F(2)=G(2)**3*B2/FAC
      F(3)=G(3)**3*B3/FAC
c...Yukawa couplings
      DO 230 I=1,3
        DO 230 J=1,3
	  T1(I,J)=0.d0
	  T2(I,J)=0.d0
	  T3(I,J)=0.d0
	  T4(I,J)=0.d0
	  T5(I,J)=0.d0
	  T6(I,J)=0.d0
	  DO 230 K=1,3
	    T1(I,J)=T1(I,J)+YU(I,K)*STPSQ(K)*YUD(K,J)
	    T2(I,J)=T2(I,J)+YUD(I,K)*STPSU(K)*YU(K,J)
	    T3(I,J)=T3(I,J)+YDD(I,K)*STPSD(K)*YD(K,J)
	    T4(I,J)=T4(I,J)+YD(I,K)*STPSQ(K)*YDD(K,J)
	    T5(I,J)=T5(I,J)+YE(I,K)*STPSL(K)*YED(K,J)
	    T6(I,J)=T6(I,J)+YED(I,K)*STPSE(K)*YE(K,J)
230   CONTINUE	
      CALL MPROD2X(X,180,T1,YU)
      CALL MPROD2X(X,181,YU,T2)
      CALL MPROD2X(X,182,YU,T3)
      CALL MPROD2X(X,183,T4,YD)
      CALL MPROD2X(X,184,YD,T3)
      CALL MPROD2X(X,185,YD,T2)
      CALL MPROD2X(X,186,T5,YE)
      CALL MPROD2X(X,187,YE,T6)
      CALL MPROD3X(X,188,YU,YUD,YU)
      CALL MPROD3X(X,189,YU,YDD,YD)
      CALL MPROD3X(X,190,YD,YUD,YU)

      K=0
      DO 250 I=1,3
        DO 250 J=1,3
      F(4+K)=(YU(I,J)
     &        *(-3./5.d0*G(1)**2
     &  	 *(17./12.d0+3./4.d0*STPSHL
     &  	   -(1./36.d0*STPSQ(J)+4./9.d0*STPSU(I)
     &  	     +1./4.d0*STPSHL)*STPSB)
     &  	-G(2)**2*(9./4.d0+9./4.d0*STPSHL
     &  		  -3./4.d0*(STPSQ(J)+STPSHL)*STPSW)
     &  	-G(3)**2*(8.d0-4./3.d0*(STPSQ(J)+STPSU(I))*STPSG)
     &  	       +((SINB**2*STPHL+COSB**2)*3.d0*TRACE(X,1)
     &  	  +COSB**2*(STPHL-1)*3.d0*TRACE(X,2)
     &  	  +COSB**2*(STPHL-1)*TRACE(X,3)))
     &        +3./2.d0*(SINB**2*STPHL+COSB**2*STPHH)*X(188,I,J)
     &        +1./2.d0*(SINB**2*STPSHL+COSB**2*STPSHH)
     &         *(2.d0*X(180,I,J)+X(181,I,J))
     &        +1./2.d0*(COSB**2*(STPHL-4*(STPHL-STPHH))+SINB**2*STPHH)
     &         *X(189,I,J)
     &        +1./2.d0*(COSB**2*STPSHL+SINB**2*STPSHH)*X(182,I,J)
     &        )/FAC
      F(13+K)=(YD(I,J)
     &         *(-3./5.d0*G(1)**2
     &  	  *(5./12.d0+3./4.d0*STPSHL
     &  	    -(1./36.d0*STPSQ(J)+1./9.d0*STPSD(I)
     &  	      +1./4.d0*STPSHL)*STPSB)
     &  	 -G(2)**2*(9./4.d0+9./4.d0*STPSHL
     &  		   -3./4.d0*(STPSQ(J)+STPSHL)*STPSW)
     &  	 -G(3)**2*(8.d0-4./3.d0*(STPSQ(J)+STPSD(I))*STPSG)
     &  	 +(SINB**2*(STPHL-1)*3.d0*TRACE(X,1)
     &  	  +(COSB**2*STPHL+SINB**2)*3.d0*TRACE(X,2)
     &  	  +(COSB**2*STPHL+SINB**2)*TRACE(X,3)))
     &         +3./2.d0*(COSB**2*STPHL+SINB**2*STPHH)*X(6,I,J)
     &         +1./2.d0*(COSB**2*STPSHL+SINB**2*STPSHH)
     &         *(2.d0*X(183,I,J)+X(184,I,J))
     &         +1./2.d0*(SINB**2*(STPHL-4*(STPHL-STPHH))+COSB**2*STPHH)
     &  	*X(190,I,J)
     &         +1./2.d0*(SINB**2*STPSHL+COSB**2*STPSHH)*X(185,I,J)
     &         )/FAC
      F(22+K)=(YE(I,J)
     &         *(-3./5.d0*G(1)**2
     &  	  *(15./4.d0+3./4.d0*STPSHL
     &  	    -(1./4.d0*STPSL(J)+1.d0*STPSE(I)
     &  	      +1./4.d0*STPSHL)*STPSB)
     &  	 -G(2)**2*(9./4.d0+9./4.d0*STPSHL
     &  		   -3./4.d0*(STPSL(J)+STPSHL)*STPSW)
     &  	 +(SINB**2*(STPHL-1)*3.d0*TRACE(X,1)
     &  	  +(COSB**2*STPHL+SINB**2)*3.d0*TRACE(X,2)
     &  	  +(COSB**2*STPHL+SINB**2)*TRACE(X,3)))
     &  	+3./2.d0*(COSB**2*STPHL+SINB**2*STPHH)*X(8,I,J)
     &  	+1./2.d0*(COSB**2*STPSHL+SINB**2*STPSHH)
     &  	*(2.d0*X(186,I,J)+X(187,I,J))
     &         )/FAC
          K=K+1
250   CONTINUE
c...compute convenient quantity
      S=G(63)-G(64)+TR3X3(M2Q)-TR3X3(M2L)-2.d0*TR3X3(M2U)
     &  +TR3X3(M2D)+TR3X3(M2E)
      
c...gaugino masses      
      F(31)=2.d0*B1*G(1)**2*G(31)/FAC
      F(32)=2.d0*B2*G(2)**2*G(32)/FAC
      F(33)=2.d0*B3*G(3)**2*G(33)/FAC
c...higgs mixing parameters      
      F(61)=G(61)*(3.d0*TRACE(X,1)+3.d0*TRACE(X,2)+TRACE(X,3)
     &             -3.d0*G(2)**2-3./5.d0*G(1)**2)/FAC
      F(62)=(G(62)*(3.d0*TRACE(X,1)+3.d0*TRACE(X,2)+TRACE(X,3)
     &              -3.d0*G(2)**2-3./5.d0*G(1)**2)
     &       +G(61)*(6.d0*TRACE(X,11)+6.d0*TRACE(X,16)+2.d0*TRACE(X,17)
     &               +6.d0*G(2)**2*G(32)+6./5.d0*G(1)**2*G(31)))/FAC
c...higgs mass^2
      F(63)=(6.d0*(G(63)*TRACE(X,1)+TRACE(X,31)+TRACE(X,34)+TRACE(X,24))
     &       -6.d0*G(2)**2*G(32)**2-6./5.d0*G(1)**2*G(31)**2
     &       +3./5.d0*S*G(1)**2)/FAC
      F(64)=(6.d0*(G(64)*TRACE(X,2)+TRACE(X,25))+6.d0*TRACE(X,26)
     &       +2.d0*(G(64)*TRACE(X,3)+TRACE(X,27))+2.d0*TRACE(X,28)
     &       +6.d0*TRACE(X,29)+2.d0*TRACE(X,30)
     &       -6.d0*G(2)**2*G(32)**2-6./5.d0*G(1)**2*G(31)**2
     &       -3./5.d0*S*G(1)**2)/FAC

      K=0
      DO 300 I=1,3
        DO 300 J=1,3
c...soft trilinear scalar couplings
      F(34+K)=(TU(I,J)*3.d0*TRACE(X,1)+5.d0*X(9,I,J)+X(10,I,J)
     &         -TU(I,J)*(16./3.d0*G(3)**2+3.*G(2)**2+13./15.d0*G(1)**2)
     &         +YU(I,J)*6.d0*TRACE(X,11)+4.d0*X(12,I,J)+2.d0*X(13,I,J)
     &         +YU(I,J)*(32./3.d0*G(3)**2*G(33)+6.*G(2)**2*G(32)
     &                   +26./15.d0*G(1)**2*G(31)))/FAC
      F(43+K)=(TD(I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         +5.d0*X(14,I,J)+X(15,I,J)
     &         -TD(I,J)*(16./3.d0*G(3)**2+3.*G(2)**2+7./15.d0*G(1)**2)
     &         +YD(I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))
     &         +4.d0*X(18,I,J)+2.d0*X(19,I,J)
     &         +YD(I,J)*(32./3.d0*G(3)**2*G(33)+6.d0*G(2)**2*G(32)
     &                   +14./15.d0*G(1)**2*G(31)))/FAC
      F(52+K)=(TE(I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))+5.d0*X(20,I,J)
     &         -TE(I,J)*(3.d0*G(2)**2+9./5.d0*G(1)**2)
     &         +YE(I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))+4.*X(21,I,J)
     &         +YE(I,J)*(6.*G(2)**2*G(32)+18./5.d0*G(1)**2*G(31)))/FAC
c...the rest of SSB mass^2
      F(65+K)=(X(22,I,J)+2.*G(63)*X(1,I,J)+X(25,I,J)+2.*G(64)*X(2,I,J)
     &         +X(31,I,J)+X(32,I,J)+2.d0*X(23,I,J)+2.d0*X(26,I,J)
     &         +2.d0*X(24,I,J)+2.d0*X(29,I,J)
     &         +I3(I,J)*(-32./3.d0*G(3)**2*G(33)**2-6.*G(2)**2*G(32)**2
     &                   -2./15.d0*G(1)**2*G(31)**2
     &                   +G(1)**2*S/5.d0))/FAC
      F(74+K)=(X(27,I,J)+2.*G(64)*X(3,I,J)+2.*X(28,I,J)+X(33,I,J)
     &         +2.d0*X(30,I,J)
     &         +I3(I,J)*(-6.*G(2)**2*G(32)**2-6./5.d0*G(1)**2*G(31)**2
     &                   -3./5.d0*G(1)**2*S))/FAC
      F(83+K)=(2.*X(34,I,J)+4.*G(63)*X(46,I,J)+4.*X(35,I,J)
     &         +2.d0*X(36,I,J)+4.d0*X(37,I,J)
     &         +I3(I,J)*(-32./3.d0*G(3)**2*G(33)**2
     &                   -32./15.d0*G(1)**2*G(31)**2
     &                   -4./5.d0*G(1)**2*S))/FAC
      F(92+K)=(2.*X(38,I,J)+4.*G(64)*X(47,I,J)+4.*X(39,I,J)
     &         +2.*X(40,I,J)+4.*X(41,I,J)
     &         +I3(I,J)*(-32./3.d0*G(3)**2*G(33)**2
     &                   -8./15.d0*G(1)**2*G(31)**2
     &                   +2./5.*G(1)**2*S))/FAC
      F(101+K)=(2.*X(42,I,J)+4.*G(64)*X(48,I,J)+4.*X(43,I,J)
     &          +2.*X(44,I,J)+4.*X(45,I,J)
     &          +I3(I,J)*(-24./5.d0*G(1)**2*G(31)**2
     &                    +6./5.d0*G(1)**2*S))/FAC
	  K=K+1
300   CONTINUE

c...VEVs
      F(110)=-G(110)*(3.*TRACE(X,1)-3./4.d0*(G(2)**2+G(1)**2/5.d0))/FAC
      F(111)=-G(111)*(3.*TRACE(X,2)+TRACE(X,3)
     &  	      -3./4.d0*(G(2)**2+G(1)**2/5.d0))/FAC

c--------- neutrino sector ------------------------------
      IF (IRHN.GT.0) THEN
        IF(Q.LT.MIN(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
	  DO I=112,147
	    F(I)=0.d0     ! all Majorana RHN are decoupled
	  ENDDO
        ELSE
	  CALL MPROD2X(X,200,YND,YN)
	  CALL MPROD2X(X,201,YN,YND)
	  CALL MPROD2X(X,215,YND,TN)
	  CALL MPROD2X(X,225,TN,TND)
	  CALL MPROD2X(X,229,TND,TN)
	  CALL MPROD3X(X,202,YE,YND,YN)
	  CALL MPROD3X(X,203,YN,YND,YN)
	  CALL MPROD3X(X,204,YN,YED,YE)
	  CALL MPROD3X(X,205,YN,YND,MRHN)
	  CALL MPROD3X(X,216,YN,YND,TN)
	  CALL MPROD3X(X,217,TN,YND,YN)
	  CALL MPROD3X(X,218,YN,YED,TE)
	  CALL MPROD3X(X,219,TN,YED,YE)
	  CALL MPROD3X(X,220,YE,YND,TN)
	  CALL MPROD3X(X,221,TE,YND,YN)
	  CALL MPROD3X(X,222,M2N,YN,YND)
	  CALL MPROD3X(X,223,YN,YND,M2N)
	  CALL MPROD3X(X,224,YN,M2L,YND)
	  CALL MPROD3X(X,226,M2L,YND,YN)
	  CALL MPROD3X(X,227,YND,YN,M2L)
	  CALL MPROD3X(X,228,YND,M2N,YN)
	  CALL MPROD3X(X,230,YND,M2L,YN)
	  
	  K=0
          DO 320 I=1,3
            DO 320 J=1,3
c...contributions to up-quark and charged lepton Yukawas
      F(4+K)=F(4+K) + YU(I,J)*TRACE(X,200)/FAC
      F(22+K)=F(22+K) + X(202,I,J)/FAC
c...neutrino Yukawa matrix
      F(112+K)=(3.d0*X(203,I,J)+X(204,I,J)
     &          +YN(I,J)*(TRACE(X,200)+3.d0*TRACE(X,1)
     &                    -3./5.d0*G(1)**2-3.d0*G(2)**2))/FAC
c...Majorana mass matrix
      F(121+K)=2.d0*(X(205,I,J)+X(205,J,I))/FAC
c...neutrino soft trilinear scalar couplings
      F(130+K)=(TN(I,J)*(3.d0*TRACE(X,1)+TRACE(X,200)
     &                   -3./5.d0*G(1)**2-3.d0*G(2)**2)
     &         +2.d0*YN(I,J)*(3.d0*TRACE(X,11)+TRACE(X,215)
     &                       -3./5.d0*G(1)**2*G(31)-3.d0*G(2)**2*G(32))
     &         +4.d0*X(216,I,J)+5.d0*X(217,I,J)+2.d0*X(218,I,J)
     &         +X(219,I,J))/FAC
c...contributions to the other soft trilinear scalar couplings
      F(34+K)=F(34+K)
     &        +(TU(I,J)*TRACE(X,200)+2.d0*YU(I,J)*TRACE(X,215))/FAC
      F(52+K)=F(52+K)+(2.d0*X(220,I,J)+X(221,I,J))/FAC
c...neutrino mass^2
      F(139+K)=(2.d0*(X(222,I,J)+X(223,I,J))
     &          +4.d0*(X(224,I,J)+G(63)*X(201,I,J)+X(225,I,J)))/FAC
c...contributions to slepton doublet mass^2
      F(74+K)= F(74+K)
     &         +(X(226,I,J)+X(227,I,J)
     &           +2.d0*(X(228,I,J)+G(63)*X(200,I,J)+X(229,I,J)))/FAC
	      K=K+1
320       CONTINUE        
C...contribution to up-higgs mass^2
        F(63)=F(63)+(2.d0*(TRACE(X,230)+TRACE(X,228)+G(63)*TRACE(X,200)
     &                     +TRACE(X,229)))/FAC
        ENDIF
c...effective dim-5 neutrino operator
        IF(LND5ON.AND.
     &     Q.LE.MAX(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
          CALL MPROD3X(X,231,KA,YED,YE)
          CALL MPROD3X(X,232,KA,YND,YN)
          CALL MPROD4X(X,233,YED,YE,YED,YE)
          CALL MPROD4X(X,234,YND,YN,YND,YN)
          CALL MPROD4X(X,235,YDD,YD,YDD,YD)
          CALL MPROD4X(X,236,YUD,YU,YUD,YU)

	  IF (Q.LT.MSUSY) THEN     !  use SM RGEs
	    K=0
            DO I=1,3
              DO J=1,3
      F(148+K)=(-3./2.d0*(X(231,J,I)+X(231,I,J))*COSB**2
     &          +1./2.d0*(X(232,J,I)+X(232,I,J))*SINB**2
     &          +KA(I,J)*(2.d0*TRACE(X,3)*COSB**2
     &                    +2.d0*TRACE(X,200)*SINB**2
     &                    +6.d0*TRACE(X,1)*SINB**2
     &                    +6.d0*TRACE(X,3)*COSB**2
     &                    -3.d0*G(2)**2+G(157)))/FAC
	        K=K+1
	      ENDDO
	    ENDDO
c...SM quartic higgs coupling
      F(157)=(6.d0*G(157)**2-3.d0*G(157)*(3.d0*G(2)**2+3./5.d0*G(1)**2)
     &        +3.d0*G(2)**4+3./2.d0*(3./5.d0*G(1)**2+G(2)**2)**2
     &        +4.d0*G(157)*(TRACE(X,3)+TRACE(X,200)+3.d0*TRACE(X,2)
     &                      +3.d0*TRACE(X,1))
     &        -8.d0*(TRACE(X,233)+TRACE(X,234)+3.d0*TRACE(X,235)
     &               +3.d0*TRACE(X,236)))/FAC
	  ELSE                     !  use MSSM RGEs
	    K=0
            DO I=1,3
              DO J=1,3
      F(148+K)=(X(231,J,I)+X(231,I,J)+X(232,J,I)+X(232,I,J)
     &         +KA(I,J)*(2.d0*TRACE(X,200)+6.d0*TRACE(X,1)
     &                   -6./5.d0*G(1)**2-6.d0*G(2)**2))/FAC
	        K=K+1
	      ENDDO
	    ENDDO
	    F(157)=0.d0
	  ENDIF
	ELSE
	  DO I=148,157
	    F(I)=0.d0
	  ENDDO
	ENDIF
      ENDIF
C
C  2-loop part for gauge and Yukawa couplings
C
      IF (ITWOLP.GE.1) THEN
c...evaluate auxiliary matrices
        CALL MPROD4X(X,49,YU,YUD,YU,YUD)
        CALL MPROD4X(X,50,YU,YDD,YD,YUD)
        CALL MPROD4X(X,54,YD,YDD,YD,YDD)
        CALL MPROD4X(X,55,YE,YED,YE,YED)
        CALL MPROD5X(X,51,YU,YUD,YU,YUD,YU)
        CALL MPROD5X(X,52,YU,YDD,YD,YDD,YD)
        CALL MPROD5X(X,53,YU,YDD,YD,YUD,YU)
        CALL MPROD5X(X,56,YD,YDD,YD,YDD,YD)
        CALL MPROD5X(X,57,YD,YUD,YU,YUD,YU)
        CALL MPROD5X(X,58,YD,YUD,YU,YDD,YD)
        CALL MPROD5X(X,59,YE,YED,YE,YED,YE)

        IF (Q.LT.MSUSY) THEN	  ! use SM RGEs below M_SUSY
c...Gauge couplings
      F2(1)=G(1)**3
     &      *(199./50.d0*G(1)**2+27./10.d0*G(2)**2+44./5.d0*G(3)**2
     &        -17./10.d0*SINB**2*TRACE(X,1)-1./2.d0*COSB**2*TRACE(X,2)
     &        -3./2.d0*COSB**2*TRACE(X,3))
      F2(2)=G(2)**3
     &      *(9./10.d0*G(1)**2+35./6.d0*G(2)**2+12.d0*G(3)**2
     &        -3./2.d0*SINB**2*TRACE(X,1)-3./2.d0*COSB**2*TRACE(X,2)
     &        -1./2.d0*COSB**2*TRACE(X,3))
      F2(3)=G(3)**3
     &      *(11./10.d0*G(1)**2+9./2.d0*G(2)**2-26.d0*G(3)**2
     &        -2.d0*SINB**2*TRACE(X,1)-2.d0*COSB**2*TRACE(X,2))
c...Yukawa couplings
          K=0
	  DO I=1,3
            DO J=1,3
              F2(4+K)=0.d0	 ! FOR NOW a la standard ISAJET
              F2(13+K)=0.d0
	      F2(22+K)=0.d0
	      K=K+1
	    ENDDO
	  ENDDO	
        ELSE                       ! use MSSM RGEs above M_SUSY
c...Gauge couplings
      F2(1)=G(1)**3
     &      *(199./25.d0*G(1)**2+27./5.d0*G(2)**2+88./5.d0*G(3)**2
     &        -26./5.d0*TRACE(X,1)-14./5.d0*TRACE(X,2)
     &        -18./5.d0*TRACE(X,3))
      F2(2)=G(2)**3
     &      *(9./5.d0*G(1)**2+25.d0*G(2)**2+24.d0*G(3)**2
     &        -6.d0*TRACE(X,1)-6.d0*TRACE(X,2)-2.d0*TRACE(X,3))
      F2(3)=G(3)**3
     &      *(11./5.d0*G(1)**2+9.d0*G(2)**2+14.d0*G(3)**2
     &        -4.d0*TRACE(X,1)-4.d0*TRACE(X,2))
c...Yukawa couplings
          K=0
	  DO 325 I=1,3
            DO 325 J=1,3
      F2(4+K)=YU(I,J)*(-3.d0)*(3.d0*TRACE(X,49)+TRACE(X,50))
     &        -X(5,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &        -9.d0*X(4,I,J)*TRACE(X,1)
     &        -4.d0*X(51,I,J)-2.d0*X(52,I,J)-2.d0*X(53,I,J)
     &        +YU(I,J)*(16.d0*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &        +(6.d0*G(2)**2+2./5.d0*G(1)**2)*X(4,I,J)
     &        +2./5.d0*G(1)**2*X(5,I,J)
     &        +YU(I,J)*(-16./9.d0*G(3)**4+8.d0*G(3)**2*G(2)**2
     &          	+136./45.d0*G(3)**2*G(1)**2
     &          	+15./2.d0*G(2)**4+G(2)**2*G(1)**2
     &          	+2743./450.d0*G(1)**4)
      F2(13+K)=YD(I,J)*(-3.d0)*(3.d0*TRACE(X,54)+TRACE(X,50)
     &          		+TRACE(X,55))
     &         -3.d0*X(7,I,J)*TRACE(X,1)
     &         -3.d0*X(6,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &         -4.d0*X(56,I,J)-2.d0*X(57,I,J)-2.d0*X(58,I,J)
     &         +YD(I,J)*((16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +4./5.d0*G(1)**2*X(7,I,J)
     &         +(6.d0*G(2)**2+4./5.d0*G(1)**2)*X(6,I,J)
     &         +YD(I,J)*(-16./9.d0*G(3)**4+8.d0*G(3)**2*G(2)**2
     &          	 +8./9.d0*G(3)**2*G(1)**2
     &          	 +15./2.d0*G(2)**4+G(2)**2*G(1)**2
     &          	 +287./90.d0*G(1)**4)
      F2(22+K)=YE(I,J)*(-3.d0)*(3.d0*TRACE(X,54)+TRACE(X,50)
     &          		+TRACE(X,55))
     &         -3.d0*X(8,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &         -4.d0*X(59,I,J)
     &         +YE(I,J)*((16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +6.d0*G(2)**2*X(8,I,J)
     &         +YE(I,J)*(15./2.d0*G(2)**4+9./5.d0*G(2)**2*G(1)**2
     &  		 +27./2.d0*G(1)**4)
	      K=K+1
325       continue
	ENDIF
c--------- neutrino sector ------------------------------
        IF(IRHN.GT.0.AND.
     &     Q.GE.MIN(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
          CALL MPROD4X(X,209,YN,YED,YE,YND)
          CALL MPROD4X(X,210,YN,YND,YN,YND)
          CALL MPROD5X(X,206,YN,YED,YE,YED,YE)
          CALL MPROD5X(X,207,YN,YED,YE,YND,YN)
          CALL MPROD5X(X,208,YN,YND,YN,YND,YN)
          CALL MPROD5X(X,211,YE,YND,YN,YED,YE)
          CALL MPROD5X(X,212,YE,YND,YN,YND,YN)
          CALL MPROD5X(X,213,YN,YED,YE,YND,MRHN)
          CALL MPROD5X(X,214,YN,YND,YN,YND,MRHN)

c...contributions to gauge couplings
          F2(1)=F2(1)-6./5.d0*G(1)**3*TRACE(X,200)
          F2(2)=F2(2)-2.d0*G(2)**3*TRACE(X,200)
      
	  K=0
          DO 330 I=1,3
            DO 330 J=1,3
c...neutrino Yukawa matrix
      F2(112+K)=-2.d0*X(206,I,J)-2.d0*X(207,I,J)-4.d0*X(208,I,J)
     &          -X(204,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &          -X(203,I,J)*(3.d0*TRACE(X,201)+9.d0*TRACE(X,1))
     &          +6./5.d0*G(1)**2*(X(204,I,J)+X(203,I,J))
     &          +6.d0*G(2)**2*X(203,I,J)
     &          -YN(I,J)*(TRACE(X,209)+3.d0*TRACE(X,210)
     &                    +3.d0*TRACE(X,50)+9.d0*TRACE(X,49)
     &                    +(4./5.d0*G(1)**2+16.d0*G(3)**2)*TRACE(X,1)
     &                    +207./50.d0*G(1)**4+9./5.d0*G(1)**2*G(2)**2
     &                    +15./2.d0*G(2)**4)
c...contributions to the other Yukawas
      F2(4+K)=F2(4+K)-3.d0*X(4,I,J)*TRACE(X,201)
     &        -YU(I,J)*(TRACE(X,209)+3.d0*TRACE(X,210))
      F2(13+K)=F2(13+K)-YD(I,J)*TRACE(X,209)-X(7,I,J)*TRACE(X,201)
      F2(22+K)=F2(22+K)-2.d0*X(211,I,J)-2.d0*X(212,I,J)
     &         -X(202,I,J)*(TRACE(X,201)+3.d0*TRACE(X,1))
     &         -YE(I,J)*TRACE(X,209)
	      K=K+1
330       continue        
	  DO I=112,120
	    F(I)=F(I)+F2(I)/FAC**2
	  ENDDO
	ENDIF
	DO I=1,30
	  F(I)=F(I)+F2(I)/FAC**2
	ENDDO
      ENDIF
C
C  2-loop part for SSB terms
C
      IF (ITWOLP.GE.2) THEN
c...evaluate auxiliary matrices
        CALL MPROD2X(X,112,TUD,YU)
        CALL MPROD2X(X,113,TDD,YD)
        CALL MPROD2X(X,114,TED,YE)
        CALL MPROD2X(X,145,TU,YUD)
        CALL MPROD2X(X,146,TD,YDD)
        CALL MPROD2X(X,147,TE,YED)
        CALL MPROD2X(X,148,YU,TUD)
        CALL MPROD2X(X,149,YD,TDD)
        CALL MPROD2X(X,150,YE,TED)
        CALL MPROD4X(X,63,TU,YUD,YU,YUD)
        CALL MPROD4X(X,64,TU,YDD,YD,YUD)
        CALL MPROD4X(X,65,TD,YUD,YU,YDD)
        CALL MPROD4X(X,75,TD,YDD,YD,YDD)
        CALL MPROD4X(X,76,TE,YED,YE,YED)
        CALL MPROD4X(X,92,TUD,TU,YUD,YU)
        CALL MPROD4X(X,93,TUD,YU,YUD,TU)
        CALL MPROD4X(X,94,TDD,TD,YUD,YU)
        CALL MPROD4X(X,95,YDD,YD,TUD,TU)
        CALL MPROD4X(X,96,TDD,YD,YUD,TU)
        CALL MPROD4X(X,97,YDD,TD,TUD,YU)
        CALL MPROD4X(X,102,TDD,TD,YDD,YD)
        CALL MPROD4X(X,103,TDD,YD,YDD,TD)
        CALL MPROD4X(X,104,TED,TE,YED,YE)
        CALL MPROD4X(X,105,TED,YE,YED,TE)
        CALL MPROD4X(X,115,YUD,YU,TUD,TU)
        CALL MPROD4X(X,116,YUD,TU,TUD,YU)
        CALL MPROD4X(X,117,YDD,YD,TDD,TD)
        CALL MPROD4X(X,118,YDD,TD,TDD,YD)
        CALL MPROD4X(X,122,YED,YE,TED,TE)
        CALL MPROD4X(X,123,YED,TE,TED,YE)
        CALL MPROD4X(X,125,YUD,YU,YUD,YU)
        CALL MPROD4X(X,126,YDD,YD,YDD,YD)
        CALL MPROD4X(X,127,YED,YE,YED,YE)
        CALL MPROD4X(X,137,TU,TUD,YU,YUD)
        CALL MPROD4X(X,138,YU,YUD,TU,TUD)
        CALL MPROD4X(X,139,TU,YUD,YU,TUD)
        CALL MPROD4X(X,140,YU,TUD,TU,YUD)
        CALL MPROD4X(X,141,TU,TDD,YD,YUD)
        CALL MPROD4X(X,142,YU,YDD,TD,TUD)
        CALL MPROD4X(X,143,TU,YDD,YD,TUD)
        CALL MPROD4X(X,144,YU,TDD,TD,YUD)
        CALL MPROD4X(X,157,YD,YUD,YU,YDD)
        CALL MPROD4X(X,162,TD,TDD,YD,YDD)
        CALL MPROD4X(X,163,YD,YDD,TD,TDD)
        CALL MPROD4X(X,164,TD,YDD,YD,TDD)
        CALL MPROD4X(X,165,YD,TDD,TD,YDD)
        CALL MPROD4X(X,166,TD,TUD,YU,YDD)
        CALL MPROD4X(X,167,YD,YUD,TU,TDD)
        CALL MPROD4X(X,168,TD,YUD,YU,TDD)
        CALL MPROD4X(X,169,YD,TUD,TU,YDD)
        CALL MPROD4X(X,175,TE,TED,YE,YED)
        CALL MPROD4X(X,176,YE,YED,TE,TED)
        CALL MPROD4X(X,177,TE,YED,YE,TED)
        CALL MPROD4X(X,178,YE,TED,TE,YED)
        CALL MPROD5X(X,60,TU,YUD,YU,YUD,YU)
        CALL MPROD5X(X,61,TU,YDD,YD,YDD,YD)
        CALL MPROD5X(X,62,TU,YDD,YD,YUD,YU)
        CALL MPROD5X(X,66,YU,YUD,YU,YUD,TU)
        CALL MPROD5X(X,67,YU,YUD,TU,YUD,YU)
        CALL MPROD5X(X,68,YU,YDD,YD,YDD,TD)
        CALL MPROD5X(X,69,YU,YDD,TD,YDD,YD)
        CALL MPROD5X(X,70,YU,YDD,YD,YUD,TU)
        CALL MPROD5X(X,71,YU,YDD,TD,YUD,YU)
        CALL MPROD5X(X,72,TD,YDD,YD,YDD,YD)
        CALL MPROD5X(X,73,TD,YUD,YU,YUD,YU)
        CALL MPROD5X(X,74,TD,YUD,YU,YDD,YD)
        CALL MPROD5X(X,77,YD,YDD,YD,YDD,TD)
        CALL MPROD5X(X,78,YD,YDD,TD,YDD,YD)
        CALL MPROD5X(X,79,YD,YUD,TU,YUD,YU)
        CALL MPROD5X(X,80,YD,YUD,YU,YUD,TU)
        CALL MPROD5X(X,81,YD,YUD,TU,YDD,YD)
        CALL MPROD5X(X,82,YD,YUD,YU,YDD,TD)
        CALL MPROD5X(X,83,TE,YED,YE,YED,YE)
        CALL MPROD5X(X,84,YE,YED,YE,YED,TE)
        CALL MPROD5X(X,85,YE,YED,TE,YED,YE)
        CALL MPROD5X(X,86,M2Q,YUD,YU,YUD,YU)
        CALL MPROD5X(X,87,YUD,M2U,YU,YUD,YU)
        CALL MPROD5X(X,88,M2Q,YUD,YU,YDD,YD)
        CALL MPROD5X(X,89,YUD,M2U,YU,YDD,YD)
        CALL MPROD5X(X,90,YUD,YU,M2Q,YDD,YD)
        CALL MPROD5X(X,91,YUD,YU,YDD,M2D,YD)
        CALL MPROD5X(X,98,M2Q,YDD,YD,YDD,YD)
        CALL MPROD5X(X,99,YDD,M2D,YD,YDD,YD)
        CALL MPROD5X(X,100,M2L,YED,YE,YED,YE)
        CALL MPROD5X(X,101,YED,M2E,YE,YED,YE)
        CALL MPROD5X(X,106,YUD,YU,M2Q,YUD,YU)
        CALL MPROD5X(X,107,YUD,YU,YUD,M2U,YU)
        CALL MPROD5X(X,108,YUD,YU,YUD,YU,M2Q)
        CALL MPROD5X(X,109,YDD,YD,M2Q,YDD,YD)
        CALL MPROD5X(X,110,YDD,YD,YDD,M2D,YD)
        CALL MPROD5X(X,111,YDD,YD,YDD,YD,M2Q)
        CALL MPROD5X(X,119,YED,YE,M2L,YED,YE)
        CALL MPROD5X(X,120,YED,YE,YED,M2E,YE)
        CALL MPROD5X(X,121,YED,YE,YED,YE,M2L)
        CALL MPROD5X(X,124,M2U,YU,YUD,YU,YUD)
        CALL MPROD5X(X,128,YU,M2Q,YUD,YU,YUD)
        CALL MPROD5X(X,129,YU,YUD,M2U,YU,YUD)
        CALL MPROD5X(X,130,YU,YUD,YU,M2Q,YUD)
        CALL MPROD5X(X,131,YU,YUD,YU,YUD,M2U)
        CALL MPROD5X(X,132,M2U,YU,YDD,YD,YUD)
        CALL MPROD5X(X,133,YU,M2Q,YDD,YD,YUD)
        CALL MPROD5X(X,134,YU,YDD,M2D,YD,YUD)
        CALL MPROD5X(X,135,YU,YDD,YD,M2Q,YUD)
        CALL MPROD5X(X,136,YU,YDD,YD,YUD,M2U)
        CALL MPROD5X(X,151,M2D,YD,YDD,YD,YDD)
        CALL MPROD5X(X,152,YD,M2Q,YDD,YD,YDD)
        CALL MPROD5X(X,153,YD,YDD,M2D,YD,YDD)
        CALL MPROD5X(X,154,YD,YDD,YD,M2Q,YDD)
        CALL MPROD5X(X,155,YD,YDD,YD,YDD,M2D)
        CALL MPROD5X(X,156,M2D,YD,YUD,YU,YDD)
        CALL MPROD5X(X,158,YD,M2Q,YUD,YU,YDD)
        CALL MPROD5X(X,159,YD,YUD,M2U,YU,YDD)
        CALL MPROD5X(X,160,YD,YUD,YU,M2Q,YDD)
        CALL MPROD5X(X,161,YD,YUD,YU,YDD,M2D)
        CALL MPROD5X(X,170,M2E,YE,YED,YE,YED)
        CALL MPROD5X(X,171,YE,M2L,YED,YE,YED)
        CALL MPROD5X(X,172,YE,YED,M2E,YE,YED)
        CALL MPROD5X(X,173,YE,YED,YE,M2L,YED)
        CALL MPROD5X(X,174,YE,YED,YE,YED,M2E)
c...compute some convenient quantities
	SP=-(3.d0*G(63)*TRACE(X,1)+TRACE(X,31))+4.d0*TRACE(X,34)
     &     +3.d0*G(64)*TRACE(X,2)-TRACE(X,25)-2.d0*TRACE(X,26)
     &     +G(64)*TRACE(X,3)+TRACE(X,27)-2.d0*TRACE(X,28)
     &     +(3./2.d0*G(2)**2+3./10.d0*G(1)**2)
     &      *(G(63)-G(64)-TR3X3(M2L))
     &     +(8./3.d0*G(3)**2+3./2.d0*G(2)**2+G(1)**2/30.d0)*TR3X3(M2Q)
     &     -(16./3.d0*G(3)**2+16./15.d0*G(1)**2)*TR3X3(M2U)
     &     +(8./3.d0*G(3)**2+2./15.d0*G(1)**2)*TR3X3(M2D)
     &     +6./5.d0*G(1)**2*TR3X3(M2E)
        SIG1=G(1)**2/5.d0*(3.d0*(G(63)+G(64))+TR3X3(M2Q)
     &			   +3.d0*TR3X3(M2L)+8.d0*TR3X3(M2U)
     &			   +2.d0*TR3X3(M2D)+6.d0*TR3X3(M2E))
	SIG2=G(2)**2*(G(63)+G(64)+3.d0*TR3X3(M2Q)+TR3X3(M2L))
	SIG3=G(3)**2*(2.d0*TR3X3(M2Q)+TR3X3(M2U)+TR3X3(M2D))
c...gaugino masses      
      F2(31)=2.d0*G(1)**2
     &	     *(199./25.d0*G(1)**2*(G(31)+G(31))
     &	       +27./5.d0*G(2)**2*(G(31)+G(32))
     &	       +88./5.d0*G(3)**2*(G(31)+G(33))
     &	       +26./5.d0*(TRACE(X,11)-G(31)*TRACE(X,1))
     &	       +14./5.d0*(TRACE(X,16)-G(31)*TRACE(X,2))
     &	       +18./5.d0*(TRACE(X,17)-G(31)*TRACE(X,3)))
      F2(32)=2.d0*G(2)**2
     &	     *(9./5.d0*G(1)**2*(G(32)+G(31))
     &	       +25.d0*G(2)**2*(G(32)+G(32))
     &	       +24.d0*G(3)**2*(G(32)+G(33))
     &	       +6.d0*(TRACE(X,11)-G(32)*TRACE(X,1))
     &	       +6.d0*(TRACE(X,16)-G(32)*TRACE(X,2))
     &	       +2.d0*(TRACE(X,17)-G(32)*TRACE(X,3)))
      F2(33)=2.d0*G(3)**2
     &	     *(11./5.d0*G(1)**2*(G(33)+G(31))
     &	       +9.d0*G(2)**2*(G(33)+G(32))
     &	       +14.d0*G(3)**2*(G(33)+G(33))
     &	       +4.d0*(TRACE(X,11)-G(33)*TRACE(X,1))
     &	       +4.d0*(TRACE(X,16)-G(33)*TRACE(X,2)))
c...higgs mixing parameters      
      F2(61)=G(61)*(-3.*(3.d0*TRACE(X,49)+3.d0*TRACE(X,54)
     &			 +2.d0*TRACE(X,50)+TRACE(X,55))
     &  	    +(16.d0*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &  	    +(16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &  	    +6./5.d0*G(1)**2*TRACE(X,3)
     &  	    +15./2.d0*G(2)**4+9./5.d0*G(1)**2*G(2)**2
     &  	    +207./50.d0*G(1)**4)
      F2(62)=G(62)*(-3.*(3.d0*TRACE(X,49)+3.d0*TRACE(X,54)
     &	        	 +2.d0*TRACE(X,50)+TRACE(X,55))
     &              +(16.d0*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &              +(16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &              +6./5.d0*G(1)**2*TRACE(X,3)
     &              +15./2.d0*G(2)**4+9./5.d0*G(1)**2*G(2)**2
     &              +207./50.d0*G(1)**4)
     &      +G(61)*(-12.*(3.d0*TRACE(X,63)+3.d0*TRACE(X,75)
     &      		  +TRACE(X,64)+TRACE(X,65)+TRACE(X,76))
     &      	    +(32.d0*G(3)**2+8./5.d0*G(1)**2)*TRACE(X,11)
     &      	    +(32.d0*G(3)**2-4./5.d0*G(1)**2)*TRACE(X,16)
     &      	    +12./5.d0*G(1)**2*TRACE(X,17)
     &      	    -(32.d0*G(3)**2*G(33)+8./5.d0*G(1)**2*G(31))
     &      	     *TRACE(X,1)
     &      	    -(32.d0*G(3)**2*G(33)-4./5.d0*G(1)**2*G(31))
     &      	     *TRACE(X,2)
     &      	    -12./5.d0*G(1)**2*G(31)*TRACE(X,3)
     &      	    -30.d0*G(2)**4*G(32)
     &      	    -18./5.d0*G(1)**2*G(2)**2*(G(31)+G(32))
     &      	    -414./25.d0*G(1)**4*G(31))
c...higgs mass^2
      F2(63)=-6.*(6.*(G(63)*TRACE(X,49)+TRACE(X,86))+6.*TRACE(X,87)
     &       	  +(G(63)+G(64))*TRACE(X,50)+TRACE(X,88)+TRACE(X,89)
     &       	  +TRACE(X,90)+TRACE(X,91)+6.*TRACE(X,92)
     &       	  +6.*TRACE(X,93)+TRACE(X,94)+TRACE(X,95)+TRACE(X,96)
     &       	  +TRACE(X,97))
     &       +(32.*G(3)**2+8./5.d0*G(1)**2)
     &        *(G(63)*TRACE(X,1)+TRACE(X,31)+TRACE(X,34)+TRACE(X,24))
     &       +32.*G(3)**2*(2.*G(33)**2*TRACE(X,1)-G(33)*TRACE(X,11)
     &       		   -G(33)*TRACE(X,112))
     &       +8./5.d0*G(1)**2*(2.*G(31)**2*TRACE(X,1)-G(31)*TRACE(X,11)
     &       		       -G(31)*TRACE(X,112))
     &       +6./5.d0*G(1)**2*SP
     &       +33.*G(2)**4*G(32)**2
     &       +18./5.d0*G(2)**2*G(1)**2*(G(32)**2+G(31)**2+G(31)*G(32))
     &       +621./25.d0*G(1)**4*G(31)**2
     &       +3.*G(2)**2*SIG2+3./5.d0*G(1)**2*SIG1
      F2(64)=-6.*(6.*(G(64)*TRACE(X,54)+TRACE(X,98))+6.*TRACE(X,99)
     &       	  +(G(63)+G(64))*TRACE(X,50)+TRACE(X,88)+TRACE(X,89)
     &       	  +TRACE(X,90)+TRACE(X,91)+2.*G(64)*TRACE(X,55)
     &       	  +2.*TRACE(X,100)+2.*TRACE(X,101)+6.*TRACE(X,102)
     &       	  +6.*TRACE(X,103)+TRACE(X,95)+TRACE(X,94)+TRACE(X,97)
     &       	  +TRACE(X,96)+2.*TRACE(X,104)+2.*TRACE(X,105))
     &       +(32.*G(3)**2-4./5.d0*G(1)**2)
     &        *(G(64)*TRACE(X,2)+TRACE(X,25)+TRACE(X,26)+TRACE(X,29))
     &       +32.*G(3)**2*(2.*G(33)**2*TRACE(X,2)-G(33)*TRACE(X,16)
     &       		   -G(33)*TRACE(X,113))
     &       -4./5.d0*G(1)**2*(2.*G(31)**2*TRACE(X,2)-G(31)*TRACE(X,16)
     &       		       -G(31)*TRACE(X,113))
     &       +12./5.d0*G(1)**2*(G(64)*TRACE(X,3)+TRACE(X,27)
     &                          +TRACE(X,28)+TRACE(X,45)
     &                          +2.*G(31)**2*TRACE(X,3)
     &                          -G(31)*TRACE(X,114)-G(31)*TRACE(X,17))
     &       -6./5.d0*G(1)**2*SP
     &       +33.*G(2)**4*G(32)**2
     &       +18./5.d0*G(2)**2*G(1)**2*(G(32)**2+G(31)**2+G(31)*G(32))
     &       +621./25.d0*G(1)**4*G(31)**2
     &       +3.*G(2)**2*SIG2+3./5.d0*G(1)**2*SIG1

        K=0
	DO 400 I=1,3
          DO 400 J=1,3
C...soft trilinear scalar couplings
      F2(34+K)=TU(I,J)*(-3.)*(3.*TRACE(X,49)+TRACE(X,50))
     &         -X(10,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -15.d0*X(9,I,J)*TRACE(X,1)
     &         -6.d0*X(60,I,J)-2.d0*X(61,I,J)-4.d0*X(62,I,J)
     &         +TU(I,J)*(16.*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &         +12.d0*G(2)**2*X(9,I,J)
     &         +2./5.d0*G(1)**2*X(10,I,J)
     &         +TU(I,J)*(-16./9.d0*G(3)**4+8.*G(3)**2*G(2)**2
     &                   +136./45.d0*G(3)**2*G(1)**2+15./2.d0*G(2)**4
     &                   +G(2)**2*G(1)**2+2743./450.d0*G(1)**4)
     &         +YU(I,J)*(-6.)*(6.*TRACE(X,63)+TRACE(X,64)+TRACE(X,65))
     &         -18.d0*X(4,I,J)*TRACE(X,11)
     &         -X(5,I,J)*(6.d0*TRACE(X,16)+2.d0*TRACE(X,17))
     &         -12.d0*X(12,I,J)*TRACE(X,1)
     &         -X(13,I,J)*(6.d0*TRACE(X,2)+2.d0*TRACE(X,3))
     &         -6.*X(66,I,J)-8.*X(67,I,J)-4.*X(68,I,J)-4.*X(69,I,J)
     &         -2.*X(70,I,J)-4.*X(71,I,J)
     &         +YU(I,J)*(32.*G(3)**2+8./5.d0*G(1)**2)*TRACE(X,11)
     &         +(6.*G(2)**2+6./5.d0*G(1)**2)*X(12,I,J)
     &         +4./5.d0*G(1)**2*X(13,I,J)
     &         -YU(I,J)*(32.*G(3)**2*G(33)+8./5.d0*G(1)**2*G(31))
     &          *TRACE(X,1)
     &         -(12.*G(2)**2*G(32)+4./5.d0*G(1)**2*G(31))*X(4,I,J)
     &         -4./5.d0*G(1)**2*G(31)*X(5,I,J)
     &         +YU(I,J)*(64./9.d0*G(3)**4*G(33)
     &                   -16.*G(3)**2*G(2)**2*(G(33)+G(32))
     &                   -272./45.d0*G(3)**2*G(1)**2*(G(33)+G(31))
     &                   -30.*G(2)**4*G(32)
     &                   -2.*G(2)**2*G(1)**2*(G(32)+G(31))
     &                   -5486./225.d0*G(1)**4*G(31))
      F2(43+K)=TD(I,J)*(-3.)*(3.*TRACE(X,54)+TRACE(X,50)+TRACE(X,55))
     &         -3.*X(15,I,J)*TRACE(X,1)
     &         -5.*X(14,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(72,I,J)-2.*X(73,I,J)-4.*X(74,I,J)
     &         +TD(I,J)*((16.*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +4./5.d0*G(1)**2*X(15,I,J)
     &         +(12.*G(2)**2+6./5.d0*G(1)**2)*X(14,I,J)
     &         +TD(I,J)*(-16./9.d0*G(3)**4+8.*G(3)**2*G(2)**2
     &                   +8./9.d0*G(3)**2*G(1)**2+15./2.d0*G(2)**4
     &                   +G(2)**2*G(1)**2+287./90.d0*G(1)**4)
     &         +YD(I,J)*(-6.)*(6.*TRACE(X,75)+TRACE(X,64)+TRACE(X,65)
     &                         +2.*TRACE(X,76))
     &         -6.*X(7,I,J)*TRACE(X,11)
     &         -6.*X(6,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &         -6.*X(19,I,J)*TRACE(X,1)
     &         -4.*X(18,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(77,I,J)-8.*X(78,I,J)-4.*X(79,I,J)-4.*X(80,I,J)
     &         -4.*X(81,I,J)-2.*X(82,I,J)
     &         +YD(I,J)*((32.*G(3)**2-4./5.d0*G(1)**2)*TRACE(X,16)
     &                   +12./5.d0*G(1)**2*TRACE(X,17))
     &         +8./5.d0*G(1)**2*X(19,I,J)
     &         +(6.*G(2)**2+6./5.d0*G(1)**2)*X(18,I,J)
     &         -YD(I,J)*((32.*G(3)**2*G(33)-4./5.d0*G(1)**2*G(31))
     &                    *TRACE(X,2)
     &                   +12./5.d0*G(1)**2*G(31)*TRACE(X,3))
     &         -(12.*G(2)**2*G(32)+8./5.d0*G(1)**2*G(31))*X(6,I,J)
     &         -8./5.d0*G(1)**2*G(31)*X(7,I,J)
     &         +YD(I,J)*(64./9.d0*G(3)**4*G(33)
     &                   -16.*G(3)**2*G(2)**2*(G(33)+G(32))
     &                   -16./9.d0*G(3)**2*G(1)**2*(G(33)+G(31))
     &                   -30.*G(2)**4*G(32)
     &                   -2.*G(2)**2*G(1)**2*(G(32)+G(31))
     &                   -574./45.d0*G(1)**4*G(31))
      F2(52+K)=TE(I,J)*(-3.)*(3.*TRACE(X,54)+TRACE(X,50)+TRACE(X,55))
     &         -5.*X(20,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(83,I,J)
     &         +TE(I,J)*((16.*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +(12.*G(2)**2-6./5.d0*G(1)**2)*X(20,I,J)
     &         +TE(I,J)*(15./2.d0*G(2)**4+9./5.d0*G(2)**2*G(1)**2
     &                   +27./2.d0*G(1)**4)
     &         +YE(I,J)*(-6.)*(6.*TRACE(X,75)+TRACE(X,64)+TRACE(X,65)
     &                         +2.*TRACE(X,76))
     &         -4.*X(21,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(8,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &         -6.*X(84,I,J)-8.*X(85,I,J)
     &         +YE(I,J)*((32.*G(3)**2-4./5.d0*G(1)**2)*TRACE(X,16)
     &                   +12./5.d0*G(1)**2*TRACE(X,17))
     &         +(6.*G(2)**2+6./5.d0*G(1)**2)*X(21,I,J)
     &         -YE(I,J)*((32.*G(3)**2*G(33)-4./5.d0*G(1)**2*G(31))
     &                    *TRACE(X,2)
     &                   +12./5.d0*G(1)**2*G(31)*TRACE(X,3))
     &         -12.*G(2)**2*G(32)*X(8,I,J)
     &         +YE(I,J)*(-30.*G(2)**4*G(32)
     &                   -18./5.d0*G(2)**2*G(1)**2*(G(31)+G(32))
     &                   -54.*G(1)**4*G(31))
c...Squark doublet mass^2
      F2(65+K)=-(2.*X(86,I,J)+8.*G(63)*X(125,I,J))-4.*X(87,I,J)
     &         -4.*X(106,I,J)-4.*X(107,I,J)-2.*X(108,I,J)
     &         -(2.*X(98,I,J)+8.*G(64)*X(126,I,J))-4.*X(99,I,J)
     &         -4.*X(109,I,J)-4.*X(110,I,J)-2.*X(111,I,J)
     &         -(X(22,I,J)+4.*G(63)*X(1,I,J)+2.*X(23,I,J)+X(31,I,J))
     &          *3.*TRACE(X,1)
     &         -(X(25,I,J)+4.*G(64)*X(2,I,J)+2.*X(26,I,J)+X(32,I,J))
     &          *(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(1,I,J)*(TRACE(X,22)+TRACE(X,23))
     &         -X(2,I,J)*(6.*TRACE(X,25)+6.*TRACE(X,26)+2.*TRACE(X,27)
     &          	  +2.*TRACE(X,28))
     &         -4.*(X(115,I,J)+X(92,I,J)+X(116,I,J)+X(93,I,J))
     &         -4.*(X(117,I,J)+X(102,I,J)+X(118,I,J)+X(103,I,J))
     &         -X(24,I,J)*6.*TRACE(X,1)-X(1,I,J)*6.*TRACE(X,24)
     &         -X(112,I,J)*6.*TRACE(X,11)-X(11,I,J)*6.*TRACE(X,112)
     &         -X(29,I,J)*(6.*TRACE(X,2)+2.*TRACE(X,3))
     &         -X(2,I,J)*(6.*TRACE(X,29)+2.*TRACE(X,30))
     &         -X(113,I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))
     &         -X(16,I,J)*(6.*TRACE(X,113)+2.*TRACE(X,114))
     &         +2./5.d0*G(1)**2*(2.*X(22,I,J)+4.*G(63)*X(1,I,J)
     &          	        +4.*X(23,I,J)+2.*X(31,I,J)+4.*X(24,I,J)
     &          	        -4.*G(31)*X(112,I,J)-4.*G(31)*X(11,I,J)
     &          	        +8.*G(31)**2*X(1,I,J)
     &                          +X(25,I,J)+2.*G(64)*X(2,I,J)
     &                          +2.*X(26,I,J)+X(32,I,J)+2.*X(29,I,J)
     &          	        -2.*G(31)*X(113,I,J)-2.*G(31)*X(16,I,J)
     &          	        +4.*G(31)**2*X(2,I,J))
     &         +I3(I,J)*(2./5.d0*G(1)**2*SP
     &          	 -128./3.*G(3)**4*G(33)**2
     &          	 +32.*G(3)**2*G(2)**2
     &          	  *(G(33)**2+G(32)**2+G(32)*G(33))
     &          	 +32./45.d0*G(3)**2*G(1)**2
     &          	  *(G(33)**2+G(31)**2+G(31)*G(33))
     &          	 +33.*G(2)**4*G(32)**2
     &          	 +2./5.d0*G(2)**2*G(1)**2
     &          	  *(G(32)**2+G(31)**2+G(31)*G(32))
     &          	 +199./75.d0*G(1)**4*G(31)**2
     &          	 +16./3.d0*G(3)**2*SIG3+3.*G(2)**2*SIG2
     &          	 +1./15.d0*G(1)**2*SIG1)
c...Slepton doublet mass^2
      F2(74+K)=-(2.*X(100,I,J)+8.*G(64)*X(127,I,J))-4.*X(101,I,J)
     &         -4.*X(119,I,J)-4.*X(120,I,J)-2.*X(121,I,J)
     &         -(X(27,I,J)+4.*G(64)*X(3,I,J)+2.*X(28,I,J)+X(33,I,J))
     &          *(3.*TRACE(X,2)+TRACE(X,3))
     &         -X(3,I,J)*(6.*X(25,I,J)+6.*X(26,I,J)+2.*X(27,I,J)
     &          	  +2.*X(28,I,J))
     &         -4.*(X(122,I,J)+X(104,I,J)+X(123,I,J)+X(105,I,J))
     &         -X(30,I,J)*(6.*TRACE(X,2)+2.*TRACE(X,3))
     &         -X(3,I,J)*(6.*TRACE(X,29)+2.*TRACE(X,30))
     &         -X(114,I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))
     &         -X(17,I,J)*(6.*TRACE(X,113)+2.*TRACE(X,114))
     &         +6./5.d0*G(1)**2*(X(27,I,J)+2.*G(64)*X(3,I,J)
     &          	        +2.*X(28,I,J)+X(33,I,J)+2.*X(30,I,J)
     &          	        -2.*G(31)*X(114,I,J)-2.*G(31)*X(17,I,J)
     &          	        +4.*G(31)**2*X(3,I,J))
     &         +I3(I,J)*(-6./5.d0*G(1)**2*SP
     &          	 +33.*G(2)**4*G(32)**2
     &          	 +18./5.d0*G(2)**2*G(1)**2
     &          	  *(G(32)**2+G(31)**2+G(31)*G(32))
     &          	 +621./25.d0*G(1)**4*G(31)**2
     &          	 +3.*G(2)**2*SIG2+3./5.d0*G(1)**2*SIG1)
c...Up squark singlet mass^2
      F2(83+K)=-(2.*X(124,I,J)+8.*G(63)*X(49,I,J))-4.*X(128,I,J)
     &         -4.*X(129,I,J)-4.*X(130,I,J)-2.*X(131,I,J)
     &         -(2.*X(132,I,J)+4.*(G(63)+G(64))*X(50,I,J))
     &         -4.*X(133,I,J)-4.*X(134,I,J)-4.*X(135,I,J)
     &         -2.*X(136,I,J)
     &         -(X(34,I,J)+4.*G(63)*X(46,I,J)+2.*X(35,I,J)+X(36,I,J))
     &          *6.*TRACE(X,1)
     &         -12.*X(46,I,J)*(TRACE(X,22)+TRACE(X,23))
     &         -4.*(X(137,I,J)+X(138,I,J)+X(139,I,J)+X(140,I,J))
     &         -4.*(X(141,I,J)+X(142,I,J)+X(143,I,J)+X(144,I,J))
     &         -12.*(X(37,I,J)*TRACE(X,1)+X(46,I,J)*TRACE(X,24)
     &               +X(145,I,J)*TRACE(X,112)+X(148,I,J)*TRACE(X,11))
     &         +(6.*G(2)**2-2./5.d0*G(1)**2)
     &          *(X(34,I,J)+2.*G(63)*X(46,I,J)+2.*X(35,I,J)+X(36,I,J)
     &            +2.*X(37,I,J))
     &         +12.*G(2)**2*(2.*G(32)**2*X(46,I,J)-G(32)*X(145,I,J)
     &          	     -G(32)*X(148,I,J))
     &         -4./5.d0*G(1)**2*(2.*G(31)**2*X(46,I,J)-G(31)*X(145,I,J)
     &          	         -G(31)*X(148,I,J))
     &         +I3(I,J)*(-8./5.d0*G(1)**2*SP
     &          	 -128./3.d0*G(3)**4*G(33)**2
     &          	 +512./45.d0*G(3)**2*G(1)**2
     &          	  *(G(33)**2+G(31)**2+G(31)*G(33))
     &          	 +3424./75.d0*G(1)**4*G(31)**2
     &          	 +16./3.d0*G(3)**2*SIG3
     &          	 +16./15.d0*G(1)**2*SIG1)
c...Down squark singlet mass^2
      F2(92+K)=-(2.*X(151,I,J)+8.*G(64)*X(54,I,J))-4.*X(152,I,J)
     &         -4.*X(153,I,J)-4.*X(154,I,J)-2.*X(155,I,J)
     &         -(2.*X(156,I,J)+4.*(G(63)+G(64))*X(157,I,J))
     &         -4.*X(158,I,J)-4.*X(159,I,J)-4.*X(160,I,J)
     &         -2.*X(161,I,J)
     &         -(X(38,I,J)+4.*G(64)*X(47,I,J)+2.*X(39,I,J)+X(40,I,J))
     &          *(6.*TRACE(X,2)+2.*TRACE(X,3))
     &         -4.*X(47,I,J)*(3.*TRACE(X,25)+3.*TRACE(X,26)
     &          	      +TRACE(X,27)+TRACE(X,28))
     &         -4.*(X(162,I,J)+X(163,I,J)+X(164,I,J)+X(165,I,J))
     &         -4.*(X(166,I,J)+X(167,I,J)+X(168,I,J)+X(169,I,J))
     &         -4.*X(41,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -4.*X(47,I,J)*(3.*TRACE(X,29)+TRACE(X,30))
     &         -4.*X(146,I,J)*(3.*TRACE(X,113)+TRACE(X,114))
     &         -4.*X(149,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &         +(6.*G(2)**2+2./5.d0*G(1)**2)
     &          *(X(38,I,J)+2.*G(64)*X(47,I,J)+2.*X(39,I,J)
     &            +X(40,I,J)+2.*X(41,I,J))
     &         +12.*G(2)**2*(2.*G(32)**2*X(47,I,J)-G(32)*X(146,I,J)
     &          	     -G(32)*X(149,I,J))
     &         +4./5.d0*G(1)**2*(2.*G(31)**2*X(47,I,J)-G(31)*X(146,I,J)
     &          	         -G(31)*X(149,I,J))
     &         +I3(I,J)*(4./5.d0*G(1)**2*SP
     &          	 -128./3.d0*G(3)**4*G(33)**2
     &          	 +128./45.d0*G(3)**2*G(1)**2
     &          	  *(G(33)**2+G(31)**2+G(31)*G(33))
     &          	 +808./75.d0*G(1)**4*G(31)**2
     &          	 +16./3.d0*G(3)**2*SIG3
     &          	 +4./15.d0*G(1)**2*SIG1)
c...Charged slepton singlet mass^2
      F2(101+K)=-(2.*X(170,I,J)+8.*G(64)*X(55,I,J))-4.*X(171,I,J)
     &          -4.*X(172,I,J)-4.*X(173,I,J)-2.*X(174,I,J)
     &          -(X(42,I,J)+4.*G(64)*X(48,I,J)+2.*X(43,I,J)+X(44,I,J))
     &           *(6.*TRACE(X,2)+2.*TRACE(X,3))
     &          -4.*X(48,I,J)*(3.*TRACE(X,25)+3.*TRACE(X,26)
     &          	       +TRACE(X,27)+TRACE(X,44))
     &          -4.*(X(175,I,J)+X(176,I,J)+X(177,I,J)+X(178,I,J))
     &          -4.*X(45,I,J)*(2.*TRACE(X,2)+TRACE(X,3))
     &          -4.*X(48,I,J)*(3.*TRACE(X,29)+TRACE(X,30))
     &          -4.*X(147,I,J)*(3.*TRACE(X,113)+TRACE(X,114))
     &          -4.*X(150,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &          +(6.*G(2)**2-6./5.d0*G(1)**2)
     &           *(X(42,I,J)+2.*G(64)*X(48,I,J)+2.*X(43,I,J)+X(44,I,J)
     &             +2.*X(45,I,J))
     &          +12.*G(2)**2*(2.*G(32)**2*X(48,I,J)-G(32)*X(147,I,J)
     &          	      -G(32)*X(150,I,J))
     &          -12./5.d0*G(1)**2*(2.*G(31)**2*X(48,I,J)
     &          		   -G(31)*X(147,I,J)-G(31)*X(150,I,J))
     &          +I3(I,J)*(12./5.d0*G(1)**2*SP
     &          	  +2808./25.d0*G(1)**4*G(31)**2
     &          	  +12./5.d0*G(1)**2*SIG1)
	    K=K+1
400     CONTINUE
c...VEVs
        F2(110)=0.d0
	F2(111)=0.d0
c--------- neutrino sector ------------------------------
        IF (IRHN.GT.0) THEN
          IF(
     &Q.GE.MIN(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
	    K=0
            DO I=1,3
              DO J=1,3
c...contribution to gaugino masses
            F2(31)=F2(31)+2.d0*G(1)**2 
     &                    *6./5.d0*(TRACE(X,215)-G(31)*TRACE(X,200))
            F2(32)=F2(32)+2.d0*G(2)**2
     &                    *2.d0*(TRACE(X,215)-G(32)*TRACE(X,200))
c...Majorana mass matrix
      F2(121+K)=-2.d0*X(213,J,I)-2.d0*X(214,J,I)
     &          -2.d0*X(213,I,J)-2.d0*X(214,I,J)
     &          -(X(205,J,I)+X(205,I,J))
     &           *(6.d0*TRACE(X,1)+2.d0*TRACE(X,201)
     &             -6./5.d0*G(1)**2-6.d0*G(2)**2)
	        K=K+1
             ENDDO
	    ENDDO
	    DO I=121,129
	      F(I)=F(I)+F2(I)/FAC**2
	    ENDDO
	  ENDIF
c...effective dim-5 neutrino operator
          IF (LND5ON.AND.Q.GE.DBLE(MSUSY).AND.
     &   Q.LE.MAX(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
            CALL MPROD5X(X,237,KA,YED,YE,YED,YE)
            CALL MPROD5X(X,238,KA,YND,YN,YND,YN)
	    K=0
            DO I=1,3
              DO J=1,3
      F2(148+K)=KA(I,J)*(-6.d0*TRACE(X,50)-18.d0*TRACE(X,49)
     &                   -2.d0*TRACE(X,209)-6.d0*TRACE(X,210)
     &                   +(8./5.d0*G(1)**2+32.d0*G(3)**2)*TRACE(X,1)
     &                   +207./25.d0*G(1)**4+18./5.d0*G(1)**2*G(2)**2
     &                   +15.d0*G(2)**4)
     &           -2.d0*X(237,J,I)-2.d0*X(238,J,I)
     &           -2.d0*X(237,I,J)-2.d0*X(238,I,J)
     &           -(X(232,J,I)+X(232,I,J))*(TRACE(X,201)+3.*TRACE(X,1))
     &           -(X(231,J,I)+X(231,I,J))
     &            *(-6./5.d0*G(1)**2+TRACE(X,3)+3.*TRACE(X,2))            
	        K=K+1
              ENDDO
	    ENDDO
	    DO I=148,156
	      F(I)=F(I)+F2(I)/FAC**2
	    ENDDO
	  ENDIF
	ENDIF
	DO I=31,111
	  F(I)=F(I)+F2(I)/FAC**2
	ENDDO
      ENDIF
      
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE GAMMASM(T,CSM,F)            
C-----------------------------------------------------------------
C     Eq. 24 of Anlauf
      IMPLICIT NONE
      REAL*8 T,CSM(9),F(9)
      COMMON /G3/ G3         
      REAL*8 G3
      SAVE /G3/         
      
      REAL*8 PI,FAC
      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)
C
      F(1)=FAC*(20./3.d0*CSM(1)-8.d0*CSM(2)+6.d0*CSM(4)
     $         +4.d0*CSM(5)+2.d0*CSM(6))
      F(2)=FAC*(CSM(1)+2./3.d0*CSM(2)+2.d0*CSM(4)
     $         +3./2.d0*CSM(5)+CSM(6)+1./2.d0*CSM(7))
      F(3)=FAC*(-2.d0*CSM(1)+4./3.d0*CSM(2)+16./3.d0*CSM(3)
     $         -CSM(4)+CSM(6)+2.d0*CSM(7)+4./3.d0*CSM(9))
      F(4)=FAC*(2./3.d0*CSM(4)-113./36.d0*CSM(5)-2.*CSM(6)
     $         -113./36.d0*CSM(7))
      F(5)=FAC*(2.d0*CSM(4)+137./18.d0*CSM(5)+2.d0*CSM(6)
     $         +89./18.d0*CSM(7))
      F(6)=FAC*(-2.d0*CSM(4)-113./36.d0*CSM(5)+2./3.d0*CSM(6)
     $         -113./36.d0*CSM(7))
      F(7)=FAC*(-2.d0*CSM(4)-4./3.d0*CSM(5)-2.d0*CSM(6)
     $         +4./3.d0*CSM(7))
      F(8)=FAC*(9./4.d0*CSM(5)+9./4.d0*CSM(7))
      F(9)=0.d0
      
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE GAMMAHP(T,CH2,F)            
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,CH2(19),F(19)
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)     
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX
      SAVE /CHRGN/
      COMMON/SCALEQ/Q
      REAL*8 Q
      SAVE /SCALEQ/
      COMMON /G3/G3
      REAL*8 G3
      SAVE /G3/         
      INTEGER NS,NF,NFT
      REAL*8 PI,FAC,BB
      
      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)

      NS=0
      IF(MSQU(1).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(2).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(3).LT.Q) THEN
       NS=NS+10
      ENDIF             
      NF=6
      NFT=3-2*3
      BB=-11./2.d0+1./3.d0*NF+1./12.d0*NS
C
      F(1)=FAC*(20./3.d0*CH2(1)-8.d0*CH2(2)+6.d0*CH2(4)
     $         +4.d0*CH2(5)+2.d0*CH2(6))
      F(2)=FAC*(CH2(1)+2./3.d0*CH2(2)+2.d0*CH2(4)
     $         +3./2.d0*CH2(5)+CH2(6)+1./2.d0*CH2(7)
     $         -3./2.d0*CH2(10)-119./54.d0*CH2(11)
     $         +(70./27.d0+3./2.d0*NF)*CH2(12)
     $         +(3.d0+35./27.d0*NF)*CH2(13)
     $         +(-7./3.d0-3./2.d0*NF)*CH2(14)
     $         +(-2.d0-119./54.d0*NF)*CH2(15)
     $         -1./2.d0*CH2(17)+1./2.d0*CH2(19))
      F(3)=FAC*(-2.d0*CH2(1)+4./3.d0*CH2(2)+16./3.d0*CH2(3)
     $         -CH2(4)+CH2(6)+2.d0*CH2(7)+4./3.d0*CH2(9)
     $         +224./27.d0*CH2(11)+232./27.d0*CH2(12)
     $         +(8./27.d0*NF+4.d0*NFT)*CH2(13)
     $         -16./3.d0*CH2(14)
     $         +(8./27.d0*NF-4.d0*NFT)*CH2(15)
     $         +3.d0*CH2(16)+CH2(17)
     $         -3.d0*CH2(18)-CH2(19))
      F(4)=FAC*(2./3.d0*CH2(4)-113./36.d0*CH2(5)-2.d0*CH2(6)
     $         -113./36.d0*CH2(7)
     $         +2./3.d0*CH2(11)+4./3.d0*CH2(12)
     $         +2.*NF/3.d0*CH2(13)+2.*NF/3.d0*CH2(15))    
      F(5)=FAC*(2.d0*CH2(4)+137./18.d0*CH2(5)+2.d0*CH2(6)
     $         +89./18.d0*CH2(7)
     $         -4./3.d0*CH2(11)-8./3.d0*CH2(12)
     $         -4.*NF/3.d0*CH2(13)-4.*NF/3.d0*CH2(15))    
      F(6)=FAC*(-2.d0*CH2(4)-113./36.d0*CH2(5)+2./3.d0*CH2(6)
     $         -113./36.d0*CH2(7)
     $         +2./3.d0*CH2(11)+4./3.d0*CH2(12)
     $         +2.*NF/3.d0*CH2(13)+2.*NF/3.d0*CH2(15))    
      F(7)=FAC*(-2.d0*CH2(4)-4./3.d0*CH2(5)-2.d0*CH2(6)
     $         +4./3.d0*CH2(7))
      F(8)=FAC*(9./4.d0*CH2(5)+9./4.d0*CH2(7)
     $         -4.d0*CH2(10)-2.d0*CH2(11)-6.d0*CH2(12)
     $         +2./3.d0*(3-NF)*CH2(13)-6.d0*CH2(14)    
     $         +2./3.d0*(-3-NF)*CH2(15))     
      F(9)=0.d0
      F(10)=FAC*(CH2(10))
      F(11)=FAC*(-3.d0*CH2(10)-8.d0*CH2(11))
      F(12)=FAC*(-1./9.d0*CH2(11)-11./9.d0*CH2(12)
     $           +(3.d0-NF/9.d0)*CH2(13)-NF/9.d0*CH2(15)    
     $          -1./36.d0*RMIX(1)+1./36.d0*RMIX(2))   
      F(13)=FAC*(1./3.d0*CH2(11)+11./3.d0*CH2(12)
     $           +(NF/3.d0-1.d0)*CH2(13)+NF/3.d0*CH2(15)     
     $          +1./12.d0*RMIX(1)-1./12.d0*RMIX(2))   
      F(14)=FAC*(-1./9.d0*CH2(11)-2./9.d0*CH2(12)
     $           -NF/9.d0*CH2(13)+CH2(14)-NF/9.d0*CH2(15) 
     $          -1./36.d0*RMIX(1)+1./36.d0*RMIX(2))   
      F(15)=FAC*(1./3.d0*CH2(11)+2./3.d0*CH2(12)
     $           +NF/3.d0*CH2(13)-3.d0*CH2(14)
     $           +(NF/3.d0-8.d0)*CH2(15) 
     $          +1./12.d0*RMIX(1)-1./12.d0*RMIX(2))   
      F(16)=FAC*((1.d0-2.d0*BB)*CH2(16)-7.d0*CH2(18)-6.d0*CH2(19))    
      F(17)=FAC*(-3.d0*CH2(16)+(-8.d0-2.d0*BB)*CH2(17)
     $           -3.d0*CH2(18)+2.d0*CH2(19))    
      F(18)=FAC*(-7./3.d0*CH2(16)-2.d0*CH2(17)
     $           +(-19./3.d0-2.d0*BB)*CH2(18))   
      F(19)=FAC*(-CH2(16)+2./3.d0*CH2(17)
     $          +(8./3.d0-2.d0*BB)*CH2(19))
      
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE GAMMAC1(T,CPOA,F)            
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,CPOA(12),F(12)
      COMMON /G3/ G3
      REAL*8 G3
      SAVE /G3/         
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX             
      SAVE /CHRGN/
      COMMON /SCALEQ/Q
      REAL*8 Q
      SAVE /SCALEQ/
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
C
      INTEGER NS,NF,NFT
      REAL*8 PI,FAC,BB

      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)
      
      NS=0
      IF(MSQU(1).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(2).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(3).LT.Q) THEN
       NS=NS+10
      ENDIF             
      IF(Q.GT.MT) THEN 
       NF=6
       NFT=3-2*3
      ELSE
       NF=5
       NFT=3-2*2
      ENDIF
      BB=-11./2.d0+1./3.d0*DBLE(NF)+1./12.d0*DBLE(NS)
C
      F(1)=FAC*(20./3.d0*CPOA(1)-8.d0*CPOA(2)+6.d0*CPOA(4)
     $         +4.d0*CPOA(5)+2.d0*CPOA(6))
      F(2)=FAC*(CPOA(1)+2./3.d0*CPOA(2)+2.d0*CPOA(4)
     $         +3./2.d0*CPOA(5)+CPOA(6)+1./2.d0*CPOA(7))
      F(3)=FAC*(-2.d0*CPOA(1)+4./3.d0*CPOA(2)+16./3.d0*CPOA(3)
     $         -CPOA(4)+CPOA(6)+2.d0*CPOA(7)+4./3.d0*CPOA(9)
     $         +3.d0*CPOA(11))
      F(4)=FAC*(2./3.d0*CPOA(4)-113./36.d0*CPOA(5)-2.d0*CPOA(6)
     $         -113./36.d0*CPOA(7))    
      F(5)=FAC*(2.d0*CPOA(4)+137./18.d0*CPOA(5)+2.d0*CPOA(6)
     $         +89./18.*CPOA(7))
      F(6)=FAC*(-2.d0*CPOA(4)-113./36.d0*CPOA(5)+2./3.d0*CPOA(6)
     $         -113./36.d0*CPOA(7))
      F(7)=FAC*(-2.d0*CPOA(4)-4./3.d0*CPOA(5)-2.d0*CPOA(6)
     $          +4./3.d0*CPOA(7))
      F(8)=FAC*(9./4.d0*CPOA(5)+9./4.d0*CPOA(7)
     $         +2.d0*CPOA(12))
      F(9)=0.d0
      F(10)=FAC*(-2.d0*BB*CPOA(10))
      F(11)=FAC*(16./3.d0-2.d0*BB)*CPOA(11)
      F(12)=FAC*(-2.d0*BB*CPOA(12))   
      
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE GAMMAC2(T,CPOB,F)            
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,CPOB(17),F(17)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/         
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX             
      SAVE /CHRGN/
      COMMON /SCALEQ/Q
      REAL*8 Q
      SAVE /SCALEQ/
      COMMON /BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      INTEGER EN,NS,NF,NFT
      REAL*8 PI,FAC

      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)
      
      EN=3
      NS=0
      IF(MSQU(1).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(2).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(3).LT.Q) THEN
       NS=NS+10
      ENDIF             
      IF(Q.GT.MT) THEN 
       NF=6
       NFT=3-2*3
      ELSE
       NF=5
       NFT=3-2*2
      ENDIF
C
      F(1)=FAC*(20./3.d0*CPOB(1)-8.d0*CPOB(2)+6.d0*CPOB(4)
     $         +4.d0*CPOB(5)+2.d0*CPOB(6))
      F(2)=FAC*(CPOB(1)+2./3.d0*CPOB(2)+2.*CPOB(4)
     $         +3./2.d0*CPOB(5)+CPOB(6)+1./2.d0*CPOB(7)
     $         +(EN**2-2)/8.d0/EN*CPOB(10)+1./4.d0*CPOB(11)
     $         +(EN**2-2)/8.d0/EN*CPOB(12)
     $         -((EN**2-2)/16.d0/EN+(EN**2+2)/72.d0/EN)*CPOB(13)
     $         -((EN**2-2)/16.d0/EN-(EN**2+2)/72.d0/EN)*CPOB(14)
     $         +1./4.d0*CPOB(15)-1./8.d0*CPOB(16)-1./8.d0*CPOB(17))
      F(3)=FAC*(-2.d0*CPOB(1)+4./3.d0*CPOB(2)+16./3.d0*CPOB(3)
     $         -CPOB(4)+CPOB(6)+2.*CPOB(7)+4./3.d0*CPOB(9)
     $         -2.*(EN**2-1)/4.d0/EN*CPOB(10)
     $         -2.*(EN**2-1)/4.d0/EN*CPOB(12)
     $         +10./18.*(EN**2-1)/2./EN*CPOB(13)
     $         +8./18.*(EN**2-1)/2./EN*CPOB(14))
      F(4)=FAC*(2./3.d0*CPOB(4)-113./36.d0*CPOB(5)-2.d0*CPOB(6)
     $          -113./36.d0*CPOB(7))    
      F(5)=FAC*(2.d0*CPOB(4)+137./18.d0*CPOB(5)+2.d0*CPOB(6)
     $         +89./18.d0*CPOB(7))
      F(6)=FAC*(-2.d0*CPOB(4)-113./36.d0*CPOB(5)+2./3.d0*CPOB(6)
     $          -113./36.d0*CPOB(7))
      F(7)=FAC*(-2.*CPOB(4)-4./3.d0*CPOB(5)-2.*CPOB(6)
     $          +4./3.d0*CPOB(7))
      F(8)=FAC*(9./4.d0*CPOB(5)+9./4.d0*CPOB(7))
      F(9)=0.
      F(10)=FAC*(1./2.d0*CPOB(10)+9./2.d0*CPOB(12)-1./2.d0*CPOB(13)
     $           -4.d0*CPOB(14)+3./2.d0*CPOB(16)-3./2.d0*CPOB(17))
      F(11)=FAC*(-3./2.d0*CPOB(10)-4.d0*CPOB(11)-3./2.d0*CPOB(12)
     $          +3./2.d0*CPOB(13)-1./2.d0*CPOB(16)+1./2.d0*CPOB(17))
      F(12)=FAC*(-4.d0*CPOB(12))
      F(13)=FAC*(-25./6.d0*CPOB(13)+1./6.d0*CPOB(14)
     $          +1./3.d0*SMIX(2)+2./3.d0*SMIX(3) 
     $          +NF/3.d0*SMIX(4)+NF/3.d0*SMIX(6)
     $          +1./12.d0*ABMIX(1)-1./12.d0*ABMIX(2))   
      F(14)=FAC*(-25./6.d0*CPOB(14)+1./6.d0*CPOB(13)   
     $          -1./3.d0*SMIX(2)-2./3.d0*SMIX(3) 
     $          -NF/3.d0*SMIX(4)-NF/3.d0*SMIX(6)   
     $          -1./12.d0*ABMIX(1)+1./12.d0*ABMIX(2))   
      F(15)=FAC*(-4.d0*CPOB(15))
      F(16)=FAC*(13./18.d0*CPOB(13)-13./18.d0*CPOB(14)
     $          -2.d0*CPOB(16)-2.d0*CPOB(17)
     $          -1./9.d0*SMIX(2)-2./9.d0*SMIX(3) 
     $          -NF/9.d0*SMIX(4)-NF/9.d0*SMIX(6)
     $          -1./36.d0*ABMIX(1)+1./36.d0*ABMIX(2))   
         F(17)=FAC*(-13./18.d0*CPOB(13)+13./18.d0*CPOB(14)
     $          -2.d0*CPOB(16)-2.d0*CPOB(17)
     $          +1./9.d0*SMIX(2)+2./9.d0*SMIX(3) 
     $          +NF/9.d0*SMIX(4)+NF/9.d0*SMIX(6)   
     $          +1./36.d0*ABMIX(1)-1./36.d0*ABMIX(2))   
      
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE GAMMAWB1(T,CI1,F)            
C-----------------------------------------------------------------
C
C     Right hand side of 1-loop RGEs for evolution of 
C     Wilson coefficients CI1 below M_W scale
C          dCI1_i/dT = F_i(CI1) ~ \gamma_ji * CI1_j
C
      IMPLICIT NONE
      REAL*8 T,CI1(8),F(8)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/
      INTEGER NF,NU,ND
      REAL*8 PI,FAC
      
      NF=5
      NU=2
      ND=3
C     
      PI=4.*ATAN(1.d0)
      FAC=G3**2/(16.d0*PI**2)
C
      F(1)=FAC*(-2.d0*CI1(1)+6.d0*CI1(2))       
      F(2)=FAC*(6.d0*CI1(1)-2.d0*CI1(2))
      F(3)=FAC*(-2./9.d0*CI1(2)-22./9.d0*CI1(3)
     $         -2./9.d0*NF*CI1(4)-2./9.d0*NF*CI1(6))
      F(4)=FAC*(2./3.d0*CI1(2)+22./3.d0*CI1(3)
     $          +(-2.d0+2./3.d0*NF)*CI1(4)+2./3.d0*NF*CI1(6))    
      F(5)=FAC*(-2./9.d0*CI1(2)-4./9.d0*CI1(3)
     $         -2./9.d0*NF*CI1(4)+2.d0*CI1(5)-2./9.d0*NF*CI1(6))
      F(6)=FAC*(2./3.d0*CI1(2)+4./3.d0*CI1(3)+2./3.d0*NF*CI1(4)
     $          -6.d0*CI1(5)+(-16.d0+2./3.d0*NF)*CI1(6))
      F(7)=FAC*(416./81.d0*CI1(2)-464./81.d0*CI1(3)
     $          +(416./81.d0*NU-232./81.d0*ND)*CI1(4)
     $          +32./9.d0*CI1(5)
     $          +(-448./81.d0*NU+200./81.d0*ND)*CI1(6)
     $          +32./3.d0*CI1(7)-32./9.d0*CI1(8))
      F(8)=FAC*(3.d0*CI1(1)+70./27.d0*CI1(2)
     $          +(140./27.d0+3.d0*NF)*CI1(3)
     $          +(6.d0+70./27.d0*NF)*CI1(4)
     $          +(-14./3.d0-3.d0*NF)*CI1(5)
     $          +(-4.d0-119./27.d0*NF)*CI1(6)
     $          +28./3.d0*CI1(8))        
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE GAMMAWB2(T,CI2,F)            
C-----------------------------------------------------------------
C
C     Right hand side of 2-loop RGEs for evolution of 
C     Wilson coefficients CI2 below M_W scale
C          dCI2_i/dT = F_i(CI2) ~ \gamma_ji * CI2_j
C
C
      IMPLICIT NONE
      REAL*8 T,CI2(8),F(8)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/
      REAL*8 PI,FAC1,FAC2
      INTEGER EN,NF,NU,ND,PAR
      
      EN=3
      NF=5
      NU=2
      ND=3
      par=0
C     
      PI=4.d0*ATAN(1.d0)
      FAC1=G3**2/(16.d0*PI**2)
      FAC2=G3**4/(16.d0*PI**2)**2
C
      F(1)=FAC1*(-2.d0*CI2(1)+6.d0*CI2(2))       
     $    +FAC2*((-22./3.d0-57./EN**2/2.d0-2./EN*NF/3.d0)*CI2(1)
     $          +(-19.*EN/6.d0+39./EN+2./3.d0*NF)*CI2(2))
      F(2)=FAC1*(6.d0*CI2(1)-2.d0*CI2(2))
     $    +FAC2*((-19.*EN/6.d0+39.d0/EN+2.*NF/3.d0)*CI2(1)
     $           +(-22./3.d0-57./EN**2/2.d0-3./EN*NF/3.d0)*CI2(2))
      F(3)=FAC1*(-2./9.d0*CI2(2)-22./9.d0*CI2(3)
     $         -2./9.d0*NF*CI2(4)-2./9.d0*NF*CI2(6))
     $    +FAC2*((3.d0*EN-2./3.d0/EN)*CI2(1)
     $          +(-32./27.d0+86./27.d0/EN**2)*CI2(2)
     $          +(-262./27.d0-1195./54.d0/EN**2+3.*EN*NF
     $            -10./3.d0/EN*NF)*CI2(3)
     $          +(17.*EN/6.d0+113./3.d0/EN-2./27.d0*EN*NF
     $            +74./27.d0/EN**2*NF)*CI2(4)
     $          +(-3.d0*EN*NF+20./3.d0/EN*NF)*CI2(5)
     $          +(-56./27.d0*NF-178./27.d0/EN**2*NF)*CI2(6))
      F(4)=FAC1*(2./3.d0*CI2(2)+22./3.d0*CI2(3)
     $          +(-2.d0+2./3.d0*NF)*CI2(4)+2./3.d0*NF*CI2(6))    
     $    +FAC2*(-7./3.d0*CI2(1)
     $          +(-176.*EN/27.d0-230./27.d0/EN)*CI2(2)
     $          +(-533.*EN/54.d0+593./27.d0/EN
     $            +1./3.d0*NF)*CI2(3)
     $          +(-12.d0-57./2.d0/EN**2+110./27.d0*NF
     $            -182./27.d0/EN*NF)*CI2(4)
     $          -11./3.d0*NF*CI2(5)
     $          +(-16.*EN/27.d0*NF+250./27.d0/EN*NF)*CI2(6))
      F(5)=FAC1*(-2./9.d0*CI2(2)-4./9.d0*CI2(3)
     $         -2./9.d0*NF*CI2(4)+2.d0*CI2(5)-2./9.d0*NF*CI2(6))
     $    +FAC2*((-3.d0*EN+16./3.d0/EN)*CI2(1)
     $          +(-122./27.d0-94./27.d0/EN**2)*CI2(2)
     $          +(-244./27.d0-188./27.d0/EN**2-3.d0*EN*NF
     $            +10./3.d0/EN*NF)*CI2(3)
     $          +(-6.d0*EN+32./3.d0/EN-56./27.d0*NF
     $            +2./27.d0/EN**2*NF)*CI2(4)
     $          +(137./6.d0+15./2.d0/EN**2+3.d0*EN*NF
     $            -20./3.d0/EN*NF)*CI2(5)
     $          +(-71./2.d0*EN-18.d0/EN+178./27.d0*NF
     $            +74./27.d0/EN**2*NF)*CI2(6))
      F(6)=FAC1*(2./3.d0*CI2(2)+4./3.d0*CI2(3)+2./3.d0*NF*CI2(4)
     $          -6.d0*CI2(5)+(-16.d0+2./3.d0*NF)*CI2(6))
     $    +FAC2*(-7./3.d0*EN*CI2(1)
     $          +(86.*EN/27.d0+130./27.d0/EN)*CI2(2)
     $          +(172.*EN/27.d0+260./27.d0/EN-NF/3.)*CI2(3)
     $          +(-14./3.d0+74.*EN/27.d0*NF
     $            -20./27.d0/EN*NF)*CI2(4)
     $          +(-100.*EN/3.d0+3./EN+11./3.d0*NF)*CI2(5)
     $          +(-203.*EN**2/6.d0+479./6.d0+15./2.d0/EN**2
     $            +200.*EN/27.d0*NF)*CI2(6))
      F(7)=FAC1*(416./81.d0*CI2(2)-464./81.d0*CI2(3)
     $          +(416./81.d0*NU-232./81.d0*ND)*CI2(4)
     $          +32./9.d0*CI2(5)
     $          +(-448./81.d0*NU+200./81.d0*ND)*CI2(6)
     $          +32./3.d0*CI2(7)-32./9.d0*CI2(8))
     $    +FAC2*(4688./27.d0*CI2(7)-2192./81.d0*CI2(8)
     $           +PAR*CI2(2))
      F(8)=FAC1*(3.d0*CI2(1)+70./27.d0*CI2(2)
     $          +(140./27.d0+3.d0*NF)*CI2(3)
     $          +(6.d0+70./27.d0*NF)*CI2(4)
     $          +(-14./3.d0-3.d0*NF)*CI2(5)
     $          +(-4.d0-119./27.d0*NF)*CI2(6)
     $          +28./3.d0*CI2(8))        
     $    +FAC2*(4063./27.d0*CI2(8)
     $           -PAR*CI2(2))
      
      RETURN
      END
C---------------------------------------------------------------------- 
      SUBROUTINE WILSON(GK,HK,
     &                  CSM,CH1,CH2,CC111,CC112,CC113,CC121,CC122,CC123,
     &                  CC211,CC212,CC213,CC221,CC222,CC223,
     &                  DSSM,DSH1,DSH2,DSC111,DSC112,DSC113,DSC121,
     &                  DSC122,DSC123,DSC211,DSC212,DSC213,DSC221,
     &                  DSC222,DSC223)
C----------------------------------------------------------------------
C
C    Calculates Wilson coefficients and contributions from decoupling
C    of heavy particles at the scale of their mass.
C
      IMPLICIT NONE
      REAL*8 GK(2,3),HK(2,3)
      REAL*8 CSM(9),CH1(9),CH2(19),
     $       CC111(12),CC112(12),CC113(12),
     $       CC121(12),CC122(12),CC123(12),
     $       CC211(17),CC212(17),CC213(17),
     $       CC221(17),CC222(17),CC223(17),      
     $       DSSM(9),DSH1(9),DSH2(9),
     $       DSC111(9),DSC112(9),DSC113(9),
     $       DSC121(9),DSC122(9),DSC123(9),
     $       DSC211(9),DSC212(9),DSC213(9),
     $       DSC221(9),DSC222(9),DSC223(9)
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX
      SAVE /CHRGN/        
      REAL*8 X,FN1,FN2,FN3,FN4,SCALE,XWA,XWB,XCB,YPS
	
      CSM(1)=-1./2.d0
      CSM(2)=-1./2.d0
      CSM(3)= 1.d0
      CSM(4)= 11./18.d0
      CSM(5)=-8./9.d0
      CSM(6)= 11./18.d0
      CSM(7)= 1./2.d0
      CSM(8)=-9./4.d0
      CSM(9)= 3./2.d0
c
      X=(MT/MW)**2
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      DSSM(1)=-X*FN4+1./2.d0
      DSSM(2)=-X/2.*(FN3+FN4)+1./2.d0	     
      DSSM(3)=-2.*DSSM(2)
      DSSM(4)= (X+2.d0)/3.*(2.*FN2+FN3+2.*FN4)-11./18.d0	
      DSSM(5)= 2.*(X+2.d0)/3.*(FN2-FN3-2.*FN4)+8./9.d0   
      DSSM(6)= (X+2.d0)/3.*(2.*FN2+FN3+2.*FN4)-11./18.d0	
      DSSM(7)= (X-2.d0)*FN4-1./2.d0   
      DSSM(8)=-3.*((X+2.d0)*(1./2.*FN3+FN4-LOG(X)/6./(X-1.d0))
     $  	   -3./4.d0)	    
      DSSM(9)=-3.*(7./2.d0-2.*X*FN3-5.*X*FN4)
C
      CH1(1)= 1./2.d0
      CH1(2)= 1./2.d0
      CH1(3)=-1.d0
      CH1(4)= 11./18.d0/TANB**2
      CH1(5)=-8./9.d0/TANB**2
      CH1(6)= 11./18.d0/TANB**2
      CH1(7)= 1./2.d0/TANB**2
      CH1(8)=-9./4.d0/TANB**2
      CH1(9)= 3./2.d0/TANB**2	  
c
      X=(MT/MHPLUS)**2
      CH2(1)= 1./2.d0*X
      CH2(2)=-1./2.d0*X
      CH2(3)= X
      CH2(4)=-1./9.d0*X/TANB**2
      CH2(5)= 7./18.d0*X/TANB**2
      CH2(6)=-1./9.d0*X/TANB**2
      CH2(7)= 1./2.d0*X/TANB**2
      CH2(8)= 3./4.d0*X/TANB**2
      CH2(9)= 3./2.d0*X/TANB**2     
      CH2(10)=0.d0
      CH2(12)=0.d0
      CH2(13)=0.d0
      CH2(14)=0.d0
      CH2(15)=0.d0    
      CH2(16)=0.d0
      CH2(18)=0.d0
      CH2(19)=0.d0
c
      X=(MT/MHPLUS)**2
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      DSH1(1)= X*FN4-1./2.d0
      DSH1(2)= X/2.*(FN3+FN4)-1./2.d0	     
      DSH1(3)=-2.*DSH1(2)
      DSH1(4)= (X/3.d0*(2.d0*FN2+FN3+2*FN4)-11./18.d0)/TANB**2        
      DSH1(5)= (2.*X/3.d0*(FN2-FN3-2*FN4)+8./9.d0)/TANB**2   
      DSH1(6)= (X/3.d0*(2.*FN2+FN3+2*FN4)-11./18.d0)/TANB**2	    
      DSH1(7)= (X*FN4-1./2.d0)/TANB**2   
      DSH1(8)=-3.*(X*(1./2.*FN3+FN4-LOG(X)/6./(X-1.d0))
     $  	   -3./4.d0)/TANB**2	    
      DSH1(9)=3.*DSH1(7)
c
      X=(MT/MHPLUS)**2
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      DSH2(1)=X*(FN4-1./2.d0)
      DSH2(2)=X/2.d0*(FN3+FN4+1.d0)	   
     $        +(MT/MHPLUS)**2*LOG(MT/MHPLUS)
      DSH2(3)=-2.*DSH2(2)
c              -2.*(MT/MHPLUS)**2*LOG(MT/MHPLUS)
      DSH2(4)=(X/3.d0*(2*FN2+FN3+2*FN4+1./3.d0))/TANB**2	
      DSH2(5)=(2*X/3.d0*(FN2-FN3-2*FN4)-7./18.d0*X)/TANB**2   
      DSH2(6)=(X/3.d0*(2*FN2+FN3+2*FN4+1./3.d0))/TANB**2	
      DSH2(7)=X*(FN4-1./2.d0)/TANB**2	
      DSH2(8)=-3.*(X*(1./2.d0*FN3+FN4)-LOG(X)/6.d0/(X-1.d0)
     $  	   +1./4.d0*X)/TANB**2        
      DSH2(9)=3.*DSH2(7)		

      IF(MCHA(1).GT.MW) THEN
       SCALE=MCHA(1)
      ELSE
       SCALE=MW
      ENDIF
C
      XWA=(MW/MSQU(1))**2
      XCB=MCHA(1)/MBQ
      CC111(1)=-GK(1,1)*HK(1,1)*XCB*XWA
      CC111(2)= 0.d0
      CC111(3)= 3.*GK(1,1)*HK(1,1)*XCB*XWA
      CC111(4)= 5./18.d0*GK(1,1)*GK(1,1)*XWA
      CC111(5)=-2./9.d0*GK(1,1)*GK(1,1)*XWA
      CC111(6)= 5./18.d0*GK(1,1)*GK(1,1)*XWA
      CC111(7)= 0.d0
      CC111(8)=-3./2.d0*GK(1,1)*GK(1,1)*XWA
      CC111(9)=-3.*GK(1,1)*GK(1,1)*XWA
c
      X=(MCHA(1)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
        DSC111(1)=-GK(1,1)*HK(1,1)*XCB*XWA*(2.*FN4-1.d0)
        DSC111(2)= 0.d0
        DSC111(3)=-3.*GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(1))) 
        DSC111(4)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
        DSC111(5)=2./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
        DSC111(6)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
        DSC111(7)=0.d0
        DSC111(8)=6.*GK(1,1)*GK(1,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(1)))
        DSC111(9)=-6.*GK(1,1)*GK(1,1)*XWA*(FN4-1./2.d0)
      ELSE
        DSC111(1)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC111(1)
        DSC111(2)=0.d0
        DSC111(3)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(1)))-CC111(3) 
        DSC111(4)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC111(4) 
        DSC111(5)=GK(1,1)*GK(1,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC111(5)
        DSC111(6)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC111(6) 
        DSC111(7)=0.d0
        DSC111(8)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(1)))-CC111(8)
        DSC111(9)=GK(1,1)*GK(1,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC111(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(1))**2
      XCB=MCHA(1)/MBQ
      CC211(1)=-GK(1,1)*HK(1,1)*XCB*XWB
      CC211(2)= 0.d0
      CC211(3)=-3.*GK(1,1)*HK(1,1)*XCB*XWB
      CC211(4)=-5./18.d0*GK(1,1)*GK(1,1)*XWB
      CC211(5)= 11./9.d0*GK(1,1)*GK(1,1)*XWB
      CC211(6)=-5./18.d0*GK(1,1)*GK(1,1)*XWB
      CC211(7)= 0.d0
      CC211(8)= 9./2.d0*GK(1,1)*GK(1,1)*XWB
      CC211(9)=-3.*GK(1,1)*GK(1,1)*XWB      
      CC211(10)=-2.*GK(1,1)*HK(1,1)*XCB*XWB	 
      CC211(11)=0.d0
      CC211(12)=2.*GK(1,1)*GK(1,1)*XWB      
      CC211(13)=2.*GK(1,1)*GK(1,1)*XWB      
      CC211(14)=0.d0
      CC211(15)=0.d0   
      CC211(16)=0.d0
      CC211(17)=0.d0
C
      X=(MCHA(1)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC211(1)=-GK(1,1)*HK(1,1)*XCB*XWA*2.*FN4-CC211(1)
       DSC211(2)= 0.d0
       DSC211(3)=-3.*GK(1,1)*HK(1,1)*XCB*XWA*
     $  	    (FN3+FN4)-CC211(3)
       DSC211(4)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC211(4)
       DSC211(5)=2./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC211(5)
       DSC211(6)=1./3.d0*GK(1,1)*GK(1,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC211(6)
       DSC211(7)=0.d0
       DSC211(8)= 6.*GK(1,1)*GK(1,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC211(8)
       DSC211(9)=-6.*GK(1,1)*GK(1,1)*XWA*FN4-CC211(9)		       
      ELSE
       DSC211(1)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC211(1)
       DSC211(2)=0.d0
       DSC211(3)=GK(1,1)*HK(1,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC211(3) 
       DSC211(4)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC211(4) 
       DSC211(5)=GK(1,1)*GK(1,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC211(5)
       DSC211(6)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC211(6) 
       DSC211(7)=0.d0
       DSC211(8)=GK(1,1)*GK(1,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC211(8)
       DSC211(9)=GK(1,1)*GK(1,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC211(9)
      ENDIF
C
      XWA=(MW/MSQU(2))**2
      XCB=MCHA(1)/MBQ
      CC112(1)=-GK(1,2)*HK(1,2)*XCB*XWA
      CC112(2)= 0.d0
      CC112(3)= 3.*GK(1,2)*HK(1,2)*XCB*XWA
      CC112(4)= 5./18.d0*GK(1,2)*GK(1,2)*XWA
      CC112(5)=-2./9.d0*GK(1,2)*GK(1,2)*XWA
      CC112(6)= 5./18.d0*GK(1,2)*GK(1,2)*XWA
      CC112(7)= 0.d0
      CC112(8)=-3./2.d0*GK(1,2)*GK(1,2)*XWA
      CC112(9)=-3.*GK(1,2)*GK(1,2)*XWA         
c
      X=(MCHA(1)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC112(1)=-GK(1,2)*HK(1,2)*XCB*XWA*(2.*FN4-1.d0)
       DSC112(2)= 0.d0
       DSC112(3)=-3.*GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(2)))
       DSC112(4)=1./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC112(5)=2./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC112(6)=1./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC112(7)=0.d0
       DSC112(8)=6.*GK(1,2)*GK(1,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(2)))
       DSC112(9)=-6.*GK(1,2)*GK(1,2)*XWA*(FN4-1./2.d0)
      ELSE
       DSC112(1)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC112(1)
       DSC112(2)=0.d0
       DSC112(3)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(2)))-CC112(3) 
       DSC112(4)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC112(4) 
       DSC112(5)=GK(1,2)*GK(1,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC112(5)
       DSC112(6)=GK(1,1)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC112(6) 
       DSC112(7)=0.d0
       DSC112(8)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(2)))-CC112(8)
       DSC112(9)=GK(1,2)*GK(1,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC112(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(1))**2
      XCB=MCHA(1)/MBQ
      CC212(1)=-GK(1,2)*HK(1,2)*XCB*XWB
      CC212(2)= 0.d0
      CC212(3)=-3.*GK(1,2)*HK(1,2)*XCB*XWB
      CC212(4)=-5./18.d0*GK(1,2)*GK(1,2)*XWB
      CC212(5)= 11./9.d0*GK(1,2)*GK(1,2)*XWB
      CC212(6)=-5./18.d0*GK(1,2)*GK(1,2)*XWB
      CC212(7)= 0.d0
      CC212(8)= 9./2.d0*GK(1,2)*GK(1,2)*XWB
      CC212(9)=-3.*GK(1,2)*GK(1,2)*XWB      
      CC212(10)=-2.*GK(1,2)*HK(1,2)*XCB*XWB	 
      CC212(11)=0.d0
      CC212(12)=2.*GK(1,2)*GK(1,2)*XWB      
      CC212(13)=2.*GK(1,2)*GK(1,2)*XWB      
      CC212(14)=0.d0
      CC212(15)=0.d0	
      CC212(16)=0.d0
      CC212(17)=0.d0
C
      X=(MCHA(1)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC212(1)=-GK(1,2)*HK(1,2)*XCB*XWA*2.*FN4-CC212(1)
       DSC212(2)=0.d0
       DSC212(3)=-3.*GK(1,2)*HK(1,2)*XCB*XWA*
     $  	    (FN3+FN4)-CC212(3)
       DSC212(4)=1./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC212(4)
       DSC212(5)=2./3.d0*GK(1,2)*GK(1,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC212(5)
       DSC212(6)=1./3.*GK(1,2)*GK(1,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC212(6)
       DSC212(7)=0.d0
       DSC212(8)=6.*GK(1,2)*GK(1,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC212(8)
       DSC212(9)=-6.*GK(1,2)*GK(1,2)*XWA*FN4-CC212(9)	       
      ELSE
       DSC212(1)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC212(1)
       DSC212(2)=0.d0
       DSC212(3)=GK(1,2)*HK(1,2)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC212(3) 
       DSC212(4)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC212(4) 
       DSC212(5)=GK(1,2)*GK(1,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC212(5)
       DSC212(6)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC212(6) 
       DSC212(7)=0.d0
       DSC212(8)=GK(1,2)*GK(1,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC212(8)
       DSC212(9)=GK(1,2)*GK(1,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC212(9)
      ENDIF
c
      XWA=(MW/MSQU(3))**2
      XCB=MCHA(1)/MBQ
      CC113(1)=-HK(1,3)*XCB*XWA
      CC113(2)= 0.d0
      CC113(3)= 3.*HK(1,3)*XCB*XWA
      CC113(4)= 5./18.d0*GK(1,3)*XWA
      CC113(5)=-2./9.d0*GK(1,3)*XWA
      CC113(6)= 5./18.d0*GK(1,3)*XWA
      CC113(7)= 0.d0
      CC113(8)=-3./2.d0*GK(1,3)*XWA
      CC113(9)=-3.*GK(1,3)*XWA         
c
      X=(MCHA(1)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC113(1)=-HK(1,3)*XCB*XWA*(2.*FN4-1.d0)
       DSC113(2)=0.d0
       DSC113(3)=-3.*HK(1,3)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(3)))
       DSC113(4)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC113(5)=2./3.d0*GK(1,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC113(6)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC113(7)=0.d0
       DSC113(8)=6.*GK(1,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(3)))
       DSC113(9)=-6.*GK(1,3)*XWA*(FN4-1./2.d0)
      ELSE
       DSC113(1)=HK(1,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC113(1)
       DSC113(2)=0.d0
       DSC113(3)=HK(1,3)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(3)))-CC113(3) 
       DSC113(4)=GK(1,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC113(4) 
       DSC113(5)=GK(1,3)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC113(5)
       DSC113(6)=GK(1,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC113(6) 
       DSC113(7)=0.d0
       DSC113(8)=GK(1,3)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(3)))-CC113(8)
       DSC113(9)=GK(1,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC113(9)
      ENDIF
C       								     
      XWB=(MW/MCHA(1))**2
      XCB=MCHA(1)/MBQ
      CC213(1)=-HK(1,3)*XCB*XWB
      CC213(2)= 0.d0
      CC213(3)=-3.*HK(1,3)*XCB*XWB
      CC213(4)=-5./18.d0*GK(1,3)*XWB
      CC213(5)= 11./9.d0*GK(1,3)*XWB
      CC213(6)=-5./18.d0*GK(1,3)*XWB
      CC213(7)= 0.d0
      CC213(8)= 9./2.d0*GK(1,3)*XWB
      CC213(9)=-3.*GK(1,3)*XWB      
      CC213(10)=-2.*HK(1,3)*XCB*XWB	 
      CC213(11)=0.d0
      CC213(12)=2.*GK(1,3)*XWB      
      CC213(13)=2.*GK(1,3)*XWB      
      CC213(14)=0.d0
      CC213(15)=0.d0	
      CC213(16)=0.d0
      CC213(17)=0.d0
C
      X=(MCHA(1)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC213(1)=-HK(1,3)*XCB*XWA*2.*FN4-CC213(1)
       DSC213(2)=0.d0
       DSC213(3)=-3.*HK(1,3)*XCB*XWA*
     $  	    (FN3+FN4)-CC213(3)
       DSC213(4)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC213(4)
       DSC213(5)=2./3.d0*GK(1,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC213(5)
       DSC213(6)=1./3.d0*GK(1,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC213(6)
       DSC213(7)=0.d0
       DSC213(8)=6.*GK(1,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC213(8)
       DSC213(9)=-6.*GK(1,3)*XWA*FN4-CC213(9)	       
      ELSE
       DSC213(1)=HK(1,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC213(1)
       DSC213(2)=0.d0
       DSC213(3)=HK(1,3)*XCB*XWA*
     $  	 (-3./2.d0+YPS)-CC213(3) 
       DSC213(4)=GK(1,3)*XWA*
     $  	 (1./6.d0-1./20.d0*YPS)-CC213(4) 
       DSC213(5)=GK(1,3)*XWA*
     $  	 (-1./6.d0+1./30.d0*YPS)-CC213(5)
       DSC213(6)=GK(1,3)*XWA*
     $  	 (1./6.d0-1./20.d0*YPS)-CC213(6) 
       DSC213(7)=0.
       DSC213(8)=GK(1,3)*XWA*
     $  	 (1.d0-3./4.d0*YPS)-CC213(8)
       DSC213(9)=GK(1,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC213(9)
      ENDIF
c
      IF(MCHA(2).GT.MW) THEN
       SCALE=MCHA(2)
      ELSE
       SCALE=MW
      ENDIF
C
      XWA=(MW/MSQU(1))**2
      XCB=MCHA(2)/MBQ
      CC121(1)=-GK(2,1)*HK(2,1)*XCB*XWA
      CC121(2)= 0.d0
      CC121(3)= 3.*GK(2,1)*HK(2,1)*XCB*XWA
      CC121(4)= 5./18.d0*GK(2,1)*GK(2,1)*XWA
      CC121(5)=-2./9.d0*GK(2,1)*GK(2,1)*XWA
      CC121(6)= 5./18.d0*GK(2,1)*GK(2,1)*XWA
      CC121(7)= 0.d0
      CC121(8)=-3./2.d0*GK(2,1)*GK(2,1)*XWA
      CC121(9)=-3.*GK(2,1)*GK(2,1)*XWA         
c
      X=(MCHA(2)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC121(1)=-GK(2,1)*HK(2,1)*XCB*XWA*(2.*FN4-1.d0)
       DSC121(2)=0.
       DSC121(3)=-3.*GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(1)))
       DSC121(4)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC121(5)=2./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC121(6)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC121(7)=0.d0
       DSC121(8)=6.*GK(2,1)*GK(2,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(1)))
       DSC121(9)=-6.*GK(2,1)*GK(2,1)*XWA*(FN4-1./2.d0)
      ELSE
       DSC121(1)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC121(1)
       DSC121(2)=0.d0
       DSC121(3)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(1)))-CC121(3) 
       DSC121(4)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC121(4) 
       DSC121(5)=GK(2,1)*GK(2,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC121(5)
       DSC121(6)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC121(6) 
       DSC121(7)=0.d0
       DSC121(8)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(1)))-CC121(8)
       DSC121(9)=GK(2,1)*GK(2,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC121(9)
      ENDIF
C       							       
      XWB=(MW/MCHA(2))**2
      XCB=MCHA(2)/MBQ
      CC221(1)=-GK(2,1)*HK(2,1)*XCB*XWB
      CC221(2)= 0.d0
      CC221(3)=-3.*GK(2,1)*HK(2,1)*XCB*XWB
      CC221(4)=-5./18.d0*GK(2,1)*GK(2,1)*XWB
      CC221(5)= 11./9.d0*GK(2,1)*GK(2,1)*XWB
      CC221(6)=-5./18.d0*GK(2,1)*GK(2,1)*XWB
      CC221(7)= 0.d0
      CC221(8)= 9./2.d0*GK(2,1)*GK(2,1)*XWB
      CC221(9)=-3.*GK(2,1)*GK(2,1)*XWB      
      CC221(10)=-2.*GK(2,1)*HK(2,1)*XCB*XWB	 
      CC221(11)=0.d0
      CC221(12)=2.*GK(2,1)*GK(2,1)*XWB      
      CC221(13)=2.*GK(2,1)*GK(2,1)*XWB      
      CC221(14)=0.d0
      CC221(15)=0.d0	
      CC221(16)=0.d0
      CC221(17)=0.d0
C
      X=(MCHA(2)/MSQU(1))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC221(1)=-GK(2,1)*HK(2,1)*XCB*XWA*2.*FN4-CC221(1)
       DSC221(2)=0.
       DSC221(3)=-3.*GK(2,1)*HK(2,1)*XCB*XWA*
     $  	    (FN3+FN4)-CC221(3)
       DSC221(4)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC221(4)
       DSC221(5)=2./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC221(5)
       DSC221(6)=1./3.d0*GK(2,1)*GK(2,1)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC221(6)
       DSC221(7)=0.d0
       DSC221(8)=6.*GK(2,1)*GK(2,1)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC221(8)
       DSC221(9)=-6.*GK(2,1)*GK(2,1)*XWA*FN4-CC221(9)		       
      ELSE
       DSC221(1)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC221(1)
       DSC221(2)=0.d0
       DSC221(3)=GK(2,1)*HK(2,1)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC221(3) 
       DSC221(4)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC221(4) 
       DSC221(5)=GK(2,1)*GK(2,1)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC221(5)
       DSC221(6)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC221(6) 
       DSC221(7)=0.d0
       DSC221(8)=GK(2,1)*GK(2,1)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC221(8)
       DSC221(9)=GK(2,1)*GK(2,1)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC221(9)
      ENDIF
C
      XWA=(MW/MSQU(2))**2
      XCB=MCHA(2)/MBQ
      CC122(1)=-GK(2,2)*HK(2,2)*XCB*XWA
      CC122(2)= 0.d0
      CC122(3)= 3.*GK(2,2)*HK(2,2)*XCB*XWA
      CC122(4)= 5./18.d0*GK(2,2)*GK(2,2)*XWA
      CC122(5)=-2./9.d0*GK(2,2)*GK(2,2)*XWA
      CC122(6)= 5./18.d0*GK(2,2)*GK(2,2)*XWA
      CC122(7)= 0.d0
      CC122(8)=-3./2.d0*GK(2,2)*GK(2,2)*XWA
      CC122(9)=-3.*GK(2,2)*GK(2,2)*XWA         
c
      X=(MCHA(2)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC122(1)=-GK(2,2)*HK(2,2)*XCB*XWA*(2.*FN4-1.d0)
       DSC122(2)=0.d0
       DSC122(3)=-3.*GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (FN3+FN4+1.+2.*LOG(SCALE/MSQU(2)))
       DSC122(4)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC122(5)=2./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC122(6)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC122(7)=0.d0
       DSC122(8)=6.*GK(2,2)*GK(2,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(2)))
       DSC122(9)=-6.*GK(2,2)*GK(2,2)*XWA*(FN4-1./2.d0)
      ELSE
       DSC122(1)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC122(1)
       DSC122(2)=0.d0
       DSC122(3)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	 (-3./2.d0+YPS
     $  	  -6.*LOG(SCALE/MSQU(2)))-CC122(3) 
       DSC122(4)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC122(4) 
       DSC122(5)=GK(2,2)*GK(2,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC122(5)
       DSC122(6)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC122(6) 
       DSC122(7)=0.d0
       DSC122(8)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(2)))-CC122(8)
       DSC122(9)=GK(2,2)*GK(2,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC122(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(2))**2
      XCB=MCHA(2)/MBQ
      CC222(1)=-GK(2,2)*HK(2,2)*XCB*XWB
      CC222(2)= 0.d0
      CC222(3)=-3.*GK(2,2)*HK(2,2)*XCB*XWB
      CC222(4)=-5./18.d0*GK(2,2)*GK(2,2)*XWB
      CC222(5)= 11./9.d0*GK(2,2)*GK(2,2)*XWB
      CC222(6)=-5./18.d0*GK(2,2)*GK(2,2)*XWB
      CC222(7)= 0.d0
      CC222(8)= 9./2.d0*GK(2,2)*GK(2,2)*XWB
      CC222(9)=-3.*GK(2,2)*GK(2,2)*XWB      
      CC222(10)=-2.*GK(2,2)*HK(2,2)*XCB*XWB	 
      CC222(11)=0.d0
      CC222(12)=2.*GK(2,2)*GK(2,2)*XWB      
      CC222(13)=2.*GK(2,2)*GK(2,2)*XWB      
      CC222(14)=0.d0
      CC222(15)=0.d0	
      CC222(16)=0.d0
      CC222(17)=0.d0
C
      X=(MCHA(2)/MSQU(2))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC222(1)=-GK(2,2)*HK(2,2)*XCB*XWA*2.*FN4-CC222(1)
       DSC222(2)=0.d0
       DSC222(3)=-3.*GK(2,2)*HK(2,2)*XCB*XWA*
     $  	    (FN3+FN4)-CC222(3)
       DSC222(4)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC222(4)
       DSC222(5)=2./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC222(5)
       DSC222(6)=1./3.d0*GK(2,2)*GK(2,2)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC222(6)
       DSC222(7)=0.d0
       DSC222(8)=6.*GK(2,2)*GK(2,2)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC222(8)
       DSC222(9)=-6.*GK(2,2)*GK(2,2)*XWA*FN4-CC222(9)	       
      ELSE
       DSC222(1)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC222(1)
       DSC222(2)=0.d0
       DSC222(3)=GK(2,2)*HK(2,2)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC222(3) 
       DSC222(4)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.-1./20.*YPS)-CC222(4) 
       DSC222(5)=GK(2,2)*GK(2,2)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC222(5)
       DSC222(6)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC222(6) 
       DSC222(7)=0.d0
       DSC222(8)=GK(2,2)*GK(2,2)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC222(8)
       DSC222(9)=GK(2,2)*GK(2,2)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC222(9)
      ENDIF
c
      XWA=(MW/MSQU(3))**2
      XCB=MCHA(2)/MBQ
      CC123(1)=-HK(2,3)*XCB*XWA
      CC123(2)= 0.d0
      CC123(3)= 3.*HK(2,3)*XCB*XWA
      CC123(4)= 5./18.d0*GK(2,3)*XWA
      CC123(5)=-2./9.d0*GK(2,3)*XWA
      CC123(6)= 5./18.d0*GK(2,3)*XWA
      CC123(7)= 0.d0
      CC123(8)=-3./2.d0*GK(2,3)*XWA
      CC123(9)=-3.*GK(2,3)*XWA         
c
      X=(MCHA(2)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4)
      IF(ABS(YPS).GT.0.05) THEN
       DSC123(1)=-HK(2,3)*XCB*XWA*(2.*FN4-1.d0)
       DSC123(2)=0.d0
       DSC123(3)=-3.*HK(2,3)*XCB*XWA*
     $  	   (FN3+FN4+1.d0+2.*LOG(SCALE/MSQU(3)))
       DSC123(4)=1./3.d0*GK(2,3)*XWA*
     $  	   (4.*FN2-FN3+2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC123(5)=2./3.d0*GK(2,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0)+1./3.d0)
       DSC123(6)=1./3.d0*GK(2,3)*XWA*
     $  	  (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0)-5./6.d0)
       DSC123(7)=0.d0
       DSC123(8)=6.*GK(2,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0)+1./4.d0+
     $  	    2./3.d0*LOG(SCALE/MSQU(3)))
       DSC123(9)=-6.*GK(2,3)*XWA*(FN4-1./2.d0)
      ELSE
       DSC123(1)=HK(2,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC123(1)
       DSC123(2)=0.d0
       DSC123(3)=HK(2,3)*XCB*XWA*
     $  	   (-3./2.d0+YPS
     $  	    -6.*LOG(SCALE/MSQU(3)))-CC123(3) 
       DSC123(4)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC123(4) 
       DSC123(5)=GK(2,3)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC123(5)
       DSC123(6)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC123(6) 
       DSC123(7)=0.
       DSC123(8)=GK(2,3)*XWA*
     $  	   (1.d0-3./4.d0*YPS
     $  	   + 4.*LOG(SCALE/MSQU(3)))-CC123(8)
       DSC123(9)=GK(2,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC123(9)
      ENDIF
C       							      
      XWB=(MW/MCHA(2))**2
      XCB=MCHA(2)/MBQ
      CC223(1)=-HK(2,3)*XCB*XWB
      CC223(2)=0.d0
      CC223(3)=-3.*HK(2,3)*XCB*XWB
      CC223(4)=-5./18.d0*GK(2,3)*XWB
      CC223(5)=11./9.d0*GK(2,3)*XWB
      CC223(6)=-5./18.d0*GK(2,3)*XWB
      CC223(7)=0.d0
      CC223(8)=9./2.d0*GK(2,3)*XWB
      CC223(9)=-3.*GK(2,3)*XWB      
      CC223(10)=-2.*HK(2,3)*XCB*XWB	 
      CC223(11)=0.d0
      CC223(12)=2.*GK(2,3)*XWB      
      CC223(13)=2.*GK(2,3)*XWB      
      CC223(14)=0.d0
      CC223(15)=0.d0	
      CC223(16)=0.d0
      CC223(17)=0.d0
C
      X=(MCHA(2)/MSQU(3))**2
      YPS=X-1.d0
      CALL FUNS(X,FN1,FN2,FN3,FN4) 
      IF(ABS(YPS).GT.0.05) THEN
       DSC223(1)=-HK(2,3)*XCB*XWA*2.*FN4-CC223(1)
       DSC223(2)=0.d0
       DSC223(3)=-3.*HK(2,3)*XCB*XWA*
     $  	    (FN3+FN4)-CC223(3)
       DSC223(4)=1./3.d0*GK(2,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC223(4)
       DSC223(5)=2./3.d0*GK(2,3)*XWA*
     $  	   (2.*FN2+FN3+2.*FN4-LOG(X)/(X-1.d0))-CC223(5)
       DSC223(6)=1./3.d0*GK(2,3)*XWA*
     $  	   (4.*FN2-FN3-2.*FN4+LOG(X)/(X-1.d0))-CC223(6)
       DSC223(7)=0.d0
       DSC223(8)=6.*GK(2,3)*XWA*
     $  	   (1./2.d0*FN3+FN4-LOG(X)/6./(X-1.d0))-CC223(8)
       DSC223(9)=-6.*GK(2,3)*XWA*FN4-CC223(9)	       
      ELSE
       DSC223(1)=HK(2,3)*XCB*XWA*
     $  	   (-1./3.d0+1./6.d0*YPS)-CC223(1)
       DSC223(2)=0.d0
       DSC223(3)=HK(2,3)*XCB*XWA*
     $  	   (-3./2.d0+YPS)-CC223(3) 
       DSC223(4)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC223(4) 
       DSC223(5)=GK(2,3)*XWA*
     $  	   (-1./6.d0+1./30.d0*YPS)-CC223(5)
       DSC223(6)=GK(2,3)*XWA*
     $  	   (1./6.d0-1./20.d0*YPS)-CC223(6) 
       DSC223(7)=0.
       DSC223(8)=GK(2,3)*XWA*
     $  	   (1.d0-3./4.d0*YPS)-CC223(8)
       DSC223(9)=GK(2,3)*XWA*
     $  	 (-1.d0+1./2.d0*YPS)-CC223(9)
      ENDIF

      RETURN
      END
C---------------------------------------------------------------------- 
      SUBROUTINE CHARGINO(GK,HK)
C----------------------------------------------------------------------
C     Computes chargino-quark-squark interaction matrices as given by
C     Eq(50) of Anlauf.
C
      IMPLICIT NONE
      REAL*8 GK(2,3),HK(2,3)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX        
      SAVE /CHRGN/        
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
C
      REAL*8 BETA,SB,CB,ALT,ALB,UM(2,2),VM(2,2),TM(2,2)
      
      BETA=ATAN(TANB)
      SB=SIN(BETA)
      CB=COS(BETA)       
      ALT=MTQ/(SQRT(2.)*MW*SB)
      ALB=MBQ/(SQRT(2.)*MW*CB)
c
c  Convert chargino mixing matrices to H-K
c
      VM(1,1) =-SIN(GAMMAR)                       *SIGN(1.d0,AMW1SS)
      VM(1,2) = COS(GAMMAR)                       *SIGN(1.d0,AMW1SS)
      VM(2,1) =-COS(GAMMAR)*SIGN(1.d0,TAN(GAMMAR))*SIGN(1.d0,AMW2SS)
      VM(2,2) =-SIN(GAMMAR)*SIGN(1.d0,TAN(GAMMAR))*SIGN(1.d0,AMW2SS) 
      UM(1,1) =-SIN(GAMMAL)                       
      UM(1,2) = COS(GAMMAL)                       
      UM(2,1) =-COS(GAMMAL)*SIGN(1.d0,TAN(GAMMAL))
      UM(2,2) =-SIN(GAMMAL)*SIGN(1.d0,TAN(GAMMAL))

      MCHA(1) = ABS(AMW1SS)
      MCHA(2) = ABS(AMW2SS)
c
c  Convert stop mixing matrix to H-K
c
      TM(1,1) = COS(THETAT) 
      TM(1,2) =-SIN(THETAT)
      TM(2,1) = SIN(THETAT)
      TM(2,2) = COS(THETAT)
      
      MSQU(1)=MSTP1
      MSQU(2)=MSTP2
      MSQU(3)=1./2.d0*(MSCHL+MSCHR) 
C....
      GK(1,1)= VM(1,1)*TM(1,1)-ALT*VM(1,2)*TM(1,2)
      GK(1,2)= VM(1,1)*TM(2,1)-ALT*VM(1,2)*TM(2,2)
      GK(1,3)= (ALT*VM(1,2))**2-(GK(1,1))**2-(GK(1,2))**2
      GK(2,1)= VM(2,1)*TM(1,1)-ALT*VM(2,2)*TM(1,2)
      GK(2,2)= VM(2,1)*TM(2,1)-ALT*VM(2,2)*TM(2,2)
      GK(2,3)= (ALT*VM(2,2))**2-(GK(2,1))**2-(GK(2,2))**2
      
      HK(1,1)= ALB*UM(1,2)*TM(1,1)
      HK(1,2)= ALB*UM(1,2)*TM(2,1)
      HK(1,3)=-GK(1,1)*HK(1,1)-GK(1,2)*HK(1,2)
      HK(2,1)= ALB*UM(2,2)*TM(1,1)
      HK(2,2)= ALB*UM(2,2)*TM(2,1)
      HK(2,3)=-GK(2,1)*HK(2,1)-GK(2,2)*HK(2,2)
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE TACTIV(IQROT,FTMT,G)
C----------------------------------------------------------------------
C
C    Activates 3rd generation up quark in Yukawa matrix and adjusts
C    RGE-evolving vector G.
C
C    Created: 5/24/07 by Azar Mustafayev 
C
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IQROT
      REAL FTMT
      REAL*8 G(157)
#include "sslun.inc"
      REAL*8 YU(3,3),YUH(3,3),xYU(3,3),WORK(3),A(3,3),W(3),Z(3,3),TEMP,
     &       ZMAX,I3(3,3),VL(3,3),VR(3,3)
c  A(3,3)  - axiliary matrix to be diagonalized
c  W(3)    - vector of eigenvalues for AR
c  Z(3,3)  - matrix of eigenvalues of AR as columns
      INTEGER I,J,K,L,IERR
      DATA I3/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix
      
      IF (IQROT.EQ.0.OR.IQROT.EQ.1) THEN
        G(12)=DBLE(FTMT)
      ELSE 
c...assemble up Yukawa matrix
        CALL VEC2MAT(G,4,YU,-1)
        DO I=1,3
          DO J=1,3
	    YUH(I,J)=YU(J,I)
          ENDDO
        ENDDO
 
c...Compute left rotation matrix VL by diagonalizing  YUK*YUKH     
        DO I=1,3
        DO J=1,3
	  A(I,J)=0.d0
          DO K=1,3
	    A(I,J)=A(I,J)+YU(I,K)*YUH(K,J)
          ENDDO
        ENDDO
        ENDDO
        CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
        IF (IERR.NE.0) THEN
          WRITE(LOUT,*) 'EISRS1 ERROR IN TACTIV, IERR=',IERR
          STOP99
        ENDIF
c...Restore generation structure in W and Z
        DO I=1,3
	  DO J=I,3
	    IF (ABS(Z(I,I)).LT.ABS(Z(I,J)).AND.J.NE.I) THEN
	      DO K=1,3
		TEMP=Z(K,J)	! swap columns
		Z(K,J)=Z(K,I)
		Z(K,I)=TEMP
	      ENDDO
	      TEMP=W(J)        ! swap eigenvalues
	      W(J)=W(I)
	      W(I)=TEMP
	    ENDIF
	  ENDDO
        ENDDO
c...save left rotation matrix
        DO I=1,3
          DO J=1,3
	    VL(I,J)=Z(J,I)
          ENDDO
        ENDDO
c...check if VL has all diagonal elements positive
        DO I=1,3
          IF (VL(I,I).LT.0.d0) THEN      
	    DO J=1,3
	      VL(I,J)=-VL(I,J)   ! sign change in the row
	    ENDDO
          ENDIF
        ENDDO

c...Compute right rotation matrix VR by diagonalizing  YUKH*YUK     
        DO I=1,3
          DO J=1,3
	    A(I,J)=0.d0
            DO K=1,3
	      A(I,J)=A(I,J)+YUH(I,K)*YU(K,J)
            ENDDO
          ENDDO
        ENDDO
        CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
        IF (IERR.NE.0) THEN
          WRITE(LOUT,*) 'EISRS1 ERROR IN TACTIV, IERR=',IERR
          STOP99
        ENDIF
c...Restore generation structure in W and Z
        DO I=1,3
	  DO J=I,3
	    IF (ABS(Z(I,I)).LT.ABS(Z(I,J)).AND.J.NE.I) THEN
	      DO K=1,3
		TEMP=Z(K,J)	! swap columns
		Z(K,J)=Z(K,I)
		Z(K,I)=TEMP
	      ENDDO
	      TEMP=W(J)        ! swap eigenvalues
	      W(J)=W(I)
	      W(I)=TEMP
	    ENDIF
	  ENDDO
        ENDDO
c...save right rotation matrix
        DO I=1,3
          DO J=1,3
	    VR(I,J)=Z(J,I)
          ENDDO
        ENDDO
c...check if VR has all diagonal elements positive
        DO I=1,3
          IF (VR(I,I).LT.0.d0) THEN      
	    DO J=1,3
	      VR(I,J)=-VR(I,J)   ! sign change in the row
	    ENDDO
          ENDIF
        ENDDO

c...build diagonal Yukawa matrix
        DO I=1,3
          DO J=1,3
	    xYU(I,J)=SQRT(DABS(W(I)))*I3(I,J)
          ENDDO
        ENDDO

c...make the adjustment 
        xYU(3,3)=xYU(3,3)+DBLE(FTMT)
c...rotate back to original basis
        DO I=1,3
        DO J=1,3
          YU(I,J)=0.d0
          DO K=1,3
            DO L=1,3
              YU(I,J)=YU(I,J)+VL(K,I)*xYU(K,L)*VR(L,J)
            ENDDO
          ENDDO
        ENDDO
        ENDDO

        CALL MAT2VEC(G,4,YU,-1)
      ENDIF
      
      
      RETURN
      END
C---------------------------------------------------------------------- 
      SUBROUTINE GLUNENO(G,GAM,AGE,AHE)
C----------------------------------------------------------------------
c    Computes 6x6 down squark diagonalization matrix and 
C    neutralino-(d)quark-(d)squark coupling matrices.
C
c    Ref: H.Anlauf  hep-ph/9406286;
c         S.Bertolini et al  NPB353, 591 (1991).
c
c    Created: H.Baer and M.Brhlik
C 
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C   Input:
C     G(i) - vector of RGE evolved parameters
C   Output:
C     GAM       - 6x6 down squark rotation matrix
C     AGE, AHE  - 4x6x6 neutralino-quark-squark matrices 
c                 defined by Eq(C.9) of Bertolini et al.
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv            
      IMPLICIT NONE
      REAL*8 G(157)
      REAL*8 GAM(6,6),AGE(4,6,6),AHE(4,6,6)
      COMMON /GLNN/ MCH0(4),MDSB(6)
      REAL*8 MCH0,MDSB
      SAVE /GLNN/
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON /BSGSM/MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
      REAL*8 DM(6,6),ESQ(3,3),ESD(3,3),EMD(3,3),EAD(3,3),UNI(3,3),
     &       EDY(3,3),XORK(6),ANEU(4,4),BNEU(4,4)
c  ANEU(4,4) - neutralino mixing matrix on (photino,zino,higgs1,higgs2) basis
c              as defined by Haber-Kane
c  BNEU(4,4) - ANEU up to sign change in rows with negative mass eigenvalues
c
      REAL*8 BETA,SB,CB,S2B,C2B,VEV,SW,CW
      REAL*8 YDdiag(3,3),VDL(3,3),VDR(3,3),TD(3,3),M2Q(3,3),M2D(3,3)
      INTEGER I,J,K,L,IERR
      DATA UNI/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix
      
      BETA=ATAN(TANB)
      SB=SIN(BETA)
      CB=COS(BETA)       
      S2B=SIN(2.*BETA)
      C2B=COS(2.*BETA) 
      VEV=SQRT(V**2+VP**2)          
      SW=SQRT(XW)      
      CW=SQRT(1.d0-XW)
c
c  Convert neutralino mixing matrix to H-K
c
      DO I=1,4
        ANEU(I,1)= ZMIXSS(4,I)
        ANEU(I,2)= ZMIXSS(3,I)
        ANEU(I,3)=-ZMIXSS(2,I)
        ANEU(I,4)=-ZMIXSS(1,I)
	MCH0(I)=ABS(AMZISS(I))
      ENDDO
      DO I=1,4
        DO J=1,4
          BNEU(I,J)=ANEU(I,J)*SIGN(1.d0,-AMZISS(I))
        ENDDO
      ENDDO
c
c  6x6 d-squark mass matrix
c
      CALL YUKDIAG(G,13,YDdiag,VDL,VDR)
      K=0
      DO I=1,3
        DO J=1,3
	  TD(J,I)=G(43+K)
	  M2Q(I,J)=G(65+K)
	  M2D(I,J)=G(92+K)
	  K=K+1
        ENDDO
      ENDDO
c...Rotate matrices to superCKM basis      
      DO I=1,3
        DO J=1,3
          EMD(I,J)=YDdiag(I,J)*VEV*CB
	  EAD(I,J)=0.d0
	  ESQ(I,J)=0.d0
	  ESD(I,J)=0.d0
	  DO K=1,3
	    DO L=1,3
	      EAD(I,J)=EAD(I,J)+VDR(I,K)*TD(L,K)*VDL(J,L)*CB*VEV
              ESQ(I,J)=ESQ(I,J)+VDL(I,K)*M2Q(K,L)*VDL(J,L)
              ESD(I,J)=ESD(I,J)+VDR(I,K)*M2D(K,L)*VDR(J,L)
	    ENDDO
	  ENDDO
	ENDDO
      ENDDO
c...Coupling lambda given by Eq(47) of Anlauf
      EDY(1,1)=EMD(1,1)/(SQRT(2.d0)*CB*MW)
      EDY(1,2)=0.d0
      EDY(1,3)=0.d0
      EDY(2,1)=0.d0
      EDY(2,2)=EMD(2,2)/(SQRT(2.d0)*CB*MW)
      EDY(2,3)=0.d0
      EDY(3,1)=0.d0
      EDY(3,2)=0.d0
      EDY(3,3)=EMD(3,3)/(SQRT(2.d0)*CB*MW)
C...Construct 6x6 down squark mass^2 matrix
      DO I=1,3
      DO J=1,3
        DM(I,J)=0.d0
        DO K=1,3
          DM(I,J)=DM(I,J)+ESQ(I,K)*UNI(K,J)+EMD(K,I)*EMD(K,J)-
     $      MZ**2*(1./2.d0-1./3.d0*XW)*C2B*UNI(I,K)*UNI(K,J)
        ENDDO
c	DM(I,J+3)=EAD(I,J)
c	DO K=1,3
c	  DM(I,J+3)=DM(I,J+3)-EMD(I,K)*MU*TANB*UNI(K,J)
c	ENDDO
c	DM(I+3,J)=EAD(J,I)
c	DO K=1,3
c	  DM(I+3,J)=DM(I+3,J)-MU*TANB*UNI(I,K)*EMD(J,K)
c	ENDDO
c sign in trilinear and mu fixed
	DM(I,J+3)=-EAD(I,J)
	DO K=1,3
	  DM(I,J+3)=DM(I,J+3)+EMD(I,K)*MU*TANB*UNI(K,J)
	ENDDO
	DM(I+3,J)=-EAD(J,I)
	DO K=1,3
	  DM(I+3,J)=DM(I+3,J)+MU*TANB*UNI(I,K)*EMD(J,K)
	ENDDO
        DM(I+3,J+3)=0.d0
        DO K=1,3
          DM(I+3,J+3)=DM(I+3,J+3)+ESD(I,K)*UNI(K,J)+EMD(K,I)*EMD(K,J)
     $  	     +MZ**2*1./3.d0*XW*C2B*UNI(I,K)*UNI(K,J)
        ENDDO
      ENDDO
      ENDDO
c...Diagonalization
      CALL EIGSYS(6,6,DM,MDSB,GAM,IERR,XORK)
      IF (IERR.NE.0) THEN
        PRINT*, 'EISRS1 ERROR IN GLUNENO, IERR=',IERR
        STOP99
      END IF
       
      DO I=1,6
        MDSB(I)=SQRT(abs(MDSB(I)))
      ENDDO
c
c  Neutralino-quark-squark coupling matrices
c
      DO I=1,4
      DO J=1,3
      DO K=1,6
        AGE(I,J,K)=((-1./3.d0)*SW*(ANEU(I,1)*CW+ANEU(I,2)*SW)
     $              +1./CW*(-1./2.d0+1./3.d0*XW)
     $               *(-ANEU(I,1)*SW+ANEU(I,2)*CW))*GAM(J,K)
        AGE(I,J+3,K)=((-1./3.d0)*SW*(BNEU(I,1)*CW+BNEU(I,2)*SW)
     $                +1./CW*(1./3.*XW)
     $                 *(-BNEU(I,1)*SW+BNEU(I,2)*CW))*GAM(J+3,K)
        AHE(I,J,K)=0.d0
        AHE(I,J+3,K)=0.d0
        DO L=1,3
          AHE(I,J,K)=AHE(I,J,K)+BNEU(I,3)*EDY(J,L)*GAM(L,K)
          AHE(I,J+3,K)=AHE(I,J,K)+ANEU(I,3)*EDY(J,L)*GAM(L+3,K)
        ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE EIGSYS(NM,N,AR,WR,ZR,IERR,WORK)
C----------------------------------------------------------------------
C     Computes all eigenvalues and eigenvectors of real symmetric
C     NM x NM matrix AR.
C     This is a DOUBLE PRECISION version of the EISRS1 subroutine (F224)
C     from CERN program library.
C
C     Note: Full input matrix AR is preserved.
C     Note: Convention  ZR^T AR ZR = diagonal  is used.
C     Note: Eigenvalues are arranged in non-decreasing order.
C
C     Created: 02/23/07 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER NM,N,IERR
      REAL*8 AR(N,N),WR(N),ZR(N,N),WORK(N)
C     WR - vector containing eigenvalues
C     ZR - matrix containing eigenvectors as columns
C     IERR - error parameter, if non-zero the computation has failed
C
      CALL TRDIAG(NM,N,AR,WR,WORK,ZR)
      CALL TQLEIG(NM,N,WR,WORK,ZR,IERR)
c      CALL IMTQLEIG(NM,N,WR,WORK,ZR,IERR) 
      
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE TRDIAG(NM,N,A,D,E,Z)
C----------------------------------------------------------------------
C     Reduces real symmetric matrix A to trigiagonal form.
C     This is a DOUBLE PRECISION version of the TRED2 subroutine (F220)
C     from CERN program library.
C
C     Created: 02/23/07 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      REAL*8 F,G,H,HH,SCALE
      
      DO 100 I = 1, N
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
100   CONTINUE
      IF (N .EQ. 1) GO TO 320
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.D0
         SCALE = 0.D0
         IF (L .LT. 2) GO TO 130
         DO 120 K = 1, L
           SCALE = SCALE + DABS(Z(I,K))
120      CONTINUE
         IF (SCALE .NE. 0.D0) GO TO 140
130      E(I) = Z(I,L)
         GO TO 290
140      DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
150      CONTINUE
         F = Z(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.D0
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.D0
            DO 180 K = 1, J
              G = G + Z(J,K) * Z(I,K)
180         CONTINUE
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
              G = G + Z(K,J) * Z(I,K)
200         CONTINUE
220         E(J) = G / H
            F = F + E(J) * Z(I,J)
240      CONTINUE
         HH = F / (H + H)
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
260      CONTINUE
         DO 280 K = 1, L
           Z(I,K) = SCALE * Z(I,K)
280      CONTINUE
290        D(I) = H
300   CONTINUE
320   D(1) = 0.D0
      E(1) = 0.D0
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.D0) GO TO 380
         DO 360 J = 1, L
            G = 0.D0
            DO 340 K = 1, L
              G = G + Z(I,K) * Z(K,J)
340         CONTINUE
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
360      CONTINUE
380      D(I) = Z(I,I)
         Z(I,I) = 1.D0
         IF (L .LT. 1) GO TO 500
         DO 400 J = 1, L
            Z(I,J) = 0.D0
            Z(J,I) = 0.D0
400      CONTINUE
500   CONTINUE
      
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE TQLEIG(NM,N,D,E,Z,IERR)
C----------------------------------------------------------------------
C     Computes eigenvalues and eigenvectors of the symmetric tridiagonal 
C     matrix using QL-algorithm.
C     
C     This is a DOUBLE PRECISION version of the TQL2 subroutine (F220)
C     from CERN program library.
C
C     Created: 02/23/07 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      
      MACHEP=2.D0**(-23)
      MACHEP=2.D0**(-47)
      IERR = 0
      
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
        E(I-1) = E(I)
100   CONTINUE
      F = 0.D0
      B = 0.D0
      E(N) = 0.D0
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))
         IF (B .LT. H) B = H
         DO 110 M = L, N
            IF (DABS(E(M)) .LE. B) GO TO 120
110      CONTINUE
120      IF (M .EQ. L) GO TO 220
130      IF (J .EQ. 30) GO TO 1000
         J = J + 1
         P = (D(L+1) - D(L)) / (2.D0 * E(L))
         R = DSQRT(P*P+1.D0)
         H = D(L) - E(L) / (P + DSIGN(R,P))
         DO 140 I = L, N
           D(I) = D(I) - H
140      CONTINUE
         F = F + H
         P = D(M)
         C = 1.D0
         S = 0.D0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150
            C = E(I) / P
            R = DSQRT(C*C+1.D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.D0 / R
            GO TO 160
150         C = P / E(I)
            R = DSQRT(C*C+1.D0)
            E(I+1) = S * E(I) * R
            S = 1.D0 / R
            C = C * S
160         P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
180         CONTINUE
200      CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (DABS(E(L)) .GT. B) GO TO 130
220      D(L) = D(L) + F
240   CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
260      CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
280      CONTINUE
300   CONTINUE
      GO TO 1001
1000  IERR = L

1001  RETURN
      END
C---------------------------------------------------------------------
      SUBROUTINE SUG2BSG            
C----------------------------------------------------------------------
C     SUGRA - BSG interface.
C     Fills BSG common blocks from SUGRA ones.
C
C    Note: some of the parameters are transferred to the main IsaBSG
C          routine directly
C
C     Created: 02/22/07 by Azar Mustafayev
C     Modified: 6/12/07 by Azar Mustafayev - compatibility with ISAJET 7.75
C
      IMPLICIT NONE
c************************
c   ISAJET common blocks
c************************
#include "sugxin.inc"
#include "sugpas.inc"
#include "sugnu.inc"
#include "sugmg.inc"
#include "sssm.inc"
#include "sspar.inc"
c+CDE,SSNU
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RSIGT,RSIGB,RSIGL
c       GISA(i) - values of RGE parameters at MZ in DRbar:
C     GISA( 1) = g_1        GISA( 2) = g_2        GISA( 3) = g_3
C     GISA( 4) = y_tau      GISA( 5) = y_b        GISA( 6) = y_t
C     GISA( 7) = M_1        GISA( 8) = M_2        GISA( 9) = M_3
C     GISA(10) = A_tau      GISA(11) = A_b        GISA(12) = A_t
C     GISA(13) = M_hd^2     GISA(14) = M_hu^2     GISA(15) = M_er^2
C     GISA(16) = M_el^2     GISA(17) = M_dnr^2    GISA(18) = M_upr^2
C     GISA(19) = M_upl^2    GISA(20) = M_taur^2   GISA(21) = M_taul^2
C     GISA(22) = M_btr^2    GISA(23) = M_tpr^2    GISA(24) = M_tpl^2
C     GISA(25) = mu         GISA(26) = B          GISA(27) = Y_N
C     GISA(28) = M_nr       GISA(29) = A_n        GISA(30) = vdq
C     GISA(31) = vuq
c
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
c     RSIGT,RSIGB,RSIGL - radiative corrections to top, bottom and tau
c                         Yukawas at MSUSY
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RSIGT,RSIGB,RSIGL
      SAVE /BSG/
c************************
c   ISABSG common blocks
c************************
      COMMON /BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,BSGXW,S12,S23,S13,BSGAEM,
     &               BSGTHW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,BSGXW,S12,S23,S13,BSGAEM,BSGTHW
c  NB: The following name changes has to be made in order not to interfere
c      with ISAJET names, but only in this subroutine. In all other
c      subroutines original names are restored.
c       XW     -> BSGXW
c       ALFAEM -> BSGAEM
c       SN2THW -> BSGTHW
      SAVE /BSGSM/
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON /GLNN/ MCH0(4),MDSB(6)
      REAL*8 MCH0,MDSB
      SAVE /GLNN/
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),BSGRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER BSGRHN
      SAVE /BSGDEC/
      COMMON /BSGSUG/ TANB,BSGV,BSGVP,BSGQ,BSGMU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,BSGMZI(4),
     &                BSGZMIX(4,4),BSGMW1,BSGMW2,BSGGAML,BSGGAMR,
     &                BSGTHET,BSGMTQ,BSGMBQ,
     &                MSTLQ,MSTRQ,BSGMGUT,BSGNGUT,BSGFT,BSGNRIN(21),
     &                BSGGMIN(14),BSGGUN,BSGNUS(20),BSGAMIN(7),
     &                BSGENU,BSGFN(3),BSGMHU,BSGMHD,
     &                BSGNUHM,BSG3UN,BSG5ON
      REAL*8 TANB,BSGV,BSGVP,BSGQ,BSGMU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,BSGMZI,BSGZMIX,BSGMW1,BSGMW2,
     &       BSGGAML,BSGGAMR,BSGTHET,BSGMTQ,BSGMBQ,
     &       MSTLQ,MSTRQ,BSGMGUT,BSGNGUT,BSGFT,BSGNRIN,BSGGMIN,BSGGUN,
     &       BSGNUS,BSGAMIN,BSGENU,BSGFN,BSGMHU,BSGMHD
      INTEGER BSG3UN,BSGNUHM
      LOGICAL BSG5ON
c  NB: The following name changes has to be made in order not to interfere
c      with ISAJET names, but only in this subroutine. In all other
c      subroutines original names are restored.
c   V      -> BSGV               XRHNIN -> BSGNRIN
c   VP     -> BSGVP              XGMIN  -> BSGMIN
c   MU     -> BSGMU              IAL3UN -> BSG3UN
c   MSUSY  -> BSGQ               GGUTSS -> BSGGUN
c   MTQ    -> BSGMTQ             XNUSUG -> BSGNUS
c   MBQ    -> BSGMBQ             XAMIN  -> BSGAMIN
c   MGUT   -> BSGMGUT            AMZISS -> BSGMZI
c   FNGUT  -> BSGNGUT            ZMIXSS -> BSGZMIX
c   FTMT   -> BSGFT              AMWiSS -> BSGMWi
c   THETAT -> BSGTHET            EPSNU  -> BSGENU
c   GAMMAL -> BSGGAML            FTRHLD -> BSGFN
c   GAMMAR -> BSGGAMR            LND5ON -> BSG5ON
c   IRHN   -> BSGRHN             INUHM  -> BSGNUHM
c   MHUSMG -> BSGMHU             MHDSMG -> BSGMHD
      SAVE /BSGSUG/
      INTEGER I,J
      
      
      M1 = GISA(7)
      M2 = GISA(8)
      M3 = GISA(9)
      ABOT = GISA(11)
      ATOP = GISA(12)
      ATAU = GISA(10)
      
      MZ = DBLE(AMZ)
      MW = DBLE(AMW)
      MB = DBLE(AMBT)
      MC = 1.5d0
      MS = 0.2d0
      MT = DBLE(AMTP)
      MTAU = DBLE(AMTAU)
      
      S12 = 0.22d0
      S23 = 0.04d0
      S13 = 0.003d0
      
      BSGXW = DBLE(XW)
      BSGAEM = DBLE(ALFAEM)
      BSGTHW = DBLE(SN2THW)
      
      TANB = DBLE(XTANB)
      BSGV = DBLE(V)            ! v_u
      BSGVP = DBLE(VP)          ! v_d 
      BSGMU = DBLE(MU)
      
      BSGQ = DBLE(MSUSY)
      BSGMGUT = DBLE(MGUT)
      
      MSTP1 = DBLE(MSS(12))
      MSTP2 = DBLE(MSS(13))
      BSGTHET = DBLE(THETAT)
      
      MSCHL = DBLE(MSS(8))
      MSCHR = DBLE(MSS(9))
      MSUPL = DBLE(MSS(2))
      MSEL  = DBLE(MSS(17))
      MSW1  = DBLE(MSS(27))
      MGLU  = DBLE(MSS(1))
      MHPLUS= DBLE(MSS(32))
      MHA0  = DBLE(MSS(31))
      
      BSGMZI(1)=DBLE(AMZ1SS)
      BSGMZI(2)=DBLE(AMZ2SS)
      BSGMZI(3)=DBLE(AMZ3SS)
      BSGMZI(4)=DBLE(AMZ4SS)
      DO I=1,4
        DO J=1,4
          BSGZMIX(I,J)=DBLE(ZMIXSS(I,J))
        ENDDO
      ENDDO
      BSGGAML = DBLE(GAMMAL)
      BSGGAMR = DBLE(GAMMAR)
      BSGMW1 = DBLE(AMW1SS)
      BSGMW2 = DBLE(AMW2SS)
      
      BSGMTQ = DBLE(MTQ)
      BSGMBQ = DBLE(MBQ)
      MSTLQ = DBLE(AMTLSS)
      MSTRQ = DBLE(AMTRSS)
      
      BSGFT = DBLE(FTMT)
      BSGNGUT = DBLE(FNGUT)
      
c      DO I=1,21
c        BSGNRIN(I) = DBLE(XRHNIN(I))
c      ENDDO
c      BSGENU = EPSNU
c      DO I=1,3
c        BSGFN(I) = FTRHLD(I)
c      ENDDO
c      BSGRHN = IRHN
c      BSG5ON = LND5ON
      BSGENU = 1.d0
      IF (XNRIN(2).LT.1.E19) THEN 
        BSGRHN=1
      ELSE
        BSGRHN=0
      ENDIF
      BSG5ON = .FALSE.
      BSGFN(1)=0.d0
      BSGFN(2)=0.d0
      BSGFN(3)=XNRIN(2)
      
      
      DO I=1,14
        BSGGMIN(I) = DBLE(XGMIN(i))
      ENDDO
      DO I=1,20
        BSGNUS(I) = DBLE(XNUSUG(i))
      ENDDO
      DO I=1,7
        BSGAMIN(I) = DBLE(XAMIN(i))
      ENDDO
      
      BSG3UN = IAL3UN
      BSGGUN = DBLE(GGUTSS)
      
      BSGNUHM = INUHM
      BSGMHU = DBLE(MHUSMG)
      BSGMHD = DBLE(MHDSMG)

      DO I=1,3
        MSQDEC(I)=MSQISA(I)
	MSLDEC(I)=MSLISA(I)
	MSUDEC(I)=MSUISA(I)
	MSDDEC(I)=MSDISA(I)
	MSEDEC(I)=MSEISA(I)
	MRNDEC(I)=MRNISA(I)
      ENDDO
      
      RETURN
      END
C---------------------------------------------------------------------
      SUBROUTINE BSGGUT(YUGUT,YDGUT,YEGUT,YNGUT,G,IMODEL,M0,MHF,A0)            
C----------------------------------------------------------------------
C
C     Sets GUT-scale boundary conditions for G(157).
C
c    Created: 03/13/07 by Azar Mustafayev
c    Modified: 6/12/07 by Azar Mustafayev - compatible with ISAJET 7.75
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 G(157),YUGUT(3,3),YDGUT(3,3),YEGUT(3,3),YNGUT(3,3)
      INTEGER IMODEL
      REAL M0,MHF,A0
c
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),IRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER IRHN
      SAVE /BSGDEC/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
c     
      REAL*8 DDILOG
      REAL*8 A(3),XLAMGM,XMESGM,XN5GM,XLM,THRF,THRG,BLHAT,BBHAT,BTHAT,
     &       XC,DY,DG(31),PI
      INTEGER I,J,K

      PI=4.d0*DATAN(1.d0)

C...Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLAMGM=M0
        XMESGM=MHF
        XN5GM=A0
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
        XC=2*THRF*XLAMGM**2
        DY=SQRT(3./5.)*G(1)*XGMIN(11)
      ENDIF

c...gauge couplings
      G(1)=G(1)
      G(2)=G(2)
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) THEN
        G(3)=(G(1)+G(2))/2.d0
      ELSE
        G(3)=G(3)
      ENDIF
c...Yukawa couplings      
      CALL MAT2VEC(G,4,YUGUT,-1)
      CALL MAT2VEC(G,13,YDGUT,-1)
      CALL MAT2VEC(G,22,YEGUT,-1)
c...gaugino masses        
      IF (IMODEL.EQ.1) THEN
        DO J=1,3
          IF (XNUSUG(J).LT.1.d19) THEN
C       Set possible non-universal boundary conditions
            G(J+30)=XNUSUG(J)
          ELSE
	    G(J+30)= DBLE(MHF)
	  END IF
        ENDDO
      ELSEIF (IMODEL.EQ.2) THEN
	DO J=1,3
          G(J+30)=XGMIN(11+J)*XGMIN(8)*THRG*(G(J)/4./PI)**2*XLAMGM
        ENDDO
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
       G(31)= 33.d0*DBLE(MHF)*G(1)**2/5./16./PI**2
       G(32)= DBLE(MHF)*G(2)**2/16./PI**2
       G(33)=-3.d0*DBLE(MHF)*G(3)**2/16./PI**2
      ENDIF

C...Compute anomaly mediated functions
      IF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        BLHAT=G(30)*(-9*G(1)**2/5.-3*G(2)**2+3*G(21)**2+4*G(30)**2)
        BBHAT=G(21)*(-7*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             G(12)**2+6*G(21)**2+G(30)**2)
        BTHAT=G(12)*(-13*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             6*G(12)**2+G(21)**2)
      ENDIF

c...soft trilinar couplings      
      IF (IMODEL.EQ.1) THEN
        DO J=1,3
          IF (XNUSUG(J+3).LT.1.d19) THEN
C       Set possible non-universal boundary conditions
            A(J)=XNUSUG(J+3)
          ELSE
	    A(J)= DBLE(A0)
          ENDIF
        ENDDO
      ELSE IF (IMODEL.EQ.2) THEN
        DO J=1,3
          A(J)=0.
	ENDDO
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        A(3)=-BLHAT*MHF/G(30)/16./PI**2
        A(2)=-BBHAT*MHF/G(21)/16./PI**2
        A(1)=-BTHAT*MHF/G(12)/16./PI**2
      ENDIF
      DO I=1,9
        G(33+i) = G(3+i)*A(1)
        G(42+i) = G(12+i)*A(2)
        G(51+i) = G(21+i)*A(3)
      ENDDO
c...higgs mixing parameters      
      G(61)= G(61)
      G(62)= G(62)
c...higgs  mass^2
      IF (IMODEL.EQ.1) THEN
        DO J=1,2
          IF (XNUSUG(J+6).LT.1.E19) THEN
C       Set possible non-universal boundary conditions
C-FP    ANSI/gfortran fix
            G(J+62)=SIGN(1.D0,XNUSUG(J+6))*(XNUSUG(J+6))**2
          ELSE
            G(J+62)=DBLE(M0**2)
          ENDIF
        ENDDO
        IF (INUHM.EQ.1) THEN
          G(63) = MHUSMG
          G(64) = MHDSMG
        ENDIF
      ELSEIF (IMODEL.EQ.2) THEN
        G(64)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     ,             XGMIN(12)*(G(1)/4./PI)**4)+XGMIN(9)-DY
        G(63)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     ,             XGMIN(12)*(G(1)/4./PI)**4)+XGMIN(10)+DY
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(64)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(21)*BBHAT+G(30)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2+XAMIN(6)*M0**2
        G(63)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(12)*BTHAT)*
     ,        MHF**2/(16*PI**2)**2+XAMIN(7)*M0**2
      ENDIF
c
c   Set the rest of SSB mass^2      
c
c...Reset mass^2 elements
      DO I=65,147
  	IF(I.LE.109.OR.I.GE.121) G(I)=0.d0
      ENDDO

c...diagonal elements of m^2_Q
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(13).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(65)=SIGN(1.D0,XNUSUG(13))*(XNUSUG(13))**2
        ELSE
          G(65)= DBLE(M0**2)
        ENDIF
        G(69)=G(65)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(18).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(73)=SIGN(1.D0,XNUSUG(18))*(XNUSUG(18))**2
        ELSE
          G(73)= DBLE(M0**2)
        ENDIF
      ELSEIF (IMODEL.EQ.2) THEN
        G(65)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.75*XGMIN(13)*
     ,        (G(2)/4./PI)**4+.6*XGMIN(12)*(G(1)/4./PI)**4/36.)+DY/3.
        G(69)=G(65)
        G(73)=G(65)
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(65)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4)*
     ,         MHF**2/(16*PI**2)**2+XAMIN(1)*M0**2
        G(69)=G(65)
        G(73)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4+G(21)*BBHAT+
     ,         G(12)*BTHAT)*MHF**2/(16*PI**2)**2+XAMIN(1)*M0**2
      ENDIF  
c...diagonal elements of m^2_L
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(10).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(74)=SIGN(1.D0,XNUSUG(10))*(XNUSUG(10))**2
        ELSE
          G(74)= DBLE(M0**2)
        END IF
        G(78)=G(74)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(15).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(82)=SIGN(1.D0,XNUSUG(15))*(XNUSUG(15))**2
        ELSE
          G(82)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
       G(74)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     ,        XGMIN(12)*(G(1)/4./PI)**4)-DY
        G(78)=G(74) 
        G(82)=G(74) 
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(74)=(-99*G(1)**4/50.-3*G(2)**4/2.)*MHF**2/(16*PI**2)**2
     ,+XAMIN(4)*M0**2
        G(78)=G(74) 
        G(82)=(-99*G(1)**4/50.-3*G(2)**4/2.+G(30)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2+XAMIN(4)*M0**2
      ENDIF  
c...diagonal elements of m^2_u
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(12).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(83)=SIGN(1.D0,XNUSUG(12))*(XNUSUG(12))**2
        ELSE
          G(83)= DBLE(M0**2)
        END IF
        G(87)=G(83)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(17).LT.1.E19) THEN
C-FP      ANSI/gfortran fix
          G(91)=SIGN(1.D0,XNUSUG(17))*(XNUSUG(17))**2
        ELSE
          G(91)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
        G(83)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.6*4*XGMIN(12)*
     ,        (G(1)/4./PI)**4/9.)-4*DY/3.
        G(87)=G(83)
        G(91)=G(83)
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(83)=(-88*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2+
     ,XAMIN(3)*M0**2
        G(87)=G(83)
        G(91)=(-88*G(1)**4/25.+8*G(3)**4+2*G(12)*BTHAT)*
     , MHF**2/(16*PI**2)**2+XAMIN(3)*M0**2
      ENDIF  
c...diagonal elements of m^2_d
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(11).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(92)=SIGN(1.D0,XNUSUG(11))*(XNUSUG(11))**2
        ELSE
          G(92)= DBLE(M0**2)
        END IF
        G(96)=G(92)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(16).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(100)=SIGN(1.D0,XNUSUG(16))*(XNUSUG(16))**2
        ELSE
          G(100)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
        G(92)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.6*XGMIN(12)*
     , (G(1)/4./PI)**4/9.)+2*DY/3.
        G(96)=G(92)  
        G(100)=G(92) 
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(92)=(-22*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2+
     ,XAMIN(2)*M0**2
        G(96)=G(92)  
        G(100)=(-22*G(1)**4/25.+8*G(3)**4+2*G(21)*BBHAT)*
     , MHF**2/(16*PI**2)**2+XAMIN(2)*M0**2
      ENDIF  
c...diagonal elements of m^2_e
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(9).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(101)=SIGN(1.D0,XNUSUG(9))*(XNUSUG(9))**2
        ELSE
          G(101)= DBLE(M0**2)
        END IF
        G(105)=G(101)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(14).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(109)=SIGN(1.D0,XNUSUG(16))*(XNUSUG(14))**2
        ELSE
          G(109)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
        G(101)=XC*(.6*XGMIN(12)*(G(1)/4./PI)**4)+2*DY
        G(105)=G(101) 
        G(109)=G(101) 
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(101)=(-198*G(1)**4/25.)*MHF**2/(16*PI**2)**2+XAMIN(5)*M0**2
        G(105)=G(101) 
        G(109)=(-198*G(1)**4/25.+2*G(30)*BLHAT)*MHF**2/(16*PI**2)**2+
     ,XAMIN(5)*M0**2
      ENDIF  

c...additional contributions for mixed moduli-AMSB model
      IF (IMODEL.EQ.9) THEN
        DO I=1,3
          DG(I)=G(I)
        ENDDO
	DG(4)=G(30)
	DG(5)=G(21)
	DG(6)=G(12)
	DG(7)=G(31)
	DG(8)=G(32)
	DG(9)=G(33)
	DG(9)=G(33)
	DG(10)=G(60)/G(30)
	DG(11)=G(51)/G(21)
	DG(12)=G(42)/G(12)
	DG(13)=G(64)
	DG(14)=G(63)
	DG(16)=G(101)
	DG(17)=G(74)
	DG(18)=G(92)
	DG(19)=G(83)
	DG(20)=G(109)
	DG(21)=G(82)
	DG(22)=G(100)
	DG(23)=G(91)
	DG(24)=G(73)
	
        CALL MMAMSB(M0,MHF,DG)
        
	G(31)=DG(7)
	G(32)=DG(8)
	G(33)=DG(9)
	G(33)=DG(9)
	G(60)=DG(10)*G(30)
	G(51)=DG(11)*G(21)
	G(42)=DG(12)*G(12)
	G(64)=DG(13)
	G(63)=DG(14)
	G(101)=DG(16)
	G(74)=DG(17)
	G(92)=DG(18)
	G(83)=DG(19)
	G(109)=DG(20)
	G(82)=DG(21)
	G(100)=DG(22)
	G(91)=DG(23)
	G(73)=DG(24)
      ENDIF

c...higgs VEVs
      G(110) = G(110)
      G(111) = G(111)
C
C  neutrino sector
C
      IF (IRHN.GT.0) THEN
c...Yukawa matrix
        IF (IRHN.EQ.2) THEN     ! impose YN-YU unification with degree EPSNU
          DO I=1,3
	    DO J=1,3
	      YNGUT(I,J)=YUGUT(I,J)*EPSNU
	    ENDDO
	  ENDDO  
        ELSEIF (IRHN.EQ.3) THEN
          K=0
          DO I=1,3
            DO J=1,3
	      YNGUT(J,I)=XRHNIN(10+K)
	      K=K+1
            ENDDO
          ENDDO
        ENDIF
	CALL MAT2VEC(G,112,YNGUT,-1)
c...RHN majorana mass matrix
        DO I=1,9
          G(120+I) = XRHNIN(I) 
        ENDDO
c...soft trilinear coupling
        DO I=0,8
          G(130+I) = G(112+I)*XRHNIN(19)
        ENDDO
c...diagonal elements of SSB mass^2
C-FP    ANSI/gfortran fix
        G(139) = SIGN(1.D0,XRHNIN(11))*(XRHNIN(20))**2
        G(143) = G(139)             ! degernerate 1st and 2nd generations
        G(147) = SIGN(1.D0,XRHNIN(12))*(XRHNIN(21))**2
c...dim-5 related terms
        DO I=148,157
	  G(I)=0.d0
	ENDDO
      ENDIF
      
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE YUKDIAG(G,INIT,YDIAG,VL,VR)
C----------------------------------------------------------------------
C
C  Assembles Yukawa matrix in weak eigenbasis from RGE running vector G(i)
C  starting from i=INIT and diagonalizes it.
C
C   NOTE: In our notation fermion masses are proportional to transposed
C         Yukawas in gauge eigenbasis, i.e.  m ~ VR Y^T VL^dagger,
C         while SURG111 follows m ~ Y convention.
C         One has to transpose Yukasas when passing
C         into or from those subroutines.
C
C   NOTE: VL and VR contain eigenvectors as rows, while Z - as columns. 
c
c  Ref: Misak, Pokorski & Rosiek  hep-ph/9703442
c
c  Created: 01/19/06 by Azar Mustafayev.
C 
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C   Input:
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
c   Output:
c    VL,VR - Left and right rotation matrices: YDIAG = VL^* YUK VR^T
C    YDIAG - diagonal Yukawa matrix
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv            
      IMPLICIT NONE
      REAL*8 G(157),YDIAG(3,3),VL(3,3),VR(3,3)
      INTEGER INIT
#include "sslun.inc"
      REAL*8 YUK(3,3),YUKH(3,3),A(3,3),W(3),Z(3,3),WORK(3),I3(3,3)
c  YUK(3,3)  - Yukawa matrix in weak eigenbasis
c  YUKH(3,3) - hermitian conjugated YUK
c  A(3,3)  - axiliary matrix to be diagonalized
c  W(3)    - vector of eigenvalues for AR
c  Z(3,3)  - matrix of eigenvalues of AR as columns
      INTEGER I,J,K,IERR
      DATA I3/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix
      
c...assemble Yukawa matrix in weak eigenbasis
      CALL VEC2MAT(G,INIT,YUK,-1)
      YUKH(1,1)=YUK(1,1)
      YUKH(2,1)=YUK(1,2)
      YUKH(3,1)=YUK(1,3)
      YUKH(1,2)=YUK(2,1)
      YUKH(2,2)=YUK(2,2)
      YUKH(3,2)=YUK(2,3)
      YUKH(1,3)=YUK(3,1)
      YUKH(2,3)=YUK(3,2)
      YUKH(3,3)=YUK(3,3)

c...Compute left rotation matrix VL by diagonalizing  YUK*YUKH     
      CALL MPROD2(A,YUK,YUKH)
      CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
      IF (IERR.NE.0) THEN
        WRITE(LOUT,*) 'EISRS1 ERROR IN YUKDIAG, IERR=',IERR
        STOP99
      END IF
      DO I=1,3
        DO J=1,3
	  VL(I,J)=Z(J,I)
        ENDDO
      ENDDO
c...check if VL has all diagonal elements positive
      DO I=1,3
       IF (VL(I,I).LT.0.d0) THEN      
	 DO J=1,3
	  VL(I,J)=-VL(I,J)   ! sign change in the row
	 ENDDO
       ENDIF
      ENDDO
c...Compute right rotation matrix VR by diagonalizing  YUKH*YUK     
      CALL MPROD2(A,YUKH,YUK)
      CALL EIGSYS(3,3,A,W,Z,IERR,WORK)
      IF (IERR.NE.0) THEN
        WRITE(LOUT,*) 'EISRS1 ERROR IN YUKDIAG, IERR=',IERR
        STOP99
      END IF
      DO I=1,3
        DO J=1,3
	  VR(I,J)=Z(J,I)
        ENDDO
      ENDDO
c...check if VR has all diagonal elements positive
      DO I=1,3
       IF (VR(I,I).LT.0.d0) THEN      
	 DO J=1,3
	  VR(I,J)=-VR(I,J)   ! sign change in the row
	 ENDDO
       ENDIF
      ENDDO
c...Build diagonal Yukawa matrix
      DO I=1,3
        DO J=1,3
	  YDIAG(I,J)=SQRT(W(I))*I3(I,J)
        ENDDO
      ENDDO

      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE AYUKDIAG(G,INIT,YDIAG,VL,VR)
C----------------------------------------------------------------------
C
C  This is inversed of YUKDIAG.
c  It rotates diagonal Yukawa matrix to weak eigenbasis and 
c  "diassembles" it filling corresponding elements of G(i).
c
C   NOTE: In our notation fermion masses are proportional to transposed
C         Yukawas in gauge eigenbasis, i.e.  m ~ VR Y^T VL^dagger,
C         while SURG111 follows m ~ Y convention.
C         One has to transpose Yukasas when passing
C         into or from those subroutines.
C
c  Ref: Misak, Pokorski & Rosiek  hep-ph/9703442
c
c  Created: 02/13/07 by Azar Mustafayev.
C 
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C   Input:
c    VL,VR - Left and right rotation matrices: YDIAG=VL^* YUK VR^T
C    YDIAG - diagonal Yukawa matrix
c   Output:
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv            
      IMPLICIT NONE
      REAL*8 G(157),YDIAG(3,3),VL(3,3),VR(3,3)
      INTEGER INIT
c      
      REAL*8 YUK(3,3),VLT(3,3)

      VLT(1,1)=VL(1,1)
      VLT(2,1)=VL(1,2)
      VLT(3,1)=VL(1,3)
      VLT(1,2)=VL(2,1)
      VLT(2,2)=VL(2,2)
      VLT(3,2)=VL(2,3)
      VLT(1,3)=VL(3,1)
      VLT(2,3)=VL(3,2)
      VLT(3,3)=VL(3,3)
      
      CALL MPROD3(YUK,VLT,YDIAG,VR)
 
      CALL MAT2VEC(G,INIT,YUK,-1)

      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE MAT2VEC(G,INIT,MAT,IORD)
C----------------------------------------------------------------------
C
C   This is reciprocal to VEC2MAT.
C  "Disassembles" 3x3 matrix into vector G(i) starting from element i=INIT.
C   Matrix can be broken in row- or column-vise order.
C
C    Created: 5/03/07 by Azar Mustafayev.
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      IMPLICIT NONE
      REAL*8 G(157),MAT(3,3)
      INTEGER INIT,IORD
#include "sslun.inc"

      IF (IORD.EQ.1) THEN        ! row-vise
	G(INIT)  =MAT(1,1)
	G(INIT+1)=MAT(1,2)
	G(INIT+2)=MAT(1,3)
	G(INIT+3)=MAT(2,1)
	G(INIT+4)=MAT(2,2)
	G(INIT+5)=MAT(2,3)
	G(INIT+6)=MAT(3,1)
	G(INIT+7)=MAT(3,2)
	G(INIT+8)=MAT(3,3)
      ELSEIF (IORD.EQ.-1) THEN    ! column-vise
	G(INIT)  =MAT(1,1)
	G(INIT+1)=MAT(2,1)
	G(INIT+2)=MAT(3,1)
	G(INIT+3)=MAT(1,2)
	G(INIT+4)=MAT(2,2)
	G(INIT+5)=MAT(3,2)
	G(INIT+6)=MAT(1,3)
	G(INIT+7)=MAT(2,3)
	G(INIT+8)=MAT(3,3)
      ELSE
        WRITE(LOUT,*) 'MAT2VEC: unknown assembly order'
        STOP99
      ENDIF  

      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE VEC2MAT(G,INIT,MAT,IORD)
C----------------------------------------------------------------------
C
C    Assembles 3x3 matrix from vector G(i) starting from element i=INIT.
C    Matrix can be filled in row- or column-vise order.
C
C    Created: 5/03/07 by Azar Mustafayev.
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      IMPLICIT NONE
      REAL*8 G(157),MAT(3,3)
      INTEGER INIT,IORD
#include "sslun.inc"

      IF (IORD.EQ.1) THEN        ! row-vise
	MAT(1,1)=G(INIT)
	MAT(1,2)=G(INIT+1)
	MAT(1,3)=G(INIT+2)
	MAT(2,1)=G(INIT+3)
	MAT(2,2)=G(INIT+4)
	MAT(2,3)=G(INIT+5)
	MAT(3,1)=G(INIT+6)
	MAT(3,2)=G(INIT+7)
	MAT(3,3)=G(INIT+8)
      ELSEIF (IORD.EQ.-1) THEN    ! column-vise
	MAT(1,1)=G(INIT)
	MAT(2,1)=G(INIT+1)
	MAT(3,1)=G(INIT+2)
	MAT(1,2)=G(INIT+3)
	MAT(2,2)=G(INIT+4)
	MAT(3,2)=G(INIT+5)
	MAT(1,3)=G(INIT+6)
	MAT(2,3)=G(INIT+7)
	MAT(3,3)=G(INIT+8)
      ELSE
        WRITE(LOUT,*) 'VEC2MAT: unknown assembly order'
        STOP99
      ENDIF  

      RETURN
      END
C----------------------------------------------------------------------
      REAL*8 FUNCTION TRACE(X,I)
C----------------------------------------------------------------------
C
C     Takes the trace of 250x3x3 matrix X with respect to last two indices.
C
      IMPLICIT NONE
      REAL*8 X(250,3,3)
      INTEGER I
      
      TRACE=X(I,1,1)+X(I,2,2)+X(I,3,3)
      
      RETURN
      END      
C----------------------------------------------------------------------
      REAL*8 FUNCTION TR3X3(Y)
C----------------------------------------------------------------------
C
C     Takes the trace of 3x3 matrix Y.
C
      IMPLICIT NONE
      REAL*8 Y(3,3)
      
      TR3X3=Y(1,1)+Y(2,2)+Y(3,3)
      
      RETURN
      END      
C-----------------------------------------------------------------------
      SUBROUTINE MPROD2(Y,A,B)
C-----------------------------------------------------------------------
C
C     Multiplies two 3x3 matrices A and B and stores result in 
C     3x3 matrix Y.
C
C     Note: Uses unfolded DO-loop that greatly increases speed.
C           6*9=54 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 Y(3,3),A(3,3),B(3,3)
      
      Y(1,1)= A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      Y(1,2)= A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      Y(1,3)= A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      
      Y(2,1)= A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      Y(2,2)= A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      Y(2,3)= A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      
      Y(3,1)= A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      Y(3,2)= A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      Y(3,3)= A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      
      RETURN
      END      
C-----------------------------------------------------------------------
      SUBROUTINE MPROD3(Y,A,B,C)
C-----------------------------------------------------------------------
C
C     Multiplies three 3x3 matrices A, B and C and stores result in 
C     3x3 matrix Y.
C
C     Note: Uses unfolded DO-loop and two-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 Y(3,3),A(3,3),B(3,3),C(3,3),Z(3,3)

c...compute axiliary product matrix Z      
      Z(1,1)= A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      Z(1,2)= A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      Z(1,3)= A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      
      Z(2,1)= A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      Z(2,2)= A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      Z(2,3)= A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      
      Z(3,1)= A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      Z(3,2)= A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      Z(3,3)= A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)

c...multiply axiliary matrix Z by the third matrix C      
      Y(1,1)= Z(1,1)*C(1,1)+Z(1,2)*C(2,1)+Z(1,3)*C(3,1)
      Y(1,2)= Z(1,1)*C(1,2)+Z(1,2)*C(2,2)+Z(1,3)*C(3,2)
      Y(1,3)= Z(1,1)*C(1,3)+Z(1,2)*C(2,3)+Z(1,3)*C(3,3)
      
      Y(2,1)= Z(2,1)*C(1,1)+Z(2,2)*C(2,1)+Z(2,3)*C(3,1)
      Y(2,2)= Z(2,1)*C(1,2)+Z(2,2)*C(2,2)+Z(2,3)*C(3,2)
      Y(2,3)= Z(2,1)*C(1,3)+Z(2,2)*C(2,3)+Z(2,3)*C(3,3)
      
      Y(3,1)= Z(3,1)*C(1,1)+Z(3,2)*C(2,1)+Z(3,3)*C(3,1)
      Y(3,2)= Z(3,1)*C(1,2)+Z(3,2)*C(2,2)+Z(3,3)*C(3,2)
      Y(3,3)= Z(3,1)*C(1,3)+Z(3,2)*C(2,3)+Z(3,3)*C(3,3)
      
      RETURN
      END      
C-----------------------------------------------------------------------
      SUBROUTINE MPROD2X(X,I,A,B)
C-----------------------------------------------------------------------
C
C     Multiplies two 3x3 matrices A and B and stores result in 
C     i-th submatrix of 250x3x3 matrix X.
C
C     Note: Uses unfolded DO-loop that greatly increases speed.
C           6*9=54 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 X(250,3,3),A(3,3),B(3,3)
      INTEGER I
      
      X(I,1,1)= A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      X(I,1,2)= A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      X(I,1,3)= A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      
      X(I,2,1)= A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      X(I,2,2)= A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      X(I,2,3)= A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      
      X(I,3,1)= A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      X(I,3,2)= A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      X(I,3,3)= A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      
      RETURN
      END      
C-----------------------------------------------------------------------
      SUBROUTINE MPROD3X(X,I,A,B,C)
C-----------------------------------------------------------------------
C
C     Multiplies three 3x3 matrices A, B and C and stores result in 
C     i-th submatrix of 250x3x3 matrix X.
C
C     Note: Uses unfolded DO-loop and two-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 X(250,3,3),A(3,3),B(3,3),C(3,3),Y(3,3)
      INTEGER I

c...compute axiliary product matrix Y      
      Y(1,1)= A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      Y(1,2)= A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      Y(1,3)= A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      
      Y(2,1)= A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      Y(2,2)= A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      Y(2,3)= A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      
      Y(3,1)= A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      Y(3,2)= A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      Y(3,3)= A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)

c...multiply axiliary matrix Y by the third matrix C      
      X(I,1,1)= Y(1,1)*C(1,1)+Y(1,2)*C(2,1)+Y(1,3)*C(3,1)
      X(I,1,2)= Y(1,1)*C(1,2)+Y(1,2)*C(2,2)+Y(1,3)*C(3,2)
      X(I,1,3)= Y(1,1)*C(1,3)+Y(1,2)*C(2,3)+Y(1,3)*C(3,3)
      
      X(I,2,1)= Y(2,1)*C(1,1)+Y(2,2)*C(2,1)+Y(2,3)*C(3,1)
      X(I,2,2)= Y(2,1)*C(1,2)+Y(2,2)*C(2,2)+Y(2,3)*C(3,2)
      X(I,2,3)= Y(2,1)*C(1,3)+Y(2,2)*C(2,3)+Y(2,3)*C(3,3)
      
      X(I,3,1)= Y(3,1)*C(1,1)+Y(3,2)*C(2,1)+Y(3,3)*C(3,1)
      X(I,3,2)= Y(3,1)*C(1,2)+Y(3,2)*C(2,2)+Y(3,3)*C(3,2)
      X(I,3,3)= Y(3,1)*C(1,3)+Y(3,2)*C(2,3)+Y(3,3)*C(3,3)
      
      RETURN
      END      
C-----------------------------------------------------------------------
      SUBROUTINE MPROD4X(X,I,A,B,C,D)
C-----------------------------------------------------------------------
C
C     Multiplies four 3x3 matrices A, B, C and D and stores result in 
C     i-th submatrix of 250x3x3 matrix X.
C
C     Note: Uses unfolded DO-loop and two-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 X(250,3,3),A(3,3),B(3,3),C(3,3),D(3,3),Y(3,3)
      INTEGER I
      
c...compute axiliary product matrix Y      
      CALL MPROD3(Y,A,B,C)

c...multiply axiliary matrix Y by the forth matrix D      
      CALL MPROD2X(X,I,Y,D)
      
      RETURN
      END      
C-----------------------------------------------------------------------
      SUBROUTINE MPROD5X(X,I,A,B,C,D,F)
C-----------------------------------------------------------------------
C
C     Multiplies five 3x3 matrices A, B, C, D and F and stores result in 
C     i-th submatrix of 250x3x3 matrix X.
C
C     Note: Uses unfolded DO-loop and three-stage multiplication
c           that greatly increases speed.
C           54*2=108 floating point operations
C
C     Created: 12/13/07 by Azar Mustafayev
C      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 X(250,3,3),A(3,3),B(3,3),C(3,3),D(3,3),F(3,3),Y(3,3),Z(3,3)
      INTEGER I

c...compute axiliary product matrix Y      
      CALL MPROD3(Y,A,B,C)
      
c...multiply axiliary matrix Y by D      
      CALL MPROD2(Z,Y,D)

c...multiply axiliary matrix Z by remaining matrix F      
      CALL MPROD2X(X,I,Z,F)
      
      RETURN
      END      
C---------------------------------------------------------------------
      SUBROUTINE FUNS(X,FN1,FN2,FN3,FN4)            
C----------------------------------------------------------------------
C
C    Computes convenient functions given by Eq(75) of Anlauf.
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 FN1,FN2,FN3,FN4
C
      IF(ABS(X-1.d0).GT.1.d-2) THEN
       FN1=(X**2-5.*X-2.d0)/12.d0/(X-1.d0)**3
     $     +X*LOG(X)/2.d0/(X-1.d0)**4
       FN2=(2.*X**2+5.*X-1.d0)/12.d0/(X-1.d0)**3
     $    -X**2*LOG(X)/2.d0/(X-1.d0)**4
       FN3=(X-3.d0)/2.d0/(X-1.d0)**2
     $    +LOG(X)/(X-1.d0)**3
       FN4=(X+1.d0)/2.d0/(X-1.d0)**2
     $    -X*LOG(X)/(X-1.d0)**3           
      ELSE
       FN1=1./24.d0-(X-1.d0)/40.d0+(X-1.d0)**2/60.d0
       FN2=1./24.d0-(X-1.d0)/60.d0+(X-1.d0)**2/120.d0
       FN3=1./3.d0-(X-1.d0)/4.d0+(X-1.d0)**2/5.d0
       FN4=1./6.d0-(X-1.d0)/12.d0+(X-1.d0)**2/20.d0
      ENDIF
C
      RETURN 
      END
C----------------------------------------------------------------------
      REAL*8 FUNCTION BI(X)            
C----------------------------------------------------------------------
C     Inami-Lim B function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X
C
      BI=1./4.d0*(X/(1.d0-X)+X*LOG(X)/(X-1.d0)**2)
C 
      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION CI(X)            
C----------------------------------------------------------------------
C     Inami-Lim C function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X

      CI=1./8.d0*X*((X-6.d0)/(X-1.d0)+(3.d0*X+2.d0)/(X-1.d0)**2*LOG(X))
C 
      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION DI(X)            
C----------------------------------------------------------------------
C     Inami-Lim D function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X

      DI=-4./9.d0*LOG(X)+(-19.d0*X**3+25.d0*X**2)/36./(X-1.d0)**3
     $   +X**2*(5.d0*X**2-2.d0*X-6.d0)/18.d0/(X-1.d0)**4*LOG(X)
C 
      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION EI(X)            
C----------------------------------------------------------------------
C     Inami-Lim E function given in Table 1 of M.Ciuchini et al.
C
      IMPLICIT NONE
      REAL*8 X
C
      EI=-2./3.d0*LOG(X)
     &   +X**2*(15.d0-16.d0*X+4.*X**2)/6.d0/(X-1.d0)**4*LOG(X)
     $   +X*(18.d0-11.d0*X-X**2)/12.d0/(1.-X)**3
C 
      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION XKARE(X)            
C----------------------------------------------------------------------
C     Real part of function \kappa defined in Eq(B6) of Greub.
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 PI,GE
C
      PI=4.d0*ATAN(1.d0)
C
      IF(X.LT.4.d0) THEN
       GE=-2.d0*ATAN(SQRT(X/(4.d0-X)))**2
      ELSE
       GE=-PI**2/2.d0+2.d0*(LOG(SQRT(X)/2.d0+SQRT(X-4.d0)/2.d0))**2
      ENDIF
      XKARE=4.d0*(2.d0*GE+X)/X
C 
      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION XKAIM(X)            
C----------------------------------------------------------------------
C     Imaginary part of function \kappa defined in Eq(B6) of Greub.
C
      IMPLICIT NONE
      REAL*8 X
      REAL*8 PI,GE
C
      PI=4.d0*ATAN(1.d0)
C
      IF(X.LT.4.d0) THEN
       GE=0.d0
      ELSE
       GE=-2.d0*PI*LOG(SQRT(X)/2.d0+SQRT(X-4.d0)/2.d0)
      ENDIF
C 
      XKAIM=8.d0*GE/X
C
      RETURN 
      END
C-----------------------------------------------------------------
      SUBROUTINE BREMS(AEMB,ASMB,MS,MC,MB,Q,CI,GBREM)            
C-----------------------------------------------------------------
C
C    Computes bremsstrahlung corrections according to formulae
C    given in Appendix B of Greub.
C
      IMPLICIT NONE
      REAL*8 AEMB,ASMB,MB,MC,MS,Q,CI(8)
      REAL*8 GBREM
      EXTERNAL FNTG
      COMMON /TAUS/ FAC22,FAC27,FAC28,FAC78
      REAL*8 FAC22,FAC27,FAC28,FAC78
      SAVE /TAUS/ 
      REAL*8 PI,FACAL,FACBL,RINT,BREMSF,BREMS88
      REAL*8 DTRINT
      
      PI=4.d0*ATAN(1.d0)          
C     
      FACAL=AEMB*ASMB/768.d0/PI**5/4.d0
      FACBL=AEMB*ASMB/96.d0/PI**5
      FAC22= 2.*4./9.d0*CI(2)**2
      FAC27=-32.*2./3.d0*CI(2)*CI(7)
      FAC28=-32.*2./3.d0*(-1./3.d0)*CI(2)*CI(8)
      FAC78=-128.d0*(-1./3.d0)*CI(7)*CI(8)
C
      RINT=DTRINT(FNTG,0,64,1.d-3,0.d0,1.d0,1.d0,0.d0,1.d0,1.d0)
      BREMSF=FACAL*RINT
      BREMS88=FACBL*(-1./3.d0*CI(8))**2
     $        *(16./3.d0-4.*PI**2/3.d0+4.d0*LOG(MB/Q))
C
      GBREM=BREMSF+BREMS88
      
      RETURN
      END 
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION FNTG(X,Y)            
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,Y
      COMMON /TAUS/ FAC22,FAC27,FAC28,FAC78
      REAL*8 FAC22,FAC27,FAC28,FAC78
      SAVE /TAUS/ 
      COMMON /BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      REAL*8 XKARE,XKAIM
      
      REAL*8 T,U,S,SC,RHO,O22,O27,O28,O78
      
      RHO=(MS/MB)**2
      T= X
      U= Y
      S=T+U-1.d0+RHO
      SC=S*(MB/MC)**2
C
      O22=(XKARE(SC)**2+XKAIM(SC)**2)*(1.d0-S)
      O27=1./2.d0*S*XKARE(SC)
      O28=1./2.d0*S*XKARE(SC)
      O78=1./2.d0*S*(1.d0+1./2.d0*T*U)/T/U
      
      FNTG=FAC22*O22+FAC27*O27+FAC28*O28+FAC78*O78

      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION GES(X)            
C----------------------------------------------------------------------
C    Phase space function given by Eq(5.10) of Greub.
C
      IMPLICIT NONE
      REAL*8 X
      
      GES=1.d0-8.d0*X**2+8.d0*X**6-X**8-24.d0*X**4*LOG(X)
      
      RETURN 
      END
C---------------------------------------------------------------------- 
      REAL*8 FUNCTION FES(X)            
C----------------------------------------------------------------------
C     Pole to MSbar conversion factor given by Eq(5.7) of Greub.
C
      IMPLICIT NONE
      REAL*8 X,PI
C
      PI=4.d0*DATAN(1.d0)       
      FES=(PI**2-31./4.d0)*(1.d0-X)**2+3./2.d0
C 
      RETURN 
      END
#endif
