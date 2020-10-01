#include "PILOT.inc"
      SUBROUTINE SUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2
     $,NSTEP,IMODEL,BADMU)
C
C          Make one complete iteration of the renormalization group
C          equations from MZ to MGUT and back, setting the boundary
C          conditions on each end.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sssm.inc"
#include "sspar.inc"
#include "sugpas.inc"
#include "sugnu.inc"
#include "sugxin.inc"
#include "sugmg.inc"
#include "ssinf.inc"
C
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
C     Common blocks for A. Box RGE code
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
C
      COMMON/WKYUK/LAMTMT,LAMBMZ,LAMTAMZ
      DOUBLE PRECISION LAMTMT,LAMBMZ,LAMTAMZ
      SAVE/WKYUK/
C
C
      EXTERNAL SURG26
      DOUBLE PRECISION DDILOG,XLM
      COMPLEX*16 SSB0,SSB1
      REAL*8 G(31),W2(93),T,DT,DY,DPI,DAS
      REAL*8 BTHAT,BBHAT,BLHAT
      REAL M0,MHF,A0,TANB,SGNMU,MT,G0(31)
      INTEGER IG(31),NSTEP,IMODEL
      REAL PI,TZ,A1I,A2I,A3I,GGUT,AGUTI,SIG1,SIG2,SIG3,
     $MH1S,MH2S,MUS,MZ,TGUT,AGUT,Q,ASMT,MTMT,
     $QNEW,XLAMGM,XMESGM,XN5GM,XC,G3GUT,THRF,THRG,MBMZ,
     $M2,AM2,MSN,MG,MT1,MT2,MB1,MB2,MW1,MW2,AMU,
     $RSIGT,RSIGL,RSIGB,DAEM,ALEMDR,ALEM,MTAMZ,
     $TANBQ,SIN2BQ,SINB,COSB,XWMSB
      REAL SSRSGT,SSRSGB,SSRSGL,SUALFS,MZQ,SIGA
      REAL CF,CA,ZETA2,ZETA3,ST2LP,COTB
      INTEGER I,II
      LOGICAL BADMU
C
      DATA ZETA2/1.644934/,ZETA3/1.202057/
      DATA MZ/91.187/
C
C          Re-initialize weak scale parameters
C
      XLAMGM=M0
      XMESGM=MHF
      XN5GM=A0
      PI=4.*ATAN(1.)
      DPI=4.D0*DATAN(1.D0)
      CF=4./3.
      CA=3.
C     Here we input alpha_s^MSbar(MZ)
      ASMZ=0.1172
      MTAMZ=1.7463
C     Value of mb(MZ)^DRbar taken from PRD66, 074007 (2002).
      MBMZ=2.83
      SINB=SIN(ATAN(TANB))
      COSB=COS(ATAN(TANB))
      COTB=1./TANB
C
C     Calculate fermion masses including loop corrections
C
      M2=G0(8)
      AM2=ABS(M2)
      MSN=MSS(16)
      MG=ABS(MSS(1))
      MT1=MSS(12)
      MT2=MSS(13)
      MB1=MSS(10)
      MB2=MSS(11)
      MW1=ABS(MSS(27))
      MW2=ABS(MSS(28))
      AMU=ABS(MU)
      MUS=MU**2
      XLAM=DLOG(DBLE(HIGFRZ**2))
C      FTMT=MTMT/V
C      FBMZ=MBMZ/VP
C      FTAMZ=MTAMZ/VP
C          Be careful in using our convention vs Pierce et al.
C          cos(th)>1/sqrt(2):eigenstates same; cs-> -cs
C          cos(th)<1/sqrt(2);flip mass eigenstates; c <-> s interchange
C          Formula remains invariant under these switches
C          Use negative gaugino masses for consistency
C          Now input self-energies for consistency: RSIG*
C
C     here, add in 2 loop QCD correction to mt(DRbar) from Bednyakov
C     et al. Eq. 61.
      ST2LP=CF*(ASMTP/4./PI)**2*(-43.-12*ZETA2+CF*(-59./8.+30*ZETA2-
     ,48*LOG(2.)*ZETA2+12*ZETA3)+
     ,CA*(1093./24.-8*ZETA2+24*LOG(2.)*ZETA2-6*ZETA3))
      MTMT=MT/(1.+5*ASMTP/3./PI+ST2LP)
      FTMT=MTMT/SINB/VEV
      FBMZ=MBMZ/COSB/VEV
      FTAMZ=MTAMZ/COSB/VEV
      RSIGT=SSRSGT(MT**2)
      RSIGB=SSRSGB(MBQ**2)
      RSIGL=SSRSGL(MLQ**2)
C
C Weak Scale Yukawas used by A. Box RGE code
C
      LAMBMZ=DBLE(MBMZ/VEV)
      LAMTAMZ=DBLE(MTAMZ/VEV)
      LAMTMT=DBLE(MTMT/VEV)
C
C     Here, conversion from MSbar to DRbar is done at MZ.
C     Effect of sparticles is included by decoupling
C     beta functions in RGEs
      DAS=DBLE(ASMZ)/2.D0/DPI*(.5)
      ALEM=1./137.036
      DAEM=0.0682-ALEM/2./PI*(-7*LOG(AMW/AMZ))
      ALEMDR=ALEM/(1.-DAEM)
      XWMSB=.23113
C      XW=.2324-1.03E-7*(AMT**2-138.**2)
C     Convert XW to DRbar
      XW=1.-(1.-XWMSB)/(1.-ALEMDR/12./PI)      
      A1MZ=5*ALEMDR/3./(1.-XW)
      A2MZ=ALEMDR/XW
C      ALEM=1./128.
C      A1MZ=5*ALEM/3./(1.-XW)
C      A2MZ=ALEM/XW
      G(1)=DSQRT(4*DPI*A1MZ)
      G(2)=DSQRT(4*DPI*A2MZ)
      G(3)=DSQRT(4*DPI*ASMZ/(1.D0-DAS))
      G(4)=DBLE(FTAMZ)
      G(5)=DBLE(FBMZ)
      G(6)=G(6)
      G(25)=DBLE(MU)
      G(26)=DBLE(B)
      G(27)=0.D0
      G(28)=0.D0
      G(29)=0.D0
      G(30)=DBLE(VP)
      G(31)=DBLE(V)
C          Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
      END IF
C
C          Run back up to mgut with approximate susy spectra
C
      IF (IMODEL.EQ.1) THEN
        IF (XSUGIN(7).EQ.0.) THEN 
          MGUT=1.E19
        ELSE
          MGUT=XSUGIN(7)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
        MGUT=XMESGM
      END IF
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TGUT-TZ)/DBLE(FLOAT(NSTEP))
      Q=MZ
      DO 250 II=1,NSTEP
        T=TZ+(TGUT-TZ)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        QNEW=SNGL(MGUT*DEXP(T+DT))
        IF (Q.LE.MT.AND.QNEW.GT.MT) G(6)=DBLE(FTMT)
C       Implement sparticle threshold corrections at Q=HIGFRZ
        IF (Q.LE.HIGFRZ.AND.QNEW.GT.HIGFRZ) THEN
          G(6)=G(6)/(1.D0-DBLE(RSIGT))
          G(5)=G(5)/(1.D0-DBLE(RSIGB))
          G(4)=G(4)/(1.D0-DBLE(RSIGL))
          IF (INUHM.EQ.1) THEN
            G(13)=DBLE(MHDSQ)
            G(14)=DBLE(MHUSQ)
          END IF
        END IF
        IF (Q.LE.XNRIN(2).AND.QNEW.GT.XNRIN(2)) THEN
          G(27)=DBLE(FNMZ)
          G(28)=DBLE(G0(28))
          G(29)=DBLE(G0(29))
        END IF
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
        A1I=SNGL(4*DPI/G(1)**2)
        A2I=SNGL(4*DPI/G(2)**2)
        A3I=SNGL(4*DPI/G(3)**2)
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (A1I.LT.A2I.AND.XSUGIN(7).EQ.0.) GO TO 30
250   CONTINUE
      IF (IMODEL.EQ.1.AND.XSUGIN(7).EQ.0.) THEN
        WRITE(LOUT,*) 'SUGRGE ERROR: NO UNIFICATION FOUND'
        NOGOOD=1
        GO TO 100
      END IF
30    IF (XSUGIN(7).EQ.0.) THEN
        MGUT=QNEW
      ELSE
        MGUT=XSUGIN(7)
      END IF
      AGUT=SNGL((G(1)**2/4.D0/DPI+G(2)**2/4.D0/DPI)/2.D0)
      GGUT=SQRT(4*PI*AGUT)
      AGUTI=1./AGUT
      FTAGUT=SNGL(G(4))
      FBGUT=SNGL(G(5))
      FTGUT=SNGL(G(6))
      IF (INUHM.EQ.1) THEN
        MHDSMG=SNGL(G(13))
        MHUSMG=SNGL(G(14))
      END IF
      MUMG=SNGL(G(25))
      BMG=SNGL(G(26))
      IF (XNRIN(2).LT.1.E19.AND.XNRIN(1).EQ.0.) THEN
C     IMPOSE FN-FT UNIFICATION
        FNGUT=SNGL(G(6))
      ELSE
        FNGUT=SNGL(G(27))
      END IF
      G3GUT=SNGL(G(3))
      MGUTSS=MGUT
      AGUTSS=AGUT
      GGUTSS=GGUT
C
C          Set GUT boundary condition
C
      DO 260 I=1,3
        IF (IMODEL.EQ.1) THEN
          G(I+6)=DBLE(MHF)
          G(I+9)=DBLE(A0)
        ELSE IF (IMODEL.EQ.2) THEN
          G(I+6)=XGMIN(11+I)*XGMIN(8)*THRG*(G(I)/4.D0/DPI)**2*XLAMGM
          G(I+9)=0.D0
        END IF
      IF (XNRIN(2).LT.1.E19) THEN
        G(27)=DBLE(FNGUT)
        G(28)=DBLE(XNRIN(4))**2
        G(29)=DBLE(XNRIN(3))
      ELSE
        G(27)=0.D0
        G(28)=0.D0
        G(29)=0.D0
      END IF
260   CONTINUE
C     OVERWRITE ALFA_3 UNIFICATION TO GET ALFA_3(MZ) RIGHT
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) G(3)=DBLE(GGUT)
      IF (IMODEL.EQ.1) THEN
        DO 270 I=13,24
          G(I)=DBLE(M0)**2
270     CONTINUE
      IF (INUHM.EQ.1) THEN
        G(13)=DBLE(MHDSMG)
        G(14)=DBLE(MHUSMG)
      END IF
C          Set possible non-universal GUT scale boundary conditions
      DO 280 I=1,6
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=DBLE(XNUSUG(I))
        END IF
280   CONTINUE
      DO 281 I=7,18
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=SIGN(1.,XNUSUG(I))*DBLE(XNUSUG(I))**2
        END IF
281   CONTINUE
        IF (IDTERM.EQ.1) THEN
          MDS=(MHDSMG-MHUSMG)/4.
          M10S=(MHDSMG+MHUSMG)/2.
          G(15)=DBLE(M0**2+MDS)
          G(16)=DBLE(M0**2-3*MDS)
          G(17)=G(16)
          G(18)=G(15)
          G(19)=G(15)
          G(20)=DBLE(M0**2+MDS)
          G(21)=DBLE(M0**2-3*MDS)
          G(22)=G(21)
          G(23)=G(20)
          G(24)=G(20)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
       XC=2*THRF*XLAMGM**2
       DY=DSQRT(3.D0/5.D0)*G(1)*XGMIN(11)
       G(13)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)+DBLE(XGMIN(9))-DY
       G(14)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)+DBLE(XGMIN(10))+DY
       G(15)=XC*(.6*XGMIN(12)*(G(1)/4.D0/DPI)**4)+2*DY
       G(16)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)-DY
       G(17)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+.6D0*XGMIN(12)*
     , (G(1)/4.D0/DPI)**4/9.D0)+2*DY/3.D0
       G(18)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+
     , .6D0*4*XGMIN(12)*(G(1)/4.D0/DPI)**4/9.D0)-4*DY/3.D0
       G(19)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+.75*XGMIN(13)*
     ,(G(2)/4.D0/DPI)**4+.6*XGMIN(12)*(G(1)/4.D0/DPI)**4/36.D0)+DY/3.D0
       G(20)=G(15)
       G(21)=G(16)
       G(22)=G(17)
       G(23)=G(18)
       G(24)=G(19)
      ELSE IF (IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
       G(1)=G(1)
       G(2)=G(2)
       G(3)=G(3)
       BLHAT=G(4)*(-9*G(1)**2/5.D0-3*G(2)**2+3*G(5)**2+4*G(4)**2)
       BBHAT=G(5)*(-7*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
       BTHAT=G(6)*(-13*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             6*G(6)**2+G(5)**2)
       G(7)=33*MHF*G(1)**2/5.D0/16.D0/DPI**2
       IF (IMODEL.EQ.10) THEN
         G(7)=G(7)+XAMIN(11)*MHF
       END IF
       G(8)=MHF*G(2)**2/16.D0/DPI**2
       G(9)=-3*MHF*G(3)**2/16.D0/DPI**2
       G(10)=-BLHAT*MHF/G(4)/16.D0/DPI**2
       G(11)=-BBHAT*MHF/G(5)/16.D0/DPI**2
       G(12)=-BTHAT*MHF/G(6)/16.D0/DPI**2
       G(13)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(5)*BBHAT+
     ,G(4)*BLHAT)*MHF**2/(16*DPI**2)**2+XAMIN(6)*DBLE(M0)**2
       G(14)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(6)*BTHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(7)*DBLE(M0)**2
       G(15)=(-198*G(1)**4/25.D0)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(5)*DBLE(M0)**2
       G(16)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(4)*DBLE(M0)**2
       G(17)=(-22*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(2)*DBLE(M0)**2
       G(18)=(-88*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(3)*DBLE(M0)**2
       G(19)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
       G(20)=(-198*G(1)**4/25.D0+2*G(4)*BLHAT)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(5)*DBLE(M0)**2
       G(21)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+G(4)*BLHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(4)*DBLE(M0)**2
       G(22)=(-22*G(1)**4/25.D0+8*G(3)**4+2*G(5)*BBHAT)*
     , MHF**2/(16*DPI**2)**2+XAMIN(2)*DBLE(M0)**2
       G(23)=(-88*G(1)**4/25.D0+8*G(3)**4+2*G(6)*BTHAT)*
     , MHF**2/(16*DPI**2)**2+XAMIN(3)*DBLE(M0)**2
       G(24)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4+G(5)*BBHAT+
     ,        G(6)*BTHAT)*MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
      END IF
      IF (IMODEL.EQ.9) THEN
        CALL MMAMSB(M0,MHF,G)
      END IF
      IF (IMODEL.EQ.12) THEN
        CALL GENMIR(XAMIN(11),MHF,G,INUHM)
      END IF
      IF (IMODEL.EQ.13) THEN
        CALL NAMSB(MHF,G,INUHM)
      END IF
      DO 285 I=1,31
        IG(I)=0
285   CONTINUE
C          Check for tachyonic sleptons at GUT scale
      IF (G(15).LT.0.D0.OR.G(16).LT.0.D0) THEN
        ITACHY=2
      ELSE
        ITACHY=0
      END IF
C     SET MHDSMG, MHUSMG FOR HS FINETUNING CALCULATION
      MHDSMG=SNGL(G(13))
      MHUSMG=SNGL(G(14))
C
C          Run back down to weak scale
C
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TZ-TGUT)/DBLE(FLOAT(NSTEP))
      DO 290 II=1,NSTEP+2
        T=TGUT+(TZ-TGUT)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
C       Here, DRKSTP advances T by DT
        QNEW=SNGL(MGUT*DEXP(T))
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $    G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        CALL SUGFRZ(QNEW,G,G0,IG)
        IF (Q.GE.AMNRMJ.AND.QNEW.LT.AMNRMJ.AND.XNRIN(1).EQ.0.) THEN
          FNMZ=SNGL(G(27))
        END IF
        IF (Q.GT.HIGFRZ.AND.QNEW.LE.HIGFRZ) THEN
          G(6)=G(6)*(1.D0-DBLE(RSIGT))
          G(5)=G(5)*(1.D0-DBLE(RSIGB))
          G(4)=G(4)*(1.D0-DBLE(RSIGL))
        END IF
        IF (QNEW.LT.AMNRMJ) THEN
          G(27)=0.D0
          G(28)=0.D0
          G(29)=0.D0
        END IF
        IF (NOGOOD.NE.0) GO TO 100
        IF (QNEW.LT.MZ) GO TO 40
290   CONTINUE
40    CONTINUE
C
C          Electroweak breaking constraints; tree level
C
      VUQ=G0(31)
      VDQ=G0(30)
      TANBQ=VUQ/VDQ
      SIN2BQ=SIN(2*ATAN(TANBQ))
      MZQ=SQRT((G0(2)**2+.6*G0(1)**2)*(VUQ**2+VDQ**2)/2.)
      BADMU=.FALSE.
      IF (INUHM.NE.1) THEN
      MUS=(G0(13)-G0(14)*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          Compute loop corrections using masses from last iteration
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          If MUS<0, set it to MZ**2 and continue
      IF (MUS.LT.0.) THEN
        MUS=AMZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(G0(13)+G0(14)+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
      CALL SUGMAS(G0,0,IMODEL,SIGA)
      IF (NOGOOD.NE.0) GO TO 100
C
C           Electroweak breaking constraints; loop level
C
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
      IF (MUS.LT.0.) THEN
C        NOGOOD=2
C        GO TO 100
         MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
C
C     Once more, with feeling!
C
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
      IF (MUS.LT.0.) THEN
C        NOGOOD=2
C        GO TO 100
         BADMU=.TRUE.
         MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
      CALL SUGMAS(G0,1,IMODEL,SIGA)
      ELSE
        MUS=MU**2
        B=AMHA**2/MU/(COTB+TANB)+SIG3/MU
        CALL SUGMAS(G0,1,IMODEL,SIGA)
        CALL SUGEFF(G0,SIG1,SIG2,SIG3)
        MHDSQ=(AMHA**2+MZQ**2)*TANB**2/(TANB**2+1.)-SIG1-MUS-MZQ**2/2.
	MHUSQ=(AMHA**2+MZQ**2)/(TANB**2+1.)-MUS-MZQ**2/2.-SIG2
      END IF
C
C  Save radiative corrections to Yukawas for b->s gamma computation
C
      RTISA=RSIGT
      RBISA=RSIGB
      RLISA=RSIGL
C
C Initial value of MU passed to A. Box RGE code
C
      RGEMU=-MU !THE -VE SIGN FIXES CONVENTIONS
C

100   RETURN
      END
