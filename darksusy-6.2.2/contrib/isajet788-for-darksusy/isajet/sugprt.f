#include "PILOT.inc"
C--------------------------------------------------------------------
      SUBROUTINE SUGPRT(IMODEL,IMODIN)
C--------------------------------------------------------------------
C
C     Print SUGRA parameters and results
C     IMODEL = model type for SUGRA
C     IMODIN = input model type to control formatting
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sugxin.inc"
#include "sugmg.inc"
#include "sugpas.inc"
#include "sspar.inc"
#include "sugnu.inc"
#include "kphggs.inc"
      REAL PI,GPX,SIN2W,ALEMI,AS,TANBQ
      REAL PCEWFT,PCHSFT,DELEW,DELHS
      REAL M10,MD
      INTEGER IMODEL,J,K,IMODIN
C
C          Entry
C
      PI=4.*ATAN(1.)
      GPX=SQRT(.6)*GSS(1)
      SIN2W=GPX**2/(GSS(2)**2+GPX**2)
      ALEMI=4*PI/GSS(2)**2/SIN2W
      AS=GSS(3)**2/4./PI
      TANBQ=VUQ/VDQ
C
C          Print inputs and GUT couplings for SUGRA/AMSB models
C
      IF(IMODEL.EQ.1.OR.IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10
     $.OR.IMODEL.EQ.12.OR.IMODEL.EQ.13) THEN
        IF(IMODEL.EQ.1) THEN
          WRITE(LOUT,1000) XSUGIN(1),XSUGIN(2),XSUGIN(3),XSUGIN(4),
     $    XSUGIN(5),XSUGIN(6)
1000      FORMAT(
     $    ' M_0,  M_(1/2),  A_0,  tan(beta),  sgn(mu),  M_t ='
     $    /4F12.3,2X,F4.1,2X,F8.3)
        ELSE IF (IMODEL.EQ.7) THEN
          WRITE(LOUT,1018) XSUGIN(1),XSUGIN(2),XSUGIN(4),XSUGIN(5),
     $    XSUGIN(6)
1018      FORMAT(
     $    ' M_0,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
     $    /3F10.3,2X,F6.1,2F10.3)
        ELSE IF (IMODEL.EQ.9) THEN
          WRITE(LOUT,1019) XSUGIN(1),XSUGIN(2),XSUGIN(4),XSUGIN(5),
     $    XSUGIN(6)
1019      FORMAT(
     $    ' alpha,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
     $    /3F10.3,2X,F6.1,2F10.3)
          WRITE(LOUT,1020) XAMIN(1),XAMIN(2),XAMIN(3),XAMIN(4),XAMIN(5),
     $    XAMIN(6),XAMIN(7),XAMIN(8),XAMIN(9),XAMIN(10)
1020      FORMAT(
     $    ' Moduli nQ, nD, nU, nL, nE, nHd, nHu, L1, L2, L3 ='
     $    /10F4.1)
        ELSE IF (IMODEL.EQ.10) THEN
          WRITE(LOUT,1040) XAMIN(11),XSUGIN(2),XSUGIN(4),XSUGIN(5),
     $    XSUGIN(6)
1040      FORMAT(
     $    ' alpha,  M_(3/2),  tan(beta),  sgn(mu),  M_t ='
     $    /3F10.3,2X,F6.1,2F10.3)
        ELSE IF (IMODEL.EQ.12) THEN
          WRITE(LOUT,1041) XAMIN(11),XSUGIN(2),XAMIN(1),XAMIN(5),
     $XAMIN(4),XSUGIN(4),XSUGIN(5),XSUGIN(6)
1041  FORMAT(' alpha, m(3/2), c_m, c_m3, a_3, tan(beta), sgn(mu), m_t ='
     $/F4.1,2X,F8.1,1X,F6.2,1X,F5.2,2X,F5.2,2X,F5.2,7X,F3.1,5X,F7.2)
        ELSE IF (IMODEL.EQ.13) THEN
          WRITE(LOUT,1042) XAMIN(10),XAMIN(11),XSUGIN(2),XAMIN(1),
     $XSUGIN(4),XSUGIN(5),XSUGIN(6)
1042  FORMAT(
     $' m0(12),   m0(3),  m(3/2),    A_0,   tan(beta), sgn(mu), m_t ='
     $/F8.1,1X,F8.1,1X,F10.2,1X,F8.2,2X,F5.2,4X,F3.1,4X,F7.2)
          IF (INUHM.EQ.1) THEN
            WRITE(LOUT,*) 'mu= ',XNUSUG(19),'  mA= ',XNUSUG(20)
          ELSE
            WRITE(LOUT,*) 'mHu(B)= ',XAMIN(2),' mHd(B)= ',XAMIN(3)
          END IF
        END IF
C
C          Write out non-universal GUT scale parameters
        IF (INUHM.EQ.1) THEN
          WRITE(LOUT,1023)
          WRITE(LOUT,1021) MHDSQ,MHUSQ
          WRITE(LOUT,1022) MHDSMG,MHUSMG
1021      FORMAT(/,' M_Hd^2(Q)= ',E10.3,3X,'M_Hu^2(Q)= ',E10.3)
1022      FORMAT(' M_Hd^2(MGUT)= ',E10.3,3X,'M_Hu^2(MGUT)= ',E10.3)
1023      FORMAT(/,' NUHM model has been selected:')
        END IF
        IF (IDTERM.EQ.1) THEN
          WRITE(LOUT,1024)
          M10=SIGN(1.,M10S)*SQRT(ABS(M10S))
          MD=SIGN(1.,MDS)*SQRT(ABS(MDS))
          WRITE(LOUT,1025) M10,MD
1025      FORMAT(/,' M(10)= ',E10.3,3X,'MD= ',E10.3)
1024      FORMAT(/,' NUHM-d-TERM model has been selected:')
        END IF
        IF(XNUSUG(1).LT.1.E19.OR.XNUSUG(2).LT.1.E19.OR.XNUSUG(3)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1010) XNUSUG(1),XNUSUG(2),XNUSUG(3)
1010      FORMAT(/' M_1(GUT)= ',F8.2,'    M_2(GUT)= ',F8.2,
     $    '    M_3(GUT)= ',F8.2) 
        END IF
        IF(XNUSUG(4).LT.1.E19.OR.XNUSUG(5).LT.1.E19.OR.XNUSUG(6)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1011) XNUSUG(4),XNUSUG(5),XNUSUG(6)
1011      FORMAT(/' A_tau(GUT)= ',F8.2,'    A_b(GUT)= ',F8.2,
     $    '    A_t(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(7).LT.1.E19.OR.XNUSUG(8).LT.1.E19) THEN
          WRITE(LOUT,1012) XNUSUG(7),XNUSUG(8)
1012      FORMAT(/' M_Hd(GUT)= ',F8.2,'    M_Hu(GUT)= ',F8.2)
        END IF
        IF (XNUSUG(9).LT.1.E19.OR.XNUSUG(10).LT.1.E19) THEN
          WRITE(LOUT,1013) XNUSUG(9),XNUSUG(10)
1013      FORMAT(/' M_eR(GUT)= ',F8.2,'    M_eL(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(11).LT.1.E19.OR.XNUSUG(12).LT.1.E19.OR.XNUSUG(13)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1014) XNUSUG(11),XNUSUG(12),XNUSUG(13)
1014      FORMAT(' M_dR(GUT)= ',F8.2,'    M_uR(GUT)= ',F8.2,
     $    '    M_uL(GUT)=',F8.2)
        END IF
        IF(XNUSUG(14).LT.1.E19.OR.XNUSUG(15).LT.1.E19) THEN
          WRITE(LOUT,1015) XNUSUG(14),XNUSUG(15)
1015      FORMAT(/' M_tauR(GUT)= ',F8.2,'    M_tauL(GUT)= ',F8.2)
        END IF
        IF(XNUSUG(16).LT.1.E19.OR.XNUSUG(17).LT.1.E19.OR.XNUSUG(18)
     $  .LT.1.E19) THEN
          WRITE(LOUT,1016) XNUSUG(16),XNUSUG(17),XNUSUG(18)
1016      FORMAT(' M_bR(GUT)= ',F8.2,'    M_tR(GUT)= ',F8.2,
     $    '    M_tL(GUT)=',F8.2)
        END IF
        IF(XSUGIN(7).NE.0) THEN
          WRITE(LOUT,1026) XSUGIN(7)
1026      FORMAT(' Q_max= ',E12.4)
        ENDIF
C
C          Right-handed neutrino parameters
        IF (XNRIN(2).LT.1.E19) THEN
          WRITE(LOUT,1017) XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4),
     $    FNMZ,FNGUT
1017      FORMAT(' Right-handed neutrino parameters:'/
     $    ' M(nu_tau)=',E10.3,'   M(N_R) =',E10.3,
     $    '   A_N=',F8.2,'   M(NRSS)=',F8.2/
     $    ' FN(M_Z)  =',F8.4, '   FN(M_GUT) =',F8.4)
        END IF
C
C          Unification results
        WRITE(LOUT,1001) MGUTSS,GGUTSS,AGUTSS
1001    FORMAT(/' ISASUGRA unification:'/' M_GUT      =',E10.3,
     $  '   g_GUT          =',F5.3,3X,'   alpha_GUT =',F5.3)
        WRITE(LOUT,999) FTGUT,FBGUT,FTAGUT
999     FORMAT(' FT_GUT     =',F6.3,
     $  '       FB_GUT         =',F6.3,3X,'  FL_GUT =',F6.3)
C
C          Print inputs for GMSB models
C
      ELSE IF (IMODEL.EQ.2) THEN
        WRITE(LOUT,1002) (XGMIN(J),J=1,7)
1002    FORMAT(
     $  ' Lambda,  M_mes,  N_5,  tan(beta),  sgn(mu),  M_t,  C_grav='
     $  /2E10.3,2F10.3,2X,F6.1,F10.3,1X,E10.3)
        WRITE(LOUT,1030) (XGMIN(J),J=8,14)
1030    FORMAT(/' GMSB2 model input:'/
     $  ' Rsl,    dmH_d^2,   dmH_u^2,     d_Y,     N5_1,  N5_2,  N5_3='
     $  /F7.3,1X,E10.3,1X,E10.3,1X,E10.3,2X,3F7.3)
        WRITE(LOUT,1003) AMGVSS
1003    FORMAT(/' M(gravitino)=',E10.3)
      END IF
C
C          Weak scale couplings
C
      WRITE(LOUT,1004) ALEMI,SIN2W,AS
1004  FORMAT(/' 1/alpha_em =',F8.2,2X,
     $'   sin**2(thetaw) =',F6.4,2X,'   a_s^DRB   =  ',F5.3)
      WRITE(LOUT,1005) GSS(7),GSS(8),GSS(9)
1005  FORMAT(' M_1        =',F8.2,2X,
     $'   M_2            =',F8.2,'   M_3       =',F8.2)
      WRITE(LOUT,1006) MU,B,HIGFRZ
1006  FORMAT(' mu(Q)      =',F8.2,2X,
     $'   B(Q)           =',F8.2,'   Q         =',F8.2)
      WRITE(LOUT,1007) GSS(13),GSS(14),TANBQ
1007  FORMAT(' M_Hd^2     =',E10.3,'   M_Hu^2         =',E10.3,
     $' TANBQ     =   ',F6.3)
C
C          Print mass spectrum from ISASUGRA 
C
      WRITE(LOUT,2000) MSS(1),MSS(2),MSS(3),MSS(4),MSS(5),MSS(10),
     $MSS(11),MSS(12),MSS(13),MSS(14),MSS(17),MSS(18),MSS(16),
     $MSS(21),MSS(22),MSS(23),MSS(24),MSS(25),MSS(26),MSS(27),
     $MSS(28),MSS(29),MSS(30),MSS(31),MSS(32)
2000  FORMAT(/' ISAJET masses (with signs):'/
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
      WRITE(LOUT,2001) THETAT,THETAB,THETAL,ALFAH
2001  FORMAT(/,' theta_t=',F9.4,'   theta_b=',F9.4,
     $'   theta_l=',F9.4,'   alpha_h=',F9.4)
C
C     Write out chargino /neutralino masses/eigenvectors
C
      WRITE(LOUT,3100) AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS
3100  FORMAT(/' NEUTRALINO MASSES (SIGNED) =',4F10.3)
      DO 100 J=1,4
        WRITE(LOUT,3200) J,(ZMIXSS(K,J),K=1,4)
3200    FORMAT(' EIGENVECTOR ',I1,'       =',4F10.5)
100   CONTINUE
      WRITE(LOUT,3300) AMW1SS,AMW2SS
3300  FORMAT(/' CHARGINO MASSES (SIGNED)  =',2F10.3)
      WRITE(LOUT,3400) GAMMAL,GAMMAR
3400  FORMAT(' GAMMAL, GAMMAR             =',2F10.5/)

C
C          Print ISAJET MSSMi equivalent input
C
      WRITE(LOUT,3000)
3000  FORMAT(/' ISAJET equivalent input:')
      WRITE(LOUT,3001) MSS(1),MU,MSS(31),XSUGIN(4)
3001  FORMAT(' MSSMA: ',4F9.2)
      WRITE(LOUT,3002) SQRT(GSS(19)),SQRT(GSS(17)),SQRT(GSS(18)),
     $SQRT(GSS(16)),SQRT(GSS(15))
3002  FORMAT(' MSSMB: ',5F9.2)
      WRITE(LOUT,3003) SIGN(1.,GSS(24))*SQRT(ABS(GSS(24))),
     $SQRT(GSS(22)),SIGN(1.,GSS(23))*SQRT(ABS(GSS(23))),
     $SQRT(GSS(21)),SQRT(GSS(20)),GSS(12),GSS(11),GSS(10)
3003  FORMAT(' MSSMC: ',5F9.2,3F11.2)
      WRITE(LOUT,3004)
3004  FORMAT(' MSSMD: SAME AS MSSMB (DEFAULT)')
      WRITE(LOUT,3005) GSS(7),GSS(8)
3005  FORMAT(' MSSME: ',2F8.2/)
      PCEWFT=100.*.5*91.17**2/C5MAX
      DELEW=2*C5MAX/91.17/91.17
      PCHSFT=100.*.5*91.17**2/CHSMAX
      DELHS=2*CHSMAX/91.17/91.17
      WRITE(LOUT,3006) PCEWFT,DELEW
3006  FORMAT(' EW finetuning:         ',F8.5,'% ;','   DelEW=',F8.2)
      WRITE(LOUT,3007) PCHSFT,DELHS
3007  FORMAT(' High scale finetuning: ',F8.5,'% ;','   DelHS=',F8.2/)
      IF(C5MAX.LT.0.) 
     &WRITE(LOUT,*)'WARNING: EWFT unreliable <- degenerate spartcles'
      WRITE(LOUT,3008) 
3008  FORMAT(/' Higgs h coupling strengths')
      WRITE(LOUT,*) 'KPB,KPT,KPtau,KPW,KPZ,KPGL,KPGM='
      WRITE(LOUT,*) KPB,KPT,KPL,KPW,KPZ,KPGL,KPGM
      RETURN
      END
