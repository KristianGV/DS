#include "PILOT.inc"
#ifdef ISATOOLS_X
C  04/27/02
C
C  BR(B_s -> mu+ m-)  and BR(B_d -> tau+ tau-)
C
C  Formula (4), (18) and  CORRECT (19)  from Borzumati, Huang, Kolda 
C  and Logan (BHKL)
C    
C  Eq. (19) has a wrong sign
C
C  Warning!   ISAJET notation: sin(alpha) -> - sin(alpha)
C                              cos(alpha) ->   cos(alpha)
C
C  b-mass correction from Pierce at al., Nucl. Phys. B491, 3 (1997)
C  programmed in the ISAJET gives (in leading tan\beta)
C
C  del_mb/mb = epsg tanb + epsu y_t^2 tanb + aw tanb (c_t^2 C_0(...) +
C  s_t^2 C_0(...))
C
C  ag * C_0(...)   = epsg,   ag = -2 alphas/(3pi) \mu mg 
C  au * C_0(...)   = epsu,   au = -g^2/(16 pi^2) \mu A_t
C                            aw = g^2/(16 pi^2) \mu M_2
C
C Matching with ISAJET notation
C
C      mhl = mss(29)
C      mhh = mss(30)
C      mha = mss(31)
C      mg = mss(1)
C      mb1 = mss(10)
C      mb2 = mss(11)
C      mt1 = mss(12)
C      mt2 = mss(13)
C      msl = mss(6)
C      mcl = mss(8) 
C      mer = mss(18)
C      mtau1 = mss(21)
C      mz1 = mss(23)
C      mw1 = mss(27)
C      m2 = gss(8)
C      aab = gss(11)
C      aat = gss(12)
C      mu = mu  WARNING: gss(25) does not give the same value as mu
C      alphah = alfah
C
C-------------------------------------------------------------------------
      SUBROUTINE ISABMM(MT,TANB,MHL,MHH,MHA,MG,MB1,MB2,AAB,
     &THETAB,MT1,MT2,AAT,THETAT,MU,M2,ALPHAH,MSL,MCL,BRBS,BRBD)
      IMPLICIT NONE
      REAL MBS, TAUBS, FBS, VTB, VTS, CQ1MU, CQ2MU, GF
      REAL MBD, TAUBD, FBD, VTD, CQ1TAU, CQ2TAU
      REAL PI
      PARAMETER(PI=3.1415926)
      REAL MZ, MW, G, G2, SW, CW, SW2, VEV, MB, MBMB, MTMT, MMU
      REAL MTAU, ALPHA, E, MBMZ
      REAL AM2, AMU, BETA, SINB, COSB, SINA, COSA
      REAL FTMT,  ASMB, ASMT, ASMZ, FBMZ
      REAL COS2B, COSBA, SINBA
      REAL EPSG, EPSU, EPSUA, EPSUB, CHIFC, CHIFCBABU
      REAL SUALFS, SUALFE
      COMPLEX*16 SSB0, SSB1, ZZZ
      REAL*8 REAL8
      REAL PSGEV, C0, D0
      REAL EPSGBABU, EPSUBABU, TERM1, TERM2, TERM3, CCOUNT
      REAL C10, AG, AW, AU, DEN, NUM, TEST1, TEST2
      REAL MT, TANB, MHL, MHH, MHA, MG, MB1, MB2, AAB, 
     &       THETAB, MT1, MT2, AAT, THETAT, MU,  M2, ALPHAH,
     &       MSL, MCL, BRBS, BRBD
C
C*************************************************************************
C
      REAL8(ZZZ) = DREAL(ZZZ)
C
      MB = 4.9
      MZ = 91.188
      MW = 80.419
      MMU = 0.1056
      MTAU = 1.777
C
C B_s meson
      MBS = 5.3696
      FBS = 0.250  
      TAUBS = 1.493   !(IN UNITS OF PS)
C
C B_d meson
      MBD = 5.2794
      FBD = 0.208
      TAUBD = 1.548   !(IN UNITS OF PS)
C
      PSGEV = 1.E+13/6.582
      TAUBS = TAUBS*PSGEV
      TAUBD = TAUBD*PSGEV
C
C CKM matrix
      VTB = 0.999
      VTS = 0.039
      VTD = 0.009
C
CCC      SW2 = 0.2309
      SW2 = .2324-1.03E-7*(MT**2-138.**2)
CC      PRINT*, 'SW2=', SW2
      SW  = SQRT(SW2)
      CW  = SQRT(1. - SW2)
      VEV = 246./SQRT(2.)
      ALPHA = 1./128.
      E = SQRT(4.*PI*ALPHA)
      G = E/SW
      G2 = G
CC      GF = SQRT(2.)*G**2/8./MW**2
      GF = 1.166E-5
C
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      COS2B=COS(2*BETA)
C
C Running mass for B
C
      ASMB = SUALFS(MB**2, .36, MT, 3)
CC      MBMB = MB*(1. - 4.*ASMB/3./PI)
      MBMB=MB*(1.-4*ASMB/3./PI-12.4*(ASMB/PI)**2)
C
      ASMZ=0.118
C
      MBMZ = 2.92
C
      COSA = COS(ALPHAH)
      SINA = SIN(ALPHAH)
      COSBA = COSB*COSA - SINB*SINA
      SINBA = SINB*COSA + COSB*SINA
      AMU = ABS(MU)
      AM2 = ABS(M2)
C
      FTMT = 0.94
C
C-------------------------------------------------------------------
C
C squark/gluino contribution
C Pierce et al., second term of the Eq. (8) with t -> b
C
C In an effective Lagrangian method
C
C                               H_u
C                                
C       _    H_d                 |
C      / \                      / \
C     /   \  |                 /   \
C  ---+++++--X---          ---/+++++\----
C  Q_L     Q_L  d_R        Q_L        d_R
C
C      Diagram 1             Diagram 2
C
C In the leading tan\beta contribution Pierce et al. and effective 
C Lagrangian, with sq_L sq_R H_u vertex, give same answer. In this 
C approximaion, set 
C 
C \sin 2\theta  = 2 m_b (\mu tanb - A_b)/(...) -> 2 m_b \mu tanb /(...)
C
C----------------------------------------------------------------------
C
C According to Pierce's definition, f_1 is mostly left-handed
C  
C M2 -> -M2, MG -> -MG, sc -> -sc
C squark/gluino contribution, diagram 2 contribution
C
ccc      delmb1 = -ASMZ/3./PI*SIN(2*THETAB)*(-MG)/MBMZ*
ccc     $        (REAL8(SSB0(MZ**2,MG,MB1))-REAL8(SSB0(MZ**2,MG,MB2)))
C
      AG = -2.*ASMZ/3./PI*MU*MG
      EPSG = AG*C0(MG,MB1,MB2)
C
C---------------------------------------------------------------------
C squark/chargino contribution
C
C       H_u                       
C                                 
C        |                      _  st_L, st_R
C st_R  / \  sQ_L              / \
C      /   \                  /   \
C   ---+++++----          ---/++X++\----
C   Q_L      d_R          Q_L        d_R
C       chi1, chi2           chi1 chi2  mass insertion
C
C      Diagram 1             Diagram 2
C
C----------------------------------------------------------------------
C 
C M2 -> -M2, MG -> -MG
C squark/chargino contribution, diagram 1 (delmb2) + 2 (delmb3)
C
CC       delmb2 = FTMT**2*MU*(-AAT*TANB+MU)/16./PI**2/(MT1**2-MT2**2)*
CC     $     (REAL8(SSB0(MZ**2,AMU,MT1))-REAL8(SSB0(MZ**2,AMU,MT2)))
C
CC       delmb3 = G2**2*MU*(-M2)*TANB/16./PI**2/(AMU**2-M2**2)*
CC     $     (COS(THETAT)**2*(REAL8(SSB0(MZ**2,AM2,MT1))-
CC     $     REAL8(SSB0(MZ**2,AMU,MT1)))+
CC     $     SIN(THETAT)**2*(REAL8(SSB0(MZ**2,AM2,MT2))-
CC     $     REAL8(SSB0(MZ**2,AMU,MT2))))
C
      AU = -1./16./PI**2*MU*AAT
      EPSU = AU*C0(MU,MT1,MT2)
C
      AW = G2**2/16./PI**2*MU*M2
C
      NUM = EPSU*FTMT**2 + AG*(C0(MG,MB1,MB2) - SIN(THETAB)**2*
     &      C0(MG,MSL,MB1) - COS(THETAB)**2*C0(MG,MSL,MB2)) +
     &      AW*(COS(THETAT)**2*C0(MU,M2,MT1) + SIN(THETAT)**2*
     &      C0(MU,M2,MT2) - C0(MU,M2,MCL))
C
      DEN = (1. + AG*TANB*(SIN(THETAB)**2*C0(MG,MSL,MB1) +
     &      COS(THETAB)**2*C0(MG,MSL,MB2)) + AW*TANB*C0(MU,M2,MCL))*
     &      (1. + EPSG*TANB + AW*TANB*(COS(THETAT)**2*C0(MU,M2,MT1) + 
     &      SIN(THETAT)**2*C0(MU,M2,MT2)) + EPSU*FTMT**2*TANB)
C
      CHIFC = -NUM*TANB/DEN
C
C----------------------------------------------------------------------
C FINAL STATE WITH MUONS
C
      CQ1MU = -2.*PI/ALPHA*CHIFC*MBMB*MMU/COSB**2/SINB*
     &      (-COSBA*SINA/MHL**2 + SINBA*COSA/MHH**2)

      CQ2MU = -2.*PI/ALPHA*CHIFC*MBMB*MMU/COSB**2*(-1./MHA**2)
C
C WITH TAUS
C
      CQ1TAU = -2.*PI/ALPHA*CHIFC*MBMB*MTAU/COSB**2/SINB*
     &      (-COSBA*SINA/MHL**2 + SINBA*COSA/MHH**2)

      CQ2TAU = -2.*PI/ALPHA*CHIFC*MBMB*MTAU/COSB**2*(-1./MHA**2)
C
C------------------------------------------------------------------------
C BRANCHING RATIOS
C
      C10 = -4.221
C
C B_S -> MU+ MU-
C
      BRBS = GF**2*ALPHA**2*MBS**3*TAUBS*FBS**2/64./PI**3*(VTB*VTS)**2*
     &     SQRT(1.-4.*MMU**2/MBS**2)*
     &     ((1.-4.*MMU**2/MBS**2)*CQ1MU**2 + 
     &      (CQ2MU - 2.*MMU/MBS*C10)**2)
C
C B_D -> TAU+ TAU-
C
      BRBD = GF**2*ALPHA**2*MBD**3*TAUBD*FBD**2/64./PI**3*(VTB*VTD)**2*
     &     SQRT(1.-4.*MTAU**2/MBD**2)*
     &     ((1.-4.*MTAU**2/MBD**2)*CQ1TAU**2 + 
     &      (CQ2TAU - 2.*MTAU/MBD*C10)**2)
C
      RETURN
      END
C
C
C**********************************************************************      
C Given by Pierce, (Eq. C19)
C
C----------------------------------------------------------------------
      FUNCTION C0(M1,M2,M3)
C
      REAL C0, M1, M2, M3
      REAL*8 X, Y, Z, X2, Y2, Z2, FRACXY, FRACXZ, FRACYZ
C
      X = M1*M1
      Y = M2*M2
      Z = M3*M3
C
      FRACXY = ABS(X-Y)/(X+Y)
      FRACXZ = ABS(X-Z)/(X+Z)
      FRACYZ = ABS(Y-Z)/(Y+Z)
      IF(FRACXY.LT.1.D-4) Y = X
      IF(FRACXZ.LT.1.D-4) Z = X
      IF(FRACYZ.LT.1.D-4) Y = Z
C
      IF(X.NE.0.AND.Y.NE.0.AND.Z.NE.0)THEN
         IF(X.NE.Y.AND.X.NE.Z.AND.Y.NE.Z)THEN
            C0 = 1.D0/(Y-Z)*(Y/(X-Y)*DLOG(Y/X)-Z/(X-Z)*DLOG(Z/X))
         ELSE
            IF(X.EQ.Y.AND.X.NE.Z)THEN
               C0 = -(1.D0 - Z/(Y-Z)*DLOG(Y/Z))/(Y-Z)
            ELSE
               IF(X.EQ.Y.AND.X.EQ.Z)THEN
                  C0 = -1.D0/2.D0/Z
               END IF
            END IF
         END IF
         IF(X.NE.Y)THEN
            IF(X.EQ.Z)THEN
               C0 = -(1.D0 - Y/(X-Y)*DLOG(X/Y))/(X-Y)
            ELSE
               IF(Y.EQ.Z)THEN
                  C0 = -(1.D0 - X/(Y-X)*DLOG(Y/X))/(Y-X)
               END IF
            END IF
         END IF
      ELSE
         WRITE(*,*) 'C0 COM ARGUMENTO NULO!' 
      END IF
      RETURN
      END
C
C----------------------------------------------------------------------
      FUNCTION D0(M1,M2,M3,M4)
      REAL M1, M2, M3, M4, X, Y, Z, W, C0, FRACXY, FRACXZ, FRACXW
      INTEGER FLAG
C
      X = M1*M1
      Y = M2*M2
      Z = M3*M3
      W = M4*M4
      FRACXY = (X-Y)/(X+Y)
      FRACXZ = (X-Z)/(X+Z)
      FRACXW = (X-W)/(X+W)
C
      IF(FRACXY.LT.1.E-2) Y = X
      IF(FRACXZ.LT.1.E-2) Z = X
      IF(FRACXW.LT.1.E-2) W = X
C
      FLAG = 0
      IF(X.NE.0.AND.Y.NE.0.AND.Z.NE.0.AND.W.NE.0)THEN
         IF(X.NE.Y)THEN
            D0 = (C0(X,Z,W) - C0(Y,Z,W))/(X-Y)
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.NE.Z.AND.X.NE.W)THEN
            D0 = -(-(Y**2-Z*W)/(Y-Z)**2/(Y-W)**2*LOG(Y/W) +
     &           Z/(Y-Z)**2/(Z-W)*LOG(Z/W) + 1./(Y-Z)/(Y-W))
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.EQ.Z.AND.X.NE.W)THEN
            D0 = -(1./2./Y - 1./(Y-W) + W/(Y-W)**2*LOG(Y/W))/(Y-W)
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.NE.Z.AND.X.EQ.W)THEN
            D0 = -(1./(X-Z)*(-1./2./X - Z/X/(X-Z)) + 
     &           Z/(X-Z)**3*LOG(X/Z))
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.EQ.Z.AND.X.EQ.W)THEN
            D0 = 1./6./W**2
            FLAG = 1
         END IF
      END IF
C
      IF(FLAG.EQ.0) WRITE(*,*) 'D0: ARGUMENTOS ESTRANHOS'
      RETURN
      END
C
C**********************************************************************
C FUNCTION B0 FROM PIERCE, PAG. 12
C
      FUNCTION B0FUNC(SCALE,M1,M2)
      REAL B0FUNC, SCALE, M1, M2, MC2, ML2
C
      IF(M1.GT.M2)THEN
         MC2 = M1*M1
         ML2 = M2*M2
      ELSE
         MC2 = M2*M2
         ML2 = M1*M1
      END IF
C
      B0FUNC = -LOG(MC2/SCALE**2) + 1. + ML2/(ML2-MC2)*LOG(MC2/ML2)
C
      RETURN
      END
C
#endif
