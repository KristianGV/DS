#include "PILOT.inc"
#ifdef ISATOOLS_X
C  03/09/2012
C
C  Compute BR(B_u -> tau+nu_tau)  and R=B^MSSM/B^SM
C
C Matching with ISAJET notation
C
C      mhc = mss(30)
C      mg = mss(1)
C      mb1 = mss(10)
C      mb2 = mss(11)
C      m2 = gss(8)
C      aab = gss(11)
C
C     CALL ISABTN(TANB,MSS(32),MU,MSS(1),MSS(10),MSS(11),R,BFBTN)
C-------------------------------------------------------------------------
      SUBROUTINE ISABTN(TANB,MHC,MU,MGL,MB1,MB2,R,BFBTN)
      IMPLICIT NONE
      REAL TANB,MHC,MU,MGL,MB1,MB2
      REAL MBU,TAUBU,FBU,PSGEV,VUB,GF,ASMZ,EP0
      REAL PI
      PARAMETER(PI=3.1415926)
      REAL MTAU,BFSM,R,BFBTN
      REAL*8 H2,X,Y
C
C*************************************************************************
C
      MTAU = 1.777
C
C B_U meson
      MBU = 5.2791
      FBU = 0.216
      TAUBU = 1.643   !(IN UNITS OF PS)
      PSGEV = 1.E+13/6.582
      ASMZ=0.118
C
C CKM matrix
      VUB = 3.95E-3
      GF = 1.166E-5
C
C BRANCHING RATIOS
C
      BFSM=GF**2*MBU*MTAU**2/8.D0/PI*(1.D0-MTAU**2/MBU**2)**2*
     ,FBU**2*VUB**2*TAUBU*PSGEV
C
      X=DBLE((MB1/MGL)**2)
      Y=DBLE((MB2/MGL)**2)
      EP0=-2*ASMZ*MU/3.0/PI/MGL*SNGL(H2(X,Y))
      R=(1.0-(MBU/MHC)**2*TANB**2/(1.0-EP0*TANB))**2
      BFBTN=R*BFSM
C
      RETURN
      END
C
C
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION H2(X,Y)
C
      REAL*8 X,Y
C
      H2=X*DLOG(X)/(1.D0-X)/(X-Y)+Y*DLOG(Y)/(1.D0-Y)/(Y-X)
      RETURN
      END
C
#endif
