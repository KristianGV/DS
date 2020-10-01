#include "PILOT.inc"
C---------------------------------------------------------------
      SUBROUTINE SUGMAS(G0,ILOOP,IMODEL,SIGA)
C---------------------------------------------------------------
C
C     Compute tree level sparticle masses; output to MSS, XISAIN
C     Further tadpoles added to mA by Javier, 5/20/03
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
#include "sspar.inc"
#include "sssm.inc"
#include "sugpas.inc"
#include "sugxin.inc"
#include "sugmg.inc"
#include "sugnu.inc"
#include "ssinf.inc"
      REAL MSB1,MSB2,MST1,MST2,SIGA
      REAL G0(31)
      REAL SUGMFN,SUALFS,SSPOLE,MHP,MGLMGL,MHPS,
     $RDEL,ASMGL,DELHPS,M1S,M2S,FNB,FCN,
     $MB,FNT,MT,MW,TANB,BETA,COSB,COTB,SINB,MZ,COS2B,
     $PI,T2S,G,AL,MSSS,AT,AB,BRKT,B2S,T1S,TERM,B1S,Q,
     $FNL,MSL1,MSL2,COS2W,MAJXS,MAJX,
     $SIGST,SIGT,SIGSB,SIGB,SIGSL,SIGL,RZT,RZB,RZL
      REAL AA,BB,CC,DA,DB,DC,L1,L2,EVAL1,RL1,RL2
      REAL SIGHIG,SIGCHA
      DOUBLE PRECISION SSMQCD
      INTEGER IALLOW,ILOOP,IMODEL
C
C          Statement function
C
      SUGMFN(Q)=Q**2*(LOG(Q**2/HIGFRZ**2)-1.)
C
      MHPNEG=0
      PI=4.*ATAN(1.)
      XW=.232
      G=G0(2)
      COS2W=1.-SN2THW
      TANB=VUQ/VDQ
      MT=AMT
      MZ=AMZ
      MW=AMW
      AMTP=MT
      BETA=ATAN(TANB)
      COTB=1./TANB
      SINB=SIN(BETA)
      COSB=COS(BETA)
      SIN2B=SIN(2*BETA)
      COS2B=COS(2*BETA)
      AT=G0(12)
      AB=G0(11)
      AL=G0(10)
      MLQ=G0(4)*VDQ
      MBQ=G0(5)*VDQ
      MTQ=G0(6)*VUQ
C
C          Compute some masses from RGE solution to prepare for SSMASS,
C          which computes the rest.
C
      MSSS=G0(19)+AMUP**2+(.5-2*XW/3.)*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
C          Squark and slepton masses
      MSS(2)=SQRT(MSSS)
      MSSS=G0(18)+AMUP**2+2./3.*XW*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(3)=SQRT(MSSS)
      MSS(4)=SQRT(G0(19)+AMDN**2+(-.5+XW/3.)*MZ**2*COS2B)
      MSS(5)=SQRT(G0(17)+AMDN**2-1./3.*XW*MZ**2*COS2B)
      MSS(6)=SQRT(G0(19)+AMST**2+(-.5+XW/3.)*MZ**2*COS2B)
      MSS(7)=SQRT(G0(17)+AMST**2-1./3.*XW*MZ**2*COS2B)
      MSS(8)=SQRT(G0(19)+AMCH**2+(.5-2*XW/3.)*MZ**2*COS2B)
      MSS(9)=SQRT(G0(18)+AMCH**2+2./3.*XW*MZ**2*COS2B)
      BRKT=(.5*(G0(24)-G0(22))-COS2B*(4*MW**2-MZ**2)/12.)**2+
     $       MBQ**2*(AB-MU*TANB)**2
      TERM=.5*(G0(24)+G0(22))+MBQ**2-MZ**2*COS2B/4.
      B1S=TERM-SQRT(BRKT)
      B2S=TERM+SQRT(BRKT)
      MSS(10)=SQRT(MAX(0.,B1S))
      MSS(11)=SQRT(MAX(0.,B2S))
      THETAB=ATAN((B1S-MBQ**2+MZ**2*COS2B*(.5-XW/3.)-
     $G0(24))/MBQ/(AB-MU*TANB))
      BRKT=(.5*(G0(24)-G0(23))+COS2B*(8*MW**2-5*MZ**2)/12.)**2+
     $       MTQ**2*(AT-MU*COTB)**2
      TERM=.5*(G0(24)+G0(23))+MTQ**2+MZ**2*COS2B/4.
      T1S=TERM-SQRT(BRKT)
      IF (T1S.LE.0..OR.B1S.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      T2S=TERM+SQRT(BRKT)
      MSS(12)=SQRT(MAX(0.,T1S))
      MSS(13)=SQRT(MAX(0.,T2S))
      THETAT=ATAN((T1S-MTQ**2+MZ**2*COS2B*(-.5+2*XW/3.)-
     $G0(24))/MTQ/(AT-MU*COTB))
      MSSS=G0(16)+.5*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(14)=SQRT(MSSS)
      MSS(15)=MSS(14)
      MSSS=G0(21)+.5*MZ**2*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(16)=SQRT(MSSS)
      MSSS=G0(16)+AME**2-.5*(2*MW**2-MZ**2)*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(17)=SQRT(MSSS)
      MSSS=G0(15)+AME**2+(MW**2-MZ**2)*COS2B
      IF (MSSS.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      MSS(18)=SQRT(MSSS)
      MSS(19)=SQRT(G0(16)+AMMU**2-.5*(2*MW**2-MZ**2)*COS2B)
      MSS(20)=SQRT(G0(15)+AMMU**2+(MW**2-MZ**2)*COS2B)
      BRKT=(.5*(G0(21)-G0(20))-COS2B*(4*MW**2-3*MZ**2)/4.)**2+
     $       MLQ**2*(AL-MU*TANB)**2
      TERM=.5*(G0(21)+G0(20))+MLQ**2-MZ**2*COS2B/4.
      T1S=TERM-SQRT(BRKT)
      IF (T1S.LE.0.) THEN
        NOGOOD=1
        GO TO 100
      END IF
      T2S=TERM+SQRT(BRKT)
      MSS(21)=SQRT(MAX(0.,T1S))
      MSS(22)=SQRT(MAX(0.,T2S))
      THETAL=ATAN((T1S-MLQ**2+MZ**2*COS2B*(.5-XW)-
     $G0(21))/MLQ/(AL-MU*TANB))

C          A0 mass
      M1S=MU**2+G0(13)
      M2S=MU**2+G0(14)
      MSB1=MSS(10)
      MSB2=MSS(11)
      MST1=MSS(12)
      MST2=MSS(13)
      MSL1=MSS(21)
      MSL2=MSS(22)
      MB=AMBT
      FNT=(SUGMFN(MST2)-SUGMFN(MST1))/(MST2**2-MST1**2)
     $*AT*MTQ**2/SINB**2
      FNB=(SUGMFN(MSB2)-SUGMFN(MSB1))/(MSB2**2-MSB1**2)
     $*AB*MBQ**2/COSB**2
      FNL=(SUGMFN(MSL2)-SUGMFN(MSL1))/(MSL2**2-MSL1**2)
     $*AL*MLQ**2/COSB**2
      FCN=FNT+FNB+FNL/3.
C      DELHPS=3*G0(2)**2*MU*(COTB+TANB)/32./PI**2/MW**2*FCN
      DELHPS=6*MU*(COTB+TANB)/32./PI**2/VEV**2*FCN
      RDEL=SQRT(ABS(DELHPS))
C     Try Javier/Xerxes improved mA formula
      RZT=.5*(G0(24)-G0(23))+MZ**2*(8*COS2W-5.)*COS2B/12.
      RZB=.5*(G0(24)-G0(22))-MZ**2*(4*COS2W-1.)*COS2B/12.
      RZL=.5*(G0(21)-G0(20))-MZ**2*(4*COS2W-3.)*COS2B/4.
C      ZT=.25*(MST2**2-MST1**2)**2-(G0(6)*VEV)**2*(AT*SINB-MU*COSB)**2
C      ZB=.25*(MSB2**2-MSB1**2)**2-(G0(5)*VEV)**2*(AB*COSB-MU*SINB)**2
C      ZL=.25*(MSL2**2-MSL1**2)**2-(G0(4)*VEV)**2*(AL*COSB-MU*SINB)**2
      SIGT=3*G0(6)**2*SUGMFN(MTQ)/8./PI/PI
      SIGB=-3*G0(5)**2*SUGMFN(MBQ)/8./PI/PI
      SIGL=-1*G0(4)**2*SUGMFN(MLQ)/8./PI/PI
      SIGST=-3.*((SUGMFN(MST1)+SUGMFN(MST2))*(G0(6)**2-.25*
     $(G**2+GP**2))+(SUGMFN(MST2)-SUGMFN(MST1))/(MST2**2-MST1**2)*
     $(G0(6)**2*(AT**2-MU**2)-(G**2+GP**2)*(8*COS2W-5.)*
     $RZT/6.))/16./PI/PI
      SIGSB=3.*((SUGMFN(MSB1)+SUGMFN(MSB2))*(G0(5)**2-.25*
     $(G**2+GP**2))+(SUGMFN(MSB2)-SUGMFN(MSB1))/(MSB2**2-MSB1**2)*
     $(G0(5)**2*(AB**2-MU**2)-(G**2+GP**2)*(4*COS2W-1.)*
     $RZB/6.))/16./PI/PI
      SIGSL=1.*((SUGMFN(MSL1)+SUGMFN(MSL2))*(G0(4)**2-.25*
     $(G**2+GP**2))+(SUGMFN(MSL2)-SUGMFN(MSL1))/(MSL2**2-MSL1**2)*
     $(G0(4)**2*(AL**2-MU**2)-(G**2+GP**2)*(2*COS2W-1.5)*
     $RZL/6.))/16./PI/PI
C          Tree level mhp not needed at this point so fix if negative
      IF (INUHM.NE.1) THEN
      IF (ILOOP.EQ.0) THEN
        MHPS=M1S+M2S
        IF (MHPS.LT.0.) MHPS=0.
      ELSE
C        MHPS=B*MU*(COTB+TANB)+DELHPS
C        Use improved Javier formula
        MHPS=(TANB**2+1.)/(TANB**2-1.)*(G0(13)-G0(14)+SIGT+SIGST+
     $SIGSB+SIGB+SIGSL+SIGL)-AMZ**2
C     If MHPS neg. on this round, set to MZ to AMHL can be
C     calculated, then check MHPS on next round...
        IF (MHPS.LT.0.) THEN
c          MHPNEG=1
          MHPS=AMZ**2
        END IF
      END IF
      MHP=SQRT(MHPS)
      ELSE
      MHP=AMHA
      END IF
      MSS(31)=MHP
C
C          Initialize SUSY parameters in /SSPAR/:
C
      AMGLSS=ABS(G0(9))
      AMULSS=MSS(2)
      AMURSS=MSS(3)
      AMDLSS=MSS(4)
      AMDRSS=MSS(5)
      AMSLSS=MSS(6)
      AMSRSS=MSS(7)
      AMCLSS=MSS(8)
      AMCRSS=MSS(9)
      AMN1SS=MSS(14)
      AMN2SS=MSS(15)
      AMN3SS=MSS(16)
      AMELSS=MSS(17)
      AMERSS=MSS(18)
      AMMLSS=MSS(19)
      AMMRSS=MSS(20)
      TWOM1=-MU
      RV2V1=1./XTANB
      AMTLSS=SIGN(1.,G0(24))*SQRT(ABS(G0(24)))
      AMTRSS=SIGN(1.,G0(23))*SQRT(ABS(G0(23)))
      AMBLSS=SQRT(G0(24))
      AMBRSS=SQRT(G0(22))
      AMLLSS=SQRT(G0(21))
      AMLRSS=SQRT(G0(20))
      AMB1SS=MSS(10)
      AMB2SS=MSS(11)
      AMT1SS=MSS(12)
      AMT2SS=MSS(13)
      AML1SS=MSS(21)
      AML2SS=MSS(22)
      AAT=G0(12)
      AAB=G0(11)
      AAL=G0(10)
      AMHA=MHP
C
C          Use SSMASS to diagonalize neutralino and chargino mass
C          matrices and calculate Higgs masses.
C
      MHLNEG=0
      MHCNEG=0
      CALL SSMASS(G0(9),G0(7),G0(8),IALLOW,ILOOP,MHLNEG,MHCNEG,
     $MSQNEG,IMODEL)
c      IF(MHLNEG.EQ.1.OR.MHCNEG.EQ.1) THEN
c        NOGOOD=8
c      ENDIF
C      IF(IALLOW.NE.0.AND.ILOOP.NE.0) THEN
C        NOGOOD=5
C        GO TO 100
C      ENDIF
C
C          Save results also in MSS; re-save radiative corrected masses
C
      MSS(2)=AMULSS
      MSS(3)=AMURSS
      MSS(4)=AMDLSS
      MSS(5)=AMDRSS
      MSS(6)=AMSLSS
      MSS(7)=AMSRSS
      MSS(8)=AMCLSS
      MSS(9)=AMCRSS
      MSS(10)=AMB1SS
      MSS(11)=AMB2SS
      MSS(12)=AMT1SS
      MSS(13)=AMT2SS
      MSS(14)=AMN1SS
      MSS(15)=AMN2SS
      MSS(16)=AMN3SS
      MSS(17)=AMELSS
      MSS(18)=AMERSS
      MSS(19)=AMMLSS
      MSS(20)=AMMRSS
      MSS(21)=AML1SS
      MSS(22)=AML2SS
      MSS(23)=AMZ1SS
      MSS(24)=AMZ2SS
      MSS(25)=AMZ3SS
      MSS(26)=AMZ4SS
      MSS(27)=AMW1SS
      MSS(28)=AMW2SS
      MSS(29)=AMHL
      MSS(30)=AMHH
      MSS(31)=AMHA
      MSS(32)=AMHC
C     Azar's SGNM3 fix
C          Keep track of sign of M3; user input of Mgl>0 means M3<0
      SGNM3=-SIGN(1.,G0(9))
C
C     CALCULATE CHARGINO AND HIGGS ONE LOOP 
C     CORRECTIONS TO mA
      IF (ILOOP.EQ.0) THEN
        MHPS=M1S+M2S
        IF (MHPS.LT.0.) THEN
        MHPS=0.
        ENDIF
      ELSE
      SIGHIG= -(G**2+GP**2)/8.*AMHA**2*COS2B/2./PI**2*
     #   (SUGMFN(AMHH)-SUGMFN(AMHL))
     #   /(AMHH**2-AMHL**2)
      SIGCHA= -2*XW*(G**2+GP**2)/8.*MW**2*COS2B/PI**2*
     #   (SUGMFN(AMW1SS)-SUGMFN(AMW2SS))
     #   /(AMW1SS**2-AMW2SS**2)
      MHPS=(TANB**2+1.)/(TANB**2-1.)*(G0(13)-
     # G0(14)+SIGT+SIGST+
     # SIGSB+SIGB+SIGSL+SIGL+SIGHIG+SIGCHA)
     #     -AMZ**2
      SIGA=SIGT+SIGST+SIGSB+SIGB+SIGSL+SIGL+SIGHIG+SIGCHA
      IF (INUHM.NE.1) THEN
      IF (MHPS.GT.0) THEN
          MHP=SQRT(MHPS)
      ELSE
        MHPNEG=1
        MHPS=1.
      END IF
      MHP=SQRT(MHPS)
      ELSE
      MHP=AMHA
      END IF
      MSS(31)=MHP
      AMHA=MHP
      END IF
C          Gluino pole mass
      MGLMGL=G0(9)
      ASMGL=SUALFS(MGLMGL**2,.36,MT,3)
      XLAM=DLOG(DBLE(MGLMGL**2))
      MSS(1)=SSPOLE(MGLMGL,MGLMGL**2,ASMGL)
      AMGLSS=ABS(MSS(1))
      GSS(9)=G0(9)
C
100   RETURN
      END
