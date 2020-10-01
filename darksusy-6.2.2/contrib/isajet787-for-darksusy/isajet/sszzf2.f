#include "PILOT.inc"
      REAL FUNCTION SSZZF2(E)
C-----------------------------------------------------------------------
C          SSWZBF: ziss -> zjss f fbar
C          Baer's Z2ZFUN
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
C
      REAL E
      DOUBLE PRECISION MZ1,MZ2,SN,M,PI,F1,F2,ZWID,MF,BF
      DATA PI/3.14159265D0/
C
      MZ2=TMP(1)
      MZ1=TMP(2)
      SN=TMP(4)
      MF=TMP(6)
      M=AMZ
      ZWID=GAMZ
      BF=DSQRT(MAX(0.D0,1.D0-4*MF**2/(MZ2**2+MZ1**2-2*E*MZ2)))
C
      F1=SQRT(MAX(0.D0,E**2-MZ1**2))/
     $((MZ1**2+MZ2**2-M**2-2*MZ2*E)**2+ZWID**2*M**2)
      F2=E*(MZ1**2+MZ2**2+2*SN*MZ1*MZ2)-MZ2*(E**2+MZ1**2)-
     $BF*MZ2*(E**2-MZ1**2)/3.D0-SN*MZ1*(MZ1**2+MZ2**2-2*MF**2)
      SSZZF2=MZ2*BF*F1*F2/2.D0/PI**3
      RETURN
      END
