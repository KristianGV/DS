#include "PILOT.inc"
        REAL FUNCTION SSN1ST(SS)
C-----------------------------------------------------------------------
C          SSN1ST: l_1 -> stau_1+nu_l+nutaubar: TATA G FUNCTION
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
        REAL SS
        DOUBLE PRECISION S,M1,M2,MST1,MN1,WID
        S=SS
        M1=TMP(1)
        M2=TMP(2)
        MST1=TMP(3)
        MN1=TMP(4)
        WID=(S-MST1**2)**2/(S-M1**2)/(S-M2**2)*(S-MN1**2)**2
     $       /S**2/MN1**2
        SSN1ST=WID
        RETURN
        END
