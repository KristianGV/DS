#include "PILOT.inc"
      REAL FUNCTION SSPOLE(MGMS,MUSQ,AS)
C*********************************************************************
C* Computes the on-shell (pole) gluino mass for given running (MSbar)*
C* gluino mass, defined at scale MUSQ, and given alpha_s (AS). The   *
C* squark masses are stored in the SQUARK COMMON block.              *
C* This function needs the complex functions B0 and B1.              *
C* Contributed by M. Drees; modified by H. Baer                      *
C  B0 contributions from Pierce et al.                               *
C  included on 9/23/05 by J.Ferrandis and H. Baer                    *
C                                                                    *
C  Version 7.30: Cast COMPLEX*16 to REAL*8 in standard way. :-(      *
C*********************************************************************
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sssm.inc"
#include "sspar.inc"
#include "sstmp.inc"
#include "ssinf.inc"
      REAL MGMS,MUSQ,AS,MGSQ,FAC
      DOUBLE PRECISION DMUSQ,DFAC
      COMPLEX*16 SSB1,SSB0
      DMUSQ=MUSQ
      XLAM = LOG(DMUSQ)
      MGSQ = MGMS*MGMS
C
C          Cast COMPLEX*16 to REAL*8:
C
      DFAC = -SSB1(MGSQ,0.,AMULSS) -SSB1(MGSQ,0.,AMURSS)
     $-SSB1(MGSQ,0.,AMDLSS) -SSB1(MGSQ,0.,AMDRSS)
     $-SSB1(MGSQ,0.,AMSLSS) -SSB1(MGSQ,0.,AMSRSS)
     $-SSB1(MGSQ,0.,AMCLSS) -SSB1(MGSQ,0.,AMCRSS)-
     $ ( SSB1(MGSQ,AMTP,AMT1SS)+ SSB1(MGSQ,AMTP,AMT2SS)+
     $  SSB1(MGSQ,4.0,AMB1SS) + SSB1(MGSQ,4.0,AMB2SS) )
     $ - AMTP*SIN(2.*THETAT)*(SSB0(MGSQ,AMTP,AMT1SS)-
     $   SSB0(MGSQ,AMTP,AMT2SS))/MGMS
     $ - 4.0*SIN(2.*THETAB)*(SSB0(MGSQ,4.0,AMB1SS)-
     $   SSB0(MGSQ,4.0,AMB2SS))/MGMS
      DFAC = DFAC + 15.D0 + 9.D0*LOG(DMUSQ/MGSQ)
      FAC=DFAC
      SSPOLE = MGMS*(1.0 + .0796*AS*FAC )
      RETURN
      END
