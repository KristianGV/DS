#include "PILOT.inc"
      SUBROUTINE SIGTC2
C
C          Compute the techni-rho decay distribution cross section
C          D(SIGMA)/D(QMW**2)D(YW)D(OMEGA)
C          for the specified jet types. This is trivial but done for
C          compatibility with Drell-Yan and Higgs.
C
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "itapes.inc"
#include "const.inc"
#include "jetpar.inc"
#include "jetsig.inc"
#include "pjets.inc"
#include "wsig.inc"
#include "tcpar.inc"
C
      REAL AM12,AM22,ANGFAC,S,T,U
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT)
C
C          Angfac is (1-z**2), and is determined in terms of S,T,U.
C          Note that both rho+- and rho0 are always elastic.
      AM12=PJETS(5,1)**2
      AM22=PJETS(5,2)**2
      ANGFAC=4.*(T*U-AM12*AM22)/((S-AM12-AM22)**2-4.*AM12*AM22)
C          Differential cross section
      SIGLLQ=SIGEVT*ANGFAC*3./(8.*PI)
      RETURN
      END
