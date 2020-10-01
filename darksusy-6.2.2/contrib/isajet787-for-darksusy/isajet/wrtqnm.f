#include "PILOT.inc"
C-----------------------------------------------------------------------
      SUBROUTINE WRTQNM(LOUT)
C-----------------------------------------------------------------------
C    Write SLHA-like blocks QNUMBERS with sparticles quantum numbers
C    in LHE event file.
C
C    Ref: arXiv:0712.3311
C    Created: 07/11/13 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER LOUT
      Integer j,iPDG(31),iQ3(31),iNSS(31),iCREP(31),iAPD(31)
      CHARACTER*16 cSPN(31)
C   Particle names (in this order)
      DATA cSPN /
     $' H^0',' A^0',' H^+',
     $' dnl',' upl',' stl',' chl',' b1',' t1',
     $' el-',' nuel',' mul-',' numl',' tau1',' nutl',
     $' glss',' z1ss',' z2ss',' w1ss',' z3ss',' z4ss',' w2ss',
     $' dnr',' upr',' str',' chr',' b2',' t2',
     $' er-',' mur-',' tau2' /
C   PDG codes of the above particles
      DATA iPDG /
     &     35,     36,     37,
     &1000001,1000002,1000003,1000004,1000005,1000006,
     &1000011,1000012,1000013,1000014,1000015,1000016,
     &1000021,1000022,1000023,1000024,1000025,1000035,1000037,
     &2000001,2000002,2000003,2000004,2000005,2000006,
     &2000011,2000013,2000015 /
C   Electric charge (x3) of the above particles
      DATA iQ3 /
     &  0, 0, 3,
     & -1, 2,-1, 2,-1, 2,
     & -3, 0,-3, 0,-3, 0,
     &  0, 0, 0, 3, 0, 0, 3,
     & -1, 2,-1, 2,-1, 2,
     & -3,-3,-3 /
C   Number of spin states (2S+1) of the above particles
      DATA iNSS /
     & 1, 1, 1,
     & 1, 1, 1, 1, 1, 1,
     & 1, 1, 1, 1, 1, 1,
     & 2, 2, 2, 2, 2, 2, 2,
     & 1, 1, 1, 1, 1, 1,
     & 1, 1, 1 /
C   Color representation (1=singlet, 3=triplet, 8=octet) of the above particles
      DATA iCREP /
     & 1, 1, 1,
     & 3, 3, 3, 3, 3, 3,
     & 1, 1, 1, 1, 1, 1,
     & 8, 1, 1, 1, 1, 1, 1,
     & 3, 3, 3, 3, 3, 3,
     & 1, 1, 1 /
C   Particle/Antiparticle distinction (0=own antiparticle, 1=otherwise) of the a
      DATA iAPD /
     & 0, 0, 1,
     & 1, 1, 1, 1, 1, 1,
     & 1, 1, 1, 1, 1, 1,
     & 0, 0, 0, 1, 0, 0, 1,
     & 1, 1, 1, 1, 1, 1,
     & 1, 1, 1 /

      DO j=1,31
        WRITE(LOUT,5000) iPDG(j),cSPN(j)
	WRITE(LOUT,5001) 1,iQ3(j),'3 times electric charge'
	WRITE(LOUT,5001) 2,iNSS(j),'number of spin states (2s+1)'
	WRITE(LOUT,5001) 3,iCREP(j),
     &                  'colour rep (1: singlet, 3: triplet, 8: octet)'
	WRITE(LOUT,5001) 4,iAPD(j),
     &                  'particle/antiparticle distinction (0=own anti)'
      ENDDO

 5000 FORMAT('Block QNUMBERS',1x,I7,2x'#',1x,A)
 5001 FORMAT(1x,I2,1x,I2,2x'#',1x,A)

      RETURN
      END
