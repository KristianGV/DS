C------------------------------------------------------------
      FUNCTION EEELEL(Z)
      IMPLICIT NONE
      COMMON/BSQRK/ RS,M1,M2,AE,BE,PROPZ,E4,MZ
      COMMON/EPOL/ FLEM,FREP,FREM,FLEP
      COMMON/STUFF/ UNITS,PI
      COMMON/SFSF/QF,ALR,NC
      COMMON/SESE/ ALZ,BLZ,MZI
      COMPLEX ALZ(4),BLZ(4)
      REAL Z,SSXLAM,ESQ
      REAL RS,M1,M2,AE,BE,PROPZ,E4,MZ 
      REAL FLEM,FREP,FREM,FLEP,UNITS,PI
      REAL QF,ALR,NC
      REAL S,P,E,PHIELR,PHIELL,SIGELR,SIGELL,EEELEL
      REAL ALZS,ALZJS,MZI(4)
      INTEGER I,II
      S=RS**2
      P=SQRT(SSXLAM(S,M1**2,M2**2))/2./RS
      E=RS/2.
      ESQ=SQRT(E4)
      PHIELR=8.*QF**2/S+(2*ALR**2*(AE-BE)**2*S-8.*(AE-BE)*QF*ALR*
     ,(S-MZ**2))/PROPZ
      PHIELL=8.*QF**2/S+(2*ALR**2*(AE+BE)**2*S-8.*(AE+BE)*QF*ALR*
     ,(S-MZ**2))/PROPZ
      PHIELR=E4*(1.-Z**2)*PHIELR
      PHIELL=E4*(1.-Z**2)*PHIELL
      DO 10 I=1,4
        ALZS=ALZ(I)*CONJG(ALZ(I))
        PHIELL=PHIELL+2.*ALZS**2*S*(1.-Z**2)/
     ,(2.*E*(E-P*Z)-M1**2+ABS(MZI(I)**2))**2-
     ,8.*ESQ*(1.-Z**2)*ALZS/(2.*E*(E-P*Z)-M1**2+ABS(MZI(I)**2))*
     ,(1.+(AE-BE)**2*S*(S-MZ**2)/PROPZ)
        IF (I.LE.3) THEN
          DO 11 II=I+1,4
            ALZJS=ALZ(II)*CONJG(ALZ(II))
            PHIELL=PHIELL+4.*ALZS*ALZJS*S*(1.-Z**2)/
     ,(2.*E*(E-P*Z)-M1**2+ABS(MZI(I)**2))/
     ,(2.*E*(E-P*Z)-M1**2+ABS(MZI(II)**2))
 11       CONTINUE
        END IF
 10   CONTINUE
      SIGELR=NC*P**3*PHIELR/256./PI/E**3
      SIGELL=NC*P**3*PHIELL/256./PI/E**3
      EEELEL=(FLEM*FREP*SIGELL+FREM*FLEP*SIGELR)*UNITS
      RETURN
      END
