C------------------------------------------------------------
      FUNCTION EEERER(Z)
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
      REAL S,P,E,PHIERR,PHIERL,SIGERR,SIGERL,EEERER
      REAL BLZS,BLZJS,MZI(4)
      INTEGER I,II
      S=RS**2
      P=SQRT(SSXLAM(S,M1**2,M2**2))/2./RS
      E=RS/2.
      ESQ=SQRT(E4)
      PHIERR=8.*QF**2/S+(2*ALR**2*(AE-BE)**2*S-8.*(AE-BE)*QF*ALR*
     ,(S-MZ**2))/PROPZ
      PHIERL=8.*QF**2/S+(2*ALR**2*(AE+BE)**2*S-8.*(AE+BE)*QF*ALR*
     ,(S-MZ**2))/PROPZ
      PHIERR=E4*(1.-Z**2)*PHIERR
      PHIERL=E4*(1.-Z**2)*PHIERL
      DO 10 I=1,4
        BLZS=BLZ(I)*CONJG(BLZ(I))
        PHIERR=PHIERR+2.*BLZS**2*S*(1.-Z**2)/
     ,(2.*E*(E-P*Z)-M1**2+ABS(MZI(I)**2))**2-
     ,8.*ESQ*(1.-Z**2)*BLZS/(2.*E*(E-P*Z)-M1**2+ABS(MZI(I)**2))*
     ,(1.+(AE-BE)**2*S*(S-MZ**2)/PROPZ)
        IF (I.LE.3) THEN
          DO 11 II=I+1,4
            BLZJS=BLZ(II)*CONJG(BLZ(II))
            PHIERR=PHIERR+4.*BLZS*BLZJS*S*(1.-Z**2)/
     ,(2.*E*(E-P*Z)-M1**2+ABS(MZI(I)**2))/
     ,(2.*E*(E-P*Z)-M1**2+ABS(MZI(II)**2))
 11       CONTINUE
        END IF
 10   CONTINUE
      SIGERR=NC*P**3*PHIERR/256./PI/E**3
      SIGERL=NC*P**3*PHIERL/256./PI/E**3
      EEERER=(FLEM*FREP*SIGERL+FREM*FLEP*SIGERR)*UNITS
      RETURN
      END
