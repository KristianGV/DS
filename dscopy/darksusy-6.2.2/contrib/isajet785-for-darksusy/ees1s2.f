C------------------------------------------------------------
      FUNCTION EES1S2(Z)
      IMPLICIT NONE
      COMMON/BSQRK/ RS,M1,M2,AE,BE,PROPZ,E4,MZ
      REAL RS,M1,M2,AE,BE,PROPZ,E4,MZ 
      COMMON/S1S2/ BF,THETA
      REAL BF,THETA,SSXLAM
      COMMON/EPOL/ FLEM,FREP,FREM,FLEP
      COMMON/STUFF/ UNITS,PI
      REAL FLEM,FREP,FREM,FLEP,UNITS,PI
      COMMON/SFSF/QF,ALR,NC
      REAL Z
      REAL QF,NC,ESQ,AL,ALR
      REAL S,P,SIGLR,SIGRL,EES1S2
      S=RS**2
      P=SQRT(SSXLAM(S,M1**2,M2**2))/2./RS
      ESQ=SQRT(E4)
      AL=ESQ/4./PI
      SIGLR=NC*PI*AL**2*(AE-BE)**2*BF**2*COS(THETA)**2*SIN(THETA)**2
     ,*P**3*(1.-Z)**2/RS/PROPZ
      SIGRL=NC*PI*AL**2*(AE+BE)**2*BF**2*COS(THETA)**2*SIN(THETA)**2
     ,*P**3*(1.-Z)**2/RS/PROPZ
      EES1S2=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS
      RETURN
      END
