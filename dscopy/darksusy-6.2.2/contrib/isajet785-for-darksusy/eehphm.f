C----------------------------------------------------------
      FUNCTION EEHPHM(Z)
      IMPLICIT NONE
      COMMON/BSQRK/ RS,M1,M2,AE,BE,PROPZ,E4,MZ
      COMMON/EPOL/ FLEM,FREP,FREM,FLEP
      COMMON/STUFF/ UNITS,PI
      REAL Z,SSXLAM,EEHPHM
      REAL RS,M1,M2,AE,BE,PROPZ,E4,MZ
      REAL FLEM,FREP,FREM,FLEP,UNITS,PI
      REAL S,P,EZ,SIGLR,SIGRL
      REAL SN2THW,CS2THW,COSW,SINW
      DATA SN2THW/.232/
      S=RS**2
      P=SQRT(SSXLAM(S,M1**2,M2**2))/2./RS
      CS2THW=1.-SN2THW
      COSW=SQRT(CS2THW)
      SINW=SQRT(SN2THW)
      SIGRL=E4*P**3*(1.-Z**2)/4./PI/RS*
     ,(1./S/S+((2*SN2THW-1.)/2/COSW/SINW)**2*(AE+BE)**2/PROPZ+
     ,(2*SN2THW-1.)*(AE+BE)*(S-MZ**2)/S/COSW/SINW/PROPZ)
      SIGLR=E4*P**3*(1.-Z**2)/4./PI/RS*
     ,(1./S/S+((2*SN2THW-1.)/2/COSW/SINW)**2*(AE-BE)**2/PROPZ+
     ,(2*SN2THW-1.)*(AE-BE)*(S-MZ**2)/S/COSW/SINW/PROPZ)
      EEHPHM=(FLEM*FREP*SIGLR+FREM*FLEP*SIGRL)*UNITS
      RETURN
      END
