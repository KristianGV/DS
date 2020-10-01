C------------------------------------------------------------
      FUNCTION EEELER(Z)
      IMPLICIT NONE
      COMMON/BSQRK/ RS,M1,M2,AE,BE,PROPZ,E4,MZ
      COMMON/EPOL/ FLEM,FREP,FREM,FLEP
      COMMON/STUFF/ UNITS,PI
      COMMON/SESE/ ALZ,BLZ,MZI
      COMPLEX ALZ(4),BLZ(4)
      REAL Z,SSXLAM,ESQ
      REAL RS,M1,M2,AE,BE,PROPZ,E4,MZ 
      REAL FLEM,FREP,FREM,FLEP,UNITS,PI
      REAL QF,ALR,NC
      REAL S,P,E,SIGELR,EEELER
      REAL ALZS,BLZS,AZI,TM,AZII,BLZJS,MZI(4)
      INTEGER I,II
      S=RS**2
      P=SQRT(SSXLAM(S,M1**2,M2**2))/2./RS
      E=SQRT(P**2+M1**2)
      TM=0.
      DO 10 I=1,4
        ALZS=ALZ(I)*CONJG(ALZ(I))
        BLZS=BLZ(I)*CONJG(BLZ(I))
        AZI=(MZI(I)**2-M1**2)/2./(RS/2.)
        TM=TM+ALZS*BLZS*MZI(I)**2/(E-P*Z+AZI)**2
        IF (I.LE.3) THEN
          DO 11 II=I+1,4
            AZII=(MZI(II)**2-M1**2)/2./(RS/2.)
            TM=TM+2*ABS(MZI(I))*ABS(MZI(II))*
     ,REAL(ALZ(I)*CONJG(ALZ(II)*BLZ(I))*BLZ(II))/
     ,(E-P*Z+AZI)/(E-P*Z+AZII)
 11       CONTINUE
        END IF
 10   CONTINUE
      SIGELR=P*TM/32./PI/S/(RS/2.)
      EEELER=FLEM*FLEP*SIGELR*UNITS
      RETURN
      END
