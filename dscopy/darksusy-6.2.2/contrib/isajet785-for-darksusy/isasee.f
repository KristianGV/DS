CDECK  ID>, ISASEE.
      SUBROUTINE ISASEE(RS,PLEM,PLEP,LOUT)
C
C          Compute SIGMA(TOT) for
C          e- e+ ----> SUSY particles
c          Dump into LHA file for use by Fittino
C          See Baer et. al., PRD54 (1996) 6735 for differential sigma's
C
      IMPLICIT NONE
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C
      COMMON/BSQRK/ ROOTS,M1,M2,AE,BE,PROPZ,E4,MZ
      COMMON /EPOL/ FLEM,FLEP,FREM,FREP
      COMMON/STUFF/ UNITS,PI
      COMMON/SFSF/QF,ALR,NC
      COMMON/S1S2/ BF,THETA
      COMMON/SESE/ ALZ,BLZ,MZI
      COMMON/N1N1/ G,COSGR,SINGR,MWI
      COMMON/WIWJ/ X,Y,MSN,THY,GRTERM,COS2W
      COMMON/ZIZJ/ IZ,JZ,W,SGN,AMER,AMEL
      COMMON/EEZHGS/ APBFAC
      COMPLEX AUZ(4),BUZ(4),ADZ(4),BDZ(4),ALZ(4),BLZ(4),ANZ(4),BNZ(4)
      COMPLEX ZI,W(4,4)
      REAL SSXINT
      EXTERNAL EESFSF,EEELEL,EEERER,EES1S2,EEELER,EEEREL,EEN1N1,EEZIZJ
      EXTERNAL EEWIWI,EEWIWJ,EEZH,EEHA,EEHPHM
      REAL BF,THETA,APBFAC,BETA
      REAL PLEM,PLEP
      REAL RS,M1,M2,PROPZ,MZ,S
      REAL FLEM,FREP,FREM,FLEP,UNITS,PI,SR2,XW,ROOTS,SIG(48)
      REAL QF,ALR,NC,SGN,AMEL,AMER
      REAL XC,YC,XS,YS,X,Y,MSN,THX,THY,GRTERM
      REAL COSGL,SINGL,XM,YM
      REAL COS2W,ALEM,G,GP,TT,C,AU,BU,AD,BD,AE,BE,AN,BN,AEL,AER,ANL,
     $QU,QD,QE,QN,E4,ESQ,AL2,AL3,COSGR,SINGR
      REAL ZR(4,4),TH(4),MZI(4),MWI(2)
      CHARACTER*40 VERSN,VISAJE
      INTEGER LOUT,IDEP,IDEM,NDA,ID1,ID2,I,J,IZ,JZ,IDZI(4)
      DATA IDEP/-11/,IDEM/11/,UNITS/.389E12/,XW/.232/,NDA/2/
      DATA IDZI/1000022,1000023,1000025,1000035/
      DATA ZI/(0.,1.)/
      PI=4*ATAN(1.)
      SR2=SQRT(2.)
      MZ=AMZ
      ROOTS=RS
      S=RS**2
      BETA=ATAN(1./RV2V1)
      WRITE(LOUT,7000)
     . ' ISAJET e-e+ cross sections in SUSY Les Houches Accord 2 format'
      WRITE(LOUT,7000)
     .  ' Created by ISASEE. Last revision: H. Baer, 2014 May 31'
      VERSN=VISAJE()
      VERSN=VERSN(14:)
      WRITE(LOUT,7001)    'XSINFO', 
     ,                      'Program information'
      WRITE(LOUT,7012) 1, 'ISASUGRA/ISASUSY from ISAJET    ',
     ,                      'Spectrum Calculator'
      WRITE(LOUT,7012) 2,  VERSN, 
     ,                      'Version number'
      WRITE(LOUT,7500) IDEM,IDEP,RS,PLEM,PLEP,1,' e- e+ XS'
      WRITE(LOUT,'(A)') '#   Sigma [fb]     NDA   ID1   ID2'

C
      FLEP=(1.+PLEP)/2.
      FLEM=(1.+PLEM)/2.
      FREP=(1.-PLEP)/2.
      FREM=(1.-PLEM)/2.
C     constants
      PROPZ=(S-AMZ**2)**2+AMZ**2*GAMZ**2
      COS2W=1.-XW
      ALEM=1./128.
      G=SQRT(4*PI*ALEM/XW)
      GP=G*SQRT(XW/COS2W)
      TT=SQRT(XW/COS2W)
      C=1./TT
      AU=(C/4.-5*TT/12.)
      BU=-(C+TT)/4.
      AD=(TT/12.-C/4.)
      BD=(C+TT)/4.
      AE=(3*TT-C)/4.
      BE=(C+TT)/4.
      AN=(TT+C)/4.
      BN=-(C+TT)/4.
      AEL=2*(AE-BE)
      AER=2*(AE+BE)
      ANL=2*(AN-BN)
      QU=2./3.
      QD=-1./3.
      QE=-1.
      QN=0.
      E4=(4*PI*ALEM)**2
      ESQ=SQRT(E4)
      AL2=1./128./XW
      AL3=.118
C     neutralino couplings
      MZI(1)=AMZ1SS
      MZI(2)=AMZ2SS
      MZI(3)=AMZ3SS
      MZI(4)=AMZ4SS
      DO I=1,4
        IF (MZI(I).GT.0.) THEN
          TH(I)=0.
        ELSE
          TH(I)=1.
        END IF
        DO J=1,4
          ZR(I,J)=ZMIXSS(I,J)
        END DO
      END DO
      DO I=1,4
      AUZ(I)=ZI*((-ZI)**TH(I))*(G*ZR(3,I)+GP/3.*ZR(4,I))/SR2
      ADZ(I)=ZI*((-ZI)**TH(I))*(-G*ZR(3,I)+GP/3.*ZR(4,I))/SR2
      ALZ(I)=-ZI*((-ZI)**(TH(I)-1.))*(G*ZR(3,I)+GP*ZR(4,I))/SR2
      ANZ(I)=((-ZI)**(TH(I)-1.))*(G*ZR(3,I)+GP*ZR(4,I))/SR2
      BUZ(I)=-ZI*((ZI)**TH(I))*4*GP*ZR(4,I)/3./SR2
      BDZ(I)=ZI*((ZI)**TH(I))*2*GP*ZR(4,I)/3./SR2
      BLZ(I)=-((-ZI)**(TH(I)-1.))*SR2*GP*ZR(4,I)
      BNZ(I)=0.
        DO J=1,4
          W(I,J)=(ZI**TH(J))*((-ZI)**(TH(I)))*SQRT(G**2+GP**2)
     ,      /4.*(ZR(1,I)*ZR(1,J)-ZR(2,I)*ZR(2,J))
        END DO
      END DO
C
C----- chargino couplings
      MWI(1)=AMW1SS
      MWI(2)=AMW2SS
      COSGR=COS(GAMMAR)
      SINGR=SIN(GAMMAR)
      COSGL=COS(GAMMAL)
      SINGL=SIN(GAMMAL)
      MSN=AMN1SS
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      XC=1.-(COSGL**2+COSGR**2)/4./COS2W
      YC=(COSGR**2-COSGL**2)/4./COS2W
      XS=1.-(SINGL**2+SINGR**2)/4./COS2W
      YS=(SINGR**2-SINGL**2)/4./COS2W
C         Initialize since Fittino likes even zero cross sections
      DO I=1,48
        SIG(I)=0.
      END DO
C         Compute total cross sections
C         Squarks
      NC=3.
c-----Generation 1 
      IF (RS.GT.(2*AMULSS)) THEN
        M1=AMULSS
        M2=AMULSS
        QF=QU
        ALR=2*(AU-BU)
        SIG(1)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000002
      ID2=-ID1
      WRITE(LOUT,7600) SIG(1),NDA,ID1,ID2,' e- e+ -> uL + uLbar'
      IF (RS.GT.(2*AMURSS)) THEN
        M1=AMURSS
        M2=AMURSS
        QF=QU
        ALR=2*(AU+BU)
        SIG(2)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000002
      ID2=-ID1
      WRITE(LOUT,7600) SIG(2),NDA,ID1,ID2,' e- e+ -> uR + uRbar'
      IF (RS.GT.(2*AMDLSS)) THEN
        M1=AMDLSS
        M2=AMDLSS
        QF=QD
        ALR=2*(AD-BD)
        SIG(3)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000001
      ID2=-ID1
      WRITE(LOUT,7600) SIG(3),NDA,ID1,ID2,' e- e+ -> dL + dLbar'
      IF (RS.GT.(2*AMDRSS)) THEN
        M1=AMDRSS
        M2=AMDRSS
        QF=QD
        ALR=2*(AD+BD)
        SIG(4)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000001
      ID2=-ID1
      WRITE(LOUT,7600) SIG(4),NDA,ID1,ID2,' e- e+ -> dR + dRbar'
c-----Generation 2
      IF (RS.GT.(2*AMCLSS)) THEN
        M1=AMCLSS
        M2=AMCLSS
        QF=QU
        ALR=2*(AU-BU)
        SIG(5)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000004
      ID2=-ID1
      WRITE(LOUT,7600) SIG(5),NDA,ID1,ID2,' e- e+ -> cL + cLbar'
      IF (RS.GT.(2*AMCRSS)) THEN
        M1=AMCRSS
        M2=AMCRSS
        QF=QU
        ALR=2*(AU+BU)
        SIG(6)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000004
      ID2=-ID1
      WRITE(LOUT,7600) SIG(6),NDA,ID1,ID2,' e- e+ -> cR + cRbar'
      IF (RS.GT.(2*AMSLSS)) THEN
        M1=AMSLSS
        M2=AMSLSS
        QF=QD
        ALR=2*(AD-BD)
        SIG(7)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000003
      ID2=-ID1
      WRITE(LOUT,7600) SIG(7),NDA,ID1,ID2,' e- e+ -> sL + sLbar'
      IF (RS.GT.(2*AMSRSS)) THEN
        M1=AMSRSS
        M2=AMSRSS
        QF=QD
        ALR=2*(AD+BD)
        SIG(8)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000003
      ID2=-ID1
      WRITE(LOUT,7600) SIG(8),NDA,ID1,ID2,' e- e+ -> sR + sRbar'
c-----Generation 3
      IF (RS.GT.(2*AMT1SS)) THEN
        M1=AMT1SS
        M2=AMT1SS
        QF=QU
        ALR=2*(AU-BU)*COS(THETAT)**2+2*(AU+BU)*SIN(THETAT)**2
        SIG(9)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000006
      ID2=-ID1
      WRITE(LOUT,7600) SIG(9),NDA,ID1,ID2,' e- e+ -> t1 + t1bar'
      IF (RS.GT.(2*AMT2SS)) THEN
        M1=AMT2SS
        M2=AMT2SS
        QF=QU
        ALR=2*(AU-BU)*SIN(THETAT)**2+2*(AU+BU)*COS(THETAT)**2
        SIG(10)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000006
      ID2=-ID1
      WRITE(LOUT,7600) SIG(10),NDA,ID1,ID2,' e- e+ -> t2 + t2bar'
      IF (RS.GT.(AMT1SS+AMT2SS)) THEN
        M1=AMT1SS
        M2=AMT2SS
        THETA=THETAT
        BF=BU
        SIG(11)=SSXINT(-1.,EES1S2,1.)
        SIG(12)=SIG(11)
      END IF
      ID1=1000006
      ID2=-2000006
      WRITE(LOUT,7600) SIG(11),NDA,ID1,ID2,' e- e+ -> t1 + t2bar'
      ID1=2000006
      ID2=-1000006
      WRITE(LOUT,7600) SIG(12),NDA,ID1,ID2,' e- e+ -> t2 + t1bar'
      IF (RS.GT.(2*AMB1SS)) THEN
        M1=AMB1SS
        M2=AMB1SS
        QF=QD
        ALR=2*(AD-BD)*COS(THETAB)**2+2*(AD+BD)*SIN(THETAB)**2
        SIG(13)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000005
      ID2=-ID1
      WRITE(LOUT,7600) SIG(13),NDA,ID1,ID2,' e- e+ -> b1 + b1bar'
      IF (RS.GT.(2*AMB2SS)) THEN
        M1=AMB2SS
        M2=AMB2SS
        QF=QD
        ALR=2*(AD-BD)*SIN(THETAB)**2+2*(AD+BD)*COS(THETAB)**2
        SIG(14)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000005
      ID2=-ID1
      WRITE(LOUT,7600) SIG(14),NDA,ID1,ID2,' e- e+ -> b2 + b2bar'
      IF (RS.GT.(AMB1SS+AMB2SS)) THEN
        M1=AMB1SS
        M2=AMB2SS
        THETA=THETAB
        BF=BD
        SIG(15)=SSXINT(-1.,EES1S2,1.)
        SIG(16)=SIG(15)
      END IF
      ID1=1000005
      ID2=-2000005
      WRITE(LOUT,7600) SIG(15),NDA,ID1,ID2,' e- e+ -> b1 + b2bar'
      ID1=2000005
      ID2=-1000005
      WRITE(LOUT,7600) SIG(16),NDA,ID1,ID2,' e- e+ -> b2 + b1bar'
c-----Sleptons
      NC=1.
c-----Generation 1
      IF (RS.GT.(2*AMELSS)) THEN
        M1=AMELSS
        M2=AMELSS
        QF=QE
        ALR=2*(AE-BE)
        SIG(17)=SSXINT(-1.,EEELEL,1.)
      END IF
      ID1=1000011
      ID2=-ID1
      WRITE(LOUT,7600) SIG(17),NDA,ID1,ID2,' e- e+ -> eL + eLbar'
      IF (RS.GT.(2*AMERSS)) THEN
        M1=AMERSS
        M2=AMERSS
        QF=QE
        ALR=2*(AE+BE)
        SIG(18)=SSXINT(-1.,EEERER,1.)
      END IF
      ID1=2000011
      ID2=-ID1
      WRITE(LOUT,7600) SIG(18),NDA,ID1,ID2,' e- e+ -> eR + eRbar'
      IF (RS.GT.(AMELSS+AMERSS)) THEN
        M1=AMELSS
        M2=AMERSS
        SIG(19)=SSXINT(-1.,EEELER,1.)
      END IF
      ID1=1000011
      ID2=-2000011
      WRITE(LOUT,7600) SIG(19),NDA,ID1,ID2,' e- e+ -> eL + eRbar'
      IF (RS.GT.(AMERSS+AMELSS)) THEN
        M1=AMERSS
        M2=AMELSS
        SIG(20)=SSXINT(-1.,EEEREL,1.)
      END IF
      ID1=2000011
      ID2=-1000011
      WRITE(LOUT,7600) SIG(20),NDA,ID1,ID2,' e- e+ -> eR + eLbar'
      IF (RS.GT.(2*AMN1SS)) THEN
        M1=AMN1SS
        M2=AMN1SS
        QF=0.
        ALR=2*(AN-BN)
        SIG(21)=SSXINT(-1.,EEN1N1,1.)
      END IF
      ID1=1000012
      ID2=-ID1
      WRITE(LOUT,7600) SIG(21),NDA,ID1,ID2,' e- e+ -> neL + neLbar'
C---------Generation 2
      IF (RS.GT.(2*AMMLSS)) THEN
        M1=AMMLSS
        M2=AMMLSS
        QF=QE
        ALR=2*(AE-BE)
        SIG(22)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000013
      ID2=-ID1
      WRITE(LOUT,7600) SIG(22),NDA,ID1,ID2,' e- e+ -> muL + muLbar'
      IF (RS.GT.(2*AMMRSS)) THEN
        M1=AMMRSS
        M2=AMMRSS
        QF=QE
        ALR=2*(AE+BE)
        SIG(23)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000013
      ID2=-ID1
      WRITE(LOUT,7600) SIG(23),NDA,ID1,ID2,' e- e+ -> muR + muRbar'
      IF (RS.GT.(2*AMN2SS)) THEN
        M1=AMN2SS
        M2=AMN2SS
        QF=0.
        ALR=2*(AN-BN)
        SIG(24)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000014
      ID2=-ID1
      WRITE(LOUT,7600) SIG(24),NDA,ID1,ID2,' e- e+ -> nmuL + nmuLbar'
C---------Generation 3
      IF (RS.GT.(2*AML1SS)) THEN
        M1=AML1SS
        M2=AML1SS
        QF=QE
        ALR=2*(AE-BE)*COS(THETAL)**2+2*(AE+BE)*SIN(THETAL)**2
        SIG(25)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000015
      ID2=-ID1
      WRITE(LOUT,7600) SIG(25),NDA,ID1,ID2,' e- e+ -> tau1 + tau1bar'
      IF (RS.GT.(2*AML2SS)) THEN
        M1=AML2SS
        M2=AML2SS
        QF=QE
        ALR=2*(AE-BE)*SIN(THETAL)**2+2*(AE+BE)*COS(THETAL)**2
        SIG(26)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=2000015
      ID2=-ID1
      WRITE(LOUT,7600) SIG(26),NDA,ID1,ID2,' e- e+ -> tau2 + tau2bar'
      IF (RS.GT.(AML1SS+AML2SS)) THEN
        M1=AML1SS
        M2=AML2SS
        THETA=THETAL
        BF=BE
        SIG(27)=SSXINT(-1.,EES1S2,1.)
        SIG(28)=SIG(27)
      END IF
      ID1=1000015
      ID2=-2000015
      WRITE(LOUT,7600) SIG(27),NDA,ID1,ID2,' e- e+ -> tau1 + tau2bar'
      ID1=2000015
      ID2=-1000015
      WRITE(LOUT,7600) SIG(28),NDA,ID1,ID2,' e- e+ -> tau2 + tau1bar'
      IF (RS.GT.(2*AMN3SS)) THEN
        M1=AMN3SS
        M2=AMN3SS
        QF=0.
        ALR=2*(AN-BN)
        SIG(29)=SSXINT(-1.,EESFSF,1.)
      END IF
      ID1=1000016
      ID2=-ID1
      WRITE(LOUT,7600) SIG(29),NDA,ID1,ID2,' e- e+ -> ntauL + ntauLbar'
      I=29
      DO IZ=1,4
        DO JZ=IZ,4
          I=I+1
          IF (RS.GT.(ABS(MZI(IZ))+ABS(MZI(JZ)))) THEN
          M1=ABS(MZI(IZ))
          M2=ABS(MZI(JZ))
          AMER=AMERSS
          AMEL=AMELSS
          SGN=SIGN(1.,MZI(IZ))*SIGN(1.,MZI(JZ))
          SIG(I)=SSXINT(-1.,EEZIZJ,1.)
          END IF
          ID1=IDZI(IZ)
          ID2=IDZI(JZ)
          WRITE(LOUT,7600) SIG(I),NDA,ID1,ID2,' e- e+ -> Zi + Zj'
        END DO
      END DO
      IF (RS.GT.(2*ABS(AMW1SS))) THEN
        M1=AMW1SS
        M2=AMW1SS
        X=XC
        Y=YC
        GRTERM=SINGR**2
        SIG(40)=SSXINT(-1.,EEWIWI,1.)
      END IF
      ID1=1000024
      ID2=-ID1
      WRITE(LOUT,7600) SIG(40),NDA,ID1,ID2,' e- e+ -> W1 + W1bar'
      IF (RS.GT.(2*ABS(AMW2SS))) THEN
        M1=AMW2SS
        M2=AMW2SS
        X=XS
        Y=YS
        GRTERM=COSGR**2
        SIG(41)=SSXINT(-1.,EEWIWI,1.)
      END IF
      ID1=1000037
      ID2=-ID1
      WRITE(LOUT,7600) SIG(41),NDA,ID1,ID2,' e- e+ -> W2 + W2bar'
      IF (RS.GT.(ABS(AMW1SS)+ABS(AMW2SS))) THEN
        M1=AMW1SS
        M2=AMW2SS
        GRTERM=SINGR*COSGR
        X=(THX*SINGL*COSGL-THY*SINGR*COSGR)/2.
        Y=(THX*SINGL*COSGL+THY*SINGR*COSGR)/2.
        SIG(42)=SSXINT(-1.,EEWIWJ,1.)
        SIG(43)=SIG(42)
      END IF
      ID1=1000024
      ID2=-1000037
      WRITE(LOUT,7600) SIG(42),NDA,ID1,ID2,' e- e+ -> W1 + W2bar'
      ID1=1000037
      ID2=-1000024
      WRITE(LOUT,7600) SIG(43),NDA,ID1,ID2,' e- e+ -> W2 + W1bar'
      IF (RS.GT.(AMZ+AMHL)) THEN
        M1=AMZ
        M2=AMHL
        APBFAC=SIN(ALFAH+BETA)**2
        SIG(44)=SSXINT(-1.,EEZH,1.)
      END IF
      ID1=23
      ID2=25
      WRITE(LOUT,7600) SIG(44),NDA,ID1,ID2,' e- e+ -> Z + h'
      IF (RS.GT.(AMZ+AMHH)) THEN
        M1=AMZ
        M2=AMHH
        APBFAC=COS(ALFAH+BETA)**2
        SIG(45)=SSXINT(-1.,EEZH,1.)
      END IF
      ID1=23
      ID2=35
      WRITE(LOUT,7600) SIG(45),NDA,ID1,ID2,' e- e+ -> Z + H'
      IF (RS.GT.(AMHL+AMHA)) THEN
        M1=AMHL
        M2=AMHA
        APBFAC=COS(ALFAH+BETA)**2
        SIG(46)=SSXINT(-1.,EEHA,1.)
      END IF
      ID1=25
      ID2=36
      WRITE(LOUT,7600) SIG(46),NDA,ID1,ID2,' e- e+ -> h + A'
      IF (RS.GT.(AMHH+AMHA)) THEN
        M1=AMHH
        M2=AMHA
        APBFAC=SIN(ALFAH+BETA)**2
        SIG(47)=SSXINT(-1.,EEHA,1.)
      END IF
      ID1=35
      ID2=36
      WRITE(LOUT,7600) SIG(47),NDA,ID1,ID2,' e- e+ -> H + A'
      IF (RS.GT.(2*AMHC)) THEN
        M1=AMHC
        M2=AMHC
        SIG(48)=SSXINT(-1.,EEHPHM,1.)
      END IF
      ID1=37
      ID2=-ID1
      WRITE(LOUT,7600) SIG(48),NDA,ID1,ID2,' e- e+ -> H+ + H-'
 7000 FORMAT('# ',A)
C     Format to use for block statements
 7001 FORMAT('Block',1x,A,27x,'#',1x,A)
 7500 FORMAT('XS',1x,I3,1x,I3,3x,F6.1,1x,F5.2,1x,F5.2,2X,I2,3x,'#',1x,A)
 7600 FORMAT(E10.3,3x,I2,3x,I8,2x,I8,3x,'#',1x,A)
C     Indexed Char(12)
 7012 FORMAT(1x,I5,3x,A28,3x,'#',1x,A)
      RETURN
      END
