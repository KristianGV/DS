#include "PILOT.inc"
#ifdef ISATOOLS_X
      SUBROUTINE ISARES(IPRT)

      IMPLICIT NONE
      REAL*8 AQI(2,6),BQI(2,6),MI(2,2,6),NL(4),NLP(2),T3Q(6),EQ(6),MQ(6)
      REAL*8 FQQ(6),FQH(6),msq(2,6),cz(2),cq(2,6),rq(2,6),rq1(2,6)
      REAL*8 FNTQ(2,3),GQ(6),CSQ(2,2,6),THETAQ(6),SJ(2),ATQ(6)
      REAL*8 MH(2),DQ(6),delq(2,3),PDF(2,-5:5),mnucl(2),mrnucl(2)
      REAL*8 FNUCL(6,2),FNTG(2)

C! Two variables are added  
      REAL*8 THETA, NIXB     
      
      REAL*8 SIN2W,PI,ALPE,SW,CW,MW,MZ,E,G,G1,MNE,SA,CA,BETA,SB,
     $       CB,CAPB,SAPB,MU,MSUSY,QR,ALPS,GEVTOPBN,TQT,BD,B1D,
     $       BS,B1S,B2S,BGLU,ap,an,XMIN,XMAX,XD,Ctq5Pdf,X,
     $       SIGMA0PROT,SIGMA0NEUT,SSALFS,CMD,XJ,
     $       SIGMASPROT,SIGMASNEUT,thetaw,GPX

      INTEGER I,II,IQ,IX,JH,IN,J,kk
      
      INTEGER IPRT
      
      

      common/check/p1,p2
      real*8 p1,p2

#include "sugmg.inc"
#include "sspar.inc"
#include "sssm.inc"

c-----ISARES OUTPUT---------------------------
      
      COMMON/SUGRES/SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
       
      save/SUGRES/
c----------------------------------------------
ccccccccccccccccccc CONSTANTS ccccccccccccccccccccccccccccc
   
      PI=DACOS(-1d0)
       
      GPX=SQRT(0.6)*GSS(1)
      SIN2W=GPX**2/(GSS(2)**2+GPX**2)
      thetaw=asin(sqrt(SIN2W))
      ALPE=1d0/(4d0*PI/GSS(2)**2/SIN2W)
      
      
      SW=SQRT(SIN2W)
      CW=sqrt(1d0-SW**2)
      MW=80.419d0
      MZ=MW/CW
  
c      print *, SIN2W, ALPE,MW,MZ
       
      E=SQRT(4*PI*ALPE)
      G =E/SW
      G1=E/CW
      MNE=ABS(MSS(23))
      SA=SIN(-ALFAH)
      CA=COS(-ALFAH)
      MH(2)=MSS(29)
      MH(1)=MSS(30)
      BETA=ATAN(1d0/RV2V1)
      SB=SIN(BETA)
      CB=COS(BETA)
      CAPB=COS(-alfah+beta)
      SAPB=SIN(-alfah+beta)
      MU=GSS(25)
 
      mnucl(1)=0.938d0
C!      mnucl(2)=0.930d0
      mnucl(2)=0.939d0
      mrnucl(1)=mnucl(1)*mne/(mne+mnucl(1))
      mrnucl(2)=mnucl(2)*mne/(mne+mnucl(2))
      
      
      MSUSY=(MSS(2)+MSS(3)+MSS(4)+MSS(5)
     $       +MSS(10)+MSS(11)+MSS(12)+MSS(13))/8d0

      QR=sqrt(MSUSY**2-MSS(23)**2)
    
      ALPS=SSALFS(QR**2)
      
      GEVTOPBN= 0.3893796623*1.E+09      

c-------------------------------------------------------------c
C     update form factors to values from Hisano PRD87 (2013)035020
      FNTQ(1,1)=0.019
      FNTQ(1,2)=0.027
      FNTQ(1,3)=0.009
      
      FNTQ(2,1)=0.013
      FNTQ(2,2)=0.040
      FNTQ(2,3)=0.009
      
      FNTG(1)=1d0-FNTQ(1,1)-FNTQ(1,2)-FNTQ(1,3)
      FNTG(2)=1d0-FNTQ(2,1)-FNTQ(2,2)-FNTQ(2,3)
      
      delq(1,1)= 0.78
      delq(1,2)=-0.48
      delq(1,3)=-0.15
      
      delq(2,1)=-0.48
      delq(2,2)= 0.78
      delq(2,3)=-0.15
 
cccccccccccccccccc QUARK MASS           ccccccccccccccccccccccc
      mq(1)=AMUP
      mq(2)=AMDN
      mq(3)=AMST
      mq(4)=AMCH
      mq(5)=MBQ
      mq(6)=MTQ
c      print *,'QUARK MASSES:',(mq(kk),kk=1,6)
c-------------------------------------------------------------c

cccccccccccccccccc SQUARK MASS           ccccccccccccccccccccccc
      DO I=1,6
        msq(1,I)=mss(2*I  )
        msq(2,I)=mss(2*I+1)
      ENDDO
c-------------------------------------------------------------c

     
cccccccccccccccccc SQUARK MIXING ANGLE ccccccccccccccccccccccc
C! This part is rewritten to introduce mixing for first two generations.

      DO I=1,4       ! u d s c
        ATQ(I)=0
        IF((2*(I/2).eq.I.or.I.eq.1).and.I.ne.2) THEN
          THETAQ(I)=
     $ ATAN(-2*Mq(I)*(ATQ(I)-mu*CB/SB)/(MSQ(2,I)**2-MSQ(1,I)**2))/2  
        ELSE
          THETAQ(I)=
     $ ATAN(-2*Mq(I)*(ATQ(I)-mu*SB/CB)/(MSQ(2,I)**2-MSQ(1,I)**2))/2
        ENDIF
      ENDDO

      ATQ(5)=AAB
      ATQ(6)=AAT
      THETAQ(5)=THETAB
      THETAQ(6)=THETAT
      DO I=1,6
         THETA=THETAQ(I)
         MI(1,1,I)=  COS(THETA)
         MI(1,2,I)=  SIN(THETA)
         MI(2,1,I)= -SIN(THETA)
         MI(2,2,I)=  COS(THETA)
      ENDDO
c-------------------------------------------------------------c
cccccccccccccccc NEUTRALINO MIXING MATRIX, CZ ccccccccccccccccccccc
     
      NL(1)=  ZMIXSS(4,1)
      NL(2)=  ZMIXSS(3,1)
      NL(3)= -ZMIXSS(2,1)
      NL(4)= -ZMIXSS(1,1)

C!  Indeed NLP is not needed.           
      NLP(1) = NL(1)*CW+NL(2)*SW
      NLP(2) =-NL(1)*SW+NL(2)*CW
      CZ(2)= 1./2.*(G*NL(2)-G1*NL(1))*(NL(3)*SA+NL(4)*CA)
      CZ(1)= 1./2.*(G*NL(2)-G1*NL(1))*(NL(4)*SA-NL(3)*CA)
c      print *,'line 130:',CZ(1),CZ(2)

c------------------------------------------------------------c

ccccccccccccccccc T3Q, CQ      ccccccccccccccccccccccccccccccccccc
c        do I=1,6
c         thetaq(i)=-thetaq(i)
c        enddo

       
      DO I=1,6
        IF((2*(I/2).eq.I.or.I.eq.1).and.I.ne.2) THEN 
	  T3Q(I)= 1d0/2d0
	  EQ(I) = 2d0/3d0
	  RQ(1,I)= -sa/sb
	  RQ(2,I)= -ca/sb
	  RQ1(1,I)=-CA/SB
	  RQ1(2,I)= SA/SB
	ELSE
	  T3Q(I)= -1d0/2d0
	  EQ(I) = -1d0/3d0
	  RQ(1,I)= -ca/cb
	  RQ(2,I)=  sa/cb
	  RQ1(1,I)=-SA/CB
	  RQ1(2,I)=-CA/CB
	ENDIF

	CQ(1,I)=g/2d0/MW*RQ(1,I)
	CQ(2,I)=g/2d0/MW*RQ(2,I)
c        print *,'line 156:',CQ(1,I),CQ(2,I)
	
	SJ(1)=-CAPB
	SJ(2)= SAPB
	
	DO J=1,2
	  CSQ(J,1,I)=
     $    G*MZ/CW*SJ(J)*(T3Q(I)*cos(thetaq(I))**2
     $    -eq(i)*sw**2*cos(2.*thetaq(I)))
     $    +G*mq(I)**2/MW*RQ(J,I)
     $    -G*MQ(I)*SIN(2.*THETAQ(I))/2./MW*(mu*RQ1(J,I)-ATQ(I)*RQ(J,I))


	  CSQ(J,2,I)=
     $    G*MZ/CW*SJ(J)*(T3Q(I)*sin(thetaq(I))**2
     $    +eq(i)*sw**2*cos(2.*thetaq(I)))
     $    +G*mq(I)**2/MW*RQ(J,I)
     $    +G*MQ(I)*SIN(2.*THETAQ(I))/2./MW*(mu*RQ1(J,I)-ATQ(I)*RQ(J,I))


        ENDDO
	
	
      ENDDO
c------------------------------------------------------------c

ccccccccccccccccc AQI , BQI , FQQ, FQH , GQ  ccccccccccccccccccccccccccccccccccc
        DO IQ=1,6
        FQQ(IQ)=0d0
	FQH(IQ)=0d0
	GQ(IQ)=0d0
	DQ(IQ) =0d0
 
        DO II=1,2
 
C!        IX=3
C!        IF((2*(IQ/2).eq.IQ.or.IQ.eq.1).and.IQ.ne.2) IX=4
        NIXB=NL(3)/CB
        IF((2*(IQ/2).eq.IQ.or.IQ.eq.1).and.IQ.ne.2) NIXB=NL(4)/SB
      
      AQI(II,IQ)  =1d0/sqrt(2d0)*(
     $   -MI(II,1,IQ)*(
C!     $               E*EQ(IQ)*NLP(1)+G/CW*NLP(2)*(T3Q(IQ)-EQ(IQ)*SW**2) 
     $                G*(T3Q(IQ)*NL(2)+(SW/CW)*NL(1)/6) 
     $               +G*MQ(IQ)/2./MW*NIXB )
     $   +MI(II,2,IQ)*(
     $               E*EQ(IQ)*NLP(1)-G/CW*NLP(2)*eq(IQ)*SW**2
C*  the same as above, but more simple:   G*(SW/CW)*EQ(IQ)*NL(1)
C!     $               +G*MQ(IQ)/2./MW*NL(IX) )
     $               -G*MQ(IQ)/2./MW*NIXB
     $                )
     $                          )
     
      BQI(II,IQ)  =1./sqrt(2.)*(
     $   -MI(II,1,IQ)*(
C!     $               E*EQ(IQ)*NLP(1)+G/CW*NLP(2)*(T3Q(IQ)-EQ(IQ)*SW**2) 
     $                 G*(T3Q(IQ)*NL(2)+(SW/CW)*NL(1)/6)    
     $               -G*MQ(IQ)/2./MW*NIXB )
     $   +MI(II,2,IQ)*(
     $               -E*EQ(IQ)*NLP(1)+G/CW*NLP(2)*eq(IQ)*SW**2
C* the same as above, but more simple:  -G*(SW/CW)*EQ(IQ)*NL(1)
C!     $               -G*MQ(IQ)/2./MW*NL(IX)
     $               -G*MQ(IQ)/2./MW*NIXB
     $                )
     $                          )

        FQQ(IQ)=FQQ(IQ)+ 
     $     (-1d0/4d0)*(AQI(II,IQ)**2-BQI(II,IQ)**2)/
     $                  (msq(II,IQ)**2-(mne+mq(iq))**2)
  

c      print *,'FQQ(',II,IQ,')=', 
c     $ (-1d0/4d0)*(AQI(II,IQ)**2-BQI(II,IQ)**2)/
c     $  (msq(II,IQ)**2-(mne+mq(iq))**2)
 
         
         FQH(IQ)=FQH(IQ)+
     $    MQ(IQ)*CZ(II)*CQ(II,IQ)/MH(II)**2

c      print *,'FQH(',II,IQ,')=',MQ(IQ)*CZ(II)*CQ(II,IQ)/MH(II)**2
c     $,MQ(IQ),CZ(II),CQ(II,IQ),MH(II)
     
     
        GQ(IQ)=GQ(IQ)+ 
     $     (-1d0/4d0)*(AQI(II,IQ)**2+BQI(II,IQ)**2)/
     $                  (msq(II,IQ)**2-(mne+mq(iq))**2)**2
     
        DQ(IQ)=DQ(IQ)+
     $     (1d0/4d0)*(AQI(II,IQ)**2+BQI(II,IQ)**2)/
     $                  (msq(II,IQ)**2-(mne+mq(iq))**2)
        
       
c        if(IQ.le.3)   print *,'aq,bq',II,IQ,aqi(II,IQ),bqi(II,IQ)
c     .,DQ(IQ),mq(iq)
	ENDDO
c        if(IQ.le.3)  print *,DQ(IQ)
C!	DQ(IQ)=DQ(IQ)-G**2/4d0/MW**2*(NL(4)**2-NL(3)**2)/2d0*T3Q(IQ)
        DQ(IQ)=DQ(IQ)+G**2/4d0/MW**2*(NL(4)**2-NL(3)**2)/2d0*T3Q(IQ)
c      if(IQ.le.3)print *,'G=',G,' MW=',MW,' NL(4)=',NL(4),'NL(3)=',NL(3)	
      ENDDO
c------------------------------------------------------------c

cccccccccccccccccccccc gluon partccccccccccccccccccccccccccccccccccc      
       TQT=0d0
       BD =0d0
       BS=0d0
       B1D=0d0
       B1S=0d0
       B2S=0d0
c          print *,B1S

c      print *,((AQI(II,IQ),II=1,2),IQ=1,6)
    
       DO IQ=1,6
	DO II=1,2
          
	  DO JH=1,2
	    TQT=TQT+CZ(JH)/MH(JH)**2*
     $              CSQ(JH,II,IQ)/MSQ(II,IQ)**2/96d0/PI
          ENDDO
	  
	  IF(IQ.GT.3)BD =BD+(AQI(II,IQ)**2-BQI(II,IQ)**2)*MQ(IQ)*
     $            CMD(1,MSQ(II,IQ),MQ(IQ),MNE)/32d0/PI	  
	  BS =BS+(AQI(II,IQ)**2+BQI(II,IQ)**2)*MNE*
     $            CMD(2,MSQ(II,IQ),MQ(IQ),MNE)/32d0/PI
     
     
     	  
	  IF(IQ.GT.3)B1D=B1D+(AQI(II,IQ)**2-BQI(II,IQ)**2)*MQ(IQ)*
     $            CMD(3,MSQ(II,IQ),MQ(IQ),MNE)/12d0/PI	  
	    
         	
     
	  B1S=B1S+(AQI(II,IQ)**2+BQI(II,IQ)**2)*MNE*
     $             CMD(4,MSQ(II,IQ),MQ(IQ),MNE)/12d0/PI	
     
c       print *,'b1s',
c     $ II,IQ,CMD(4,MSQ(II,IQ),MQ(IQ),MNE) ,p1,p2

     

       
	  B2S=B2S+(AQI(II,IQ)**2+BQI(II,IQ)**2)*
     $             CMD(5,MSQ(II,IQ),MQ(IQ),MNE)/48d0/PI
          
        ENDDO
      ENDDO
      	
      
      BGLU=-TQT+BD+BS-MNE/2d0*B2S-MNE**2/4d0*(B1D+B1S)

c       print *,'BGLU=',BGLU      
c      print *,'TQT=',TQT*96*PI
       
      ap=0d0
      an=0d0
      DO IQ=1,3
        AP=AP+DQ(IQ)*delq(1,IQ)/sqrt(2d0)
        AN=AN+DQ(IQ)*delq(2,IQ)/sqrt(2d0)
c 	print *,IQ, DQ(IQ),delq(1,IQ),delq(2,IQ),ap,an,IQ,SW
      ENDDO
c------------------------------------------------------------c
      
ccccccccccccccccccccc STRUCTURE FUNCTION ccccccccccccccccccccccccccc
CC      Call SetCtq5(3)
      XMIN=1.5d-05
      XMAX=0.999d0
      XD  =(XMAX-XMIN)/10000d0
      DO IQ=-5,5,1
        PDF(1,IQ)=0d0
        PDF(2,IQ)=0d0
        DO I=1,10001
	  X=XMIN+XD*(I-1d0)
          PDF(1,IQ)=PDF(1,IQ)+X*Ctq5Pdf(IQ, X, QR,IPRT)*XD
	ENDDO
        PDF(2,IQ)=PDF(1,IQ)
      ENDDO   
      PDF(2,1) =PDF(1,2)
      PDF(2,2) =PDF(1,1)
      PDF(2,-1)=PDF(1,-2)
      PDF(2,-2)=PDF(1,-1)
c------------------------------------------------------------c

cccccccccccccccccccccc  PARTS OF FN  ccccccccccccccccccccccccccccccc      
      DO IN=1,2
      FNUCL(1,IN)=0.
      FNUCL(2,IN)=0.
      FNUCL(3,IN)=0.
      FNUCL(4,IN)=0.
      FNUCL(5,IN)=0.
      DO IQ=1,3
          FNUCL(1,IN)=FNUCL(1,IN)+ FNTQ(IN,IQ)/MQ(IQ)*
     $         (FQQ(IQ)+FQH(IQ)-mq(IQ)*mne*GQ(IQ)/2d0)

c          print *, 'FQQ', IN,IQ, FNTQ(IN,IQ), FQQ(IQ),FQH(IQ), GQ(IQ)
        ENDDO
	
c	print *,'part1',FNUCL(1,IN)
	DO  IQ=4,6
          FNUCL(2,IN)=FNUCL(2,IN)+2./27.*FNTG(IN)*FQH(IQ)/MQ(IQ)
	ENDDO
c	print *,'part2',FNUCL(2,IN)
	
	DO IQ=1,5
	  FNUCL(3,IN)=FNUCL(3,IN) -3./2.*MNE*GQ(IQ)*
     $    (PDF(IN,IQ)+PDF(IN,-IQ))
c          print *,'GQ(IQ)',GQ(IQ),IQ
        ENDDO
c	print *,'part3',FNUCL(3,IN)
	
	FNUCL(4,IN)= -8d0*PI/9d0*BGLU*FNTG(IN)
c	print *,'part4',FNUCL(4,IN),BGLU,FNTG(IN),IN
	
	FNUCL(5,IN)=  3d0/2d0*ALPS*MNE*
     $    (B2S+1d0/2d0*MNE*(B1D+B1S))*PDF(IN,0)
c 	print *,'part5',FNUCL(5,IN),IN,ALPS,
c     $  B2S,b1d,b1s,PDF(IN,0)
    
     
	
c        print *,'***** FNUCL**',IN,'******'
c        print *,(FNUCL(kk,IN),KK=1,5)
c        print *,'**************************'
	
      ENDDO
      
     
      DO IN=1,2
        FNUCL(6,IN)=0.
        DO I=1,5 
	  FNUCL(6,IN)=FNUCL(6,IN)+MNUCL(IN)*FNUCL(I,IN)
   
c	  print *, FNUCL(6,IN),MNUCL(IN),FNUCL(I,IN),I,IN
	  
	ENDDO
      ENDDO
c------------------------------------------------------------c
	 
       
cccccccccccccccccccccc  THE FINAL ANSWER ccccccccccccccccccccccccccccccc      
 
        
 
        SIGMA0PROT=4d0*MRNUCL(1)**2/PI*FNUCL(6,1)**2*GEVTOPBN
        SIGMA0NEUT=4d0*MRNUCL(2)**2/PI*FNUCL(6,2)**2*GEVTOPBN

        XJ=1d0/2d0
	
        SIGMASPROT=32d0*MRNUCL(1)**2/PI*AP**2*(XJ+1d0)*XJ*GEVTOPBN
        SIGMASNEUT=32d0*MRNUCL(2)**2/PI*AN**2*(XJ+1d0)*XJ*GEVTOPBN
	
c	print *, SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT


	RETURN
	END     
      
c------------------------------------------------------------c
c------------------------------------------------------------c
c------------------------------------------------------------c
      
ccccccccccccccccccccccccccc FUNCTION CMD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      REAL*8 FUNCTION CMD(I,MSQ,MQ,MNE)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      REAL*8 LAM
      
      common/check/p1,p2

      DSQN=MSQ**2-MNE**2
      DQN =MQ**2 -MNE**2
      DSQQ=MSQ**2-MQ**2
      
      D3M=  msq**2+mq**2-MNE**2
      D3M1= msq**2-mq**2-MNE**2
      
      R1  =DSQQ/MNE**2
      R2  =DQN /MSQ**2
      R3  =DSQN/MQ**2
      
      DEL =2d0*MNE**2*(MQ**2+MSQ**2)-MNE**4-(MSQ**2-MQ**2)**2
      ADE =abs(DEL)
      IF(ADE.lt.1E-20)PRINT *,'DEL=0(!!!)'
      LAM=0d0
      IF(DEL.GE.0)LAM=2d0*ATAN(sqrt(ADE)/D3M)/sqrt(ADE)
      IF(DEL.LT.0)LAM=LOG((D3M+sqrt(ADE))/(D3M-sqrt(ADE)))/sqrt(ADE)

      CMD=0d0
      
      IF    (I.eq.1) THEN
        
	CMD=1d0/DEL*(R2/3d0 - 2d0/3d0*R3-5d0/3d0+
     $              (2d0*msq**2-2d0/3d0*mne**2)*LAM)
c       print *,'CI1:',DEL, R2/3d0,- 2d0/3d0*R3, 5d0/3d0,
c     $                (2d0*msq**2-2d0/3d0*mne**2)*LAM
	
      ELSEIF(I.eq.2) THEN
      
        CMD=(LOG(MSQ**2/MQ**2)-D3M1*LAM)/2d0/MNE**4+
     $ (((mq**4-mq**2*msq**2)/mne**2-7./3.*mq**2-2./3.*DSQN)*LAM+
     $ R2/3.+R1+2./3.)/DEL
      
      ELSEIF(I.eq.3) THEN
      
        CMD=-3d0/DEL**2*D3M+LAM/DEL*(-1d0+6d0*mq**2*msq**2/del)
	
      ELSEIF(I.eq.4) THEN
      
        CMD=
     $   ((2d0*log(msq/mq) - D3M1*LAM)/2d0/mne**2
     $   -1d0/msq**2
     $   -mq**2*d3m1/DEL*LAM)/mne**4
     $ +((mq**2/mne**4-(1d0-mq**2/mne**2)**2/msq**2+1d0/2d0/mne**2)
     $   +3d0*mq**2/DEL*
     $   (1d0 +  R1 + (-R1*mq**2-2d0*mq**2-msq**2+mne**2)*LAM))/DEL
     
c      p1=del*lam
           
c      p2=(-mq**2*d3m1)/mne**4

      
      ELSEIF(I.eq.5) THEN
    	
         CMD=1d0/2d0/mne**4*(log(msq**2/mq**2)-D3M1*LAM)-
     $      1d0/DEL*(LAM*(2d0*DSQN+3d0*MQ**2+R1*mq**2)-3d0-R1)
     
      ELSE
      
      print *,'wrong "I" IN CMD(I,...) !!!! check I:=',I
     
      ENDIF
      
      RETURN
      END
c------------------------------------------------------------c


c------------------------------------------------------------c
C============================================================================
C Interface to CTEQ5L parton distribution
C Uses CTEQ5L from main ISAJET
C============================================================================
      DOUBLE PRECISION FUNCTION CTQ5PDF (IPARTON, X, Q, IPRT)
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
      INTEGER IPARTON,IPRT,IFL
      DOUBLE PRECISION X,Q,CTEQ5L,SUM5L,RAT5L
      IFL=IPARTON
      IF(IFL.GE.3) IFL=-IFL
      IF(IFL.EQ.-1) THEN
        SUM5L=CTEQ5L(-1,X,Q)
        RAT5L=CTEQ5L(-2,X,Q)
        CTQ5PDF=SUM5L/(1.D0+RAT5L)
      ELSEIF(IFL.EQ.-2) THEN
        SUM5L=CTEQ5L(-1,X,Q)
        RAT5L=CTEQ5L(-2,X,Q)
        CTQ5PDF=SUM5L*RAT5L/(1.D0+RAT5L)
      ELSEIF(IFL.GE.-5.AND.IFL.LE.5) THEN
        CTQ5PDF=CTEQ5L(IFL,X,Q)
      ELSE
        CTQ5PDF=0
      ENDIF
      RETURN
      END
#endif
