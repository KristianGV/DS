#include "PILOT.inc"
      SUBROUTINE BK2MVSM(GIN,GOUT)
!
!Purpose: To convert the book notation G's to Martin and Vaughn notation.
!         The same subroutine will convert back as well.
!         For SMRGE it is sufficient to convert the Yukawa matrices.
!
!With f_{u,d,e} in Baer/Tata notation and Y_{u,d,e} in Martin/Vaughn
!notation, the conversion is:
!
!     Y_{u,d,e}= {f_{u,d,e}}^T
!
      IMPLICIT NONE
!
      DOUBLE PRECISION GIN(32),GOUT(32),YU(3,3),YD(3,3),YE(3,3)
      INTEGER I,J
!
!First set GOUT to GIN
!
     	DO I=1,32
     	  GOUT(I)=GIN(I)
     	END DO
!
!Convert input into 3x3 matrices
!
      DO I=1,3
        DO J=1,3
          YU(I,J)=GIN(3+(I-1)*3+J)
          YD(I,J)=GIN(12+(I-1)*3+J)
          YE(I,J)=GIN(21+(I-1)*3+J)
        END DO
      END DO
!
!Reset relevant parts of GOUT
!
      DO I=1,3
        DO J=1,3
          GOUT(3+(I-1)*3+J)=YU(J,I)
          GOUT(12+(I-1)*3+J)=YD(J,I)
          GOUT(21+(I-1)*3+J)=YE(J,I)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE CBK2MVSM(GIN,GOUT)
!
!Purpose: To convert the book notation G's to Martin and Vaughn notation.
!         The same subroutine will convert back as well.
!         For SMRGE it is sufficient to convert the Yukawa matrices.
!
!With f_{u,d,e} in Baer/Tata notation and Y_{u,d,e} in Martin/Vaughn
!notation, the conversion is:
!
!     Y_{u,d,e}= {f_{u,d,e}}^T
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX GIN(32),GOUT(32),YU(3,3),YD(3,3),YE(3,3)
      INTEGER I,J
!
!First set GOUT to GIN
!
     	DO I=1,32
     	  GOUT(I)=GIN(I)
     	END DO
!
!Convert input into 3x3 matrices
!
      DO I=1,3
        DO J=1,3
          YU(I,J)=GIN(3+(I-1)*3+J)
          YD(I,J)=GIN(12+(I-1)*3+J)
          YE(I,J)=GIN(21+(I-1)*3+J)
        END DO
      END DO
!
!Reset relevant parts of GOUT
!
      DO I=1,3
        DO J=1,3
          GOUT(3+(I-1)*3+J)=YU(J,I)
          GOUT(12+(I-1)*3+J)=YD(J,I)
          GOUT(21+(I-1)*3+J)=YE(J,I)
        END DO
      END DO
!
      RETURN
      END
!
      FUNCTION CCON(A)
!
!Purpose: To find the conjugate of A. Again only necessary if
!         if A is complex
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A,CCON
!
      CCON=CONJG(A)
!
      RETURN
      END
!
      SUBROUTINE CDAGGER(A,DAGA)
!
!Purpose: To compute the dagger of matrix A
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),DAGA(3,3)
      INTEGER I,J
!
      DO I=1,3
        DO J=1,3
          DAGA(I,J)=CONJG(A(J,I))
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE CHDEC(Q,G,LAST)
!
!Purpose: To check if we are at the correct scale for decoupling
!         and if so, fill the vector which will be used to change
!         the threshold locations.
!         If called with LAST=1, the thresholds from the previous
!         run will be used.
!
      IMPLICIT NONE
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      DOUBLE COMPLEX MQCURR(3,3),MUCURR(3,3),MDCURR(3,3)
      DOUBLE COMPLEX MQTMP(3,3),MUTMP(3,3),MDTMP(3,3)
      DOUBLE COMPLEX MQQB(3,3),MUQB(3,3),MDQB(3,3)
      DOUBLE COMPLEX MLQB(3,3),MEQB(3,3),CMATMUL
      DOUBLE COMPLEX CID(3,3),GMASS(601),SUM,G(601),COMBR(3,3)
      DOUBLE COMPLEX VALQ(3),VALU(3),VALD(3),VALL(3),VALE(3)
      DOUBLE COMPLEX VLQ(3,3)
      DOUBLE PRECISION Q,MAX
      INTEGER I,J,K,NSQ,NSU,NSD,NSL,NSE,LAST
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!Check if user wanted to keep isajet's thresholds
!
      IF(ISADEC.EQ.1.OR.LAST.EQ.1)THEN
        DO I=1,20
          BELOW(I)=1
        END DO
!
!I can only use this method of checking for the thresholds
!if the location of the thresholds has been set in advance.
!The running will hit each threshold exactly, and EPS will be
!set in the calling routine DOWNMHIGHMZ.
!
        DO I=1,3
          THSQ(I)=1
          IF((Q-QTHQL(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHQL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
          THSU(I)=1
          IF((Q-QTHUR(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHUR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
          THSD(I)=1
          IF((Q-QTHDR(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHDR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
          THSL(I)=1
          IF((Q-QTHLL(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHLL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSL(I)=0
          THSE(I)=1
          IF((Q-QTHER(I)).LT.-ABS(EPS).OR.
     $           (ABS(Q-QTHER(I)).LT.ABS(EPS).AND.EPS.LT.0))THSE(I)=0
        END DO
!
!Otherwise, this is the way each threshold is checked
!
      ELSE
        DO I=1,3
          IF(BELOW(I).EQ.0.AND.Q.LT.DSQRT(ABS(MQVA(I)))
     $                           .AND.DBLE(MQVA(I)).GT.0.D0)THEN
            BELOW(I)=1
            NEWTH(I)=Q
            THSQ(I)=0
          END IF
          IF(BELOW(I+3).EQ.0.AND.Q.LT.DSQRT(ABS(MUPVA(I)))
     $                           .AND.DBLE(MUPVA(I)).GT.0.D0)THEN
            BELOW(I+3)=1
            NEWTH(I+3)=Q
            THSU(I)=0
          END IF
          IF(BELOW(I+6).EQ.0.AND.Q.LT.DSQRT(ABS(MDVA(I)))
     $                           .AND.DBLE(MDVA(I)).GT.0.D0)THEN
            BELOW(I+6)=1
            NEWTH(I+6)=Q
            THSD(I)=0
          END IF
          IF(BELOW(I+9).EQ.0.AND.Q.LT.DSQRT(ABS(MLVA(I)))
     $                           .AND.DBLE(MLVA(I)).GT.0.D0)THEN
            BELOW(I+9)=1
            NEWTH(I+9)=Q
            THSL(I)=0
          END IF
          IF(BELOW(I+12).EQ.0.AND.Q.LT.DSQRT(ABS(MEVA(I)))
     $                           .AND.DBLE(MEVA(I)).GT.0.D0)THEN
            BELOW(I+12)=1
            NEWTH(I+12)=Q
            THSE(I)=0
          END IF
        END DO
!
        IF(BELOW(16).EQ.0.AND.Q.LT.ABS(G(108)))THEN
          BELOW(16)=1
          NEWTH(16)=Q
        END IF
        IF(Q.LT.QNSG)BELOW(17)=1
        IF(Q.LT.QNH)BELOW(18)=1
        IF(BELOW(19).EQ.0.AND.Q.LT.ABS(G(31)))THEN
          BELOW(19)=1
          NEWTH(19)=Q
        END IF
        IF(BELOW(20).EQ.0.AND.Q.LT.ABS(G(32)))THEN
          BELOW(20)=1
          NEWTH(20)=Q
        END IF
      END IF
!
!Now I can set the new Ns and save the frozen eigenvector
!
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      DO I=1,3
        NSQ=NSQ+THSQ(I)
        NSU=NSU+THSU(I)
        NSD=NSD+THSD(I)
        NSL=NSL+THSL(I)
        NSE=NSE+THSE(I)
      END DO
!
!Find the quark basis squark mass matrices for use later
!Note that for mq we need to choose which rotation
!we are using, and we use the choice made in the input
!file.
!
      DO I=1,3
        DO J=1,3
          IF(SVLQ.EQ.1)THEN
            VLQ(I,J)=VLU(I,J)
          ELSE
            VLQ(I,J)=VLD(I,J)
          END IF
        END DO
      END DO
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15.AND.
     $   (OLDNSQ-NSQ.NE.0.OR.OLDNSU-NSU.NE.0.OR.OLDNSD-NSD.NE.0.OR.
     $    OLDNSL-NSL.NE.0.OR.OLDNSE-NSE.NE.0))THEN
        DO I=1,3
          DO J=1,3
            MQCURR(I,J)=G(62+(I-1)*3+J)
            MUCURR(I,J)=G(80+(I-1)*3+J)
            MDCURR(I,J)=G(89+(I-1)*3+J)
            MLQB(I,J)=G(71+(I-1)*3+J)
            MEQB(I,J)=G(98+(I-1)*3+J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            MQTMP(I,J)=CMATMUL(0,MQCURR,VLQ,I,J)
            MUTMP(I,J)=CMATMUL(0,MUCURR,VRU,I,J)
            MDTMP(I,J)=CMATMUL(0,MDCURR,VRD,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            MQQB(I,J)=CMATMUL(1,VLQ,MQTMP,I,J)
            MUQB(I,J)=CMATMUL(1,VRU,MUTMP,I,J)
            MDQB(I,J)=CMATMUL(1,VRD,MDTMP,I,J)
          END DO
        END DO
      END IF
!
!Find the largest off-diagonal entry for use by ORTH
!
      IF(OLDNSQ-NSQ.NE.0)THEN
        OFFMAXQVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXQVAL).LE.ABS(MQQB(I,J)))THEN
                OFFMAXQ(1)=I
                OFFMAXQ(2)=J
                OFFMAXQVAL=MQQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSU-NSU.NE.0)THEN
        OFFMAXUVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXUVAL).LE.ABS(MUQB(I,J)))THEN
                OFFMAXU(1)=I
                OFFMAXU(2)=J
                OFFMAXUVAL=MUQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSD-NSD.NE.0)THEN
        OFFMAXDVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXDVAL).LE.ABS(MDQB(I,J)))THEN
                OFFMAXD(1)=I
                OFFMAXD(2)=J
                OFFMAXDVAL=MDQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSL-NSL.NE.0)THEN
        OFFMAXLVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXLVAL).LE.ABS(MLQB(I,J)))THEN
                OFFMAXL(1)=I
                OFFMAXL(2)=J
                OFFMAXLVAL=MLQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
      IF(OLDNSE-NSE.NE.0)THEN
        OFFMAXEVAL=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            IF(I.NE.J)THEN
              IF(ABS(OFFMAXEVAL).LE.ABS(MEQB(I,J)))THEN
                OFFMAXE(1)=I
                OFFMAXE(2)=J
                OFFMAXEVAL=MEQB(I,J)
              END IF
            END IF
          END DO
        END DO
      END IF
!
!Now store the eigenvalues for use by ORTH
!
      IF(OLDNSQ.EQ.3)THEN
        DO I=1,3
          VALQ(I)=MQVA(I)
        END DO
      ELSE IF(OLDNSQ.EQ.2)THEN
        DO I=1,2
          VALQ(I)=MQVA(I)
        END DO
        VALQ(3)=MQSAV(3,4,3)
      ELSE IF(OLDNSQ.EQ.1)THEN
        VALQ(1)=MQVA(1)
        VALQ(2)=MQSAV(2,4,2)
        VALQ(3)=MQSAV(3,4,3)
      END IF
!
      IF(OLDNSU.EQ.3)THEN
        DO I=1,3
          VALU(I)=MUPVA(I)
        END DO
      ELSE IF(OLDNSU.EQ.2)THEN
        DO I=1,2
          VALU(I)=MUPVA(I)
        END DO
        VALU(3)=MUPSAV(3,4,3)
      ELSE IF(OLDNSU.EQ.1)THEN
        VALU(1)=MUPVA(1)
        VALU(2)=MUPSAV(2,4,2)
        VALU(3)=MUPSAV(3,4,3)
      END IF
!
      IF(OLDNSD.EQ.3)THEN
        DO I=1,3
          VALD(I)=MDVA(I)
        END DO
      ELSE IF(OLDNSD.EQ.2)THEN
        DO I=1,2
          VALD(I)=MDVA(I)
        END DO
        VALD(3)=MDSAV(3,4,3)
      ELSE IF(OLDNSD.EQ.1)THEN
        VALD(1)=MDVA(1)
        VALD(2)=MDSAV(2,4,2)
        VALD(3)=MDSAV(3,4,3)
      END IF
!
      IF(OLDNSL.EQ.3)THEN
        DO I=1,3
          VALL(I)=MLVA(I)
        END DO
      ELSE IF(OLDNSL.EQ.2)THEN
        DO I=1,2
          VALL(I)=MLVA(I)
        END DO
        VALL(3)=MLSAV(3,4,3)
      ELSE IF(OLDNSL.EQ.1)THEN
        VALL(1)=MLVA(1)
        VALL(2)=MLSAV(2,4,2)
        VALL(3)=MLSAV(3,4,3)
      END IF
!
      IF(OLDNSE.EQ.3)THEN
        DO I=1,3
          VALE(I)=MEVA(I)
        END DO
      ELSE IF(OLDNSE.EQ.2)THEN
        DO I=1,2
          VALE(I)=MEVA(I)
        END DO
        VALE(3)=MESAV(3,4,3)
      ELSE IF(OLDNSE.EQ.1)THEN
        VALE(1)=MEVA(1)
        VALE(2)=MESAV(2,4,2)
        VALE(3)=MESAV(3,4,3)
      END IF
!
!Only set the Rs to be different from identity if we are below
!a threshold.
!
      DO I=1,3
        DO J=1,3
          RQTOT(I,J)=CID(I,J)
          RUPTOT(I,J)=CID(I,J)
          RDTOT(I,J)=CID(I,J)
          RLTOT(I,J)=CID(I,J)
          RETOT(I,J)=CID(I,J)
        END DO
      END DO
!
!Set the rotation if we are below at least one threshold. The first
!dimension in the saved matrix is the number of active squarks in
!the region we use the rotation.
!For only 1 active squark, we need to find the compound rotation.
!MATMUL cannot be used due to the nature of the matrix with three
!indices.
!Each time I save a rotation, I need to carry out my fix using ORTH
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15)THEN
        IF(NSQ.LT.3.AND.OLDNSQ.NE.0)THEN
          IF(OLDNSQ.EQ.3)THEN
            CALL ORTH(VALQ,MQVE,OFFMAXQ,OFFMAXQVAL,VLQ,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RQSAV(2,I,K)*MQVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALQ,COMBR,OFFMAXQ,OFFMAXQVAL,VLQ,0)
          END IF
          IF(OLDNSQ.EQ.3.AND.OLDNSQ-NSQ.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RQSAV(2,I,J)=MQVE(I,J)
                IF(NSQ.LT.2)RQSAV(1,I,J)=MQVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSQ.EQ.2.AND.OLDNSQ-NSQ.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RQSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSQ.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RQTOT(I,J)=RQSAV(NSQ,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSU.LT.3.AND.OLDNSU.NE.0)THEN
          IF(OLDNSU.EQ.3)THEN
            CALL ORTH(VALU,MUPVE,OFFMAXU,OFFMAXUVAL,VRU,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RUPSAV(2,I,K)*MUPVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALU,COMBR,OFFMAXU,OFFMAXUVAL,VRU,0)
          END IF
          IF(OLDNSU.EQ.3.AND.OLDNSU-NSU.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RUPSAV(2,I,J)=MUPVE(I,J)
                IF(NSU.LT.2)RUPSAV(1,I,J)=MUPVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSU.EQ.2.AND.OLDNSU-NSU.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RUPSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSU.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RUPTOT(I,J)=RUPSAV(NSU,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSD.LT.3.AND.OLDNSD.NE.0)THEN
          IF(OLDNSD.EQ.3)THEN
            CALL ORTH(VALD,MDVE,OFFMAXD,OFFMAXDVAL,VRD,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RDSAV(2,I,K)*MDVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALD,COMBR,OFFMAXD,OFFMAXDVAL,VRD,0)
          END IF
          IF(OLDNSD.EQ.3.AND.OLDNSD-NSD.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RDSAV(2,I,J)=MDVE(I,J)
                IF(NSD.LT.2)RDSAV(1,I,J)=MDVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSD.EQ.2.AND.OLDNSD-NSD.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RDSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSD.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RDTOT(I,J)=RDSAV(NSD,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSL.LT.3.AND.OLDNSL.NE.0)THEN
          IF(OLDNSL.EQ.3)THEN
            CALL ORTH(VALL,MLVE,OFFMAXL,OFFMAXLVAL,CID,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RLSAV(2,I,K)*MLVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALL,COMBR,OFFMAXL,OFFMAXLVAL,CID,0)
          END IF
          IF(OLDNSL.EQ.3.AND.OLDNSL-NSL.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RLSAV(2,I,J)=MLVE(I,J)
                IF(NSL.LT.2)RLSAV(1,I,J)=MLVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSL.EQ.2.AND.OLDNSL-NSL.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RLSAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSL.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RLTOT(I,J)=RLSAV(NSL,I,J)
              END DO
            END DO
          END IF
        END IF
!
        IF(NSE.LT.3.AND.OLDNSE.NE.0)THEN
          IF(OLDNSE.EQ.3)THEN
            CALL ORTH(VALE,MEVE,OFFMAXE,OFFMAXEVAL,CID,0)
          ELSE
            DO I=1,3
              DO J=1,3
                SUM=(0.D0,0.D0)
                DO K=1,3
                  SUM=SUM+RESAV(2,I,K)*MEVE(K,J)
                END DO
                COMBR(I,J)=SUM
              END DO
            END DO
            CALL ORTH(VALE,COMBR,OFFMAXE,OFFMAXEVAL,CID,0)
          END IF
          IF(OLDNSE.EQ.3.AND.OLDNSE-NSE.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RESAV(2,I,J)=MEVE(I,J)
                IF(NSE.LT.2)RESAV(1,I,J)=MEVE(I,J)
              END DO
            END DO
          ELSE IF(OLDNSE.EQ.2.AND.OLDNSE-NSE.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RESAV(1,I,J)=COMBR(I,J)
              END DO
            END DO
          END IF
          IF(NSE.NE.0)THEN
            DO I=1,3
              DO J=1,3
                RETOT(I,J)=RESAV(NSE,I,J)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!Save the values of the sfermion matrices at each scale
!to use when running up. The first dimension of each array
!is the number of active sfermions above the threshold.
!The fourth entry in the second index is for the eigenvalues.
!
      DO I=1,3
        IF(OLDNSQ-NSQ.NE.0)THEN
          DO K=OLDNSQ,NSQ+1,-1
            DO J=1,3
              MQSAV(K,I,J)=G(62+(I-1)*3+J)
            END DO
            MQSAV(K,4,I)=VALQ(I)
          END DO
        END IF
        IF(OLDNSU-NSU.NE.0)THEN
          DO K=OLDNSU,NSU+1,-1
            DO J=1,3
              MUPSAV(K,I,J)=G(80+(I-1)*3+J)
            END DO
            MUPSAV(K,4,I)=VALU(I)
          END DO
        END IF
        IF(OLDNSD-NSD.NE.0)THEN
          DO K=OLDNSD,NSD+1,-1
            DO J=1,3
              MDSAV(K,I,J)=G(89+(I-1)*3+J)
            END DO
            MDSAV(K,4,I)=VALD(I)
          END DO
        END IF
        IF(OLDNSL-NSL.NE.0)THEN
          DO K=OLDNSL,NSL+1,-1
            DO J=1,3
              MLSAV(K,I,J)=G(71+(I-1)*3+J)
            END DO
            MLSAV(K,4,I)=VALL(I)
          END DO
        END IF
        IF(OLDNSE-NSE.NE.0)THEN
          DO K=OLDNSE,NSE+1,-1
            DO J=1,3
              MESAV(K,I,J)=G(98+(I-1)*3+J)
            END DO
            MESAV(K,4,I)=VALE(I)
          END DO
        END IF
      END DO
!
!Construct the squark matrices with saved entries when some squarks
!have decoupled.
!
      CALL MASSSQM(G)
!
!Rotate into the squark mass basis, remove the appropriate terms and
!rotate back.
!
      IF(NSQ-OLDNSQ.NE.0.OR.NSU-OLDNSU.NE.0.OR.NSD-OLDNSD.NE.0
     $     .OR.NSL-OLDNSL.NE.0.OR.NSE-OLDNSE.NE.0)THEN
        CALL REMSF(G)
      END IF
!
      RETURN
      END
!
      SUBROUTINE CHINT(Q,G)
!
!Purpose: To check if the sfermions must be reintroduced
!         into the running. If so, we must use the saved values
!         of the sfermion mass matrices as a boundary
!         condition on the upwards running.
!
      IMPLICIT NONE
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      DOUBLE COMPLEX CID(3,3),G(601)
      DOUBLE PRECISION Q
      INTEGER I,J,K,NSQ,NSU,NSD,NSL,NSE
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!We can use the same statements as before to set the various thetas.
!This is because we know the position of each threshold and the
!calling routine is hitting each one exactly.
!
      DO I=1,3
        THSQ(I)=1
        IF((Q-QTHQL(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHQL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
        THSU(I)=1
        IF((Q-QTHUR(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHUR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
        THSD(I)=1
        IF((Q-QTHDR(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHDR(I)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
        THSL(I)=1
        IF((Q-QTHLL(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHLL(I)).LT.ABS(EPS).AND.EPS.LT.0))THSL(I)=0
        THSE(I)=1
        IF((Q-QTHER(I)).LT.-ABS(EPS).OR.
     $         (ABS(Q-QTHER(I)).LT.ABS(EPS).AND.EPS.LT.0))THSE(I)=0
      END DO
!
!Next set the values of the number of active sfermions
!
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      DO I=1,3
        NSQ=NSQ+THSQ(I)
        NSU=NSU+THSU(I)
        NSD=NSD+THSD(I)
        NSL=NSL+THSL(I)
        NSE=NSE+THSE(I)
      END DO
!
!Only set the Rs to be different from identity if we are below
!a threshold. I do not need to reset the entries of the saved
!matrices as in CHDEC...
!
      DO I=1,3
        DO J=1,3
          RQTOT(I,J)=CID(I,J)
          RUPTOT(I,J)=CID(I,J)
          RDTOT(I,J)=CID(I,J)
          RLTOT(I,J)=CID(I,J)
          RETOT(I,J)=CID(I,J)
        END DO
      END DO
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15.AND.NSQ+NSU+NSD+NSL+NSE.GT.0)THEN
        DO I=1,3
          DO J=1,3
            IF(NSQ.LT.3.AND.NSQ.GT.0)RQTOT(I,J)=RQSAV(NSQ,I,J)
            IF(NSU.LT.3.AND.NSU.GT.0)RUPTOT(I,J)=RUPSAV(NSU,I,J)
            IF(NSD.LT.3.AND.NSD.GT.0)RDTOT(I,J)=RDSAV(NSD,I,J)
            IF(NSL.LT.3.AND.NSL.GT.0)RLTOT(I,J)=RLSAV(NSL,I,J)
            IF(NSE.LT.3.AND.NSE.GT.0)RETOT(I,J)=RESAV(NSE,I,J)
          END DO
        END DO
      END IF
!
!Now insert the values of the sfermion matrices at each scale.
!
      DO I=1,3
        DO J=1,3
          IF(OLDNSQ-NSQ.NE.0)G(62+(I-1)*3+J)=MQSAV(NSQ,I,J)
          IF(OLDNSU-NSU.NE.0)G(80+(I-1)*3+J)=MUPSAV(NSU,I,J)
          IF(OLDNSD-NSD.NE.0)G(89+(I-1)*3+J)=MDSAV(NSD,I,J)
          IF(OLDNSL-NSL.NE.0)G(71+(I-1)*3+J)=MLSAV(NSL,I,J)
          IF(OLDNSE-NSE.NE.0)G(98+(I-1)*3+J)=MESAV(NSE,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      FUNCTION CMATMUL(DAG,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B
!         If DAG=1 compute dagger of A
!         If DAG=2 compute dagger of B
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),B(3,3)
      DOUBLE COMPLEX CMATMUL,DAGA(3,3),DAGB(3,3)
      INTEGER LOOP,DAG,I,J
!
      CMATMUL=(0.D0,0.D0)
      IF(DAG.EQ.1)THEN
        CALL CDAGGER(A,DAGA)
        DO LOOP=1,3
          CMATMUL=CMATMUL+DAGA(I,LOOP)*B(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.2)THEN
        CALL CDAGGER(B,DAGB)
        DO LOOP=1,3
          CMATMUL=CMATMUL+A(I,LOOP)*DAGB(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.0)THEN
        DO LOOP=1,3
          CMATMUL=CMATMUL+A(I,LOOP)*B(LOOP,J)
        END DO
      ELSE
        WRITE(*,*)'WRONG DAG IN MATMUL'
      END IF
!
      RETURN
      END
!
      FUNCTION CMODSQ(A)
!
!Purpose: Calculates the mod squared of a value. This will only
!         be a necessary function if A is complex
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A,CMODSQ
!
      CMODSQ=CONJG(A)*A
!
      RETURN
      END
!
      FUNCTION CON(A)
!
!Purpose: To find the conjugate of A. Again only necessary if
!         if A is complex
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A,CON
!
      CON=A
!
      RETURN
      END
!
      FUNCTION CRE(STAR,A,B)
!
!Purpose: To find the real part of A*B with:
!         If STAR=1 compute conjugate of A
!         If STAR=2 compute conjugate of B
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A,B
      DOUBLE COMPLEX CCON,CRE
      INTEGER STAR
!
      CRE=0.D0
      IF(STAR.EQ.1)THEN
        CRE=DBLE(CCON(A)*B)
      ELSE IF(STAR.EQ.2)THEN
        CRE=DBLE(A*CCON(B))
      ELSE IF(STAR.EQ.0)THEN
        CRE=DBLE(A*B)
      ELSE
        WRITE(*,*)'WRONG STAR'
      END IF
!
      RETURN
      END
!
      SUBROUTINE CRGE215(T,G,F)
!
!Contains threshold RGEs for gauge and yukawas.
!
!CRGE215 IS FOR RUNNNG THRESHOLDS, NO TILDES, FOR THE FIRST TIME
!The sfermion thresholds are not all distinct. All left-squarks
!are at the same point, all up-right-squarks are at another point,
!all down-right-squarks at another, all left-sleptons another and
!all right-sleptons another.
!
!     G(  1) = g_1         G(  2) = g_2         G(  3) = g_3
!     G(  4) = FU(1,1)     G(  5) = FU(1,2)     G( 12) = FU(3,3)
!     G( 13) = FD(1,1)     G( 22) = FE(1,1)     G( 30) = FE(3,3)
!     G( 31) = mu          G( 32) = V_U         G( 33) = V_D
!
!     G( 34) = LU(1,1)     G( 43) = LD(1,1)     G( 52) = LE(1,1)
!
!     G( 61) = VEV_SM      G( 62) = LAM_SM
!
!     G( 63) = GTPQ(1,1)   G( 72) = GTPL(1,1)   G( 81) = GTPU(1,1)
!     G( 90) = GTPD(1,1)   G( 99) = GTPE(1,1)   G(108) = GTP_Hu
!     G(109) = GTP_Hd      G(110) = GTQ(1,1)    G(119) = GTL(1,1)
!     G(128) = GT_Hu       G(129) = GT_Hd       G(130) = GTSQ(1,1)
!     G(139) = GTSU(1,1)   G(148) = GTSD(1,1)   G(157) = FTUQ(1,1)
!     G(166) = FTDQ(1,1)   G(175) = FTEL(1,1)   G(184) = FTUU(1,1)
!     G(193) = FTDD(1,1)   G(202) = FTEE(1,1)   G(211) = sGTP_Hu
!     G(212) = cGTP_Hd     G(213) = sGT_Hu      G(214) = cGT_Hd
!
!     G(215) = mu(M)
!
!This is the BT version which receives G in book notation
!
      IMPLICIT NONE
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE PRECISION T
      DOUBLE COMPLEX G(215)
      DOUBLE COMPLEX F(215)
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3)
      DOUBLE COMPLEX YU(3,3),YD(3,3),YE(3,3)
      DOUBLE COMPLEX AU(3,3),AD(3,3),AE(3,3)
      DOUBLE COMPLEX HU(3,3),HD(3,3),HE(3,3),DH(3,3,3)
      DOUBLE COMPLEX MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3),DM(5,3,3)
      DOUBLE COMPLEX LU(3,3),LD(3,3),LE(3,3),LYU(3,3),LYD(3,3),LYE(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3),GTPE(3,3)
      DOUBLE COMPLEX GTQ(3,3),GTL(3,3),GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE COMPLEX FTUQ(3,3),FTDQ(3,3),FTEL(3,3)
      DOUBLE COMPLEX FTUU(3,3),FTDD(3,3),FTEE(3,3)
!
      DOUBLE COMPLEX CMATMUL,CSFMUL,CTRACE,SUM
      DOUBLE COMPLEX CMODSQ,CCON,CRE
!
!These are used in the calculation of the RGEs which contain
!thresholds
!
      DOUBLE COMPLEX DUMU1(3,3),DUMU2(3,3),DUMD1(3,3),DUMD2(3,3)
      DOUBLE COMPLEX DUME1(3,3),DUME2(3,3),DUMGRKMU1(3,3),TDUMGRKMU
      DOUBLE COMPLEX DUMLUD1(3,3),DUMLUD2(3,3),DUMLUD(3,3),TDUMLUD
!
      DOUBLE COMPLEX FUS(3,3),FDS(3,3),FES(3,3)
      DOUBLE COMPLEX FTUQD(3,3),FTUUD(3,3),FTDDD(3,3),FTDQD(3,3)
      DOUBLE COMPLEX FTELD(3,3),FTEED(3,3)
      DOUBLE COMPLEX GTPUS(3,3),GTSUS(3,3),GTQS(3,3),GTLS(3,3)
      DOUBLE COMPLEX GTPQS(3,3),GTPLS(3,3),GTSQS(3,3),GTPDS(3,3)
      DOUBLE COMPLEX GTPES(3,3),GTSDS(3,3),GTPUT(3,3),GTQT(3,3)
      DOUBLE COMPLEX GTLT(3,3),GTPQT(3,3),GTPLT(3,3),GTPDT(3,3)
      DOUBLE COMPLEX GTPET(3,3),GTSUT(3,3),GTSQT(3,3),GTSDT(3,3)
!
      DOUBLE COMPLEX FUFUD(3,3),FDFDD(3,3),FEFED(3,3)
      DOUBLE COMPLEX YUDYU(3,3),YDDYD(3,3),YEDYE(3,3)
      DOUBLE COMPLEX LULUD(3,3),LDLDD(3,3),LELED(3,3)
      DOUBLE COMPLEX LYUDLYU(3,3),LYDDLYD(3,3),LYEDLYE(3,3)
      DOUBLE COMPLEX LYUDLYU2(3,3),LYDDLYD2(3,3),LYEDLYE2(3,3)
      DOUBLE COMPLEX LYUDLYU3(3,3),LYDDLYD3(3,3),LYEDLYE3(3,3)
      DOUBLE COMPLEX TLYUDLYU,TLYDDLYD,TLYEDLYE
      DOUBLE COMPLEX TLYUDLYU2,TLYDDLYD2,TLYEDLYE2
      DOUBLE COMPLEX TLYUDLYU3,TLYDDLYD3,TLYEDLYE3
      DOUBLE COMPLEX LYUDLYULYDDLYD(3,3),LYDDLYDLYUDLYU(3,3)
      DOUBLE COMPLEX TLYUDLYULYDDLYD
      DOUBLE COMPLEX FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
      DOUBLE COMPLEX YUYUD(3,3),YDYDD(3,3),YEYED(3,3)
      DOUBLE COMPLEX YUYDD(3,3),YDYUD(3,3)
      DOUBLE COMPLEX TFUDFU,TFDDFD,TFEDFE
      DOUBLE COMPLEX TYUDYU,TYDDYD,TYEDYE
      DOUBLE COMPLEX TLUDLU,TLDDLD,TLEDLE
      DOUBLE COMPLEX FUFUDFU(3,3),FUFUDFD(3,3),FDFDDFU(3,3)
      DOUBLE COMPLEX FDFDDFD(3,3),FEFEDFE(3,3)
      DOUBLE COMPLEX LULUDLU(3,3),LULUDLD(3,3),LDLDDLU(3,3)
      DOUBLE COMPLEX LDLDDLD(3,3),LELEDLE(3,3)
      DOUBLE COMPLEX SQFTUQDFTUQ(3,3),SQFTDQDFTDQ(3,3),SQGTQTGTQS(3,3)
      DOUBLE COMPLEX SQGTQTFTUQ(3,3),SQGTQTFTDQ(3,3),SQGTPQTGTPQS(3,3)
      DOUBLE COMPLEX SQGTPQTFTUQ(3,3),SQGTPQTFTDQ(3,3),SQGTSQTGTSQS(3,3)
      DOUBLE COMPLEX SUGTPUSGTPUT(3,3),SUGTSUSGTSUT(3,3)
      DOUBLE COMPLEX SUFTUUFTUUD(3,3),SUFTUUGTPUT(3,3),SDFTDDFTDDD(3,3)
      DOUBLE COMPLEX SDFTDDGTPDT(3,3),SDGTPDSGTPDT(3,3)
      DOUBLE COMPLEX SDGTSDSGTSDT(3,3),SLFTELDFTEL(3,3),SLGTLTGTLS(3,3)
      DOUBLE COMPLEX SLGTLTFTEL(3,3),SLGTPLTGTPLS(3,3),SLGTPLTFTEL(3,3)
      DOUBLE COMPLEX SEGTPESGTPET(3,3),SEFTEEFTEED(3,3),SEFTEEGTPET(3,3)
      DOUBLE COMPLEX MGTPHUSQ,MGTPHDSQ,MGTHUSQ,MGTHDSQ
      DOUBLE COMPLEX MSGTPHUSQ,MCGTPHDSQ,MSGTHUSQ,MCGTHDSQ
      DOUBLE COMPLEX Y2,H,Y4,CHI4
!
      DOUBLE COMPLEX B1U(3,3),B1D(3,3),B1E(3,3)
      DOUBLE COMPLEX BGRKMU,BETALAM1,BETALAM2

!
!The following are for two loop and soft RGEs
!
      DOUBLE COMPLEX DUM2U1(3,3),DUM2U2(3,3)
      DOUBLE COMPLEX DUM2D1(3,3),DUM2D2(3,3)
      DOUBLE COMPLEX DUM2E1(3,3),DUM2E2(3,3)
      DOUBLE COMPLEX DUM2U(3,3),DUM2D(3,3),DUM2E(3,3)
      DOUBLE COMPLEX DUM1GRKMU(3,3),DUM2GRKMU(3,3)
!
      DOUBLE COMPLEX YUYUDYUYUD(3,3),YDYDDYDYDD(3,3)
      DOUBLE COMPLEX YEYEDYEYED(3,3),YUDYUYUDYU(3,3)
      DOUBLE COMPLEX YDDYDYDDYD(3,3),YEDYEYEDYE(3,3)
      DOUBLE COMPLEX YDDYDYUDYU(3,3),YUDYUYDDYD(3,3)
      DOUBLE COMPLEX YUYDDYDYUD(3,3)
!
      DOUBLE COMPLEX TYUDYUYUDYU,TYDDYDYDDYD,TYEDYEYEDYE
      DOUBLE COMPLEX TYUDYUYDDYD,TYDDYDYUDYU
!
      DOUBLE COMPLEX BETA1U(3,3),BETA2U(3,3),BETA1D(3,3),BETA2D(3,3)
      DOUBLE COMPLEX BETA1E(3,3),BETA2E(3,3)
      DOUBLE COMPLEX BETA2GRKMU,BETA1VU,BETA1VD,BETA2VU,BETA2VD,BETA2B
      DOUBLE COMPLEX B2YMU(3,3),B2YMD(3,3),B2YME(3,3),BETAVEV1,BETAVEV2
!
      DOUBLE COMPLEX ID(3,3),MVMU,MVMUM
      DOUBLE PRECISION B1LP(3),B1LPM(3),B2LPSM(3,3),B2LPM(3,3)
      DOUBLE PRECISION CM(3,3),CSM(3,3)
      DOUBLE PRECISION PI,Q
!
      INTEGER I,J,NG,ND,NE,NNU,NU,NSQ,NSU,NSD,NSL,NSE,NSH,NH,NSW,NSG
      INTEGER THLH,THHH,THSH,THSQ(3),THSU(3),THSD(3),THSL(3)
      INTEGER THSE(3),THSB,THSW,THGL
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $    ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $    ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $    ,ID(3,3)/(1.D0,0.D0)/
      DATA B1LPM(1)/6.6D0/,B1LPM(2)/1.D0/,B1LPM(3)/-3.D0/
      DATA B2LPM(1,1)/7.96D0/,B2LPM(1,2)/5.4D0/,B2LPM(1,3)/17.6D0/
      DATA B2LPM(2,1)/1.8D0/,B2LPM(2,2)/25.D0/,B2LPM(2,3)/24.D0/
      DATA B2LPM(3,1)/2.2D0/,B2LPM(3,2)/9.D0/,B2LPM(3,3)/14.D0/
      DATA CM(1,1)/5.2D0/,CM(1,2)/2.8D0/,CM(1,3)/3.6D0/
      DATA CM(2,1)/6.D0/,CM(2,2)/6.D0/,CM(2,3)/2.D0/
      DATA CM(3,1)/4.D0/,CM(3,2)/4.D0/,CM(3,3)/0.D0/
!
!Set all F's and betas to zero
!
      DO I=1,215
        F(I)=(0.D0,0.D0)
      END DO
      DO I=1,3
        DO J=1,3
          B1U(I,J)=(0.D0,0.D0)
          B1D(I,J)=(0.D0,0.D0)
          B1E(I,J)=(0.D0,0.D0)
          BETA1U(I,J)=(0.D0,0.D0)
          BETA2U(I,J)=(0.D0,0.D0)
          BETA1D(I,J)=(0.D0,0.D0)
          BETA2D(I,J)=(0.D0,0.D0)
          BETA1E(I,J)=(0.D0,0.D0)
          BETA2E(I,J)=(0.D0,0.D0)
          B2YMU(I,J)=(0.D0,0.D0)
          B2YMD(I,J)=(0.D0,0.D0)
          B2YME(I,J)=(0.D0,0.D0)
        END DO
      END DO
      BGRKMU=(0.D0,0.D0)
      BETA2GRKMU=(0.D0,0.D0)
      BETA1VU=(0.D0,0.D0)
      BETA1VD=(0.D0,0.D0)
      BETA2VU=(0.D0,0.D0)
      BETA2VD=(0.D0,0.D0)
      BETALAM1=(0.D0,0.D0)
      BETALAM2=(0.D0,0.D0)
      BETAVEV1=(0.D0,0.D0)
      BETAVEV2=(0.D0,0.D0)
!
      PI=4.D0*ATAN(1.D0)
      Q=SSQSTEP
      IF(Q.LT.1.D0)THEN
        WRITE(*,*)'ERROR IN Q: ',Q
        STOP
      END IF
      NG=3.D0
      NU=3
      ND=3
      NE=3
      NNU=3
!
      NSQ=3
      IF((Q-QTHQL(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHQL(1)).LT.ABS(EPS).AND.EPS.LT.0))NSQ=0
      NSU=3
      IF((Q-QTHUR(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHUR(1)).LT.ABS(EPS).AND.EPS.LT.0))NSU=0
      NSD=3
      IF((Q-QTHDR(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHDR(1)).LT.ABS(EPS).AND.EPS.LT.0))NSD=0
      NSL=3
      IF((Q-QTHLL(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHLL(1)).LT.ABS(EPS).AND.EPS.LT.0))NSL=0
      NSE=3
      IF((Q-QTHER(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHER(1)).LT.ABS(EPS).AND.EPS.LT.0))NSE=0
!
      IF ((Q-QNSH).GT.EPS.OR.
     $         (ABS(Q-QNSH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSH=2
      ELSE
        NSH=0
      END IF
      IF ((Q-QNSG).GT.EPS.OR.
     $         (ABS(Q-QNSG).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF ((Q-QNH).GT.EPS.OR.
     $         (ABS(Q-QNH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF ((Q-QTHSB).GT.EPS.OR.
     $         (ABS(Q-QTHSB).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        THSB=1
      ELSE
        THSB=0
      END IF
      IF ((Q-QTHSW).GT.EPS.OR.
     $         (ABS(Q-QTHSW).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSW=1
        THSW=1
      ELSE
        NSW=0
        THSW=0
      END IF
      THLH=1
      THHH=NH/2 !This works so long as THHH is an integer variable
      IF(THHH.NE.1)THEN !Perform check
        IF(THHH.NE.0)WRITE(*,*)'ERROR IN THHH'
      END IF
      THSH=NSH/2
      DO I=1,3
        THSQ(I)=1
        IF((Q-QTHQL(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHQL(1)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
        THSU(I)=1
        IF((Q-QTHUR(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHUR(1)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
        THSD(I)=1
        IF((Q-QTHDR(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHDR(1)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
        THSL(I)=NSL/3
        THSE(I)=NSE/3
      END DO
      THGL=NSG
!
!Convert input into 3x3 matrices.
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
!
!MV notation Yukawas
!
          YU(I,J)=G(3+(J-1)*3+I)
          YD(I,J)=G(12+(J-1)*3+I)
          YE(I,J)=G(21+(J-1)*3+I)
!
          LU(I,J)=G(33+(I-1)*3+J)
          LD(I,J)=G(42+(I-1)*3+J)
          LE(I,J)=G(51+(I-1)*3+J)
          GTPQ(I,J)=G(62+(I-1)*3+J)
          GTPL(I,J)=G(71+(I-1)*3+J)
          GTPU(I,J)=G(80+(I-1)*3+J)
          GTPD(I,J)=G(89+(I-1)*3+J)
          GTPE(I,J)=G(98+(I-1)*3+J)
          GTQ(I,J)=G(109+(I-1)*3+J)
          GTL(I,J)=G(118+(I-1)*3+J)
          GTSQ(I,J)=G(129+(I-1)*3+J)
          GTSU(I,J)=G(138+(I-1)*3+J)
          GTSD(I,J)=G(147+(I-1)*3+J)
          FTUQ(I,J)=G(156+(I-1)*3+J)
          FTDQ(I,J)=G(165+(I-1)*3+J)
          FTEL(I,J)=G(174+(I-1)*3+J)
          FTUU(I,J)=G(183+(I-1)*3+J)
          FTDD(I,J)=G(192+(I-1)*3+J)
          FTEE(I,J)=G(201+(I-1)*3+J)
        END DO
      END DO
!
      MVMU=G(31)
      MVMUM=G(215)
!
!SET THE TILDE TERMS TO THEIR SM COUNTERPARTS.
!
      DO I=1,3
        DO J=1,3
          GTPQ(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(62+(I-1)*3+J)=GTPQ(I,J)
          GTPL(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(71+(I-1)*3+J)=GTPL(I,J)
          GTPU(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(80+(I-1)*3+J)=GTPU(I,J)
          GTPD(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(89+(I-1)*3+J)=GTPD(I,J)
          GTPE(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(98+(I-1)*3+J)=GTPE(I,J)
          GTQ(I,J)=G(2)*ID(I,J)
          G(109+(I-1)*3+J)=GTQ(I,J)
          GTL(I,J)=G(2)*ID(I,J)
          G(118+(I-1)*3+J)=GTL(I,J)
          GTSQ(I,J)=G(3)*ID(I,J)
          G(129+(I-1)*3+J)=GTSQ(I,J)
          GTSU(I,J)=G(3)*ID(I,J)
          G(138+(I-1)*3+J)=GTSU(I,J)
          GTSD(I,J)=G(3)*ID(I,J)
          G(147+(I-1)*3+J)=GTSD(I,J)
          FTUQ(I,J)=YU(J,I)
          G(156+(I-1)*3+J)=FTUQ(I,J)
          FTDQ(I,J)=YD(J,I)
          G(165+(I-1)*3+J)=FTDQ(I,J)
          FTEL(I,J)=YE(J,I)
          G(174+(I-1)*3+J)=FTEL(I,J)
          FTUU(I,J)=YU(J,I)
          G(183+(I-1)*3+J)=FTUU(I,J)
          FTDD(I,J)=YD(J,I)
          G(192+(I-1)*3+J)=FTDD(I,J)
          FTEE(I,J)=YE(J,I)
          G(201+(I-1)*3+J)=FTEE(I,J)
        END DO
      END DO
      G(108)=SQRT(3.D0/5.D0)*G(1)
      G(109)=SQRT(3.D0/5.D0)*G(1)
      G(128)=G(2)
      G(129)=G(2)
      G(211)=(0.D0,0.D0)
      G(212)=(0.D0,0.D0)
      G(213)=(0.D0,0.D0)
      G(214)=(0.D0,0.D0)
!
!The separated out contributions are from PRD 49 4882 (1194),
!Castano,Piard,Ramond
!
      B1LP(1)=2.D0/5.D0*(17.D0/12.D0*DBLE(NU)+5.D0/12.D0*DBLE(ND)
     $        +5.D0/4.D0*DBLE(NE)+1.D0/4.D0*DBLE(NNU))+1.D0/30.D0
     $        *DBLE(NSQ)+4.D0/15.D0*DBLE(NSU)+1.D0/15.D0*DBLE(NSD)
     $        +1.D0/10.D0*DBLE(NSL)+1.D0/5.D0*DBLE(NSE)
     $        +1.D0/5.D0*DBLE(NSH)+1.D0/10.D0*DBLE(NH)
      B1LP(2)=-22.D0/3.D0+1.D0/2.D0*(DBLE(NU)+DBLE(ND))+1.D0/6.D0
     $        *(DBLE(NE)+DBLE(NNU))+1.D0/2.D0*DBLE(NSQ)
     $        +1.D0/6.D0*DBLE(NSL)+1.D0/3.D0*DBLE(NSH)+1.D0/6.D0
     $        *DBLE(NH)+4.D0/3.D0*DBLE(NSW)
      B1LP(3)=-11.D0+2.D0/3.D0*(DBLE(NU)+DBLE(ND))+1.D0/3.D0
     $        *DBLE(NSQ)+1.D0/6.D0*DBLE(NSU)+1.D0/6.D0*DBLE(NSD)
     $        +2.D0*DBLE(NSG)
      IF(THHH.EQ.0)THEN
        B2LPSM(1,1)=-(-NG*19.D0/15.D0-9.D0/50.D0)
        B2LPSM(1,2)=-(-NG*3.D0/5.D0-9.D0/10.D0)
        B2LPSM(1,3)=-(-NG*44.D0/15.D0)
        B2LPSM(2,1)=-(-NG*1.D0/5.D0-3.D0/10.D0)
        B2LPSM(2,2)=-(136.D0/3.D0-NG*49.D0/3.D0-13.D0/6.D0)
        B2LPSM(2,3)=-(-NG*4.D0)
        B2LPSM(3,1)=-(-NG*11.D0/30.D0)
        B2LPSM(3,2)=-(-NG*3.D0/2.D0)
        B2LPSM(3,3)=-(102.D0-NG*76.D0/3.D0)
        CSM(1,1)=1.7D0
        CSM(1,2)=.5D0
        CSM(1,3)=1.5D0
        CSM(2,1)=1.5D0
        CSM(2,2)=1.5D0
        CSM(2,3)=.5D0
        CSM(3,1)=2.D0
        CSM(3,2)=2.D0
        CSM(3,3)=0.D0
      END IF
!
!I need many variations on the 3x3 matrices.
!
      CALL CDAGGER(FTUQ,FTUQD)
      CALL CDAGGER(FTUU,FTUUD)
      CALL CDAGGER(FTDD,FTDDD)
      CALL CDAGGER(FTDQ,FTDQD)
      CALL CDAGGER(FTEL,FTELD)
      CALL CDAGGER(FTEE,FTEED)
      DO I=1,3
        DO J=1,3
          GTPUS(I,J)=CONJG(GTPU(I,J))
          GTSUS(I,J)=CONJG(GTSU(I,J))
          GTQS(I,J)=CONJG(GTQ(I,J))
          GTLS(I,J)=CONJG(GTL(I,J))
          GTPQS(I,J)=CONJG(GTPQ(I,J))
          GTPLS(I,J)=CONJG(GTPL(I,J))
          GTSQS(I,J)=CONJG(GTSQ(I,J))
          GTPDS(I,J)=CONJG(GTPD(I,J))
          GTPES(I,J)=CONJG(GTPE(I,J))
          GTSDS(I,J)=CONJG(GTSD(I,J))
          GTPUT(I,J)=GTPU(J,I)
          GTQT(I,J)=GTQ(J,I)
          GTLT(I,J)=GTL(J,I)
          GTPQT(I,J)=GTPQ(J,I)
          GTPLT(I,J)=GTPL(J,I)
          GTPDT(I,J)=GTPD(J,I)
          GTPET(I,J)=GTPE(J,I)
          GTSUT(I,J)=GTSU(J,I)
          GTSQT(I,J)=GTSQ(J,I)
          GTSDT(I,J)=GTSD(J,I)
        END DO
      END DO
!
!Now all the matrix multiples
!
      DO I=1,3
        DO J=1,3
          YUDYU(I,J)=CMATMUL(1,YU,YU,I,J)
          YDDYD(I,J)=CMATMUL(1,YD,YD,I,J)
          YEDYE(I,J)=CMATMUL(1,YE,YE,I,J)
          YUYUD(I,J)=CMATMUL(2,YU,YU,I,J)
          YDYDD(I,J)=CMATMUL(2,YD,YD,I,J)
          YEYED(I,J)=CMATMUL(2,YE,YE,I,J)
          YUYDD(I,J)=CMATMUL(2,YU,YD,I,J)
          YDYUD(I,J)=CMATMUL(2,YD,YU,I,J)
!
          FUFUD(I,J)=CMATMUL(2,FU,FU,I,J)
          FDFDD(I,J)=CMATMUL(2,FD,FD,I,J)
          FEFED(I,J)=CMATMUL(2,FE,FE,I,J)
          FUDFU(I,J)=CMATMUL(1,FU,FU,I,J)
          FDDFD(I,J)=CMATMUL(1,FD,FD,I,J)
          FEDFE(I,J)=CMATMUL(1,FE,FE,I,J)
          LULUD(I,J)=CMATMUL(2,LU,LU,I,J)
          LDLDD(I,J)=CMATMUL(2,LD,LD,I,J)
          LELED(I,J)=CMATMUL(2,LE,LE,I,J)
        END DO
      END DO
!
      IF(THHH.EQ.0)THEN
        MSGTPHUSQ=CONJG(G(211))*G(211)
        MCGTPHDSQ=CONJG(G(212))*G(212)
        MSGTHUSQ=CONJG(G(213))*G(213)
        MCGTHDSQ=CONJG(G(214))*G(214)
      END IF
      MGTPHUSQ=CONJG(G(108))*G(108)
      MGTPHDSQ=CONJG(G(109))*G(109)
      MGTHUSQ=CONJG(G(128))*G(128)
      MGTHDSQ=CONJG(G(129))*G(129)
!
      TYUDYU=CTRACE(YUDYU)
      TYDDYD=CTRACE(YDDYD)
      TYEDYE=CTRACE(YEDYE)
!
      TFUDFU=CTRACE(FUDFU)
      TFDDFD=CTRACE(FDDFD)
      TFEDFE=CTRACE(FEDFE)
      TLUDLU=CTRACE(LULUD)
      TLDDLD=CTRACE(LDLDD)
      TLEDLE=CTRACE(LELED)
!
      DO I=1,3
        DO J=1,3
          FUFUDFU(I,J)=CMATMUL(0,FUFUD,FU,I,J)
          FUFUDFD(I,J)=CMATMUL(0,FUFUD,FD,I,J)
          FDFDDFU(I,J)=CMATMUL(0,FDFDD,FU,I,J)
          FDFDDFD(I,J)=CMATMUL(0,FDFDD,FD,I,J)
          FEFEDFE(I,J)=CMATMUL(0,FEFED,FE,I,J)
          LULUDLU(I,J)=CMATMUL(0,LULUD,LU,I,J)
          LULUDLD(I,J)=CMATMUL(0,LULUD,LD,I,J)
          LDLDDLU(I,J)=CMATMUL(0,LDLDD,LU,I,J)
          LDLDDLD(I,J)=CMATMUL(0,LDLDD,LD,I,J)
          LELEDLE(I,J)=CMATMUL(0,LELED,LE,I,J)
          SQFTUQDFTUQ(I,J)=CSFMUL(THSQ,FTUQD,FTUQ,I,J)
          SQFTDQDFTDQ(I,J)=CSFMUL(THSQ,FTDQD,FTDQ,I,J)
          SQGTQTGTQS(I,J)=CSFMUL(THSQ,GTQT,GTQS,I,J)
          SQGTQTFTUQ(I,J)=CSFMUL(THSQ,GTQT,FTUQ,I,J)
          SQGTQTFTDQ(I,J)=CSFMUL(THSQ,GTQT,FTDQ,I,J)
          SQGTPQTGTPQS(I,J)=CSFMUL(THSQ,GTPQT,GTPQS,I,J)
          SQGTPQTFTUQ(I,J)=CSFMUL(THSQ,GTPQT,FTUQ,I,J)
          SQGTPQTFTDQ(I,J)=CSFMUL(THSQ,GTPQT,FTDQ,I,J)
          SQGTSQTGTSQS(I,J)=CSFMUL(THSQ,GTSQT,GTSQS,I,J)
          SUFTUUFTUUD(I,J)=CSFMUL(THSU,FTUU,FTUUD,I,J)
          SUFTUUGTPUT(I,J)=CSFMUL(THSU,FTUU,GTPUT,I,J)
          SUGTPUSGTPUT(I,J)=CSFMUL(THSU,GTPUS,GTPUT,I,J)
          SUGTSUSGTSUT(I,J)=CSFMUL(THSU,GTSUS,GTSUT,I,J)
          SDFTDDFTDDD(I,J)=CSFMUL(THSD,FTDD,FTDDD,I,J)
          SDFTDDGTPDT(I,J)=CSFMUL(THSD,FTDD,GTPDT,I,J)
          SDGTPDSGTPDT(I,J)=CSFMUL(THSD,GTPDS,GTPDT,I,J)
          SDGTSDSGTSDT(I,J)=CSFMUL(THSD,GTSDS,GTSDT,I,J)
          SLFTELDFTEL(I,J)=CSFMUL(THSL,FTELD,FTEL,I,J)
          SLGTLTGTLS(I,J)=CSFMUL(THSL,GTLT,GTLS,I,J)
          SLGTLTFTEL(I,J)=CSFMUL(THSL,GTLT,FTEL,I,J)
          SLGTPLTGTPLS(I,J)=CSFMUL(THSL,GTPLT,GTPLS,I,J)
          SLGTPLTFTEL(I,J)=CSFMUL(THSL,GTPLT,FTEL,I,J)
          SEFTEEFTEED(I,J)=CSFMUL(THSE,FTEE,FTEED,I,J)
          SEFTEEGTPET(I,J)=CSFMUL(THSE,FTEE,GTPET,I,J)
          SEGTPESGTPET(I,J)=CSFMUL(THSE,GTPES,GTPET,I,J)
        END DO
      END DO
!
!These are the two loop terms. All in MV notation
!
      IF(SW2LP.EQ.1)THEN
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=CMATMUL(0,YUYUD,YUYUD,I,J)
            YDYDDYDYDD(I,J)=CMATMUL(0,YDYDD,YDYDD,I,J)
            YEYEDYEYED(I,J)=CMATMUL(0,YEYED,YEYED,I,J)
            YUDYUYUDYU(I,J)=CMATMUL(0,YUDYU,YUDYU,I,J)
            YDDYDYDDYD(I,J)=CMATMUL(0,YDDYD,YDDYD,I,J)
            YEDYEYEDYE(I,J)=CMATMUL(0,YEDYE,YEDYE,I,J)
            YDDYDYUDYU(I,J)=CMATMUL(0,YDDYD,YUDYU,I,J)
            YUDYUYDDYD(I,J)=CMATMUL(0,YUDYU,YDDYD,I,J)
            YUYDDYDYUD(I,J)=CMATMUL(0,YUYDD,YDYUD,I,J)
          END DO
        END DO
!
        TYUDYUYUDYU=CTRACE(YUDYUYUDYU)
        TYUDYUYDDYD=CTRACE(YUDYUYDDYD)
        TYDDYDYUDYU=CTRACE(YDDYDYUDYU)
        TYDDYDYDDYD=CTRACE(YDDYDYDDYD)
        TYEDYEYEDYE=CTRACE(YEDYEYEDYE)
!
!These are SM terms for the two loop running below m_H
!I am going to use LYU for the MV notation SM Yukawa
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=CMATMUL(1,LYU,LYU,I,J)
              LYDDLYD(I,J)=CMATMUL(1,LYD,LYD,I,J)
              LYEDLYE(I,J)=CMATMUL(1,LYE,LYE,I,J)
            END DO
          END DO
          TLYUDLYU=CTRACE(LYUDLYU)
          TLYDDLYD=CTRACE(LYDDLYD)
          TLYEDLYE=CTRACE(LYEDLYE)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=CMATMUL(0,LYUDLYU,LYUDLYU,I,J)
              LYDDLYD2(I,J)=CMATMUL(0,LYDDLYD,LYDDLYD,I,J)
              LYEDLYE2(I,J)=CMATMUL(0,LYEDLYE,LYEDLYE,I,J)
              LYUDLYULYDDLYD(I,J)=CMATMUL(0,LYUDLYU,LYDDLYD,I,J)
              LYDDLYDLYUDLYU(I,J)=CMATMUL(0,LYDDLYD,LYUDLYU,I,J)
              DUMLUD1(I,J)=LYUDLYU(I,J)+LYDDLYD(I,J)
            END DO
          END DO
          TLYUDLYU2=CTRACE(LYUDLYU2)
          TLYDDLYD2=CTRACE(LYDDLYD2)
          TLYEDLYE2=CTRACE(LYEDLYE2)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=CMATMUL(0,DUMLUD1,LYDDLYD,I,J)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=CMATMUL(0,LYUDLYU2,LYUDLYU,I,J)
              LYDDLYD3(I,J)=CMATMUL(0,LYDDLYD2,LYDDLYD,I,J)
              LYEDLYE3(I,J)=CMATMUL(0,LYEDLYE2,LYEDLYE,I,J)
              DUMLUD(I,J)=CMATMUL(0,LYUDLYU,DUMLUD2,I,J)
            END DO
          END DO
          TLYUDLYU3=CTRACE(LYUDLYU3)
          TLYDDLYD3=CTRACE(LYDDLYD3)
          TLYEDLYE3=CTRACE(LYEDLYE3)
          TLYUDLYULYDDLYD=CTRACE(LYUDLYULYDDLYD)
          TDUMLUD=CTRACE(DUMLUD)
!
          Y2=3.D0*TLYUDLYU+3.D0*TLYDDLYD+TLYEDLYE
          H=3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
          Y4=(83.D0/40.D0*G(1)**2+27.D0/8.D0*G(2)**2
     $       +28.D0*G(3)**2)*TLYUDLYU+(-1.D0/40.D0*G(1)**2
     $       +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TLYDDLYD
     $       +(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TLYEDLYE
          CHI4=9.D0/4.D0*(3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
     $         -2.D0/3.D0*TLYUDLYULYDDLYD)
        END IF
      ELSE
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=(0.D0,0.D0)
            YDYDDYDYDD(I,J)=(0.D0,0.D0)
            YEYEDYEYED(I,J)=(0.D0,0.D0)
            YUDYUYUDYU(I,J)=(0.D0,0.D0)
            YDDYDYDDYD(I,J)=(0.D0,0.D0)
            YEDYEYEDYE(I,J)=(0.D0,0.D0)
            YDDYDYUDYU(I,J)=(0.D0,0.D0)
            YUDYUYDDYD(I,J)=(0.D0,0.D0)
            YUYDDYDYUD(I,J)=(0.D0,0.D0)
          END DO
        END DO
        TYUDYUYUDYU=(0.D0,0.D0)
        TYUDYUYDDYD=(0.D0,0.D0)
        TYDDYDYUDYU=(0.D0,0.D0)
        TYDDYDYDDYD=(0.D0,0.D0)
        TYEDYEYEDYE=(0.D0,0.D0)
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=(0.D0,0.D0)
              LYDDLYD(I,J)=(0.D0,0.D0)
              LYEDLYE(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU=(0.D0,0.D0)
          TLYDDLYD=(0.D0,0.D0)
          TLYEDLYE=(0.D0,0.D0)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=(0.D0,0.D0)
              LYDDLYD2(I,J)=(0.D0,0.D0)
              LYEDLYE2(I,J)=(0.D0,0.D0)
              LYUDLYULYDDLYD(I,J)=(0.D0,0.D0)
              LYDDLYDLYUDLYU(I,J)=(0.D0,0.D0)
              DUMLUD1(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU2=(0.D0,0.D0)
          TLYDDLYD2=(0.D0,0.D0)
          TLYEDLYE2=(0.D0,0.D0)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=(0.D0,0.D0)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=(0.D0,0.D0)
              LYDDLYD3(I,J)=(0.D0,0.D0)
              LYEDLYE3(I,J)=(0.D0,0.D0)
              DUMLUD(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU3=(0.D0,0.D0)
          TLYDDLYD3=(0.D0,0.D0)
          TLYEDLYE3=(0.D0,0.D0)
          TLYUDLYULYDDLYD=(0.D0,0.D0)
          TDUMLUD=(0.D0,0.D0)
!
          Y2=(0.D0,0.D0)
          H=(0.D0,0.D0)
          Y4=(0.D0,0.D0)
          CHI4=(0.D0,0.D0)
        END IF
      END IF
!
!The threshold gauge running with a change at m_H to SM
!
      DO I=1,3
        SUM=0.D0
        DO J=1,3
          IF(THHH.EQ.0)THEN
            SUM=SUM+B2LPSM(I,J)*G(J)**2
          ELSE
            SUM=SUM+B2LPM(I,J)*G(J)**2
          END IF
        END DO
        IF(THHH.EQ.0)THEN
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CSM(I,1)*TLYUDLYU+CSM(I,2)*TLYDDLYD
     $         +CSM(I,3)*TLYEDLYE)))
        ELSE
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CM(I,1)*TYUDYU+CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        END IF
      END DO
!
!Next the full Yukawas
!
      DO I=1,3
        DO J=1,3
          DUMU1(I,J)=THSH*SQFTUQDFTUQ(I,J)
     $               +4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $               +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
          DUMU2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMD1(I,J)=THSH*SQFTDQDFTDQ(I,J)
     $               +1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $               +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
          DUMD2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUME1(I,J)=THSH*SLFTELDFTEL(I,J)+THSB*SEGTPESGTPET(I,J)
          DUME2(I,J)=2.D0*THSH*SEFTEEFTEED(I,J)
     $               +3.D0*THSW*SLGTLTGTLS(I,J)+THSB*SLGTPLTGTPLS(I,J)
!
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U2(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2D1(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)+
     $                  YEYEDYEYED(I,J)
            DUM2D2(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E1(I,J)=DUM2D1(I,J)
            DUM2E2(I,J)=DUM2U2(I,J)
          ELSE
            DUM2U1(I,J)=(0.D0,0.D0)
            DUM2U2(I,J)=(0.D0,0.D0)
            DUM2D1(I,J)=(0.D0,0.D0)
            DUM2D2(I,J)=(0.D0,0.D0)
            DUM2E1(I,J)=(0.D0,0.D0)
            DUM2E2(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Here are the two loop terms
!
          IF(SW2LP.EQ.1)THEN
            DUM2U(I,J)=-3.D0*CTRACE(DUM2U1)*ID(I,J)-YDDYD(I,J)
     $                 *CTRACE(DUM2U2)-9.D0*YUDYU(I,J)*TYUDYU
     $                 -4.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)-2.D0
     $                 *YDDYDYUDYU(I,J)+(16.D0*G(3)**2+4.D0/5.D0
     $                 *G(1)**2)*TYUDYU*ID(I,J)+(6.D0*G(2)**2
     $                 +2.D0/5.D0*G(1)**2)*YUDYU(I,J)+2.D0/5.D0
     $                 *G(1)**2*YDDYD(I,J)+(-16.D0/9.D0*G(3)**4
     $                 +8.D0*G(3)**2*G(2)**2+136.D0/45.D0
     $                 *G(3)**2*G(1)**2+15.D0/2.D0*G(2)**4
     $                 +G(2)**2*G(1)**2+2743.D0/450.D0
     $                 *G(1)**4)*ID(I,J)
            DUM2D(I,J)=-3.D0*CTRACE(DUM2D1)*ID(I,J)-3.D0*YUDYU(I,J)
     $                 *TYUDYU-3.D0*YDDYD(I,J)*CTRACE(DUM2D2)
     $                 -4.D0*YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)
     $                 -2.D0*YUDYUYDDYD(I,J)+(16.D0*G(3)**2
     $                 -2.D0/5.D0*G(1)**2)*TYDDYD*ID(I,J)+6.D0/5.D0
     $                 *G(1)**2*TYEDYE*ID(I,J)+4.D0/5.D0*G(1)**2
     $                 *YUDYU(I,J)+(6.D0*G(2)**2+4.D0/5.D0
     $                 *G(1)**2)*YDDYD(I,J)+(-16.D0/9.D0*G(3)**4
     $                 +8.D0*G(3)**2*G(2)**2+8.D0/9.D0
     $                 *G(3)**2*G(1)**2+15.D0/2.D0*G(2)**4
     $                 +G(2)**2*G(1)**2+287.D0/90.D0
     $                 *G(1)**4)*ID(I,J)
            DUM2E(I,J)=-3.D0*CTRACE(DUM2E1)*ID(I,J)-3.D0*YEDYE(I,J)
     $                 *CTRACE(DUM2E2)-4.D0*YEDYEYEDYE(I,J)+(16.D0
     $                 *G(3)**2-2.D0/5.D0*G(1)**2)*TYDDYD
     $                 *ID(I,J)+6.D0/5.D0*G(1)**2*TYEDYE*ID(I,J)
     $                 +6.D0*G(2)**2*YEDYE(I,J)+(15.D0/2.D0
     $                 *G(2)**4+9.D0/5.D0*G(2)**2*G(1)**2
     $                 +27.D0/2.D0*G(1)**4)*ID(I,J)
          ELSE
            DUM2U(I,J)=(0.D0,0.D0)
            DUM2D(I,J)=(0.D0,0.D0)
            DUM2E(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the Yukawas
!
      DO I=1,3
        DO J=1,3
!
          B1U(I,J)=1.D0/2.D0*(3.D0*FUFUDFU(I,J)+FDFDDFU(I,J))
     $             +CMATMUL(0,FU,DUMU1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUMU2,FU,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(128))
     $              *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $              *CONJG(G(108))*SQGTPQTFTUQ(I,J))
     $             -4.D0/3.D0*THSB*THSH*CONJG(G(108))
     $              *SUFTUUGTPUT(I,J)
     $             +FU(I,J)*3.D0*TFUDFU
     $             +1.D0/2.D0*THSH*FU(I,J)*(3.D0*THSW*MGTHUSQ
     $              +THSB*MGTPHUSQ)
     $             -FU(I,J)*(17.D0/20.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1D(I,J)=1.D0/2.D0*(3.D0*FDFDDFD(I,J)+FUFUDFD(I,J))
     $             +CMATMUL(0,FD,DUMD1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUMD2,FD,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(129))
     $              *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $              *CONJG(G(109))*SQGTPQTFTDQ(I,J))
     $             -2.D0/3.D0*THSB*THSH*CONJG(G(109))
     $              *SDFTDDGTPDT(I,J)
     $             +FD(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FD(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FD(I,J)*(3.D0/12.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1E(I,J)=3.D0/2.D0*FEFEDFE(I,J)
     $             +CMATMUL(0,FE,DUME1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUME2,FE,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(129))
     $              *SLGTLTFTEL(I,J)+THSB
     $              *CONJG(G(109))*SLGTPLTFTEL(I,J))
     $             -2.D0*THSB*THSH*CONJG(G(109))
     $              *SEFTEEGTPET(I,J)
     $             +FE(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FE(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FE(I,J)*(9.D0/4.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2)
!
          IF(SW2LP.EQ.1)THEN
            B2YMU(I,J)=CMATMUL(0,YU,DUM2U,I,J)
            B2YMD(I,J)=CMATMUL(0,YD,DUM2D,I,J)
            B2YME(I,J)=CMATMUL(0,YE,DUM2E,I,J)
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Convert into form readable by RKSTP. The transpose in BETA2 takes
!account of the differences in notation.
!
          F(3+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I)
          F(12+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YMD(J,I)
          F(21+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YME(J,I)
        END DO
      END DO
!
!The lambdas use the same dummy matrices. I will reuse the betas
!Only find the lambdas if we are below m_H
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            B1U(I,J)=(0.D0,0.D0)
            B1D(I,J)=(0.D0,0.D0)
            B1E(I,J)=(0.D0,0.D0)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            IF(SW2LP.EQ.1)THEN
              DUM2U(I,J)=3.D0/2.D0*LYUDLYU2(I,J)
     $                  -LYUDLYULYDDLYD(I,J)
     $                  -1.D0/4.D0*LYDDLYDLYUDLYU(I,J)
     $                  +11.D0/4.D0*LYDDLYD2(I,J)
     $                  +Y2*(5.D0/4.D0*LYDDLYD(I,J)-9.D0/4.D0
     $                  *LYUDLYU(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(62)**2*ID(I,J)-2.D0*G(62)*(3.D0
     $                  *LYUDLYU(I,J)+LYDDLYD(I,J))+(221.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYUDLYU(I,J)-(17.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYDDLYD(I,J)+Y4*ID(I,J)+((7.D0/150.D0
     $                  +2.D0/3.D0*NG)*G(1)**4-9.D0/20.D0*G(1)**2
     $                  *G(2)**2+19.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2D(I,J)=3.D0/2.D0*LYDDLYD2(I,J)
     $                  -LYDDLYDLYUDLYU(I,J)
     $                  -1.D0/4.D0*LYUDLYULYDDLYD(I,J)
     $                  +11.D0/4.D0*LYUDLYU2(I,J)
     $                  +Y2*(5.D0/4.D0*LYUDLYU(I,J)-9.D0/4.D0
     $                  *LYDDLYD(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(62)**2*ID(I,J)-2.D0*G(62)*(3.D0
     $                  *LYDDLYD(I,J)+LYUDLYU(I,J))+(161.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYDDLYD(I,J)-(77.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYUDLYU(I,J)+Y4*ID(I,J)+(-(37.D0/300.D0
     $                  -4.D0/15.D0*NG)*G(1)**4-27.D0/20.D0*G(1)**2
     $                  *G(2)**2+31.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2E(I,J)=3.D0/2.D0*LYEDLYE2(I,J)
     $                  -Y2*9.D0/4.D0*LYEDLYE(I,J)-CHI4*ID(I,J)
     $                  +3.D0/2.D0*G(62)**2*ID(I,J)-6.D0*G(62)
     $                  *LYEDLYE(I,J)+(441.D0/80.D0*G(1)**2
     $                  +117.D0/16.D0*G(2)**2)*LYEDLYE(I,J)
     $                  +Y4*ID(I,J)+((21.D0/100.D0+8.D0/5.D0*NG)
     $                  *G(1)**4+27.D0/20.D0*G(1)**2*G(2)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4)*ID(I,J)
            ELSE
              DUM2U(I,J)=(0.D0,0.D0)
              DUM2D(I,J)=(0.D0,0.D0)
              DUM2E(I,J)=(0.D0,0.D0)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            B1U(I,J)=1.D0/2.D0*(3.D0*LULUDLU(I,J)+LDLDDLU(I,J)
     $                -4.D0*LDLDDLU(I,J))
     $               +CMATMUL(0,LU,DUMU1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUMU2,LU,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(213))
     $                *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $                *CONJG(G(211))*SQGTPQTFTUQ(I,J))
     $               -4.D0/3.D0*THSB*THSH*CONJG(G(211))
     $                *SUFTUUGTPUT(I,J)
     $               +LU(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LU(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LU(I,J)*(17.D0/20.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1D(I,J)=1.D0/2.D0*(3.D0*LDLDDLD(I,J)+LULUDLD(I,J)
     $                -4.D0*LULUDLD(I,J))
     $               +CMATMUL(0,LD,DUMD1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUMD2,LD,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(214))
     $                *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $                *CONJG(G(212))*SQGTPQTFTDQ(I,J))
     $               -2.D0/3.D0*THSB*THSH*CONJG(G(212))
     $                *SDFTDDGTPDT(I,J)
     $               +LD(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LD(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LD(I,J)*(3.D0/12.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1E(I,J)=3.D0/2.D0*LELEDLE(I,J)+CMATMUL(0,LE,DUME1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUME2,LE,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(214))
     $                *SLGTLTFTEL(I,J)+THSB
     $                *CONJG(G(212))*SLGTPLTFTEL(I,J))
     $               -2.D0*THSB*THSH*CONJG(G(212))
     $                *SEFTEEGTPET(I,J)
     $               +LE(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LE(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LE(I,J)*(9.D0/4.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2)
!
            IF(SW2LP.EQ.1)THEN
              BETA2U(I,J)=CMATMUL(0,LYU,DUM2U,I,J)
              BETA2D(I,J)=CMATMUL(0,LYD,DUM2D,I,J)
              BETA2E(I,J)=CMATMUL(0,LYE,DUM2E,I,J)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
!
!Convert into form readable by RKSTP.
!
            F(33+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2U(J,I)
            F(42+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2D(J,I)
            F(51+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2E(J,I)
          END DO
        END DO
      END IF
!
!Next I am going to work out the gaugino terms, \mu and M_{1,2,3}
!and the running of B in MV notation
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*SUFTUUFTUUD(I,J)+3.D0*SDFTDDFTDDD(I,J)
     $                   +SEFTEEFTEED(I,J)+3.D0*SQFTUQDFTUQ(I,J)
     $                   +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
!
          IF(SW2LP.EQ.1)THEN
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2GRKMU(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
      TDUMGRKMU=CTRACE(DUMGRKMU1)
!
      BGRKMU=1.D0/2.D0*G(31)*THSH*TDUMGRKMU
     $       +1.D0/4.D0*G(31)*THSH*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ
     $        +3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       -G(31)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMU*(-3.D0*CTRACE(DUM2GRKMU)+(16.D0*G(3)**2
     $             +4.D0/5.D0*G(1)**2)*TYUDYU+(16.D0*G(3)**2
     $             -2.D0/5.D0*G(1)**2)*TYDDYD+6.D0/5.D0*G(1)**2
     $             *TYEDYE+15.D0/2.D0*G(2)**4+9.D0/5.D0*G(1)**2
     $             *G(2)**2+207.D0/50.D0*G(1)**4)
      END IF
!
!The RKSTP compatible derivative is:
!
        F(31)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $         *BETA2GRKMU
!
!V_U and V_D - from PRD 49,4882 (1994)
!
      BETA1VU=3.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-3.D0*TYUDYU
      BETA1VD=3.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-3.D0*TYDDYD
     $        -TYEDYE
      IF(SW2LP.EQ.1)THEN
        BETA2VU=3.D0/4.D0*(3.D0*TYUDYUYUDYU+3.D0*TYUDYUYDDYD)-(19.D0
     $          /10.D0*G(1)**2+9.D0/2.D0*G(2)**2+20.D0*G(3)**2)
     $          *TYUDYU-(279.D0/800.D0+1803.D0/1600.D0*3.D0)*G(1)**4
     $          -(207.D0/32.D0+357.D0/64.D0*3.D0)*G(2)**4-(27.D0/80.D0
     $          +9.D0/80.D0*3.D0)*G(1)**2*G(2)**2
        BETA2VD=3.D0/4.D0*(3.D0*TYDDYDYDDYD+3.D0*TYDDYDYUDYU
     $          +TYEDYEYEDYE)-(2.D0/5.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $          +20.D0*G(3)**2)*TYDDYD-(9.D0/5.D0*G(1)**2+3.D0/2.D0
     $          *G(2)**2)*TYEDYE-(279.D0/800.D0+1803.D0/1600.D0*3.D0)
     $          *G(1)**4-(207.D0/32.D0+357.D0/64.D0*3.D0)*G(2)**4
     $          -(27.D0/80.D0+9.D0/80.D0*3.D0)*G(1)**2*G(2)**2
      END IF
!
      F(32)=G(32)*(1.D0/16.D0/PI**2*BETA1VU+1.D0/(16.D0*PI**2)**2
     $       *BETA2VU)
      F(33)=G(33)*(1.D0/16.D0/PI**2*BETA1VD+1.D0/(16.D0*PI**2)**2
     $       *BETA2VD)
!
      IF(THHH.EQ.0)THEN
!
!Finally we have the running of the Higgs Quartic Coupling and SM VEV.
!Programmed here is the MS-bar running. It therefore needs the MS-bar
!gauge and Yukawas.
!The gauge couplings and Yukawas needed to be converted to MS-bar using
!the Martin and Vaughn conversion in hep-ph/9308222.
!The following is after the conversion, so all Yukawas and Gauge
!couplings are still in the DR-bar scheme.
!
        BETALAM1=12*G(62)**2-(9.D0/5.D0*G(1)**2+9.D0*G(2)**2)
     $           *G(62)+9.D0/4.D0*(3.D0/25.D0*G(1)**4+2.D0/5.D0
     $           *G(1)**2*G(2)**2+G(2)**4)+4.D0*Y2*G(62)-4*H
        IF(SW2LP.EQ.1)THEN
          BETALAM2=-78.D0*G(62)**3+18.D0*(3.D0/5.D0*G(1)**2
     $            +3.D0*G(2)**2)*G(62)**2-((265.D0/8.D0-10*NG)
     $            *G(2)**4-117.D0/20.D0*G(1)**2*G(2)**2
     $            -9.D0/25.D0*(229.D0/24.D0+50.D0/9.D0*NG)*G(1)**4)
     $            *G(62)+(473.D0/8.D0-8.D0*NG)*G(2)**6-3.D0/5.D0
     $            *(121.D0/24.D0+8.D0/3.D0*NG)*G(1)**2*G(2)**4
     $            -9.D0/25.D0*(239.D0/24.D0+40.D0/9.D0*NG)
     $            *G(1)**4*G(2)**2-27.D0/125.D0*(59.D0/24.D0
     $            +40.D0/9.D0*NG)*G(1)**6+(-14.D0/5.D0*G(1)**2
     $            +18.D0*G(2)**2-128.D0*G(3)**2)*TLYUDLYU2
     $            +(34.D0/5.D0*G(1)**2+18.D0*G(2)**2-128.D0
     $            *G(3)**2)*TLYDDLYD2+(-42.D0/5.D0*G(1)**2
     $            +6.D0*G(2)**2)*TLYEDLYE2-3.D0/2.D0*G(2)**4
     $            *Y2+G(62)*((83.D0/10.D0*G(1)**2+27.D0/2.D0
     $            *G(2)**2+112.D0*G(3)**2)*TLYUDLYU+(-1.D0/10.D0
     $            *G(1)**2+27.D0/2.D0*G(2)**2+112.D0*G(3)**2)
     $            *TLYDDLYD+(93.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $            *TLYEDLYE)+3.D0/5.D0*G(1)**2*((-57.D0/10.D0
     $            *G(1)**2+21.D0*G(2)**2)*TLYUDLYU+(3.D0/2.D0
     $            *G(1)**2+9.D0*G(2)**2)*TLYDDLYD+(-15.D0/2.D0
     $            *G(1)**2+11.D0*G(2)**2)*TLYEDLYE)-24.D0
     $            *G(62)**2*Y2-G(62)*H+6.D0*G(62)
     $            *TLYUDLYULYDDLYD+20.D0*(3.D0*TLYUDLYU3+3.D0*TLYDDLYD3
     $            +TLYEDLYE3)-12.D0*TDUMLUD
        END IF
!
        F(62)=1.D0/(16.D0*PI**2)*BETALAM1+1.D0/(16.D0*PI**2)**2
     $         *BETALAM2
!
!Calculate the betas for the standard model vev.
!As with lambda this is the MS-bar running with DR-bar inputs except
!v and lambda
!
      BETAVEV1=9.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-Y2
      IF(SW2LP.EQ.1)THEN
        BETAVEV2=-3.D0/2.D0*G(62)**2-(83.D0/40.D0*G(1)**2+27.D0/8.D0
     $           *G(2)**2+28.D0*G(3)**2)*TYUDYU-(-1.D0/40.D0*G(1)**2
     $           +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TYDDYD
     $           -(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TYEDYE+CHI4
     $           -27.D0/80.D0*G(1)**2*G(2)**2-(93.D0/800.D0+1.D0/2.D0
     $           *NG)*G(1)**4+(463.D0/32.D0-5.D0/2.D0*NG)*G(2)**4
      END IF
!
      F(61)=G(61)*(1.D0/(16.D0*PI**2)*BETAVEV1+1.D0/(16.D0*PI**2)**2
     $      *BETAVEV2)
!
      END IF
!
!Finally, the MSSM mu parameter
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2GRKMU(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
      BGRKMU=MVMUM*(CTRACE(DUMGRKMU1)-3.D0*G(2)**2-3.D0/5.D0
     $           *G(1)**2)
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMUM*(-3.D0*CTRACE(DUM2GRKMU)+(16.D0*G(3)**2+4.D0
     $             /5.D0*G(1)**2)*TYUDYU+(16.D0*G(3)**2-2.D0/5.D0
     $             *G(1)**2)*TYDDYD+6.D0/5.D0*G(1)**2*TYEDYE
     $             +15.D0/2.D0*G(2)**4+9.D0/5.D0*G(1)**2*G(2)**2
     $             +207.D0/50.D0*G(1)**4)
      END IF
!
!The RKSTP compatible derivative is.
!
        F(215)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $         *BETA2GRKMU
!
      RETURN
      END
!
      SUBROUTINE CRGE601(T,GCURR,FCURR)
!
!Contains the full threshold RGEs.
!
!CRGE601 IS FOR THRESHOLDS WITH THETAS.
!
!     G(  1) = g_1         G(  2) = g_2         G(  3) = g_3
!     G(  4) = FU(1,1)     G(  5) = FU(1,2)     G( 12) = FU(3,3)
!     G( 13) = FD(1,1)     G( 22) = FE(1,1)     G( 30) = FE(3,3)
!     G( 31) = M_1         G( 32) = M_2         G( 33) = M_3
!     G( 34) = AU(1,1)     G( 43) = AD(1,1)     G( 52) = AE(1,1)
!     G( 61) = M_Hu^2+MT^2 G( 62) = M_Hd^2+MT^2 G( 63) = MQ(1,1)
!     G( 72) = ML(1,1)     G( 81) = MUP(1,1)    G( 90) = MD(1,1)
!     G( 99) = ME(1,1)     G(108) = mu          G(109) = b
!     G(110) = V_U         G(111) = V_D
!
!     G(112) = LU(1,1)     G(121) = LD(1,1)     G(130) = LE(1,1)
!
!     G(139) = GTPQ(1,1)   G(148) = GTPL(1,1)   G(157) = GTPU(1,1)
!     G(166) = GTPD(1,1)   G(175) = GTPE(1,1)   G(184) = GTP_Hu
!     G(185) = GTP_Hd      G(186) = GTQ(1,1)    G(195) = GTL(1,1)
!     G(204) = GT_Hu       G(205) = GT_Hd       G(206) = GTSQ(1,1)
!     G(215) = GTSU(1,1)   G(224) = GTSD(1,1)   G(233) = FTUQ(1,1)
!     G(242) = FTDQ(1,1)   G(251) = FTEL(1,1)   G(260) = FTUU(1,1)
!     G(269) = FTDD(1,1)   G(278) = FTEE(1,1)   G(287) = sGTP_Hu
!     G(288) = cGTP_Hd     G(289) = sGT_Hu      G(290) = cGT_Hd
!
!     G(291) = g(M)_1      G(292) = g(M)_2      G(293) = g(M)_3
!     G(294) = f(M)_U(1,1) G(303) = f(M)_D(1,1) G(312) = f(M)_E(1,1)
!     G(321) = M^(MV)_1(M) G(322) = M^(MV)_2(M) G(323) = M^(MV)_3(M)
!     G(324) = AU(M)(1,1)  G(333) = AD(M)(1,1)  G(342) = AE(M)(1,1)
!     G(351) = M(M)_Hu^2   G(352) = M(M)_Hd^2   G(353) = MQ(M)(1,1)
!     G(362) = ML(M)(1,1)  G(371) = MUP(M)(1,1) G(380) = MD(M)(1,1)
!     G(389) = ME(M)(1,1)  G(398) = mu(M)       G(399) = b(M)
!
!     G(400) = TRI_U(1,1)  G(409) = TRI_D(1,1)  G(418) = TRI_E(1,1)
!     G(427) = M_HUD       G(428) = VEV_SM      G(429) = LAM_SM
!
!     G(430) = MTSF_U(1,1) G(439) = MTSF_D(1,1) G(448) = MTSF_E(1,1)
!     G(457) = GLP         G(458) = G2L         G(459) = G3L
!     G(460) = FUHU(1,1)   G(469) = FDHD(1,1)   G(478) = FEHD(1,1)
!     G(487) = FUQ(1,1)    G(496) = FUUR(1,1)   G(505) = FDQ(1,1)
!     G(514) = FDDR(1.1)   G(523) = FEL(1,1)    G(532) = FEER(1,1)
!     G(541) = GLP*RT(C^2-S^2)      G(542) = G2L*RT(C^2-S^2)
!     G(543) = GLP*(C^2-S^2)        G(544) = G2L*(C^2-S^2)
!     G(545) = sFUQ(1,1)  G(554) = sFUUR(1,1)   G(563) = cFDQ(1,1)
!     G(572) = cFDDR(1.1) G(581) = cFEL(1,1)    G(590) = cFEER(1,1)
!     G(599) = M'_1       G(600) = M'_2         G(601) = M'_3
!
!NOTE: MQ, MU, ETC... ARE USED TO DENOTE *SQUARED* SOFT MASSES
!NOTE: BOTH G(31-33) AND G(599-601) SHOULD BE REAL. THEY ARE DEFINED
!      AS COMPLEX FOR EASE OF USE. ALSO, THE MSSM GAUGINO MASS TERMS
!      ARE MV NOTATION, ie M-iM'
!NOTE: FOR THE QUARTIC MATRICES, THE VALUES IN G ARE NOT READ FROM
!      OR WRITTEN TO BY THIS PROGRAMME. THE FINAL NUMBERS FOR THESE
!      ENTRIES IN G SHOULD BE IGNORED. IN ADDITION, IN THE CASE THAT
!      THIS CODE IS EXTENDED TO INCLUDE QUARTIC RUNNING SINCE IT IS
!      EXPECTED THAT QUARTIC RGES WILL ONLY BE DERIVABLE FOR THE
!      SQUARE OF THE TERMS ABOVE. EG FUQ*FUQ.
!
!This is the BT version which receives G in book notation
!
      IMPLICIT NONE
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE PRECISION T
      DOUBLE COMPLEX G(601),GCURR(601)
      DOUBLE COMPLEX F(601),FCURR(601)
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3)
      DOUBLE COMPLEX YU(3,3),YD(3,3),YE(3,3)
      DOUBLE COMPLEX AU(3,3),AD(3,3),AE(3,3)
      DOUBLE COMPLEX HU(3,3),HD(3,3),HE(3,3),DH(3,3,3)
      DOUBLE COMPLEX MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE COMPLEX MQM(3,3),MLM(3,3),MUPM(3,3),MDM(3,3),MEM(3,3)
      DOUBLE COMPLEX DM(5,3,3)
      DOUBLE COMPLEX LU(3,3),LD(3,3),LE(3,3)
      DOUBLE COMPLEX LYU(3,3),LYD(3,3),LYE(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3),GTPE(3,3)
      DOUBLE COMPLEX GTQ(3,3),GTL(3,3)
      DOUBLE COMPLEX GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE COMPLEX FTUQ(3,3),FTDQ(3,3),FTEL(3,3)
      DOUBLE COMPLEX FTUU(3,3),FTDD(3,3),FTEE(3,3)
!
      DOUBLE COMPLEX TRIU(3,3),TRID(3,3),TRIE(3,3)
      DOUBLE COMPLEX MTSFU(3,3),MTSFD(3,3),MTSFE(3,3)
      DOUBLE COMPLEX FUHU(3,3),FDHD(3,3),FEHD(3,3)
      DOUBLE COMPLEX FUQ(3,3),FUUR(3,3)
      DOUBLE COMPLEX FDQ(3,3),FDDR(3,3)
      DOUBLE COMPLEX FEL(3,3),FEER(3,3)
      DOUBLE COMPLEX SFUQ(3,3),CFDQ(3,3),SFUUR(3,3)
      DOUBLE COMPLEX CFDDR(3,3),CFEL(3,3),CFEER(3,3)
!
      DOUBLE COMPLEX CMATMUL,CSFMUL,CTRACE,TCSFMUL,SUM
      DOUBLE COMPLEX CMODSQ,CCON,CRE
!
!These are used in the calculation of the RGEs which contain
!thresholds
!
      DOUBLE COMPLEX DUMU1(3,3),DUMU2(3,3),DUMD1(3,3),DUMD2(3,3)
      DOUBLE COMPLEX DUME1(3,3),DUME2(3,3)
      DOUBLE COMPLEX DUMM1(3,3),TDUMM1,DUMM2(3,3),TDUMM2,DUMM3(3,3)
      DOUBLE COMPLEX TDUMM3,DUMGRKMU1(3,3),TDUMGRKMU
      DOUBLE COMPLEX DUMTUQ1(3,3),DUMTUQ2(3,3),DUMTUQ3(3,3),TDUMTUQ2
      DOUBLE COMPLEX DUMTUQ4(3,3),DUMTUQ5(3,3),DUMTDQ1(3,3),DUMTDQ2(3,3)
      DOUBLE COMPLEX TDUMTDQ2,DUMTDQ3(3,3),DUMTDQ4(3,3),DUMTDQ5(3,3)
      DOUBLE COMPLEX DUMTEL1(3,3),TDUMTEL2,DUMTEL3(3,3),DUMTEL4(3,3)
      DOUBLE COMPLEX DUMTEL5(3,3),TDUMTUU1,DUMTUU2(3,3),DUMTUU3(3,3)
      DOUBLE COMPLEX DUMTUU4(3,3),TDUMTDD1,DUMTDD2(3,3),DUMTDD3(3,3)
      DOUBLE COMPLEX DUMTDD4(3,3),TDUMTEE1,DUMTEE3(3,3)
      DOUBLE COMPLEX DUMTQ1(3,3),TDUMTQ1,DUMTQ2(3,3),DUMTQ3(3,3)
      DOUBLE COMPLEX DUMTQ4(3,3),DUMTQ5(3,3),TDUMTL1,DUMTL2(3,3)
      DOUBLE COMPLEX DUMTL3(3,3),DUMTPQ1(3,3),TDUMTPQ1,DUMTPQ2(3,3)
      DOUBLE COMPLEX DUMTPQ3(3,3),DUMTPQ4(3,3),DUMTPQ5(3,3),TDUMTPL1
      DOUBLE COMPLEX DUMTPL2(3,3),DUMTPL3(3,3),DUMTPU1(3,3),TDUMTPU2
      DOUBLE COMPLEX DUMTPU3(3,3),DUMTPD1(3,3),TDUMTPD2,DUMTPD3(3,3)
      DOUBLE COMPLEX DUMTPE1(3,3),TDUMTPE2
      DOUBLE COMPLEX DUMTSQ1(3,3),TDUMTSQ1,DUMTSQ2(3,3),DUMTSQ3(3,3)
      DOUBLE COMPLEX DUMTSQ4(3,3),DUMTSQ5(3,3),DUMTSU1(3,3),TDUMTSU2
      DOUBLE COMPLEX DUMTSU3(3,3),TDUMTSU3,DUMTSU4(3,3),DUMTSD1(3,3)
      DOUBLE COMPLEX TDUMTSD2,TDUMTSD3,DUMTSD4(3,3)
      DOUBLE COMPLEX DUMGTHU1(3,3),TDUMGTHU1,DUMGTHD1(3,3),TDUMGTHD1
      DOUBLE COMPLEX DUMGTPHU1(3,3),TDUMGTPHU1,DUMGTPHU2(3,3),TDUMGTPHU2
      DOUBLE COMPLEX DUMGTPHD1(3,3),TDUMGTPHD1,DUMGTPHD2(3,3),TDUMGTPHD2
      DOUBLE COMPLEX DUMLUD1(3,3),DUMLUD2(3,3),DUMLUD(3,3),TDUMLUD
      DOUBLE COMPLEX DUMMTSFU1(3,3),DUMMTSFU2(3,3),DUMMTSFU3(3,3)
      DOUBLE COMPLEX DUMMTSFD1(3,3),DUMMTSFD2(3,3),DUMMTSFD3(3,3)
      DOUBLE COMPLEX DUMMTSFE1(3,3),DUMMTSFE2(3,3),DUMMTSFE3(3,3)
!
      DOUBLE COMPLEX FUS(3,3),FDS(3,3),FES(3,3)
      DOUBLE COMPLEX FUT(3,3),FDT(3,3),FET(3,3)
      DOUBLE COMPLEX LUS(3,3),LDS(3,3),LES(3,3)
      DOUBLE COMPLEX LUT(3,3),LDT(3,3),LET(3,3)
      DOUBLE COMPLEX AUS(3,3),ADS(3,3),AES(3,3)
      DOUBLE COMPLEX AUT(3,3),ADT(3,3),AET(3,3)
      DOUBLE COMPLEX FTUQS(3,3),FTUQT(3,3),FTUQD(3,3),FTUUD(3,3)
      DOUBLE COMPLEX FTUUS(3,3),FTUUT(3,3),FTDQS(3,3),FTDQT(3,3)
      DOUBLE COMPLEX FTDDD(3,3),FTDDS(3,3),FTDDT(3,3),FTDQD(3,3)
      DOUBLE COMPLEX FTELD(3,3),FTELS(3,3),FTELT(3,3),FTEED(3,3)
      DOUBLE COMPLEX FTEES(3,3),FTEET(3,3),GTQD(3,3),GTLD(3,3)
      DOUBLE COMPLEX GTPQD(3,3),GTSQD(3,3),GTPUD(3,3),GTPDD(3,3)
      DOUBLE COMPLEX GTPED(3,3),GTPLD(3,3),GTSUD(3,3),GTSDD(3,3)
      DOUBLE COMPLEX GTPUS(3,3),GTSUS(3,3),GTQS(3,3),GTLS(3,3)
      DOUBLE COMPLEX GTPQS(3,3),GTPLS(3,3),GTSQS(3,3),GTPDS(3,3)
      DOUBLE COMPLEX GTPES(3,3),GTSDS(3,3),GTPUT(3,3),GTQT(3,3)
      DOUBLE COMPLEX GTLT(3,3),GTPQT(3,3),GTPLT(3,3),GTPDT(3,3)
      DOUBLE COMPLEX GTPET(3,3),GTSUT(3,3),GTSQT(3,3),GTSDT(3,3)
      DOUBLE COMPLEX FUHUD(3,3),FDHDD(3,3),FEHDD(3,3)
      DOUBLE COMPLEX MTSFUS(3,3),MTSFDS(3,3),MTSFES(3,3)
      DOUBLE COMPLEX MTSFUT(3,3),MTSFDT(3,3),MTSFET(3,3)
      DOUBLE COMPLEX FUHUS(3,3),FDHDS(3,3),FEHDS(3,3)
      DOUBLE COMPLEX TRIUS(3,3),TRIDS(3,3),TRIES(3,3)
      DOUBLE COMPLEX FUHUT(3,3),FDHDT(3,3),FEHDT(3,3)
      DOUBLE COMPLEX TRIUT(3,3),TRIDT(3,3),TRIET(3,3)
!
      DOUBLE COMPLEX FUFUD(3,3),FDFDD(3,3),FEFED(3,3)
      DOUBLE COMPLEX YUDYU(3,3),YDDYD(3,3),YEDYE(3,3)
      DOUBLE COMPLEX LULUD(3,3),LDLDD(3,3),LELED(3,3)
      DOUBLE COMPLEX LUDLU(3,3),LDDLD(3,3),LEDLE(3,3)
      DOUBLE COMPLEX LYUDLYU(3,3),LYDDLYD(3,3),LYEDLYE(3,3)
      DOUBLE COMPLEX LYUDLYU2(3,3),LYDDLYD2(3,3),LYEDLYE2(3,3)
      DOUBLE COMPLEX LYUDLYU3(3,3),LYDDLYD3(3,3),LYEDLYE3(3,3)
      DOUBLE COMPLEX TLYUDLYU,TLYDDLYD,TLYEDLYE
      DOUBLE COMPLEX TLYUDLYU2,TLYDDLYD2,TLYEDLYE2
      DOUBLE COMPLEX TLYUDLYU3,TLYDDLYD3,TLYEDLYE3
      DOUBLE COMPLEX LYUDLYULYDDLYD(3,3),LYDDLYDLYUDLYU(3,3)
      DOUBLE COMPLEX TLYUDLYULYDDLYD
      DOUBLE COMPLEX FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
      DOUBLE COMPLEX YUYUD(3,3),YDYDD(3,3),YEYED(3,3)
      DOUBLE COMPLEX YUYDD(3,3),YDYUD(3,3)
      DOUBLE COMPLEX FUGTPUS(3,3),FDGTPDS(3,3),FEGTPES(3,3)
      DOUBLE COMPLEX FUGTSUS(3,3),FDGTSDS(3,3)
      DOUBLE COMPLEX LUGTPUS(3,3),LDGTPDS(3,3),LEGTPES(3,3)
      DOUBLE COMPLEX LUGTSUS(3,3),LDGTSDS(3,3)
      DOUBLE COMPLEX FTUQFTUQD(3,3),FTUQSFTUQT(3,3),FTUQSFUT(3,3)
      DOUBLE COMPLEX FTUQSLUT(3,3)
      DOUBLE COMPLEX FTDQFTDQD(3,3),FTDQSFTDQT(3,3),FTDQSFDT(3,3)
      DOUBLE COMPLEX FTDQSLDT(3,3)
      DOUBLE COMPLEX FTELFTELD(3,3),FTELSFTELT(3,3),FTELSFET(3,3)
      DOUBLE COMPLEX FTELSLET(3,3)
      DOUBLE COMPLEX GTQSGTQT(3,3),GTQGTQD(3,3),GTPQSGTPQT(3,3)
      DOUBLE COMPLEX GTPQGTPQD(3,3),GTPUTGTPUS(3,3),GTPUDGTPU(3,3)
      DOUBLE COMPLEX GTPDDGTPD(3,3),GTPDTGTPDS(3,3),GTSQSGTSQT(3,3)
      DOUBLE COMPLEX GTSQGTSQD(3,3),GTSUTGTSUS(3,3),GTSUDGTSU(3,3)
      DOUBLE COMPLEX GTSDDGTSD(3,3),GTSDTGTSDS(3,3),GTLSGTLT(3,3)
      DOUBLE COMPLEX GTLGTLD(3,3),GTPLSGTPLT(3,3),GTPLGTPLD(3,3)
      DOUBLE COMPLEX GTPETGTPES(3,3),GTPEDGTPE(3,3)
      DOUBLE COMPLEX FDDRFDDRD(3,3),FDQDFUQ(3,3),FTUQGTPUS(3,3)
      DOUBLE COMPLEX FTDQGTPDS(3,3),FTELGTPES(3,3),FTUQFUDFTDD(3,3)
      DOUBLE COMPLEX FTDQFDDFTUU(3,3),FTDQLDDFTUU(3,3),FTUQLUDFTDD(3,3)
      DOUBLE COMPLEX FUURFUURD(3,3),FUQDFDQ(3,3)
      DOUBLE COMPLEX FUQDFUQ(3,3),FDQDFDQ(3,3)
      DOUBLE COMPLEX FUQTFUQS(3,3),FUURSFUURT(3,3)
      DOUBLE COMPLEX FDQTFDQS(3,3),FDDRSFDDRT(3,3),FELTFELS(3,3)
      DOUBLE COMPLEX FEERSFEERT(3,3)
      DOUBLE COMPLEX GTPQSFTUU(3,3),GTQSFTUU(3,3),FDDFTUU(3,3)
      DOUBLE COMPLEX LUDFTDD(3,3),LDDFTUU(3,3)
      DOUBLE COMPLEX FUDFTDD(3,3),FEERFEERD(3,3),FELDFEL(3,3)
      DOUBLE COMPLEX FTUUDFTUU(3,3),FTDDDFTDD(3,3),FTEEDFTEE(3,3)
      DOUBLE COMPLEX GTPQSFTDD(3,3),GTQSFTDD(3,3)
      DOUBLE COMPLEX GTPLSFTEE(3,3),GTLSFTEE(3,3)
      DOUBLE COMPLEX FUSFUT(3,3),FDSFDT(3,3),FESFET(3,3)
      DOUBLE COMPLEX LUSLUT(3,3),LDSLDT(3,3),LESLET(3,3)
      DOUBLE COMPLEX FUTFUS(3,3),FDTFDS(3,3),FETFES(3,3)
      DOUBLE COMPLEX LUTLUS(3,3),LDTLDS(3,3),LETLES(3,3)
      DOUBLE COMPLEX FUTFTUUS(3,3),FDTFTDDS(3,3),FETFTEES(3,3)
      DOUBLE COMPLEX LUTFTUUS(3,3),LDTFTDDS(3,3),LETFTEES(3,3)
      DOUBLE COMPLEX TFUDFU,TFDDFD,TFEDFE
      DOUBLE COMPLEX TYUDYU,TYDDYD,TYEDYE
      DOUBLE COMPLEX TLUDLU,TLDDLD,TLEDLE
      DOUBLE COMPLEX FUFUDFU(3,3),FUFUDFD(3,3),FUTFUSGTPU(3,3)
      DOUBLE COMPLEX LUTLUSGTPU(3,3)
      DOUBLE COMPLEX FUTFUSGTSU(3,3),LUTLUSGTSU(3,3)
      DOUBLE COMPLEX FDTFDSGTSD(3,3),LDTLDSGTSD(3,3)
      DOUBLE COMPLEX FDFDDFU(3,3)
      DOUBLE COMPLEX FDTFDSGTPD(3,3),LDTLDSGTPD(3,3)
      DOUBLE COMPLEX FDFDDFD(3,3),FEFEDFE(3,3)
      DOUBLE COMPLEX FETFESGTPE(3,3),LETLESGTPE(3,3)
      DOUBLE COMPLEX LULUDLU(3,3),LULUDLD(3,3),LDLDDLU(3,3)
      DOUBLE COMPLEX LDLDDLD(3,3),LELEDLE(3,3)
      DOUBLE COMPLEX FTUQFUDFU(3,3),FTDQFDDFD(3,3),FTELFEDFE(3,3)
      DOUBLE COMPLEX FTUQLUDLU(3,3),FTDQLDDLD(3,3),FTELLEDLE(3,3)
      DOUBLE COMPLEX SQFTUQDFTUQ(3,3),SQFTUQTFTUQS(3,3),SQFTUQTGTQ(3,3)
      DOUBLE COMPLEX SQFTUQTGTPQ(3,3),SQFTUQTGTSQ(3,3),SQFTDQDFTDQ(3,3)
      DOUBLE COMPLEX SQFTDQTFTDQS(3,3),SQFTDQTGTPQ(3,3),SQFTDQTGTQ(3,3)
      DOUBLE COMPLEX SQFTDQTGTSQ(3,3),SQGTQTGTQS(3,3),SQGTQTFTUQ(3,3)
      DOUBLE COMPLEX SQGTQTFTDQ(3,3),SQGTQDGTQ(3,3),SQGTPQTGTPQS(3,3)
      DOUBLE COMPLEX SQGTPQTFTUQ(3,3),SQGTPQTFTDQ(3,3),SQGTPQDGTPQ(3,3)
      DOUBLE COMPLEX SQGTSQTGTSQS(3,3),SQGTSQTFTUQ(3,3),SQGTSQTFTDQ(3,3)
      DOUBLE COMPLEX SQGTSQDGTSQ(3,3),SQIDMTSFU(3,3),SQIDMTSFD(3,3)
      DOUBLE COMPLEX SQFUHUDMTSFU(3,3),SQFDHDDMTSFD(3,3)
      DOUBLE COMPLEX SUGTPUFTUUT(3,3),SUGTPUGTPUD(3,3)
      DOUBLE COMPLEX SUGTPUSGTPUT(3,3),SUGTSUGTSUD(3,3)
      DOUBLE COMPLEX SUGTSUSGTSUT(3,3),SUGTSUFTUUT(3,3),SUFTUUFTUUD(3,3)
      DOUBLE COMPLEX SUFTUUGTPUT(3,3),SUFTUUSFTUUT(3,3),SUFTUUGTSUT(3,3)
      DOUBLE COMPLEX SUMTSFUID(3,3),SUSQIDMTSFUID(3,3)
      DOUBLE COMPLEX SUMTSFUFUQDFDQ(3,3)
      DOUBLE COMPLEX SDFTDDFTDDD(3,3),SDFTDDSFTDDT(3,3),SDFTDDGTPDT(3,3)
      DOUBLE COMPLEX SDFTDDGTSDT(3,3),SDGTPDFTDDT(3,3),SDGTPDSGTPDT(3,3)
      DOUBLE COMPLEX SDGTPDGTPDD(3,3),SDGTSDGTSDD(3,3),SDGTSDSGTSDT(3,3)
      DOUBLE COMPLEX SDGTSDFTDDT(3,3),SDMTSFDFDQDFUQ(3,3)
      DOUBLE COMPLEX SDMTSFDID(3,3),SDSQIDMTSFDID(3,3)
      DOUBLE COMPLEX SLFTELDFTEL(3,3),SLFTELTFTELS(3,3)
      DOUBLE COMPLEX SLFTELTGTL(3,3),SLFTELTGTPL(3,3),SLGTLTGTLS(3,3)
      DOUBLE COMPLEX SLGTLTFTEL(3,3),SLGTLDGTL(3,3),SLGTPLTGTPLS(3,3)
      DOUBLE COMPLEX SLGTPLDGTPL(3,3),SLGTPLTFTEL(3,3),SLFEHDDMTSFE(3,3)
      DOUBLE COMPLEX SLIDMTSFE(3,3)
      DOUBLE COMPLEX SEGTPESGTPET(3,3)
      DOUBLE COMPLEX SEGTPEGTPED(3,3),SEGTPEFTEET(3,3),SEFTEEFTEED(3,3)
      DOUBLE COMPLEX SEFTEESFTEET(3,3),SEFTEEGTPET(3,3),SEMTSFEID(3,3)
      DOUBLE COMPLEX SESLIDMTSFEID(3,3)
      DOUBLE COMPLEX TFUSFUT,TFDSFDT,TFESFET,TLUSLUT,TLDSLDT,TLESLET
      DOUBLE COMPLEX TSQGTQDGTQ,TSLGTLDGTL,TSUFTUUFTUUD,TSDFTDDFTDDD
      DOUBLE COMPLEX TSEFTEEFTEED,TSQFTUQDFTUQ,TSQFTDQDFTDQ
      DOUBLE COMPLEX TSLFTELDFTEL,TSQGTPQDGTPQ,TSLGTPLDGTPL
      DOUBLE COMPLEX TSDGTPDGTPDD,TSUGTPUGTPUD,TSEGTPEGTPED
      DOUBLE COMPLEX SUTSQFUHUDMTSFU,SDTSQFDHDDMTSFD
      DOUBLE COMPLEX SETSLFEHDDMTSFE
      DOUBLE COMPLEX SQTSUMTSFUSAUT,SUMTSFUSAUT(3,3)
      DOUBLE COMPLEX SQTSDMTSFDSADT,SDMTSFDSADT(3,3)
      DOUBLE COMPLEX SLTSEMTSFESAET,SEMTSFESAET(3,3)
      DOUBLE COMPLEX MGTPHUSQ,MGTPHDSQ,MGTHUSQ,MGTHDSQ
      DOUBLE COMPLEX MSGTPHUSQ,MCGTPHDSQ,MSGTHUSQ,MCGTHDSQ
      DOUBLE COMPLEX MMUSQ,M1PM1PSQ,M2PM2PSQ,M3PM3PSQ
      DOUBLE COMPLEX Y2,H,Y4,CHI4
!
      DOUBLE COMPLEX SUAUID(3,3),SUAUFUQDFUQ(3,3),SUSQIDAUID(3,3)
      DOUBLE COMPLEX SQIDAU(3,3),SUTSQFUHUDAU,SQFUHUDAU(3,3)
      DOUBLE COMPLEX SQFUURFUURDAU(3,3),GTPQSFUGTPUS(3,3)
      DOUBLE COMPLEX SDADFDQDFUQ(3,3),GTSQSFUGTSUS(3,3)
      DOUBLE COMPLEX SQIDAD(3,3),SQFDDRFDDRDAD(3,3),SDSQIDADID(3,3)
      DOUBLE COMPLEX SDTSQFDHDDAD,SQFDHDDAD(3,3)
      DOUBLE COMPLEX SETSLFEHDDAE,SLFEHDDAE(3,3),SDADID(3,3)
      DOUBLE COMPLEX SDADFDQDFDQ(3,3),GTPQSFDGTPDS(3,3)
      DOUBLE COMPLEX GTSQSFDGTSDS(3,3),SUAUFUQDFDQ(3,3)
      DOUBLE COMPLEX SLIDAE(3,3),SLFEERFEERDAE(3,3),SEAEID(3,3)
      DOUBLE COMPLEX SESLIDAEID(3,3),SEAEFELDFEL(3,3)
      DOUBLE COMPLEX GTPLSFEGTPES(3,3)
!
      DOUBLE COMPLEX SFUQDSFUQ(3,3),SFUURSFUURD(3,3),CFDDRCFDDRD(3,3)
      DOUBLE COMPLEX CFDQDSFUQ(3,3),CFDQDCFDQ(3,3)
      DOUBLE COMPLEX SFUQDCFDQ(3,3),CFEERCFEERD(3,3)
      DOUBLE COMPLEX CFELDCFEL(3,3)
      DOUBLE COMPLEX SUTRIUID(3,3),SUTRIUSFUQDSFUQ(3,3),SQIDTRIU(3,3)
      DOUBLE COMPLEX SUSQIDTRIUID(3,3),SUTSQFUHUDTRIU,SQFUHUDTRIU(3,3)
      DOUBLE COMPLEX SQSFUURSFUURDTRIU(3,3),SQCFDDRCFDDRDTRIU(3,3)
      DOUBLE COMPLEX GTPQSLUGTPUS(3,3),GTSQSLUGTSUS(3,3)
      DOUBLE COMPLEX SQIDTRID(3,3),SQSFUURSFUURDTRID(3,3)
      DOUBLE COMPLEX SQCFDDRCFDDRDTRID(3,3),SDSQIDTRIDID(3,3)
      DOUBLE COMPLEX SDTSQFDHDDTRID,SQFDHDDTRID(3,3)
      DOUBLE COMPLEX SETSLFEHDDTRIE,SLFEHDDTRIE(3,3),SDTRIDID(3,3)
      DOUBLE COMPLEX SDTRIDCFDQDCFDQ(3,3)
      DOUBLE COMPLEX GTPQSLDGTPDS(3,3),GTSQSLDGTSDS(3,3)
      DOUBLE COMPLEX SLIDTRIE(3,3),SLCFEERCFEERDTRIE(3,3)
      DOUBLE COMPLEX SESLIDTRIEID(3,3),SETRIEID(3,3)
      DOUBLE COMPLEX SETRIECFELDCFEL(3,3),GTPLSLEGTPES(3,3)
!
      DOUBLE COMPLEX SUTSUIDMU,SUIDMU(3,3),SUTSUFUQTFUQSMU
      DOUBLE COMPLEX SUFUQTFUQSMU(3,3),SQTSQIDMQ,SQIDMQ(3,3)
      DOUBLE COMPLEX SQTSQFUURSFUURTMQ,SQFUURSFUURTMQ(3,3)
      DOUBLE COMPLEX SDTSDIDMD,SDIDMD(3,3),SLTSLIDML,SLIDML(3,3)
      DOUBLE COMPLEX SETSEIDME,SEIDME(3,3),SQTSUAUSAUT,SUAUSAUT(3,3)
      DOUBLE COMPLEX SQTSDMTSFDSMTSFDT,SDMTSFDSMTSFDT(3,3)
      DOUBLE COMPLEX SLTSEMTSFESMTSFET,SEMTSFESMTSFET(3,3)
      DOUBLE COMPLEX SQTSQFDDRSFDDRTMQ,SQFDDRSFDDRTMQ(3,3)
      DOUBLE COMPLEX SDTSDFDQTFDQSMD,SDFDQTFDQSMD(3,3)
      DOUBLE COMPLEX SLTSLFEERSFEERTML,SLFEERSFEERTML(3,3)
      DOUBLE COMPLEX SETSEFELTFELSME,SEFELTFELSME(3,3)
      DOUBLE COMPLEX SQTSUMTSFUSMTSFUT,SUMTSFUSMTSFUT(3,3)
      DOUBLE COMPLEX SQTSDADSADT,SDADSADT(3,3),SLTSEAESAET,SEAESAET(3,3)
!
      DOUBLE COMPLEX SUTSUSFUQTSFUQSMU,SUSFUQTSFUQSMU(3,3)
      DOUBLE COMPLEX SFUQTSFUQS(3,3),SQTSQSFUURSSFUURTMQ
      DOUBLE COMPLEX SQSFUURSSFUURTMQ(3,3),SFUURSSFUURT(3,3)
      DOUBLE COMPLEX SQTSQCFDDRSCFDDRTMQ,SQCFDDRSCFDDRTMQ(3,3)
      DOUBLE COMPLEX CFDDRSCFDDRT(3,3),SDTSDCFDQTCFDQSMD
      DOUBLE COMPLEX SDCFDQTCFDQSMD(3,3),CFDQTCFDQS(3,3)
      DOUBLE COMPLEX SLTSLCFEERSCFEERTML,SLCFEERSCFEERTML(3,3)
      DOUBLE COMPLEX CFEERSCFEERT(3,3),SETSECFELTCFELSME
      DOUBLE COMPLEX SECFELTCFELSME(3,3),CFELTCFELS(3,3)
      DOUBLE COMPLEX SQTSUTRIUTRIUD,SUTRIUTRIUD(3,3),TRIUD(3,3)
      DOUBLE COMPLEX SQTSDTRIDTRIDD,SDTRIDTRIDD(3,3),TRIDD(3,3)
      DOUBLE COMPLEX SLTSETRIETRIED,SETRIETRIED(3,3),TRIED(3,3)
!
      DOUBLE COMPLEX SUTMU,SQTMQ,SDTMD,SLTML,SETME
      DOUBLE COMPLEX SQSQIDMQID(3,3)
      DOUBLE COMPLEX SUSUFUHUSMUFUHUT(3,3),SUFUHUSMU(3,3)
      DOUBLE COMPLEX SDSDFDHDSMDFDHDT(3,3),SDFDHDSMD(3,3)
      DOUBLE COMPLEX SUTRIUSTRIUT(3,3),SDTRIDSTRIDT(3,3)
      DOUBLE COMPLEX SUSUIDMUID(3,3),SQSQFUHUTMQFUHUS(3,3)
      DOUBLE COMPLEX SQFUHUTMQ(3,3),SQAUTAUS(3,3),SQMTSFUTMTSFUS(3,3)
      DOUBLE COMPLEX SQTRIUTTRIUS(3,3),FTUUTFTUUS(3,3)
      DOUBLE COMPLEX SDSDIDMDID(3,3)
      DOUBLE COMPLEX SQSQFDHDTMQFDHDS(3,3),SQFDHDTMQ(3,3),SQADTADS(3,3)
      DOUBLE COMPLEX SQMTSFDTMTSFDS(3,3),FTDDTFTDDS(3,3)
      DOUBLE COMPLEX SQTRIDTTRIDS(3,3),SLSLIDMLID(3,3)
      DOUBLE COMPLEX SEFEHDSME(3,3),SESEFEHDSMEFEHDT(3,3)
      DOUBLE COMPLEX SETRIESTRIET(3,3)
      DOUBLE COMPLEX SESEIDMEID(3,3),SLSLFEHDTMLFEHDS(3,3)
      DOUBLE COMPLEX SLFEHDTML(3,3)
      DOUBLE COMPLEX SLAETAES(3,3),SLMTSFETMTSFES(3,3),FTEETFTEES(3,3)
      DOUBLE COMPLEX SLTRIETTRIES(3,3)
!
      DOUBLE COMPLEX B1U(3,3),B1D(3,3),B1E(3,3)
      DOUBLE COMPLEX BFTUQ(3,3),BFTDQ(3,3),BFTEL(3,3)
      DOUBLE COMPLEX BFTUU(3,3),BFTDD(3,3),BFTEE(3,3)
      DOUBLE COMPLEX BGTQ(3,3),BGTL(3,3),BGTPQ(3,3),BGTPL(3,3)
      DOUBLE COMPLEX BGTPU(3,3),BGTPD(3,3),BGTPE(3,3)
      DOUBLE COMPLEX BGTSQ(3,3),BGTSU(3,3),BGTSD(3,3)
      DOUBLE COMPLEX BGRKMU,BM(3),BMP(3),BGTHU,BGTHD,BGTPHU,BGTPHD
      DOUBLE COMPLEX BSGTHU,BCGTHD,BSGTPHU,BCGTPHD,BETALAM1,BETALAM2
      DOUBLE COMPLEX BMTSFU(3,3),BMTSFD(3,3),BMTSFE(3,3)
      DOUBLE COMPLEX BMHUPMT,BMHDPMT,BMHUD
!
!The following are for two loop and soft RGEs
!
      DOUBLE COMPLEX S,SP,SIG1,SIG2,SIG3
      DOUBLE COMPLEX DUM2U11(3,3),DUM2U12(3,3),DUM2U21(3,3)
      DOUBLE COMPLEX DUM2U22(3,3),DUM2U23(3,3),DUM1D11(3,3)
      DOUBLE COMPLEX DUM1D12(3,3),DUM2D11(3,3),DUM2D12(3,3)
      DOUBLE COMPLEX DUM2D21(3,3),DUM2D22(3,3),DUM2D23(3,3)
      DOUBLE COMPLEX DUM1E11(3,3),DUM1E21(3,3),DUM2E11(3,3)
      DOUBLE COMPLEX DUM2E12(3,3),DUM2E21(3,3),DUM2E22(3,3)
      DOUBLE COMPLEX DUM2E23(3,3)
      DOUBLE COMPLEX DUM1U1(3,3),DUM1U2(3,3),DUM2U1(3,3),DUM2U2(3,3)
      DOUBLE COMPLEX DUM1D1(3,3),DUM1D2(3,3),DUM2D1(3,3),DUM2D2(3,3)
      DOUBLE COMPLEX DUM1E1(3,3),DUM1E2(3,3),DUM2E1(3,3),DUM2E2(3,3)
      DOUBLE COMPLEX DUM1U(3,3),DUM2U(3,3),DUM1D(3,3),DUM2D(3,3)
      DOUBLE COMPLEX DUM1E(3,3),DUM2E(3,3)
      DOUBLE COMPLEX DUMS(3,3),DUMSP1(3,3),DUMSP2(3,3),DUMSP3(3,3)
      DOUBLE COMPLEX DUMSIG1(3,3),DUMSIG2(3,3),DUMSIG3(3,3)
      DOUBLE COMPLEX DUM1HU1(3,3)
      DOUBLE COMPLEX DUMSP(3,3)
      DOUBLE COMPLEX DUM2HU11(3,3),DUM2HU12(3,3),DUM2HU21(3,3)
      DOUBLE COMPLEX DUM2HD11(3,3)
      DOUBLE COMPLEX DUM2HD12(3,3),DUM2HD13(3,3),DUM2HD21(3,3)
      DOUBLE COMPLEX DUM2HD31(3,3)
      DOUBLE COMPLEX DUM2HU1(3,3),DUM2HU2(3,3),DUM1HD1(3,3)
      DOUBLE COMPLEX DUM1HD2(3,3)
      DOUBLE COMPLEX DUM2HD1(3,3),DUM2HD2(3,3),DUM2HD3(3,3)
      DOUBLE COMPLEX DUM1HU(3,3),DUM1HD(3,3)
      DOUBLE COMPLEX DUM2Q31(3,3),DUM2Q41(3,3),DUM2L21(3,3)
      DOUBLE COMPLEX DUM2U31(3,3),DUM2D31(3,3)
      DOUBLE COMPLEX DUM1Q1(3,3),DUM1Q2(3,3),DUM1Q3(3,3),DUM2Q1(3,3)
      DOUBLE COMPLEX DUM2Q2(3,3),DUM2Q3(3,3),DUM2Q4(3,3),DUM2Q5(3,3)
      DOUBLE COMPLEX DUM2Q6(3,3),DUM2Q7(3,3),DUM2Q8(3,3),DUM2Q9(3,3)
      DOUBLE COMPLEX DUM2QA(3,3),DUM2QB(3,3),DUM2QC(3,3),DUM2QD(3,3)
      DOUBLE COMPLEX DUM1L1(3,3),DUM2L1(3,3),DUM2L2(3,3),DUM2L3(3,3)
      DOUBLE COMPLEX DUM2L4(3,3),DUM2L5(3,3),DUM2L6(3,3),DUM2L7(3,3)
      DOUBLE COMPLEX DUM2L8(3,3),DUM2L9(3,3),DUM2U3(3,3),DUM2U4(3,3)
      DOUBLE COMPLEX DUM2U5(3,3),DUM2D3(3,3),DUM2D4(3,3),DUM2D5(3,3)
      DOUBLE COMPLEX DUM2D6(3,3),DUM2D7(3,3),DUM2D8(3,3),DUM2D9(3,3)
      DOUBLE COMPLEX DUM2DA(3,3),DUM2E3(3,3),DUM2E4(3,3),DUM2E5(3,3)
      DOUBLE COMPLEX DUM2E6(3,3),DUM2E7(3,3),DUM2E8(3,3),DUM2E9(3,3)
      DOUBLE COMPLEX DUM1GRKMU(3,3),DUM2GRKMU(3,3),DUM1B1(3,3)
      DOUBLE COMPLEX DUM1B2(3,3)
      DOUBLE COMPLEX DUM2B1(3,3),DUM2B2(3,3)
!
      DOUBLE COMPLEX TYUDHU,TYDDHD,TYEDHE,THUYUD
      DOUBLE COMPLEX THDYDD,THEYED,THUDYU,THDDYD,THEDYE,THUDHU
      DOUBLE COMPLEX YUYUDYUYUD(3,3),YDYDDYDYDD(3,3)
      DOUBLE COMPLEX YEYEDYEYED(3,3),YUDYUYUDYU(3,3)
      DOUBLE COMPLEX YDDYDYDDYD(3,3),YEDYEYEDYE(3,3)
      DOUBLE COMPLEX YDDYDYUDYU(3,3),YUDYUYDDYD(3,3)
      DOUBLE COMPLEX YUYDDYDYUD(3,3),YDYUDYUYDD(3,3)
      DOUBLE COMPLEX HUYUD(3,3),HDYDD(3,3),HEYED(3,3)
      DOUBLE COMPLEX YUDHU(3,3),YDDHD(3,3),YEDHE(3,3)
      DOUBLE COMPLEX HUDHU(3,3),HDDHD(3,3),HEDHE(3,3)
      DOUBLE COMPLEX HUDYU(3,3),HDDYD(3,3),HEDYE(3,3)
      DOUBLE COMPLEX HUHUD(3,3),HDHDD(3,3),HEHED(3,3)
      DOUBLE COMPLEX YUHUD(3,3),YDHDD(3,3),YEHED(3,3)
      DOUBLE COMPLEX HUYDD(3,3),HDYUD(3,3),HUHDD(3,3)
      DOUBLE COMPLEX HDHUD(3,3),YDHUD(3,3),YUHDD(3,3)
      DOUBLE COMPLEX MUPMYU(3,3),MDMYD(3,3),MEMYE(3,3)
      DOUBLE COMPLEX YUMQM(3,3),YDMQM(3,3),YEMLM(3,3)
      DOUBLE COMPLEX YUDMUPMYU(3,3),YDDMDMYD(3,3),YEDMEMYE(3,3)
      DOUBLE COMPLEX YEDYEMLM(3,3),YUMQMYUD(3,3),YUYUDMUPM(3,3)
      DOUBLE COMPLEX YDMQMYDD(3,3),YDYDDMDM(3,3),YEMLMYED(3,3)
      DOUBLE COMPLEX YEYEDMEM(3,3),YUDYUMQM(3,3),YDDYDMQM(3,3)
      DOUBLE COMPLEX MQMYUDYU(3,3),MQMYDDYD(3,3),MLMYEDYE(3,3)
      DOUBLE COMPLEX YUDYUYDD(3,3),YUMQMYDD(3,3),YUYDDMDM(3,3)
      DOUBLE COMPLEX YDMQMYUD(3,3),YDYUDMUPM(3,3)
      DOUBLE COMPLEX HUYUDYUYUD(3,3),HUYDDYDYUD(3,3),HDYUDYUYDD(3,3)
      DOUBLE COMPLEX HDYDDYDYDD(3,3),HEYEDYEYED(3,3),YUDYUYUDHU(3,3)
      DOUBLE COMPLEX YUDHUYUDYU(3,3),YDDYDYDDHD(3,3),YDDHDYDDYD(3,3)
      DOUBLE COMPLEX YDDYDYUDHU(3,3),YDDHDYUDYU(3,3),YUDHUYDDYD(3,3)
      DOUBLE COMPLEX YUDYUYDDHD(3,3),YEDYEYEDHE(3,3),YEDHEYEDYE(3,3)
      DOUBLE COMPLEX HUDHUYUDYU(3,3),HUDYUYUDHU(3,3),HDDHDYDDYD(3,3)
      DOUBLE COMPLEX HDDYDYDDHD(3,3),HDDHDYUDYU(3,3),YDDYDHUDHU(3,3)
      DOUBLE COMPLEX HDDYDYUDHU(3,3),YDDHDHUDYU(3,3),HUDHUYDDYD(3,3)
      DOUBLE COMPLEX YUDYUHDDHD(3,3),HUDYUYDDHD(3,3),YUDHUHDDYD(3,3)
      DOUBLE COMPLEX HEDHEYEDYE(3,3),HEDYEYEDHE(3,3),YUDYUHUDHU(3,3)
      DOUBLE COMPLEX YDDYDHDDHD(3,3),YUDHUHUDYU(3,3),YDDHDHDDYD(3,3)
      DOUBLE COMPLEX YEDYEHEDHE(3,3),YEDHEHEDYE(3,3),HUHUDYUYUD(3,3)
      DOUBLE COMPLEX YUYUDHUHUD(3,3),HUYUDYUHUD(3,3),YUHUDHUYUD(3,3)
      DOUBLE COMPLEX HUHDDYDYUD(3,3),YUYDDHDHUD(3,3),HUYDDYDHUD(3,3)
      DOUBLE COMPLEX YUHDDHDYUD(3,3),HDHDDYDYDD(3,3),YDYDDHDHDD(3,3)
      DOUBLE COMPLEX HDYDDYDHDD(3,3),YDHDDHDYDD(3,3),HDHUDYUYDD(3,3)
      DOUBLE COMPLEX YDYUDHUHDD(3,3),HDYUDYUHDD(3,3),YDHUDHUYDD(3,3)
      DOUBLE COMPLEX HEHEDYEYED(3,3),YEYEDHEHED(3,3),HEYEDYEHED(3,3)
      DOUBLE COMPLEX YEHEDHEYED(3,3),YUDMUPMYUYUDYU(3,3)
      DOUBLE COMPLEX YUDMUPMYUYDDYD(3,3),YUDYUMQMYDDYD(3,3)
      DOUBLE COMPLEX YUDYUYDDMDMYD(3,3),YDDMDMYDYDDYD(3,3)
      DOUBLE COMPLEX YEDMEMYEYEDYE(3,3),YUDYUMQMYUDYU(3,3)
      DOUBLE COMPLEX YUDYUYUDMUPMYU(3,3),YUDYUYUDYUMQM(3,3)
      DOUBLE COMPLEX YDDYDMQMYDDYD(3,3),YDDYDYDDMDMYD(3,3)
      DOUBLE COMPLEX YDDYDYDDYDMQM(3,3),YEDYEMLMYEDYE(3,3)
      DOUBLE COMPLEX YEDYEYEDMEMYE(3,3),YEDYEYEDYEMLM(3,3)
      DOUBLE COMPLEX YUMQMYUDYUYUD(3,3),YUYUDMUPMYUYUD(3,3)
      DOUBLE COMPLEX YUYUDYUMQMYUD(3,3),YUYUDYUYUDMUPM(3,3)
      DOUBLE COMPLEX YUMQMYDDYDYUD(3,3),YUYDDMDMYDYUD(3,3)
      DOUBLE COMPLEX YUYDDYDMQMYUD(3,3),YUYDDYDYUDMUPM(3,3)
      DOUBLE COMPLEX YDMQMYDDYDYDD(3,3),YDYDDMDMYDYDD(3,3)
      DOUBLE COMPLEX YDYDDYDMQMYDD(3,3),YDYDDYDYDDMDM(3,3)
      DOUBLE COMPLEX YDMQMYUDYUYDD(3,3),YDYUDMUPMYUYDD(3,3)
      DOUBLE COMPLEX YDYUDYUMQMYDD(3,3),YDYUDYUYDDMDM(3,3)
      DOUBLE COMPLEX YEMLMYEDYEYED(3,3),YEYEDMEMYEYED(3,3)
      DOUBLE COMPLEX YEYEDYEMLMYED(3,3),YEYEDYEYEDMEM(3,3)
      DOUBLE COMPLEX TYUDYUYUDYU,TYDDYDYDDYD,TYEDYEYEDYE
      DOUBLE COMPLEX TYUDYUYDDYD,TYDDYDYUDYU
!
      DOUBLE COMPLEX BETA1U(3,3),BETA2U(3,3),BETA1D(3,3),BETA2D(3,3)
      DOUBLE COMPLEX BETA1E(3,3),BETA2E(3,3)
      DOUBLE COMPLEX BETAT2U(3,3),BETAT2D(3,3),BETAT2E(3,3)
      DOUBLE COMPLEX BETA1HU,BETA2HU,BETA1HD
      DOUBLE COMPLEX BETA2HD,BETA1MQ(3,3),BETA2MQ(3,3),BETA1ML(3,3)
      DOUBLE COMPLEX BETA2ML(3,3),BETA1MU(3,3),BETA2MU(3,3)
      DOUBLE COMPLEX BETA1MD(3,3)
      DOUBLE COMPLEX BETA2MD(3,3),BETA1ME(3,3),BETA2ME(3,3)
      DOUBLE COMPLEX BETA2GRKMU,BETA1B,BETA2B,BETA1VU,BETA1VD,BETA2VU
      DOUBLE COMPLEX BETA2VD,B2GM(3),B1YMU(3,3),B1YMD(3,3),B1YME(3,3)
      DOUBLE COMPLEX B1TRIU(3,3),B1TRID(3,3),B1TRIE(3,3)
      DOUBLE COMPLEX B2YMU(3,3),B2YMD(3,3),B2YME(3,3),BETAVEV1,BETAVEV2
      DOUBLE COMPLEX B2GRKMUM,B2M(3),B2HMU(3,3),B2HMD(3,3),B2HME(3,3)
      DOUBLE COMPLEX B2HUM,B2HDM,B2MQM(3,3),B2MLM(3,3),B2MUM(3,3)
      DOUBLE COMPLEX B2MDM(3,3),B2MEM(3,3),BETA2BM
!
      DOUBLE COMPLEX ID(3,3),MVMU,MVB
      DOUBLE PRECISION B1LP(3),B1LPM(3),B2LPSM(3,3),B2LPM(3,3)
      DOUBLE PRECISION CM(3,3),CSM(3,3),SINB,COSB
      DOUBLE PRECISION PI,Q
!
      INTEGER I,J,NG,ND,NE,NNU,NU,NSQ,NSU,NSD,NSL,NSE,NSH,NH,NSW,NSG
      INTEGER THLH,THHH,THSH,THSB,THSW,THGL
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $    ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $    ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $    ,ID(3,3)/(1.D0,0.D0)/
      DATA B1LPM(1)/6.6D0/,B1LPM(2)/1.D0/,B1LPM(3)/-3.D0/
      DATA B2LPM(1,1)/7.96D0/,B2LPM(1,2)/5.4D0/,B2LPM(1,3)/17.6D0/
      DATA B2LPM(2,1)/1.8D0/,B2LPM(2,2)/25.D0/,B2LPM(2,3)/24.D0/
      DATA B2LPM(3,1)/2.2D0/,B2LPM(3,2)/9.D0/,B2LPM(3,3)/14.D0/
      DATA CM(1,1)/5.2D0/,CM(1,2)/2.8D0/,CM(1,3)/3.6D0/
      DATA CM(2,1)/6.D0/,CM(2,2)/6.D0/,CM(2,3)/2.D0/
      DATA CM(3,1)/4.D0/,CM(3,2)/4.D0/,CM(3,3)/0.D0/
!
      Q=SSQSTEP !This means that thresholds depend on the q being
                !calculated by the main programme instead of the
                !integration subroutine
!
      DO I=1,601
        G(I)=GCURR(I)
      END DO
!
!Set all F's and betas to zero
!
      DO I=1,601
        F(I)=(0.D0,0.D0)
      END DO
      DO I=1,3
        DO J=1,3
          B1U(I,J)=(0.D0,0.D0)
          B1D(I,J)=(0.D0,0.D0)
          B1E(I,J)=(0.D0,0.D0)
          BETA1U(I,J)=(0.D0,0.D0)
          BETA2U(I,J)=(0.D0,0.D0)
          BETAT2U(I,J)=(0.D0,0.D0)
          BETA1D(I,J)=(0.D0,0.D0)
          BETA2D(I,J)=(0.D0,0.D0)
          BETA1E(I,J)=(0.D0,0.D0)
          BETA2E(I,J)=(0.D0,0.D0)
          BETA1MQ(I,J)=(0.D0,0.D0)
          BETA2MQ(I,J)=(0.D0,0.D0)
          BETA1ML(I,J)=(0.D0,0.D0)
          BETA2ML(I,J)=(0.D0,0.D0)
          BETA1MU(I,J)=(0.D0,0.D0)
          BETA2MU(I,J)=(0.D0,0.D0)
          BETA1MD(I,J)=(0.D0,0.D0)
          BETA2MD(I,J)=(0.D0,0.D0)
          BETA1ME(I,J)=(0.D0,0.D0)
          BETA2ME(I,J)=(0.D0,0.D0)
          BFTUQ(I,J)=(0.D0,0.D0)
          BFTDQ(I,J)=(0.D0,0.D0)
          BFTEL(I,J)=(0.D0,0.D0)
          BFTUU(I,J)=(0.D0,0.D0)
          BFTDD(I,J)=(0.D0,0.D0)
          BFTEE(I,J)=(0.D0,0.D0)
          BGTQ(I,J)=(0.D0,0.D0)
          BGTL(I,J)=(0.D0,0.D0)
          BGTPQ(I,J)=(0.D0,0.D0)
          BGTPL(I,J)=(0.D0,0.D0)
          BGTPU(I,J)=(0.D0,0.D0)
          BGTPD(I,J)=(0.D0,0.D0)
          BGTPE(I,J)=(0.D0,0.D0)
          BGTSQ(I,J)=(0.D0,0.D0)
          BGTSU(I,J)=(0.D0,0.D0)
          BGTSD(I,J)=(0.D0,0.D0)
          BMTSFU(I,J)=(0.D0,0.D0)
          BMTSFD(I,J)=(0.D0,0.D0)
          BMTSFE(I,J)=(0.D0,0.D0)
!          
          B1YMU(I,J)=(0.D0,0.D0)
          B1YMD(I,J)=(0.D0,0.D0)
          B1YME(I,J)=(0.D0,0.D0)
          B1TRIU(I,J)=(0.D0,0.D0)
          B1TRID(I,J)=(0.D0,0.D0)
          B1TRIE(I,J)=(0.D0,0.D0)
          B2YMU(I,J)=(0.D0,0.D0)
          B2YMD(I,J)=(0.D0,0.D0)
          B2YME(I,J)=(0.D0,0.D0)
          B2HMU(I,J)=(0.D0,0.D0)
          B2HMD(I,J)=(0.D0,0.D0)
          B2HME(I,J)=(0.D0,0.D0)
          B2MQM(I,J)=(0.D0,0.D0)
          B2MLM(I,J)=(0.D0,0.D0)
          B2MUM(I,J)=(0.D0,0.D0)
          B2MDM(I,J)=(0.D0,0.D0)
          B2MEM(I,J)=(0.D0,0.D0)
        END DO
        BM(I)=(0.D0,0.D0)
        BMP(I)=(0.D0,0.D0)
        B2GM(I)=(0.D0,0.D0)
        B2M(I)=(0.D0,0.D0)
      END DO
      BETA1HU=(0.D0,0.D0)
      BETA2HU=(0.D0,0.D0)
      BETA1HD=(0.D0,0.D0)
      BETA2HD=(0.D0,0.D0)
      B2HUM=(0.D0,0.D0)
      B2HDM=(0.D0,0.D0)
      BGRKMU=(0.D0,0.D0)
      BETA2GRKMU=(0.D0,0.D0)
      B2GRKMUM=(0.D0,0.D0)
      BETA1B=(0.D0,0.D0)
      BETA2B=(0.D0,0.D0)
      BETA2BM=(0.D0,0.D0)
      BETA1VU=(0.D0,0.D0)
      BETA1VD=(0.D0,0.D0)
      BETA2VU=(0.D0,0.D0)
      BETA2VD=(0.D0,0.D0)
      BGTHU=(0.D0,0.D0)
      BGTHD=(0.D0,0.D0)
      BGTPHU=(0.D0,0.D0)
      BGTPHD=(0.D0,0.D0)
      BSGTHU=(0.D0,0.D0)
      BCGTHD=(0.D0,0.D0)
      BSGTPHU=(0.D0,0.D0)
      BCGTPHD=(0.D0,0.D0)
      BETALAM1=(0.D0,0.D0)
      BETALAM2=(0.D0,0.D0)
      BETAVEV1=(0.D0,0.D0)
      BETAVEV2=(0.D0,0.D0)
!
      BMHUPMT=(0D0,0.D0)
      BMHDPMT=(0D0,0.D0)
      BMHUD=(0D0,0.D0)
!
      PI=4.D0*ATAN(1.D0)
      SINB=DSQRT(DBLE(TANB)**2/(1+DBLE(TANB)**2))
      COSB=SINB/TANB
      IF(Q.LT.1.D0)THEN
        WRITE(*,*)'ERROR IN Q: ',Q
        STOP
      END IF
      NG=3.D0
      NU=3
      ND=3
      NE=3
      NNU=3
      IF ((Q-QNSH).GT.ABS(EPS).OR.
     $         (ABS(Q-QNSH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSH=2
      ELSE
        NSH=0
      END IF
      IF ((Q-QNSG).GT.ABS(EPS).OR.
     $         (ABS(Q-QNSG).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF ((Q-QNH).GT.ABS(EPS).OR.
     $         (ABS(Q-QNH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF ((Q-QTHSB).GT.ABS(EPS).OR.
     $         (ABS(Q-QTHSB).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        THSB=1
      ELSE
        THSB=0
      END IF
      IF ((Q-QTHSW).GT.ABS(EPS).OR.
     $         (ABS(Q-QTHSW).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSW=1
        THSW=1
      ELSE
        NSW=0
        THSW=0
      END IF
      THLH=1
      THHH=NH/2 !This works so long as THHH is an integer variable
      IF(THHH.NE.1)THEN !Perform check
        IF(THHH.NE.0)WRITE(*,*)'ERROR IN THHH'
      END IF
      THSH=NSH/2
      THGL=NSG
!
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      DO I=1,3
        NSQ=NSQ+THSQ(I)
        NSU=NSU+THSU(I)
        NSD=NSD+THSD(I)
        NSL=NSL+THSL(I)
        NSE=NSE+THSE(I)
      END DO
!
!Since I remove entries in the sfermion mass matrix
!as they decouple, I need to continue rotating even
!when all NS_ are zero.
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15)THEN
!
!Rotate to the mass basis.
!
        CALL ROTSQ(GCURR,G)
!
      END IF
!
!Convert input into 3x3 matrices. For the Yukawas and a-parameters
!I also convert notation and define Y and H as Martin and Vaughn
!parameters.
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
          AU(I,J)=G(33+(I-1)*3+J)
          AD(I,J)=G(42+(I-1)*3+J)
          AE(I,J)=G(51+(I-1)*3+J)
          MQ(I,J)=G(62+(I-1)*3+J)
          ML(I,J)=G(71+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
          ME(I,J)=G(98+(I-1)*3+J)
!
          YU(I,J)=G(293+(J-1)*3+I)
          YD(I,J)=G(302+(J-1)*3+I)
          YE(I,J)=G(311+(J-1)*3+I)
          HU(I,J)=-G(323+(J-1)*3+I)
          HD(I,J)=-G(332+(J-1)*3+I)
          HE(I,J)=-G(341+(J-1)*3+I)
          MQM(I,J)=G(352+(I-1)*3+J)
          MLM(I,J)=G(361+(I-1)*3+J)
          MUPM(I,J)=G(370+(I-1)*3+J)
          MDM(I,J)=G(379+(I-1)*3+J)
          MEM(I,J)=G(388+(I-1)*3+J)
!
          LU(I,J)=G(111+(I-1)*3+J)
          LD(I,J)=G(120+(I-1)*3+J)
          LE(I,J)=G(129+(I-1)*3+J)
!
          GTPQ(I,J)=G(138+(I-1)*3+J)
          GTPL(I,J)=G(147+(I-1)*3+J)
          GTPU(I,J)=G(156+(I-1)*3+J)
          GTPD(I,J)=G(165+(I-1)*3+J)
          GTPE(I,J)=G(174+(I-1)*3+J)
          GTQ(I,J)=G(185+(I-1)*3+J)
          GTL(I,J)=G(194+(I-1)*3+J)
          GTSQ(I,J)=G(205+(I-1)*3+J)
          GTSU(I,J)=G(214+(I-1)*3+J)
          GTSD(I,J)=G(223+(I-1)*3+J)
          FTUQ(I,J)=G(232+(I-1)*3+J)
          FTDQ(I,J)=G(241+(I-1)*3+J)
          FTEL(I,J)=G(250+(I-1)*3+J)
          FTUU(I,J)=G(259+(I-1)*3+J)
          FTDD(I,J)=G(268+(I-1)*3+J)
          FTEE(I,J)=G(277+(I-1)*3+J)
!
          TRIU(I,J)=G(399+(I-1)*3+J)
          TRID(I,J)=G(408+(I-1)*3+J)
          TRIE(I,J)=G(417+(I-1)*3+J)
          MTSFU(I,J)=G(429+(I-1)*3+J)
          MTSFD(I,J)=G(438+(I-1)*3+J)
          MTSFE(I,J)=G(447+(I-1)*3+J)
!
!For now the quartics are not run independently
!If you wish to introduce independent quartic running
!see the note at the top of this subroutine.
!Care should be taken with quartics which are
!equal by SU(2) invariance.
!
          FUHU(I,J)=FU(I,J)
          FDHD(I,J)=FD(I,J)
          FEHD(I,J)=FE(I,J)
          FUQ(I,J)=FU(I,J)
          FDQ(I,J)=FD(I,J)
          FUUR(I,J)=FU(I,J)
          FDDR(I,J)=FD(I,J)
          FEL(I,J)=FE(I,J)
          FEER(I,J)=FE(I,J)
          IF(THHH.EQ.0)THEN
            SFUQ(I,J)=LU(I,J)
            CFDQ(I,J)=LD(I,J)
            SFUUR(I,J)=LU(I,J)
            CFDDR(I,J)=LD(I,J)
            CFEL(I,J)=LE(I,J)
            CFEER(I,J)=LE(I,J)
          ELSE
            SFUQ(I,J)=(0.D0,0.D0)
            CFDQ(I,J)=(0.D0,0.D0)
            SFUUR(I,J)=(0.D0,0.D0)
            CFDDR(I,J)=(0.D0,0.D0)
            CFEL(I,J)=(0.D0,0.D0)
            CFEER(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Set other quartics to their counterparts
!
      G(457)=DSQRT(3.D0/5.D0)*G(1)
      G(458)=G(2)
      G(459)=G(3)
      IF(THHH.EQ.0)THEN
        G(541)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))
     $                             *DSQRT(3.D0/5.D0)*G(1)
        G(542)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*G(2)
        G(543)=(COSB**2-SINB**2)*DSQRT(3.D0/5.D0)*G(1)
        G(544)=(COSB**2-SINB**2)*G(2)
      ELSE
        G(541)=(0.D0,0.D0)
        G(542)=(0.D0,0.D0)
        G(543)=(0.D0,0.D0)
        G(544)=(0.D0,0.D0)
      END IF
!
!The mu and b running is partly in MV notation, so I must convert
!NB: MVMU is the MSSM MV notation \mu.
!
      MVMU=G(398)
      MVB=-G(399)
!      
!The separated out contributions are from PRD 49 4882 (1194),
!Castano,Piard,Ramond
!
      B1LP(1)=2.D0/5.D0*(17.D0/12.D0*DBLE(NU)+5.D0/12.D0*DBLE(ND)
     $        +5.D0/4.D0*DBLE(NE)+1.D0/4.D0*DBLE(NNU))+1.D0/30.D0
     $        *DBLE(NSQ)+4.D0/15.D0*DBLE(NSU)+1.D0/15.D0*DBLE(NSD)
     $        +1.D0/10.D0*DBLE(NSL)+1.D0/5.D0*DBLE(NSE)
     $        +1.D0/5.D0*DBLE(NSH)+1.D0/10.D0*DBLE(NH)
      B1LP(2)=-22.D0/3.D0+1.D0/2.D0*(DBLE(NU)+DBLE(ND))+1.D0/6.D0
     $        *(DBLE(NE)+DBLE(NNU))+1.D0/2.D0*DBLE(NSQ)
     $        +1.D0/6.D0*DBLE(NSL)+1.D0/3.D0*DBLE(NSH)+1.D0/6.D0
     $        *DBLE(NH)+4.D0/3.D0*DBLE(NSW)
      B1LP(3)=-11.D0+2.D0/3.D0*(DBLE(NU)+DBLE(ND))+1.D0/3.D0
     $        *DBLE(NSQ)+1.D0/6.D0*DBLE(NSU)+1.D0/6.D0*DBLE(NSD)
     $        +2.D0*DBLE(NSG)
      B2LPSM(1,1)=-(-NG*19.D0/15.D0-9.D0/50.D0)
      B2LPSM(1,2)=-(-NG*3.D0/5.D0-9.D0/10.D0)
      B2LPSM(1,3)=-(-NG*44.D0/15.D0)
      B2LPSM(2,1)=-(-NG*1.D0/5.D0-3.D0/10.D0)
      B2LPSM(2,2)=-(136.D0/3.D0-NG*49.D0/3.D0-13.D0/6.D0)
      B2LPSM(2,3)=-(-NG*4.D0)
      B2LPSM(3,1)=-(-NG*11.D0/30.D0)
      B2LPSM(3,2)=-(-NG*3.D0/2.D0)
      B2LPSM(3,3)=-(102.D0-NG*76.D0/3.D0)
      CSM(1,1)=1.7D0
      CSM(1,2)=.5D0
      CSM(1,3)=1.5D0
      CSM(2,1)=1.5D0
      CSM(2,2)=1.5D0
      CSM(2,3)=.5D0
      CSM(3,1)=2.D0
      CSM(3,2)=2.D0
      CSM(3,3)=0.D0
!
!I need many variations on the 3x3 matrices.
!
      CALL CDAGGER(FTUQ,FTUQD)
      CALL CDAGGER(FTUU,FTUUD)
      CALL CDAGGER(FTDD,FTDDD)
      CALL CDAGGER(FTDQ,FTDQD)
      CALL CDAGGER(FTEL,FTELD)
      CALL CDAGGER(FTEE,FTEED)
      CALL CDAGGER(GTQ,GTQD)
      CALL CDAGGER(GTL,GTLD)
      CALL CDAGGER(GTPQ,GTPQD)
      CALL CDAGGER(GTPU,GTPUD)
      CALL CDAGGER(GTPD,GTPDD)
      CALL CDAGGER(GTPL,GTPLD)
      CALL CDAGGER(GTPE,GTPED)
      CALL CDAGGER(GTSQ,GTSQD)
      CALL CDAGGER(GTSU,GTSUD)
      CALL CDAGGER(GTSD,GTSDD)
      CALL CDAGGER(FUHU,FUHUD)
      CALL CDAGGER(FDHD,FDHDD)
      CALL CDAGGER(FEHD,FEHDD)
      CALL CDAGGER(TRIU,TRIUD)
      CALL CDAGGER(TRID,TRIDD)
      CALL CDAGGER(TRIE,TRIED)
      DO I=1,3
        DO J=1,3
          FUS(I,J)=CONJG(FU(I,J))
          FDS(I,J)=CONJG(FD(I,J))
          FES(I,J)=CONJG(FE(I,J))
          LUS(I,J)=CONJG(LU(I,J))
          LDS(I,J)=CONJG(LD(I,J))
          LES(I,J)=CONJG(LE(I,J))
          AUS(I,J)=CONJG(AU(I,J))
          ADS(I,J)=CONJG(AD(I,J))
          AES(I,J)=CONJG(AE(I,J))
          FTUQS(I,J)=CONJG(FTUQ(I,J))
          FTUUS(I,J)=CONJG(FTUU(I,J))
          FTDQS(I,J)=CONJG(FTDQ(I,J))
          FTDDS(I,J)=CONJG(FTDD(I,J))
          FTELS(I,J)=CONJG(FTEL(I,J))
          FTEES(I,J)=CONJG(FTEE(I,J))
          GTPUS(I,J)=CONJG(GTPU(I,J))
          GTSUS(I,J)=CONJG(GTSU(I,J))
          GTQS(I,J)=CONJG(GTQ(I,J))
          GTLS(I,J)=CONJG(GTL(I,J))
          GTPQS(I,J)=CONJG(GTPQ(I,J))
          GTPLS(I,J)=CONJG(GTPL(I,J))
          GTSQS(I,J)=CONJG(GTSQ(I,J))
          GTPDS(I,J)=CONJG(GTPD(I,J))
          GTPES(I,J)=CONJG(GTPE(I,J))
          GTSDS(I,J)=CONJG(GTSD(I,J))
          MTSFUS(I,J)=CONJG(MTSFU(I,J))
          MTSFDS(I,J)=CONJG(MTSFD(I,J))
          MTSFES(I,J)=CONJG(MTSFE(I,J))
          FUHUS(I,J)=CONJG(FUHU(I,J))
          FDHDS(I,J)=CONJG(FDHD(I,J))
          FEHDS(I,J)=CONJG(FEHD(I,J))
          TRIUS(I,J)=CONJG(TRIU(I,J))
          TRIDS(I,J)=CONJG(TRID(I,J))
          TRIES(I,J)=CONJG(TRIE(I,J))
          FUT(I,J)=FU(J,I)
          FDT(I,J)=FD(J,I)
          FET(I,J)=FE(J,I)
          LUT(I,J)=LU(J,I)
          LDT(I,J)=LD(J,I)
          LET(I,J)=LE(J,I)
          AUT(I,J)=AU(J,I)
          ADT(I,J)=AD(J,I)
          AET(I,J)=AE(J,I)
          FTUQT(I,J)=FTUQ(J,I)
          FTUUT(I,J)=FTUU(J,I)
          FTDQT(I,J)=FTDQ(J,I)
          FTDDT(I,J)=FTDD(J,I)
          FTELT(I,J)=FTEL(J,I)
          FTEET(I,J)=FTEE(J,I)
          GTPUT(I,J)=GTPU(J,I)
          GTQT(I,J)=GTQ(J,I)
          GTLT(I,J)=GTL(J,I)
          GTPQT(I,J)=GTPQ(J,I)
          GTPLT(I,J)=GTPL(J,I)
          GTPDT(I,J)=GTPD(J,I)
          GTPET(I,J)=GTPE(J,I)
          GTSUT(I,J)=GTSU(J,I)
          GTSQT(I,J)=GTSQ(J,I)
          GTSDT(I,J)=GTSD(J,I)
          MTSFUT(I,J)=MTSFU(J,I)
          MTSFDT(I,J)=MTSFD(J,I)
          MTSFET(I,J)=MTSFE(J,I)
          FUHUT(I,J)=FUHU(J,I)
          FDHDT(I,J)=FDHD(J,I)
          FEHDT(I,J)=FEHD(J,I)
          TRIUT(I,J)=TRIU(J,I)
          TRIDT(I,J)=TRID(J,I)
          TRIET(I,J)=TRIE(J,I)
        END DO
      END DO
!
!Now all the matrix multiples
!
      DO I=1,3
        DO J=1,3
          YUDYU(I,J)=CMATMUL(1,YU,YU,I,J)
          YDDYD(I,J)=CMATMUL(1,YD,YD,I,J)
          YEDYE(I,J)=CMATMUL(1,YE,YE,I,J)
          YUYUD(I,J)=CMATMUL(2,YU,YU,I,J)
          YDYDD(I,J)=CMATMUL(2,YD,YD,I,J)
          YEYED(I,J)=CMATMUL(2,YE,YE,I,J)
          YUYDD(I,J)=CMATMUL(2,YU,YD,I,J)
          YDYUD(I,J)=CMATMUL(2,YD,YU,I,J)
!
          HUYUD(I,J)=CMATMUL(2,HU,YU,I,J)
          HDYDD(I,J)=CMATMUL(2,HD,YD,I,J)
          HEYED(I,J)=CMATMUL(2,HE,YE,I,J)
          YUDHU(I,J)=CMATMUL(1,YU,HU,I,J)
          YDDHD(I,J)=CMATMUL(1,YD,HD,I,J)
          YEDHE(I,J)=CMATMUL(1,YE,HE,I,J)
          HUDHU(I,J)=CMATMUL(1,HU,HU,I,J)
          HDDHD(I,J)=CMATMUL(1,HD,HD,I,J)
          HEDHE(I,J)=CMATMUL(1,HE,HE,I,J)
          HUDYU(I,J)=CMATMUL(1,HU,YU,I,J)
          HDDYD(I,J)=CMATMUL(1,HD,YD,I,J)
          HEDYE(I,J)=CMATMUL(1,HE,YE,I,J)
          HUHUD(I,J)=CMATMUL(2,HU,HU,I,J)
          HDHDD(I,J)=CMATMUL(2,HD,HD,I,J)
          HEHED(I,J)=CMATMUL(2,HE,HE,I,J)
          YUHUD(I,J)=CMATMUL(2,YU,HU,I,J)
          YDHDD(I,J)=CMATMUL(2,YD,HD,I,J)
          YEHED(I,J)=CMATMUL(2,YE,HE,I,J)
          HUYDD(I,J)=CMATMUL(2,HU,YD,I,J)
          HDYUD(I,J)=CMATMUL(2,HD,YU,I,J)
          HUHDD(I,J)=CMATMUL(2,HU,HD,I,J)
          HDHUD(I,J)=CMATMUL(2,HD,HU,I,J)
          YDHUD(I,J)=CMATMUL(2,YD,HU,I,J)
          YUHDD(I,J)=CMATMUL(2,YU,HD,I,J)
!
          MUPMYU(I,J)=CMATMUL(0,MUPM,YU,I,J)
          MDMYD(I,J)=CMATMUL(0,MDM,YD,I,J)
          MEMYE(I,J)=CMATMUL(0,MEM,YE,I,J)
          YUMQM(I,J)=CMATMUL(0,YU,MQM,I,J)
          YDMQM(I,J)=CMATMUL(0,YD,MQM,I,J)
          YEMLM(I,J)=CMATMUL(0,YE,MLM,I,J)
!
          FUFUD(I,J)=CMATMUL(2,FU,FU,I,J)
          FUSFUT(I,J)=CMATMUL(0,FUS,FUT,I,J)
          LUSLUT(I,J)=CMATMUL(0,LUS,LUT,I,J)
          FUTFUS(I,J)=CMATMUL(0,FUT,FUS,I,J)
          LUTLUS(I,J)=CMATMUL(0,LUT,LUS,I,J)
          FUTFTUUS(I,J)=CMATMUL(0,FUT,FTUUS,I,J)
          LUTFTUUS(I,J)=CMATMUL(0,LUT,FTUUS,I,J)
          FUGTPUS(I,J)=CMATMUL(0,FU,GTPUS,I,J)
          FUGTSUS(I,J)=CMATMUL(0,FU,GTSUS,I,J)
          LUGTPUS(I,J)=CMATMUL(0,LU,GTPUS,I,J)
          LUGTSUS(I,J)=CMATMUL(0,LU,GTSUS,I,J)
          FDFDD(I,J)=CMATMUL(2,FD,FD,I,J)
          FDTFDS(I,J)=CMATMUL(0,FDT,FDS,I,J)
          LDTLDS(I,J)=CMATMUL(0,LDT,LDS,I,J)
          FDTFTDDS(I,J)=CMATMUL(0,FDT,FTDDS,I,J)
          LDTFTDDS(I,J)=CMATMUL(0,LDT,FTDDS,I,J)
          FDSFDT(I,J)=CMATMUL(0,FDS,FDT,I,J)
          LDSLDT(I,J)=CMATMUL(0,LDS,LDT,I,J)
          FDGTPDS(I,J)=CMATMUL(0,FD,GTPDS,I,J)
          FDGTSDS(I,J)=CMATMUL(0,FD,GTSDS,I,J)
          LDGTPDS(I,J)=CMATMUL(0,LD,GTPDS,I,J)
          LDGTSDS(I,J)=CMATMUL(0,LD,GTSDS,I,J)
          FEFED(I,J)=CMATMUL(2,FE,FE,I,J)
          FETFES(I,J)=CMATMUL(0,FET,FES,I,J)
          LETLES(I,J)=CMATMUL(0,LET,LES,I,J)
          FETFTEES(I,J)=CMATMUL(0,FET,FTEES,I,J)
          LETFTEES(I,J)=CMATMUL(0,LET,FTEES,I,J)
          FESFET(I,J)=CMATMUL(0,FES,FET,I,J)
          LESLET(I,J)=CMATMUL(0,LES,LET,I,J)
          FEGTPES(I,J)=CMATMUL(0,FE,GTPES,I,J)
          LEGTPES(I,J)=CMATMUL(0,LE,GTPES,I,J)
          FUDFU(I,J)=CMATMUL(1,FU,FU,I,J)
          FDDFD(I,J)=CMATMUL(1,FD,FD,I,J)
          FEDFE(I,J)=CMATMUL(1,FE,FE,I,J)
          LUDLU(I,J)=CMATMUL(1,LU,LU,I,J)
          LDDLD(I,J)=CMATMUL(1,LD,LD,I,J)
          LEDLE(I,J)=CMATMUL(1,LE,LE,I,J)
          LULUD(I,J)=CMATMUL(2,LU,LU,I,J)
          LDLDD(I,J)=CMATMUL(2,LD,LD,I,J)
          LELED(I,J)=CMATMUL(2,LE,LE,I,J)
!
          FTUQFTUQD(I,J)=CMATMUL(2,FTUQ,FTUQ,I,J)
          FTUQSFTUQT(I,J)=CMATMUL(0,FTUQS,FTUQT,I,J)
          FTUQSFUT(I,J)=CMATMUL(0,FTUQS,FUT,I,J)
          FTUQSLUT(I,J)=CMATMUL(0,FTUQS,LUT,I,J)
          FTDQFTDQD(I,J)=CMATMUL(2,FTDQ,FTDQ,I,J)
          FTDQSFTDQT(I,J)=CMATMUL(0,FTDQS,FTDQT,I,J)
          FTDQSFDT(I,J)=CMATMUL(0,FTDQS,FDT,I,J)
          FTDQSLDT(I,J)=CMATMUL(0,FTDQS,LDT,I,J)
          FTELFTELD(I,J)=CMATMUL(2,FTEL,FTEL,I,J)
          FTELSFTELT(I,J)=CMATMUL(0,FTELS,FTELT,I,J)
          FTELSFET(I,J)=CMATMUL(0,FTELS,FET,I,J)
          FTELSLET(I,J)=CMATMUL(0,FTELS,LET,I,J)
          GTQSGTQT(I,J)=CMATMUL(0,GTQS,GTQT,I,J)
          GTQGTQD(I,J)=CMATMUL(2,GTQ,GTQ,I,J)
          GTPQSGTPQT(I,J)=CMATMUL(0,GTPQS,GTPQT,I,J)
          GTPQGTPQD(I,J)=CMATMUL(2,GTPQ,GTPQ,I,J)
          GTPUTGTPUS(I,J)=CMATMUL(0,GTPUT,GTPUS,I,J)
          GTPUDGTPU(I,J)=CMATMUL(1,GTPU,GTPU,I,J)
          GTPDDGTPD(I,J)=CMATMUL(1,GTPD,GTPD,I,J)
          GTPDTGTPDS(I,J)=CMATMUL(0,GTPDT,GTPDS,I,J)
          GTSQSGTSQT(I,J)=CMATMUL(0,GTSQS,GTSQT,I,J)
          GTSQGTSQD(I,J)=CMATMUL(2,GTSQ,GTSQD,I,J)
          GTSUTGTSUS(I,J)=CMATMUL(0,GTSUT,GTSUS,I,J)
          GTSUDGTSU(I,J)=CMATMUL(1,GTSU,GTSU,I,J)
          GTSDDGTSD(I,J)=CMATMUL(1,GTSD,GTSD,I,J)
          GTSDTGTSDS(I,J)=CMATMUL(0,GTSDT,GTSDS,I,J)
          GTLSGTLT(I,J)=CMATMUL(0,GTLS,GTLT,I,J)
          GTLGTLD(I,J)=CMATMUL(2,GTL,GTL,I,J)
          GTPLSGTPLT(I,J)=CMATMUL(0,GTPLS,GTPLT,I,J)
          GTPLGTPLD(I,J)=CMATMUL(2,GTPL,GTPL,I,J)
          GTPETGTPES(I,J)=CMATMUL(0,GTPET,GTPES,I,J)
          GTPEDGTPE(I,J)=CMATMUL(1,GTPE,GTPE,I,J)
          FUURFUURD(I,J)=CMATMUL(2,FUUR,FUUR,I,J)
          FUQDFUQ(I,J)=CMATMUL(1,FUQ,FUQ,I,J)
          FUQDFDQ(I,J)=CMATMUL(1,FUQ,FDQ,I,J)
          FDQDFUQ(I,J)=CMATMUL(1,FDQ,FUQ,I,J)
          FDQDFDQ(I,J)=CMATMUL(1,FDQ,FDQ,I,J)
          FDDRFDDRD(I,J)=CMATMUL(2,FDDR,FDDR,I,J)
          FEERFEERD(I,J)=CMATMUL(2,FEER,FEER,I,J)
          FELDFEL(I,J)=CMATMUL(1,FEL,FEL,I,J)
          FTUQGTPUS(I,J)=CMATMUL(0,FTUQ,GTPUS,I,J)
          FTDQGTPDS(I,J)=CMATMUL(0,FTDQ,GTPDS,I,J)
          FTELGTPES(I,J)=CMATMUL(0,FTEL,GTPES,I,J)
          GTPQSFTUU(I,J)=CMATMUL(0,GTPQS,FTUU,I,J)
          GTPQSFTDD(I,J)=CMATMUL(0,GTPQS,FTDD,I,J)
          GTPLSFTEE(I,J)=CMATMUL(0,GTPLS,FTEE,I,J)
          GTQSFTUU(I,J)=CMATMUL(0,GTQS,FTUU,I,J)
          GTQSFTDD(I,J)=CMATMUL(0,GTQS,FTDD,I,J)
          GTLSFTEE(I,J)=CMATMUL(0,GTLS,FTEE,I,J)
          FDDFTUU(I,J)=CMATMUL(1,FD,FTUU,I,J)
          LDDFTUU(I,J)=CMATMUL(1,LD,FTUU,I,J)
          FUDFTDD(I,J)=CMATMUL(1,FU,FTDD,I,J)
          LUDFTDD(I,J)=CMATMUL(1,LU,FTDD,I,J)
          FTUUDFTUU(I,J)=CMATMUL(0,FTUUD,FTUU,I,J)
          FTDDDFTDD(I,J)=CMATMUL(0,FTDDD,FTDD,I,J)
          FTEEDFTEE(I,J)=CMATMUL(0,FTEED,FTEE,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          FUQTFUQS(I,J)=FUQDFUQ(J,I)
          FUURSFUURT(I,J)=FUURFUURD(J,I)
          FDQTFDQS(I,J)=FDQDFDQ(J,I)
          FDDRSFDDRT(I,J)=FDDRFDDRD(J,I)
          FEERSFEERT(I,J)=FEERFEERD(J,I)
          FELTFELS(I,J)=FELDFEL(J,I)
          FTUUTFTUUS(I,J)=FTUUDFTUU(J,I)
          FTDDTFTDDS(I,J)=FTDDDFTDD(J,I)
          FTEETFTEES(I,J)=FTEEDFTEE(J,I)
        END DO
      END DO
!
      IF(THHH.EQ.0)THEN
        MSGTPHUSQ=CONJG(G(287))*G(287)
        MCGTPHDSQ=CONJG(G(288))*G(288)
        MSGTHUSQ=CONJG(G(289))*G(289)
        MCGTHDSQ=CONJG(G(290))*G(290)
      END IF
      MGTPHUSQ=CONJG(G(184))*G(184)
      MGTPHDSQ=CONJG(G(185))*G(185)
      MGTHUSQ=CONJG(G(204))*G(204)
      MGTHDSQ=CONJG(G(205))*G(205)
      MMUSQ=CONJG(G(108))*G(108)
      M1PM1PSQ=G(31)**2+G(599)**2
      M2PM2PSQ=G(32)**2+G(600)**2
      M3PM3PSQ=G(33)**2+G(601)**2
!
      TYUDYU=CTRACE(YUDYU)
      TYDDYD=CTRACE(YDDYD)
      TYEDYE=CTRACE(YEDYE)
      TYUDHU=CTRACE(YUDHU)
      TYDDHD=CTRACE(YDDHD)
      TYEDHE=CTRACE(YEDHE)
      THUYUD=CTRACE(HUYUD)
!
      TFUDFU=CTRACE(FUDFU)
      TFDDFD=CTRACE(FDDFD)
      TFEDFE=CTRACE(FEDFE)
      TLUDLU=CTRACE(LULUD)
      TLDDLD=CTRACE(LDLDD)
      TLEDLE=CTRACE(LELED)
!
      DO I=1,3
        DO J=1,3
          YUDMUPMYU(I,J)=CMATMUL(1,YU,MUPMYU,I,J)
          YDDMDMYD(I,J)=CMATMUL(1,YD,MDMYD,I,J)
          YEDMEMYE(I,J)=CMATMUL(1,YE,MEMYE,I,J)
          YEDYEMLM(I,J)=CMATMUL(1,YE,YEMLM,I,J)
          YUMQMYUD(I,J)=CMATMUL(2,YUMQM,YU,I,J)
          YUYUDMUPM(I,J)=CMATMUL(0,YUYUD,MUPM,I,J)
          YDMQMYDD(I,J)=CMATMUL(2,YDMQM,YD,I,J)
          YDYDDMDM(I,J)=CMATMUL(0,YDYDD,MDM,I,J)
          YEMLMYED(I,J)=CMATMUL(2,YEMLM,YE,I,J)
          YEYEDMEM(I,J)=CMATMUL(0,YEYED,MEM,I,J)
!
          FUFUDFU(I,J)=CMATMUL(0,FUFUD,FU,I,J)
          FUFUDFD(I,J)=CMATMUL(0,FUFUD,FD,I,J)
          FUTFUSGTPU(I,J)=CMATMUL(0,FUTFUS,GTPU,I,J)
          LUTLUSGTPU(I,J)=CMATMUL(0,LUTLUS,GTPU,I,J)
          FUTFUSGTSU(I,J)=CMATMUL(0,FUTFUS,GTSU,I,J)
          LUTLUSGTSU(I,J)=CMATMUL(0,LUTLUS,GTSU,I,J)
          FDFDDFU(I,J)=CMATMUL(0,FDFDD,FU,I,J)
          FDFDDFD(I,J)=CMATMUL(0,FDFDD,FD,I,J)
          FDTFDSGTPD(I,J)=CMATMUL(0,FDTFDS,GTPD,I,J)
          LDTLDSGTPD(I,J)=CMATMUL(0,LDTLDS,GTPD,I,J)
          FDTFDSGTSD(I,J)=CMATMUL(0,FDTFDS,GTSD,I,J)
          LDTLDSGTSD(I,J)=CMATMUL(0,LDTLDS,GTSD,I,J)
          FEFEDFE(I,J)=CMATMUL(0,FEFED,FE,I,J)
          FETFESGTPE(I,J)=CMATMUL(0,FETFES,GTPE,I,J)
          LETLESGTPE(I,J)=CMATMUL(0,LETLES,GTPE,I,J)
          LULUDLU(I,J)=CMATMUL(0,LULUD,LU,I,J)
          LULUDLD(I,J)=CMATMUL(0,LULUD,LD,I,J)
          LDLDDLU(I,J)=CMATMUL(0,LDLDD,LU,I,J)
          LDLDDLD(I,J)=CMATMUL(0,LDLDD,LD,I,J)
          LELEDLE(I,J)=CMATMUL(0,LELED,LE,I,J)
          FTUQFUDFU(I,J)=CMATMUL(0,FTUQ,FUDFU,I,J)
          FTUQLUDLU(I,J)=CMATMUL(0,FTUQ,LUDLU,I,J)
          FTDQFDDFD(I,J)=CMATMUL(0,FTDQ,FDDFD,I,J)
          FTDQLDDLD(I,J)=CMATMUL(0,FTDQ,LDDLD,I,J)
          FTELFEDFE(I,J)=CMATMUL(0,FTEL,FEDFE,I,J)
          FTELLEDLE(I,J)=CMATMUL(0,FTEL,LEDLE,I,J)
          FTUQFUDFTDD(I,J)=CMATMUL(0,FTUQ,FUDFTDD,I,J)
          FTUQLUDFTDD(I,J)=CMATMUL(0,FTUQ,LUDFTDD,I,J)
          FTDQFDDFTUU(I,J)=CMATMUL(0,FTDQ,FDDFTUU,I,J)
          FTDQLDDFTUU(I,J)=CMATMUL(0,FTDQ,LDDFTUU,I,J)
          GTPQSFUGTPUS(I,J)=CMATMUL(0,GTPQS,FUGTPUS,I,J)
          GTPQSLUGTPUS(I,J)=CMATMUL(0,GTPQS,LUGTPUS,I,J)
          GTSQSFUGTSUS(I,J)=CMATMUL(0,GTSQS,FUGTSUS,I,J)
          GTSQSLUGTSUS(I,J)=CMATMUL(0,GTSQS,LUGTSUS,I,J)
          GTPQSFDGTPDS(I,J)=CMATMUL(0,GTPQS,FDGTPDS,I,J)
          GTPQSLDGTPDS(I,J)=CMATMUL(0,GTPQS,LDGTPDS,I,J)
          GTSQSFDGTSDS(I,J)=CMATMUL(0,GTSQS,FDGTSDS,I,J)
          GTSQSLDGTSDS(I,J)=CMATMUL(0,GTSQS,LDGTSDS,I,J)
          GTPLSFEGTPES(I,J)=CMATMUL(0,GTPLS,FEGTPES,I,J)
          GTPLSLEGTPES(I,J)=CMATMUL(0,GTPLS,LEGTPES,I,J)
          SQFTUQDFTUQ(I,J)=CSFMUL(THSQ,FTUQD,FTUQ,I,J)
          SQFTUQTFTUQS(I,J)=CSFMUL(THSQ,FTUQT,FTUQS,I,J)
          SQFTUQTGTQ(I,J)=CSFMUL(THSQ,FTUQT,GTQ,I,J)
          SQFTUQTGTPQ(I,J)=CSFMUL(THSQ,FTUQT,GTPQ,I,J)
          SQFTUQTGTSQ(I,J)=CSFMUL(THSQ,FTUQT,GTSQ,I,J)
          SQFTDQDFTDQ(I,J)=CSFMUL(THSQ,FTDQD,FTDQ,I,J)
          SQFTDQTFTDQS(I,J)=CSFMUL(THSQ,FTDQT,FTDQS,I,J)
          SQFTDQTGTQ(I,J)=CSFMUL(THSQ,FTDQT,GTQ,I,J)
          SQFTDQTGTPQ(I,J)=CSFMUL(THSQ,FTDQT,GTPQ,I,J)
          SQFTDQTGTSQ(I,J)=CSFMUL(THSQ,FTDQT,GTSQ,I,J)
          SQGTQTGTQS(I,J)=CSFMUL(THSQ,GTQT,GTQS,I,J)
          SQGTQTFTUQ(I,J)=CSFMUL(THSQ,GTQT,FTUQ,I,J)
          SQGTQTFTDQ(I,J)=CSFMUL(THSQ,GTQT,FTDQ,I,J)
          SQGTQDGTQ(I,J)=CSFMUL(THSQ,GTQD,GTQ,I,J)
          SQGTPQTGTPQS(I,J)=CSFMUL(THSQ,GTPQT,GTPQS,I,J)
          SQGTPQTFTUQ(I,J)=CSFMUL(THSQ,GTPQT,FTUQ,I,J)
          SQGTPQTFTDQ(I,J)=CSFMUL(THSQ,GTPQT,FTDQ,I,J)
          SQGTPQDGTPQ(I,J)=CSFMUL(THSQ,GTPQD,GTPQ,I,J)
          SQGTSQTGTSQS(I,J)=CSFMUL(THSQ,GTSQT,GTSQS,I,J)
          SQGTSQTFTUQ(I,J)=CSFMUL(THSQ,GTSQT,FTUQ,I,J)
          SQGTSQTFTDQ(I,J)=CSFMUL(THSQ,GTSQT,FTDQ,I,J)
          SQGTSQDGTSQ(I,J)=CSFMUL(THSQ,GTSQD,GTSQ,I,J)
          SDMTSFDSADT(I,J)=CSFMUL(THSD,MTSFDS,ADT,I,J)
          SUFTUUFTUUD(I,J)=CSFMUL(THSU,FTUU,FTUUD,I,J)
          SUFTUUSFTUUT(I,J)=CSFMUL(THSU,FTUUS,FTUUT,I,J)
          SUFTUUGTPUT(I,J)=CSFMUL(THSU,FTUU,GTPUT,I,J)
          SUFTUUGTSUT(I,J)=CSFMUL(THSU,FTUU,GTSUT,I,J)
          SUGTPUFTUUT(I,J)=CSFMUL(THSU,GTPU,FTUUT,I,J)
          SUGTPUSGTPUT(I,J)=CSFMUL(THSU,GTPUS,GTPUT,I,J)
          SUGTPUGTPUD(I,J)=CSFMUL(THSU,GTPU,GTPUD,I,J)
          SUGTSUGTSUD(I,J)=CSFMUL(THSU,GTSU,GTSUD,I,J)
          SUGTSUSGTSUT(I,J)=CSFMUL(THSU,GTSUS,GTSUT,I,J)
          SUGTSUFTUUT(I,J)=CSFMUL(THSU,GTSU,FTUUT,I,J)
          SUMTSFUSAUT(I,J)=CSFMUL(THSU,MTSFUS,AUT,I,J)
          SDFTDDFTDDD(I,J)=CSFMUL(THSD,FTDD,FTDDD,I,J)
          SDFTDDSFTDDT(I,J)=CSFMUL(THSD,FTDDS,FTDDT,I,J)
          SDFTDDGTPDT(I,J)=CSFMUL(THSD,FTDD,GTPDT,I,J)
          SDFTDDGTSDT(I,J)=CSFMUL(THSD,FTDD,GTSDT,I,J)
          SDGTPDFTDDT(I,J)=CSFMUL(THSD,GTPD,FTDDT,I,J)
          SDGTPDSGTPDT(I,J)=CSFMUL(THSD,GTPDS,GTPDT,I,J)
          SDGTPDGTPDD(I,J)=CSFMUL(THSD,GTPD,GTPDD,I,J)
          SDGTSDFTDDT(I,J)=CSFMUL(THSD,GTSD,FTDDT,I,J)
          SDGTSDGTSDD(I,J)=CSFMUL(THSD,GTSD,GTSDD,I,J)
          SDGTSDSGTSDT(I,J)=CSFMUL(THSD,GTSDS,GTSDT,I,J)
          SLFTELDFTEL(I,J)=CSFMUL(THSL,FTELD,FTEL,I,J)
          SLFTELTFTELS(I,J)=CSFMUL(THSL,FTELT,FTELS,I,J)
          SLFTELTGTL(I,J)=CSFMUL(THSL,FTELT,GTL,I,J)
          SLFTELTGTPL(I,J)=CSFMUL(THSL,FTELT,GTPL,I,J)
          SLGTLTGTLS(I,J)=CSFMUL(THSL,GTLT,GTLS,I,J)
          SLGTLTFTEL(I,J)=CSFMUL(THSL,GTLT,FTEL,I,J)
          SLGTLDGTL(I,J)=CSFMUL(THSL,GTLD,GTL,I,J)
          SLGTPLTGTPLS(I,J)=CSFMUL(THSL,GTPLT,GTPLS,I,J)
          SLGTPLDGTPL(I,J)=CSFMUL(THSL,GTPLD,GTPL,I,J)
          SLGTPLTFTEL(I,J)=CSFMUL(THSL,GTPLT,FTEL,I,J)
          SEMTSFESAET(I,J)=CSFMUL(THSE,MTSFES,AET,I,J)
          SEFTEEFTEED(I,J)=CSFMUL(THSE,FTEE,FTEED,I,J)
          SEFTEESFTEET(I,J)=CSFMUL(THSE,FTEES,FTEET,I,J)
          SEFTEEGTPET(I,J)=CSFMUL(THSE,FTEE,GTPET,I,J)
          SEGTPEFTEET(I,J)=CSFMUL(THSE,GTPE,FTEET,I,J)
          SEGTPESGTPET(I,J)=CSFMUL(THSE,GTPES,GTPET,I,J)
          SEGTPEGTPED(I,J)=CSFMUL(THSE,GTPE,GTPED,I,J)
        END DO
      END DO
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            SFUQDSFUQ(I,J)=CMATMUL(1,SFUQ,SFUQ,I,J)
            SFUURSFUURD(I,J)=CMATMUL(2,SFUUR,SFUUR,I,J)
            SFUQDCFDQ(I,J)=CMATMUL(1,SFUQ,CFDQ,I,J)
            CFDDRCFDDRD(I,J)=CMATMUL(2,CFDDR,CFDDR,I,J)
            CFDQDCFDQ(I,J)=CMATMUL(1,CFDQ,CFDQ,I,J)
            CFDQDSFUQ(I,J)=CMATMUL(1,CFDQ,SFUQ,I,J)
            CFEERCFEERD(I,J)=CMATMUL(2,CFEER,CFEER,I,J)
            CFELDCFEL(I,J)=CMATMUL(1,CFEL,CFEL,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SFUQTSFUQS(I,J)=SFUQDSFUQ(J,I)
            SFUURSSFUURT(I,J)=SFUURSFUURD(J,I)
            CFDDRSCFDDRT(I,J)=CFDDRCFDDRD(J,I)
            CFDQTCFDQS(I,J)=CFDQDCFDQ(J,I)
            CFEERSCFEERT(I,J)=CFEERCFEERD(J,I)
            CFELTCFELS(I,J)=CFELDCFEL(J,I)
          END DO
        END DO
      END IF
!
!These are the two loop terms. All in MV notation
!
      IF(SW2LP.EQ.1)THEN
!
        THDYDD=CTRACE(HDYDD)
        THEYED=CTRACE(HEYED)
        THUDYU=CTRACE(HUDYU)
        THDDYD=CTRACE(HDDYD)
        THEDYE=CTRACE(HEDYE)
        THUDHU=CTRACE(HUDHU)
!
        DO I=1,3
          DO J=1,3
            YUDYUMQM(I,J)=CMATMUL(0,YUDYU,MQM,I,J)
            YDDYDMQM(I,J)=CMATMUL(0,YDDYD,MQM,I,J)
            MQMYUDYU(I,J)=CMATMUL(0,MQM,YUDYU,I,J)
            MQMYDDYD(I,J)=CMATMUL(0,MQM,YDDYD,I,J)
            MLMYEDYE(I,J)=CMATMUL(0,MLM,YEDYE,I,J)
            YUDYUYDD(I,J)=CMATMUL(2,YUDYU,YD,I,J)
            YUMQMYDD(I,J)=CMATMUL(2,YUMQM,YD,I,J)
            YUYDDMDM(I,J)=CMATMUL(0,YUYDD,MDM,I,J)
            YDMQMYUD(I,J)=CMATMUL(2,YDMQM,YU,I,J)
            YDYUDMUPM(I,J)=CMATMUL(0,YDYUD,MUPM,I,J)
            YDMQMYUD(I,J)=CMATMUL(2,YDMQM,YU,I,J)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=CMATMUL(0,YUYUD,YUYUD,I,J)
            YDYDDYDYDD(I,J)=CMATMUL(0,YDYDD,YDYDD,I,J)
            YEYEDYEYED(I,J)=CMATMUL(0,YEYED,YEYED,I,J)
            YUDYUYUDYU(I,J)=CMATMUL(0,YUDYU,YUDYU,I,J)
            YDDYDYDDYD(I,J)=CMATMUL(0,YDDYD,YDDYD,I,J)
            YEDYEYEDYE(I,J)=CMATMUL(0,YEDYE,YEDYE,I,J)
            YDDYDYUDYU(I,J)=CMATMUL(0,YDDYD,YUDYU,I,J)
            YUDYUYDDYD(I,J)=CMATMUL(0,YUDYU,YDDYD,I,J)
            YUYDDYDYUD(I,J)=CMATMUL(0,YUYDD,YDYUD,I,J)
            YDYUDYUYDD(I,J)=CMATMUL(0,YDYUD,YUYDD,I,J)
!
            HUYUDYUYUD(I,J)=CMATMUL(0,HUYUD,YUYUD,I,J)
            HUYDDYDYUD(I,J)=CMATMUL(0,HUYDD,YDYUD,I,J)
            HDYUDYUYDD(I,J)=CMATMUL(0,HDYUD,YUYDD,I,J)
            HDYDDYDYDD(I,J)=CMATMUL(0,HDYDD,YDYDD,I,J)
            HEYEDYEYED(I,J)=CMATMUL(0,HEYED,YEYED,I,J)
            YUDYUYUDHU(I,J)=CMATMUL(0,YUDYU,YUDHU,I,J)
            YUDHUYUDYU(I,J)=CMATMUL(0,YUDHU,YUDYU,I,J)
            YDDYDYDDHD(I,J)=CMATMUL(0,YDDYD,YDDHD,I,J)
            YDDHDYDDYD(I,J)=CMATMUL(0,YDDHD,YDDYD,I,J)
            YDDYDYUDHU(I,J)=CMATMUL(0,YDDYD,YUDHU,I,J)
            YDDHDYUDYU(I,J)=CMATMUL(0,YDDHD,YUDYU,I,J)
            YUDHUYDDYD(I,J)=CMATMUL(0,YUDHU,YDDYD,I,J)
            YUDYUYDDHD(I,J)=CMATMUL(0,YUDYU,YDDHD,I,J)
            YEDYEYEDHE(I,J)=CMATMUL(0,YEDYE,YEDHE,I,J)
            YEDHEYEDYE(I,J)=CMATMUL(0,YEDHE,YEDYE,I,J)
            HUDHUYUDYU(I,J)=CMATMUL(0,HUDHU,YUDYU,I,J)
            HUDYUYUDHU(I,J)=CMATMUL(0,HUDYU,YUDHU,I,J)
            HDDHDYDDYD(I,J)=CMATMUL(0,HDDHD,YDDYD,I,J)
            HDDYDYDDHD(I,J)=CMATMUL(0,HDDYD,YDDHD,I,J)
            HDDHDYUDYU(I,J)=CMATMUL(0,HDDHD,YUDYU,I,J)
            YDDYDHUDHU(I,J)=CMATMUL(0,YDDYD,HUDHU,I,J)
            HDDYDYUDHU(I,J)=CMATMUL(0,HDDYD,YUDHU,I,J)
            YDDHDHUDYU(I,J)=CMATMUL(0,YDDHD,HUDYU,I,J)
            HUDHUYDDYD(I,J)=CMATMUL(0,HUDHU,YDDYD,I,J)
            YUDYUHDDHD(I,J)=CMATMUL(0,YUDYU,HDDHD,I,J)
            HUDYUYDDHD(I,J)=CMATMUL(0,HUDYU,YDDHD,I,J)
            YUDHUHDDYD(I,J)=CMATMUL(0,YUDHU,HDDYD,I,J)
            HEDHEYEDYE(I,J)=CMATMUL(0,HEDHE,YEDYE,I,J)
            HEDYEYEDHE(I,J)=CMATMUL(0,HEDYE,YEDHE,I,J)
            YUDYUHUDHU(I,J)=CMATMUL(0,YUDYU,HUDHU,I,J)
            YDDYDHDDHD(I,J)=CMATMUL(0,YDDYD,HDDHD,I,J)
            YUDHUHUDYU(I,J)=CMATMUL(0,YUDHU,HUDYU,I,J)
            YDDHDHDDYD(I,J)=CMATMUL(0,YDDHD,HDDYD,I,J)
            YEDYEHEDHE(I,J)=CMATMUL(0,YEDYE,HEDHE,I,J)
            YEDHEHEDYE(I,J)=CMATMUL(0,YEDHE,HEDYE,I,J)
            HUHUDYUYUD(I,J)=CMATMUL(0,HUHUD,YUYUD,I,J)
            YUYUDHUHUD(I,J)=CMATMUL(0,YUYUD,HUHUD,I,J)
            HUYUDYUHUD(I,J)=CMATMUL(0,HUYUD,YUHUD,I,J)
            YUHUDHUYUD(I,J)=CMATMUL(0,YUHUD,HUYUD,I,J)
            HUHDDYDYUD(I,J)=CMATMUL(0,HUHDD,YDYUD,I,J)
            YUYDDHDHUD(I,J)=CMATMUL(0,YUYDD,HDHUD,I,J)
            HUYDDYDHUD(I,J)=CMATMUL(0,HUYDD,YDHUD,I,J)
            YUHDDHDYUD(I,J)=CMATMUL(0,YUHDD,HDYUD,I,J)
            HDHDDYDYDD(I,J)=CMATMUL(0,HDHDD,YDYDD,I,J)
            YDYDDHDHDD(I,J)=CMATMUL(0,YDYDD,HDHDD,I,J)
            HDYDDYDHDD(I,J)=CMATMUL(0,HDYDD,YDHDD,I,J)
            YDHDDHDYDD(I,J)=CMATMUL(0,YDHDD,HDYDD,I,J)
            HDHUDYUYDD(I,J)=CMATMUL(0,HDHUD,YUYDD,I,J)
            YDYUDHUHDD(I,J)=CMATMUL(0,YDYUD,HUHDD,I,J)
            HDYUDYUHDD(I,J)=CMATMUL(0,HDYUD,YUHDD,I,J)
            YDHUDHUYDD(I,J)=CMATMUL(0,YDHUD,HUYDD,I,J)
            HEHEDYEYED(I,J)=CMATMUL(0,HEHED,YEYED,I,J)
            YEYEDHEHED(I,J)=CMATMUL(0,YEYED,HEHED,I,J)
            HEYEDYEHED(I,J)=CMATMUL(0,HEYED,YEHED,I,J)
            YEHEDHEYED(I,J)=CMATMUL(0,YEHED,HEYED,I,J)
!
            YUDMUPMYUYUDYU(I,J)=CMATMUL(0,YUDMUPMYU,YUDYU,I,J)
            YUDMUPMYUYDDYD(I,J)=CMATMUL(0,YUDMUPMYU,YDDYD,I,J)
            YUDYUMQMYDDYD(I,J)=CMATMUL(0,YUDYUMQM,YDDYD,I,J)
            YUDYUYDDMDMYD(I,J)=CMATMUL(0,YUDYUYDD,MDMYD,I,J)
            YDDMDMYDYDDYD(I,J)=CMATMUL(0,YDDMDMYD,YDDYD,I,J)
            YEDMEMYEYEDYE(I,J)=CMATMUL(0,YEDMEMYE,YEDYE,I,J)
            YUDYUMQMYUDYU(I,J)=CMATMUL(0,YUDYUMQM,YUDYU,I,J)
            YUDYUYUDMUPMYU(I,J)=CMATMUL(0,YUDYU,YUDMUPMYU,I,J)
            YUDYUYUDYUMQM(I,J)=CMATMUL(0,YUDYU,YUDYUMQM,I,J)
            YDDYDMQMYDDYD(I,J)=CMATMUL(0,YDDYDMQM,YDDYD,I,J)
            YDDYDYDDMDMYD(I,J)=CMATMUL(0,YDDYD,YDDMDMYD,I,J)
            YDDYDYDDYDMQM(I,J)=CMATMUL(0,YDDYD,YDDYDMQM,I,J)
            YEDYEMLMYEDYE(I,J)=CMATMUL(0,YEDYEMLM,YEDYE,I,J)
            YEDYEYEDMEMYE(I,J)=CMATMUL(0,YEDYE,YEDMEMYE,I,J)
            YEDYEYEDYEMLM(I,J)=CMATMUL(0,YEDYE,YEDYEMLM,I,J)
            YUMQMYUDYUYUD(I,J)=CMATMUL(0,YUMQMYUD,YUYUD,I,J)
            YUYUDMUPMYUYUD(I,J)=CMATMUL(0,YUYUDMUPM,YUYUD,I,J)
            YUYUDYUMQMYUD(I,J)=CMATMUL(0,YUYUD,YUMQMYUD,I,J)
            YUYUDYUYUDMUPM(I,J)=CMATMUL(0,YUYUD,YUYUDMUPM,I,J)
            YUMQMYDDYDYUD(I,J)=CMATMUL(0,YUMQMYDD,YDYUD,I,J)
            YUYDDMDMYDYUD(I,J)=CMATMUL(0,YUYDDMDM,YDYUD,I,J)
            YUYDDYDMQMYUD(I,J)=CMATMUL(0,YUYDD,YDMQMYUD,I,J)
            YUYDDYDYUDMUPM(I,J)=CMATMUL(0,YUYDD,YDYUDMUPM,I,J)
            YDMQMYDDYDYDD(I,J)=CMATMUL(0,YDMQMYDD,YDYDD,I,J)
            YDYDDMDMYDYDD(I,J)=CMATMUL(0,YDYDDMDM,YDYDD,I,J)
            YDYDDYDMQMYDD(I,J)=CMATMUL(0,YDYDD,YDMQMYDD,I,J)
            YDYDDYDYDDMDM(I,J)=CMATMUL(0,YDYDD,YDYDDMDM,I,J)
            YDMQMYUDYUYDD(I,J)=CMATMUL(0,YDMQMYUD,YUYDD,I,J)
            YDYUDMUPMYUYDD(I,J)=CMATMUL(0,YDYUDMUPM,YUYDD,I,J)
            YDYUDYUMQMYDD(I,J)=CMATMUL(0,YDYUD,YUMQMYDD,I,J)
            YDYUDYUYDDMDM(I,J)=CMATMUL(0,YDYUD,YUYDDMDM,I,J)
            YEMLMYEDYEYED(I,J)=CMATMUL(0,YEMLMYED,YEYED,I,J)
            YEYEDMEMYEYED(I,J)=CMATMUL(0,YEYEDMEM,YEYED,I,J)
            YEYEDYEMLMYED(I,J)=CMATMUL(0,YEYED,YEMLMYED,I,J)
            YEYEDYEYEDMEM(I,J)=CMATMUL(0,YEYED,YEYEDMEM,I,J)
          END DO
        END DO
!
        TYUDYUYUDYU=CTRACE(YUDYUYUDYU)
        TYUDYUYDDYD=CTRACE(YUDYUYDDYD)
        TYDDYDYUDYU=CTRACE(YDDYDYUDYU)
        TYDDYDYDDYD=CTRACE(YDDYDYDDYD)
        TYEDYEYEDYE=CTRACE(YEDYEYEDYE)
!
!These are SM terms for the two loop running below m_H
!I am going to use LYU for the MV notation SM Yukawa
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=CMATMUL(1,LYU,LYU,I,J)
              LYDDLYD(I,J)=CMATMUL(1,LYD,LYD,I,J)
              LYEDLYE(I,J)=CMATMUL(1,LYE,LYE,I,J)
            END DO
          END DO
          TLYUDLYU=CTRACE(LYUDLYU)
          TLYDDLYD=CTRACE(LYDDLYD)
          TLYEDLYE=CTRACE(LYEDLYE)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=CMATMUL(0,LYUDLYU,LYUDLYU,I,J)
              LYDDLYD2(I,J)=CMATMUL(0,LYDDLYD,LYDDLYD,I,J)
              LYEDLYE2(I,J)=CMATMUL(0,LYEDLYE,LYEDLYE,I,J)
              LYUDLYULYDDLYD(I,J)=CMATMUL(0,LYUDLYU,LYDDLYD,I,J)
              LYDDLYDLYUDLYU(I,J)=CMATMUL(0,LYDDLYD,LYUDLYU,I,J)
              DUMLUD1(I,J)=LYUDLYU(I,J)+LYDDLYD(I,J)
            END DO
          END DO
          TLYUDLYU2=CTRACE(LYUDLYU2)
          TLYDDLYD2=CTRACE(LYDDLYD2)
          TLYEDLYE2=CTRACE(LYEDLYE2)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=CMATMUL(0,DUMLUD1,LYDDLYD,I,J)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=CMATMUL(0,LYUDLYU2,LYUDLYU,I,J)
              LYDDLYD3(I,J)=CMATMUL(0,LYDDLYD2,LYDDLYD,I,J)
              LYEDLYE3(I,J)=CMATMUL(0,LYEDLYE2,LYEDLYE,I,J)
              DUMLUD(I,J)=CMATMUL(0,LYUDLYU,DUMLUD2,I,J)
            END DO
          END DO
          TLYUDLYU3=CTRACE(LYUDLYU3)
          TLYDDLYD3=CTRACE(LYDDLYD3)
          TLYEDLYE3=CTRACE(LYEDLYE3)
          TLYUDLYULYDDLYD=CTRACE(LYUDLYULYDDLYD)
          TDUMLUD=CTRACE(DUMLUD)
!
          Y2=3.D0*TLYUDLYU+3.D0*TLYDDLYD+TLYEDLYE
          H=3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
          Y4=(83.D0/40.D0*G(1)**2+27.D0/8.D0*G(2)**2
     $       +28.D0*G(3)**2)*TLYUDLYU+(-1.D0/40.D0*G(1)**2
     $       +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TLYDDLYD
     $       +(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TLYEDLYE
          CHI4=9.D0/4.D0*(3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
     $         -2.D0/3.D0*TLYUDLYULYDDLYD)
        END IF
      ELSE
        THDYDD=(0.D0,0.D0)
        THEYED=(0.D0,0.D0)
        THUDYU=(0.D0,0.D0)
        THDDYD=(0.D0,0.D0)
        THEDYE=(0.D0,0.D0)
        THUDHU=(0.D0,0.D0)
!
        DO I=1,3
          DO J=1,3
            YUDYUMQM(I,J)=(0.D0,0.D0)
            YDDYDMQM(I,J)=(0.D0,0.D0)
            MQMYUDYU(I,J)=(0.D0,0.D0)
            MQMYDDYD(I,J)=(0.D0,0.D0)
            MLMYEDYE(I,J)=(0.D0,0.D0)
            YUDYUYDD(I,J)=(0.D0,0.D0)
            YUMQMYDD(I,J)=(0.D0,0.D0)
            YUYDDMDM(I,J)=(0.D0,0.D0)
            YDMQMYUD(I,J)=(0.D0,0.D0)
            YDYUDMUPM(I,J)=(0.D0,0.D0)
            YDMQMYUD(I,J)=(0.D0,0.D0)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=(0.D0,0.D0)
            YDYDDYDYDD(I,J)=(0.D0,0.D0)
            YEYEDYEYED(I,J)=(0.D0,0.D0)
            YUDYUYUDYU(I,J)=(0.D0,0.D0)
            YDDYDYDDYD(I,J)=(0.D0,0.D0)
            YEDYEYEDYE(I,J)=(0.D0,0.D0)
            YDDYDYUDYU(I,J)=(0.D0,0.D0)
            YUDYUYDDYD(I,J)=(0.D0,0.D0)
            YUYDDYDYUD(I,J)=(0.D0,0.D0)
            YDYUDYUYDD(I,J)=(0.D0,0.D0)
!
            HUYUDYUYUD(I,J)=(0.D0,0.D0)
            HUYDDYDYUD(I,J)=(0.D0,0.D0)
            HDYUDYUYDD(I,J)=(0.D0,0.D0)
            HDYDDYDYDD(I,J)=(0.D0,0.D0)
            HEYEDYEYED(I,J)=(0.D0,0.D0)
            YUDYUYUDHU(I,J)=(0.D0,0.D0)
            YUDHUYUDYU(I,J)=(0.D0,0.D0)
            YDDYDYDDHD(I,J)=(0.D0,0.D0)
            YDDHDYDDYD(I,J)=(0.D0,0.D0)
            YDDYDYUDHU(I,J)=(0.D0,0.D0)
            YDDHDYUDYU(I,J)=(0.D0,0.D0)
            YUDHUYDDYD(I,J)=(0.D0,0.D0)
            YUDYUYDDHD(I,J)=(0.D0,0.D0)
            YEDYEYEDHE(I,J)=(0.D0,0.D0)
            YEDHEYEDYE(I,J)=(0.D0,0.D0)
            HUDHUYUDYU(I,J)=(0.D0,0.D0)
            HUDYUYUDHU(I,J)=(0.D0,0.D0)
            HDDHDYDDYD(I,J)=(0.D0,0.D0)
            HDDYDYDDHD(I,J)=(0.D0,0.D0)
            HDDHDYUDYU(I,J)=(0.D0,0.D0)
            YDDYDHUDHU(I,J)=(0.D0,0.D0)
            HDDYDYUDHU(I,J)=(0.D0,0.D0)
            YDDHDHUDYU(I,J)=(0.D0,0.D0)
            HUDHUYDDYD(I,J)=(0.D0,0.D0)
            YUDYUHDDHD(I,J)=(0.D0,0.D0)
            HUDYUYDDHD(I,J)=(0.D0,0.D0)
            YUDHUHDDYD(I,J)=(0.D0,0.D0)
            HEDHEYEDYE(I,J)=(0.D0,0.D0)
            HEDYEYEDHE(I,J)=(0.D0,0.D0)
            YUDYUHUDHU(I,J)=(0.D0,0.D0)
            YDDYDHDDHD(I,J)=(0.D0,0.D0)
            YUDHUHUDYU(I,J)=(0.D0,0.D0)
            YDDHDHDDYD(I,J)=(0.D0,0.D0)
            YEDYEHEDHE(I,J)=(0.D0,0.D0)
            YEDHEHEDYE(I,J)=(0.D0,0.D0)
            HUHUDYUYUD(I,J)=(0.D0,0.D0)
            YUYUDHUHUD(I,J)=(0.D0,0.D0)
            HUYUDYUHUD(I,J)=(0.D0,0.D0)
            YUHUDHUYUD(I,J)=(0.D0,0.D0)
            HUHDDYDYUD(I,J)=(0.D0,0.D0)
            YUYDDHDHUD(I,J)=(0.D0,0.D0)
            HUYDDYDHUD(I,J)=(0.D0,0.D0)
            YUHDDHDYUD(I,J)=(0.D0,0.D0)
            HDHDDYDYDD(I,J)=(0.D0,0.D0)
            YDYDDHDHDD(I,J)=(0.D0,0.D0)
            HDYDDYDHDD(I,J)=(0.D0,0.D0)
            YDHDDHDYDD(I,J)=(0.D0,0.D0)
            HDHUDYUYDD(I,J)=(0.D0,0.D0)
            YDYUDHUHDD(I,J)=(0.D0,0.D0)
            HDYUDYUHDD(I,J)=(0.D0,0.D0)
            YDHUDHUYDD(I,J)=(0.D0,0.D0)
            HEHEDYEYED(I,J)=(0.D0,0.D0)
            YEYEDHEHED(I,J)=(0.D0,0.D0)
            HEYEDYEHED(I,J)=(0.D0,0.D0)
            YEHEDHEYED(I,J)=(0.D0,0.D0)
!
            YUDMUPMYUYUDYU(I,J)=(0.D0,0.D0)
            YUDMUPMYUYDDYD(I,J)=(0.D0,0.D0)
            YUDYUMQMYDDYD(I,J)=(0.D0,0.D0)
            YUDYUYDDMDMYD(I,J)=(0.D0,0.D0)
            YDDMDMYDYDDYD(I,J)=(0.D0,0.D0)
            YEDMEMYEYEDYE(I,J)=(0.D0,0.D0)
            YUDYUMQMYUDYU(I,J)=(0.D0,0.D0)
            YUDYUYUDMUPMYU(I,J)=(0.D0,0.D0)
            YUDYUYUDYUMQM(I,J)=(0.D0,0.D0)
            YDDYDMQMYDDYD(I,J)=(0.D0,0.D0)
            YDDYDYDDMDMYD(I,J)=(0.D0,0.D0)
            YDDYDYDDYDMQM(I,J)=(0.D0,0.D0)
            YEDYEMLMYEDYE(I,J)=(0.D0,0.D0)
            YEDYEYEDMEMYE(I,J)=(0.D0,0.D0)
            YEDYEYEDYEMLM(I,J)=(0.D0,0.D0)
            YUMQMYUDYUYUD(I,J)=(0.D0,0.D0)
            YUYUDMUPMYUYUD(I,J)=(0.D0,0.D0)
            YUYUDYUMQMYUD(I,J)=(0.D0,0.D0)
            YUYUDYUYUDMUPM(I,J)=(0.D0,0.D0)
            YUMQMYDDYDYUD(I,J)=(0.D0,0.D0)
            YUYDDMDMYDYUD(I,J)=(0.D0,0.D0)
            YUYDDYDMQMYUD(I,J)=(0.D0,0.D0)
            YUYDDYDYUDMUPM(I,J)=(0.D0,0.D0)
            YDMQMYDDYDYDD(I,J)=(0.D0,0.D0)
            YDYDDMDMYDYDD(I,J)=(0.D0,0.D0)
            YDYDDYDMQMYDD(I,J)=(0.D0,0.D0)
            YDYDDYDYDDMDM(I,J)=(0.D0,0.D0)
            YDMQMYUDYUYDD(I,J)=(0.D0,0.D0)
            YDYUDMUPMYUYDD(I,J)=(0.D0,0.D0)
            YDYUDYUMQMYDD(I,J)=(0.D0,0.D0)
            YDYUDYUYDDMDM(I,J)=(0.D0,0.D0)
            YEMLMYEDYEYED(I,J)=(0.D0,0.D0)
            YEYEDMEMYEYED(I,J)=(0.D0,0.D0)
            YEYEDYEMLMYED(I,J)=(0.D0,0.D0)
            YEYEDYEYEDMEM(I,J)=(0.D0,0.D0)
          END DO
        END DO
!
        TYUDYUYUDYU=(0.D0,0.D0)
        TYUDYUYDDYD=(0.D0,0.D0)
        TYDDYDYUDYU=(0.D0,0.D0)
        TYDDYDYDDYD=(0.D0,0.D0)
        TYEDYEYEDYE=(0.D0,0.D0)
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=(0.D0,0.D0)
              LYD(I,J)=(0.D0,0.D0)
              LYE(I,J)=(0.D0,0.D0)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=(0.D0,0.D0)
              LYDDLYD(I,J)=(0.D0,0.D0)
              LYEDLYE(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU=(0.D0,0.D0)
          TLYDDLYD=(0.D0,0.D0)
          TLYEDLYE=(0.D0,0.D0)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=(0.D0,0.D0)
              LYDDLYD2(I,J)=(0.D0,0.D0)
              LYEDLYE2(I,J)=(0.D0,0.D0)
              LYUDLYULYDDLYD(I,J)=(0.D0,0.D0)
              LYDDLYDLYUDLYU(I,J)=(0.D0,0.D0)
              DUMLUD1(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU2=(0.D0,0.D0)
          TLYDDLYD2=(0.D0,0.D0)
          TLYEDLYE2=(0.D0,0.D0)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=(0.D0,0.D0)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=(0.D0,0.D0)
              LYDDLYD3(I,J)=(0.D0,0.D0)
              LYEDLYE3(I,J)=(0.D0,0.D0)
              DUMLUD(I,J)=(0.D0,0.D0)
            END DO
          END DO
          TLYUDLYU3=(0.D0,0.D0)
          TLYDDLYD3=(0.D0,0.D0)
          TLYEDLYE3=(0.D0,0.D0)
          TLYUDLYULYDDLYD=(0.D0,0.D0)
          TDUMLUD=(0.D0,0.D0)
!
          Y2=(0.D0,0.D0)
          H=(0.D0,0.D0)
          Y4=(0.D0,0.D0)
          CHI4=(0.D0,0.D0)
        END IF
      END IF
!
!Now I calculate the runnng of the MSSM gauge, yukawas and soft terms
!which will be kept intact for two loop running of fs, SUSY couplings
!and tildas.
!
      DO I=1,3
        SUM=(0.D0,0.D0)
        DO J=1,3
          SUM=SUM+B2LPM(I,J)*G(290+J)**2
        END DO
        F(290+I)=G(290+I)**3/16.D0/PI**2*(B1LPM(I)
     $           +DBLE(SW2LP)/16.D0/PI**2*(SUM-(CM(I,1)*TYUDYU
     $           +CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        B2GM(I)=G(290+I)**3*DBLE(SW2LP)*(SUM-(CM(I,1)*TYUDYU
     $          +CM(I,2)*TYDDYD+CM(I,3)*TYEDYE))
      END DO
!
      DO I=1,3
        DO J=1,3
          DUM1U1(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1D1(I,J)=DUM1U1(I,J)
          DUM1E1(I,J)=DUM1U1(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U2(I,J)=DUM1D1(I,J)
            DUM2D1(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)+
     $                  YEYEDYEYED(I,J)
            DUM2D2(I,J)=DUM1D1(I,J)
            DUM2E1(I,J)=DUM2D1(I,J)
            DUM2E2(I,J)=DUM2U2(I,J)
          ELSE
            DUM2U1(I,J)=(0.D0,0.D0)
            DUM2U2(I,J)=(0.D0,0.D0)
            DUM2D1(I,J)=(0.D0,0.D0)
            DUM2D2(I,J)=(0.D0,0.D0)
            DUM2E1(I,J)=(0.D0,0.D0)
            DUM2E2(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          DUM1U(I,J)=3.D0*TYUDYU*ID(I,J)+3.D0*YUDYU(I,J)+YDDYD(I,J)
     $               +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2
     $               -13.D0/15.D0*G(291)**2)*ID(I,J)              
          DUM1D(I,J)=CTRACE(DUM1D1)*ID(I,J)+3.D0*YDDYD(I,J)+YUDYU(I,J)
     $               +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2
     $               -7.D0/15.D0*G(291)**2)*ID(I,J)
          DUM1E(I,J)=CTRACE(DUM1E1)*ID(I,J)+3.D0*YEDYE(I,J)+(-3.D0
     $               *G(292)**2-9.D0/5.D0*G(291)**2)*ID(I,J)
!
!Here are the two loop terms
!
          IF(SW2LP.EQ.1)THEN
            DUM2U(I,J)=-3.D0*CTRACE(DUM2U1)*ID(I,J)-YDDYD(I,J)
     $                 *CTRACE(DUM2U2)-9.D0*YUDYU(I,J)*TYUDYU
     $                 -4.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)-2.D0
     $                 *YDDYDYUDYU(I,J)+(16.D0*G(293)**2+4.D0/5.D0
     $                 *G(291)**2)*TYUDYU*ID(I,J)+(6.D0*G(292)**2
     $                 +2.D0/5.D0*G(291)**2)*YUDYU(I,J)+2.D0/5.D0
     $                 *G(291)**2*YDDYD(I,J)+(-16.D0/9.D0*G(293)**4
     $                 +8.D0*G(293)**2*G(292)**2+136.D0/45.D0
     $                 *G(293)**2*G(291)**2+15.D0/2.D0*G(292)**4
     $                 +G(292)**2*G(291)**2+2743.D0/450.D0
     $                 *G(291)**4)*ID(I,J)
            DUM2D(I,J)=-3.D0*CTRACE(DUM2D1)*ID(I,J)-3.D0*YUDYU(I,J)
     $                 *TYUDYU-3.D0*YDDYD(I,J)*CTRACE(DUM2D2)
     $                 -4.D0*YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)
     $                 -2.D0*YUDYUYDDYD(I,J)+(16.D0*G(293)**2
     $                 -2.D0/5.D0*G(291)**2)*TYDDYD*ID(I,J)+6.D0/5.D0
     $                 *G(291)**2*TYEDYE*ID(I,J)+4.D0/5.D0*G(291)**2
     $                 *YUDYU(I,J)+(6.D0*G(292)**2+4.D0/5.D0
     $                 *G(291)**2)*YDDYD(I,J)+(-16.D0/9.D0*G(293)**4
     $                 +8.D0*G(293)**2*G(292)**2+8.D0/9.D0
     $                 *G(293)**2*G(291)**2+15.D0/2.D0*G(292)**4
     $                 +G(292)**2*G(291)**2+287.D0/90.D0
     $                 *G(291)**4)*ID(I,J)
            DUM2E(I,J)=-3.D0*CTRACE(DUM2E1)*ID(I,J)-3.D0*YEDYE(I,J)
     $                 *CTRACE(DUM2E2)-4.D0*YEDYEYEDYE(I,J)+(16.D0
     $                 *G(293)**2-2.D0/5.D0*G(291)**2)*TYDDYD
     $                 *ID(I,J)+6.D0/5.D0*G(291)**2*TYEDYE*ID(I,J)
     $                 +6.D0*G(292)**2*YEDYE(I,J)+(15.D0/2.D0
     $                 *G(292)**4+9.D0/5.D0*G(292)**2*G(291)**2
     $                 +27.D0/2.D0*G(291)**4)*ID(I,J)
          ELSE
            DUM2U(I,J)=(0.D0,0.D0)
            DUM2D(I,J)=(0.D0,0.D0)
            DUM2E(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the Yukawas
!
      DO I=1,3
        DO J=1,3
          B1YMU(I,J)=CMATMUL(0,YU,DUM1U,I,J)
          B1YMD(I,J)=CMATMUL(0,YD,DUM1D,I,J)
          B1YME(I,J)=CMATMUL(0,YE,DUM1E,I,J)
          IF(SW2LP.EQ.1)THEN
            B2YMU(I,J)=CMATMUL(0,YU,DUM2U,I,J)
            B2YMD(I,J)=CMATMUL(0,YD,DUM2D,I,J)
            B2YME(I,J)=CMATMUL(0,YE,DUM2E,I,J)
          END IF
        END DO
      END DO
!
!Convert to differentials
!
      DO I=1,3
        DO J=1,3
        F(293+(I-1)*3+J)=1.D0/16.D0/PI**2*B1YMU(J,I)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I)
        F(302+(I-1)*3+J)=1.D0/16.D0/PI**2*B1YMD(J,I)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMD(J,I)
        F(311+(I-1)*3+J)=1.D0/16.D0/PI**2*B1YME(J,I)
     $                   +1.D0/(16.D0*PI**2)**2*B2YME(J,I)
        END DO
      END DO
!
!Similarly, the MSSM Gaugino mass parameters are:
!
      DO I=1,3
        SUM=(0.D0,0.D0)
        DO J=1,3
          SUM=SUM+B2LPM(I,J)*G(290+J)**2*(G(320+I)+G(320+J))
        END DO
        IF(SW2LP.EQ.1)THEN
          B2M(I)=2.D0*G(290+I)**2*(SUM+CM(I,1)*(TYUDHU-G(320+I)
     $           *TYUDYU)+CM(I,2)*(TYDDHD-G(320+I)*TYDDYD)+CM(I,3)
     $           *(TYEDHE-G(320+I)*TYEDYE))
        END IF
        F(320+I)=2.D0*G(290+I)**2/16.D0/PI**2*B1LPM(I)*G(320+I)
     $          +1.D0/(16.D0*PI**2)**2*B2M(I)
      END DO
!
!And here is the MSSM \mu
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1B1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1B2(I,J)=6.D0*HUYUD(I,J)+6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2B1(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                  +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
            DUM2B2(I,J)=3.D0*HUYUDYUYUD(I,J)+3.D0*HDYDDYDYDD(I,J)
     $                  +HUYDDYDYUD(I,J)+HDYUDYUYDD(I,J)
     $                  +HEYEDYEYED(I,J)
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2B1(I,J)=(0.D0,0.D0)
            DUM2B2(I,J)=(0.D0,0.D0)
            DUM2GRKMU(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
        BGRKMU=MVMU*(CTRACE(DUMGRKMU1)-3.D0*G(292)**2-3.D0/5.D0
     $           *G(291)**2)
        BETA1B=MVB*(CTRACE(DUM1B1)-3.D0*G(292)**2-3.D0/5.D0
     $         *G(291)**2)+MVMU*(CTRACE(DUM1B2)+6.D0*G(292)**2
     $         *G(322)+6.D0/5.D0*G(291)**2*G(321))
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMU*(-3.D0*CTRACE(DUM2GRKMU)+(16.D0*G(293)**2+4.D0
     $             /5.D0*G(291)**2)*TYUDYU+(16.D0*G(293)**2-2.D0/5.D0
     $             *G(291)**2)*TYDDYD+6.D0/5.D0*G(291)**2*TYEDYE
     $             +15.D0/2.D0*G(292)**4+9.D0/5.D0*G(291)**2*G(292)**2
     $             +207.D0/50.D0*G(291)**4)
        B2GRKMUM=BETA2GRKMU
        BETA2B=MVB*(-3.D0*CTRACE(DUM2B1)+(16.D0*G(293)**2+4.D0/5.D0
     $         *G(291)**2)*TYUDYU+(16.D0*G(293)**2-2.D0/5.D0*G(291)**2)
     $         *TYDDYD+6.D0/5.D0*G(291)**2*TYEDYE+15.D0/2.D0*G(292)**4
     $         +9.D0/5.D0*G(291)**2*G(292)**2+207.D0/50.D0*G(291)**4)
     $         +MVMU*(-12.D0*CTRACE(DUM2B2)+(32.D0*G(293)**2
     $         +8.D0/5.D0*G(291)**2)*THUYUD+(32.D0*G(293)**2-4.D0/5.D0
     $         *G(291)**2)*THDYDD+12.D0/5.D0*G(291)**2*THEYED-(32.D0
     $         *G(293)**2*G(323)+8.D0/5.D0*G(291)**2*G(321))*TYUDYU
     $         -(32.D0*G(293)**2*G(323)-4.D0/5.D0*G(291)**2*G(321))
     $         *TYDDYD-12.D0/5.D0*G(291)**2*G(321)*TYEDYE-30.D0
     $         *G(292)**4*G(322)-18.D0/5.D0*G(291)**2*G(292)**2
     $         *(G(321)+G(322))-414.D0/25.D0*G(291)**4*G(321))
         BETA2BM=BETA2B
      END IF
!
!The RKSTP compatible derivative is (with a minus sign
!to convert the running to BT notation)
!
        F(398)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $         *BETA2GRKMU
        F(399)=-1.D0/16.D0/PI**2*BETA1B-1.D0/(16.D0*PI**2)**2*BETA2B
!
!Now the MSSM trilinear couplings - it should be ok to use the same
!dummy matrices, but the betas must be reset
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=(0.D0,0.D0)
          BETA2U(I,J)=(0.D0,0.D0)
          BETA1D(I,J)=(0.D0,0.D0)
          BETA2D(I,J)=(0.D0,0.D0)
          BETA1E(I,J)=(0.D0,0.D0)
          BETA2E(I,J)=(0.D0,0.D0)
        END DO
      END DO
!
!DUM2U11 is a dummy matrix for the calculation of DUM2U1
!
      DO I=1,3
        DO J=1,3
          DUM1D11(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1D12(I,J)=6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
          DUM1E11(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1E21(I,J)=6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2U11(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U12(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2U21(I,J)=6.D0*HUYUDYUYUD(I,J)+HUYDDYDYUD(I,J)
     $                   +HDYUDYUYDD(I,J)
            DUM2U22(I,J)=6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
            DUM2U23(I,J)=6.D0*YDYDD(I,J)+2.D0*YEYED(I,J)
            DUM2D11(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)
     $                   +YEYEDYEYED(I,J)
            DUM2D12(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2D21(I,J)=6.D0*HDYDDYDYDD(I,J)+HUYDDYDYUD(I,J)
     $                   +HDYUDYUYDD(I,J)+2.D0*HEYEDYEYED(I,J)
            DUM2D22(I,J)=3.D0*HDYDD(I,J)+HEYED(I,J)
            DUM2D23(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E11(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)
     $                   +YEYEDYEYED(I,J)
            DUM2E12(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E21(I,J)=6.D0*HDYDDYDYDD(I,J)+HUYDDYDYUD(I,J)
     $                   +HDYUDYUYDD(I,J)+2.D0*HEYEDYEYED(I,J)
            DUM2E22(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E23(I,J)=3.D0*HDYDD(I,J)+HEYED(I,J)
          ELSE
            DUM2U11(I,J)=(0.D0,0.D0)
            DUM2U12(I,J)=(0.D0,0.D0)
            DUM2U21(I,J)=(0.D0,0.D0)
            DUM2U22(I,J)=(0.D0,0.D0)
            DUM2U23(I,J)=(0.D0,0.D0)
            DUM2D11(I,J)=(0.D0,0.D0)
            DUM2D12(I,J)=(0.D0,0.D0)
            DUM2D21(I,J)=(0.D0,0.D0)
            DUM2D22(I,J)=(0.D0,0.D0)
            DUM2D23(I,J)=(0.D0,0.D0)
            DUM2E11(I,J)=(0.D0,0.D0)
            DUM2E12(I,J)=(0.D0,0.D0)
            DUM2E21(I,J)=(0.D0,0.D0)
            DUM2E22(I,J)=(0.D0,0.D0)
            DUM2E23(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!DUM1U1 is a dummy matrix for the calculation of BETA1U
!
      DO I=1,3
        DO J=1,3
          DUM1U1(I,J)=3.D0*TYUDYU*ID(I,J)+5.D0*YUDYU(I,J)+YDDYD(I,J)
     $                +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2
     $                -13.D0/15.D0*G(291)**2)*ID(I,J)
          DUM1U2(I,J)=6.D0*THUYUD*ID(I,J)+4.D0*YUDHU(I,J)+2.D0
     $                *YDDHD(I,J)+(32.D0/3.D0*G(293)**2*G(323)+6.D0
     $                *G(292)**2*G(322)+26.D0/15.D0*G(291)**2*G(321))
     $                *ID(I,J)
          DUM1D1(I,J)=CTRACE(DUM1D11)*ID(I,J)+5.D0*YDDYD(I,J)+YUDYU(I,J)
     $                +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2-7.D0/15.D0
     $                *G(291)**2)*ID(I,J)
          DUM1D2(I,J)=CTRACE(DUM1D12)*ID(I,J)+4.D0*YDDHD(I,J)+2.D0
     $                *YUDHU(I,J)+(32.D0/3.D0*G(293)**2*G(323)+6.D0
     $                *G(292)**2*G(322)+14.D0/15.D0*G(291)**2*G(321))
     $                *ID(I,J)
          DUM1E1(I,J)=CTRACE(DUM1E11)*ID(I,J)+5.D0*YEDYE(I,J)+(-3.D0
     $                *G(292)**2-9.D0/5.D0*G(291)**2)*ID(I,J)
          DUM1E2(I,J)=CTRACE(DUM1E21)*ID(I,J)+4.D0*YEDHE(I,J)+(6.D0
     $                *G(292)**2*G(322)+18.D0/5.D0*G(291)**2*G(321))
     $                *ID(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=-3.D0*CTRACE(DUM2U11)*ID(I,J)-YDDYD(I,J)
     $                  *CTRACE(DUM2U12)-15.D0*YUDYU(I,J)*TYUDYU
     $                  -6.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)
     $                  -4.D0*YDDYDYUDYU(I,J)+(16.D0*G(293)**2+4.D0/5.D0
     $                  *G(291)**2)*TYUDYU*ID(I,J)+12.D0*G(292)**2
     $                  *YUDYU(I,J)+2.D0/5.D0*G(291)**2*YDDYD(I,J)
     $                  +(-16.D0/9.D0*G(293)**4+8.D0*G(293)**2
     $                  *G(292)**2+136.D0/45.D0*G(293)**2*G(291)**2
     $                  +15.D0/2.D0*G(292)**4+G(292)**2*G(291)**2
     $                  +2743.D0/450.D0*G(291)**4)*ID(I,J)
            DUM2U2(I,J)=-6.D0*CTRACE(DUM2U21)*ID(I,J)-18.D0
     $                  *YUDYU(I,J)*THUYUD-YDDYD(I,J)*CTRACE(DUM2U22)
     $                  -12.D0*YUDHU(I,J)*TYUDYU-YDDHD(I,J)
     $                  *CTRACE(DUM2U23)-6.D0*YUDYUYUDHU(I,J)
     $                  -8.D0*YUDHUYUDYU(I,J)-4.D0*YDDYDYDDHD(I,J)-4.D0
     $                  *YDDHDYDDYD(I,J)-2.D0*YDDYDYUDHU(I,J)-4.D0
     $                  *YDDHDYUDYU(I,J)+(32.D0*G(293)**2+8.D0/5.D0
     $                  *G(291)**2)*THUYUD*ID(I,J)+(6.D0*G(292)**2
     $                  +6.D0/5.D0*G(291)**2)*YUDHU(I,J)+4.D0/5.D0
     $                  *G(291)**2*YDDHD(I,J)-(32.D0*G(293)**2*G(323)
     $                  +8.D0/5.D0*G(291)**2*G(321))*TYUDYU*ID(I,J)
     $                  -(12.D0*G(292)**2*G(322)+4.D0/5.D0*G(291)**2
     $                  *G(321))*YUDYU(I,J)-4.D0/5.D0*G(291)**2
     $                  *G(321)*YDDYD(I,J)+(64.D0/9.D0*G(293)**4
     $                  *G(323)-16.D0*G(293)**2*G(292)**2*(G(323)
     $                  +G(322))-272.D0/45.D0*G(293)**2*G(291)**2
     $                  *(G(323)+G(321))-30.D0*G(292)**4*G(322)-2.D0
     $                  *G(292)**2*G(291)**2*(G(322)+G(321))
     $                  -5486.D0/225.D0*G(291)**4*G(321))*ID(I,J)      
            DUM2D1(I,J)=-3.D0*CTRACE(DUM2D11)*ID(I,J)-3.D0*YUDYU(I,J)
     $                  *TYUDYU-5.D0*YDDYD(I,J)*CTRACE(DUM2D12)-6.D0
     $                  *YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)-4.D0
     $                  *YUDYUYDDYD(I,J)+(16.D0*G(293)**2-2.D0/5.D0
     $                  *G(291)**2)*TYDDYD*ID(I,J)+6.D0/5.D0*G(291)**2
     $                  *TYEDYE*ID(I,J)+4.D0/5.D0*G(291)**2*YUDYU(I,J)
     $                  +(12.D0*G(292)**2+6.D0/5.D0*G(291)**2)
     $                  *YDDYD(I,J)+(-16.D0/9.D0*G(293)**4+8.D0
     $                  *G(293)**2*G(292)**2+8.D0/9.D0*G(293)**2
     $                  *G(291)**2+15.D0/2.D0*G(292)**4+G(292)**2
     $                  *G(291)**2+287.D0/90.D0*G(291)**4)*ID(I,J)
            DUM2D2(I,J)=-6.D0*CTRACE(DUM2D21)*ID(I,J)-6.D0*YUDYU(I,J)
     $                  *THUYUD-6.D0*YDDYD(I,J)*CTRACE(DUM2D22)-6.D0
     $                  *YUDHU(I,J)*TYUDYU-4.D0*YDDHD(I,J)
     $                  *CTRACE(DUM2D23)-6.D0*YDDYDYDDHD(I,J)-8.D0
     $                  *YDDHDYDDYD(I,J)-4.D0*YUDHUYUDYU(I,J)-4.D0
     $                  *YUDYUYUDHU(I,J)-4.D0*YUDHUYDDYD(I,J)-2.D0
     $                  *YUDYUYDDHD(I,J)+(32.D0*G(293)**2-4.D0/5.D0
     $                  *G(291)**2)*THDYDD*ID(I,J)+12.D0/5.D0*G(291)**2
     $                  *THEYED*ID(I,J)+8.D0/5.D0*G(291)**2*YUDHU(I,J)
     $                  +(6.D0*G(292)**2+6.D0/5.D0*G(291)**2)*YDDHD(I,J)
     $                  -(32.D0*G(293)**2*G(323)-4.D0/5.D0*G(291)**2
     $                  *G(321))*TYDDYD*ID(I,J)-12.D0/5.D0*G(291)**2
     $                  *G(321)*TYEDYE*ID(I,J)-(12.D0*G(292)**2*G(322)
     $                  +8.D0/5.D0*G(291)**2*G(321))*YDDYD(I,J)
     $                  -8.D0/5.D0*G(291)**2*G(321)*YUDYU(I,J)
     $                  +(64.D0/9.D0*G(293)**4*G(323)-16.D0*G(293)**2
     $                  *G(292)**2*(G(323)+G(322))-16.D0/9.D0
     $                  *G(293)**2*G(291)**2*(G(323)+G(321))-30.D0
     $                  *G(292)**4*G(322)-2.D0*G(292)**2*G(291)**2
     $                  *(G(322)+G(321))-574.D0/45.D0*G(291)**4*G(321))
     $                  *ID(I,J)
            DUM2E1(I,J)=-3.D0*CTRACE(DUM2E11)*ID(I,J)-5.D0*YEDYE(I,J)
     $                  *CTRACE(DUM2E12)-6.D0*YEDYEYEDYE(I,J)+(16.D0
     $                  *G(293)**2-2.D0/5.D0*G(291)**2)*TYDDYD*ID(I,J)
     $                  +6.D0/5.D0*G(291)**2*TYEDYE*ID(I,J)+(12.D0
     $                  *G(292)**2-6.D0/5.D0*G(291)**2)*YEDYE(I,J)
     $                  +(15.D0/2.D0*G(292)**4+9.D0/5.D0*G(292)**2
     $                  *G(291)**2+27.D0/2.D0*G(291)**4)*ID(I,J)
            DUM2E2(I,J)=-6.D0*CTRACE(DUM2E21)*ID(I,J)-4.D0*YEDHE(I,J)
     $                  *CTRACE(DUM2E22)-6.D0*YEDYE(I,J)*CTRACE(DUM2E23)
     $                  -6.D0*YEDYEYEDHE(I,J)-8.D0*YEDHEYEDYE(I,J)
     $                  +(32.D0*G(293)**2-4.D0/5.D0*G(291)**2)*THDYDD
     $                  *ID(I,J)+12.D0/5.D0*G(291)**2*THEYED*ID(I,J)
     $                  +(6.D0*G(292)**2+6.D0/5.D0*G(291)**2)*YEDHE(I,J)
     $                  -(32.D0*G(293)**2*G(323)-4.D0/5.D0*G(291)**2
     $                  *G(321))*TYDDYD*ID(I,J)-12.D0/5.D0*G(291)**2
     $                  *G(321)*TYEDYE*ID(I,J)-12.D0*G(292)**2*G(322)
     $                  *YEDYE(I,J)+(-30.D0*G(292)**4*G(322)-18.D0/5.D0
     $                  *G(292)**2*G(291)**2*(G(321)+G(322))-54.D0
     $                  *G(291)**4*G(321))*ID(I,J)
          ELSE
            DUM2U1(I,J)=(0.D0,0.D0)
            DUM2U2(I,J)=(0.D0,0.D0)
            DUM2D1(I,J)=(0.D0,0.D0)
            DUM2D2(I,J)=(0.D0,0.D0)
            DUM2E1(I,J)=(0.D0,0.D0)
            DUM2E2(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the trilinears
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=CMATMUL(0,HU,DUM1U1,I,J)+CMATMUL(0,YU,DUM1U2,I,J)
          BETA1D(I,J)=CMATMUL(0,HD,DUM1D1,I,J)+CMATMUL(0,YD,DUM1D2,I,J)
          BETA1E(I,J)=CMATMUL(0,HE,DUM1E1,I,J)+CMATMUL(0,YE,DUM1E2,I,J)
          IF(SW2LP.EQ.1)THEN
            BETA2U(I,J)=CMATMUL(0,HU,DUM2U1,I,J)
     $                  +CMATMUL(0,YU,DUM2U2,I,J)
            BETA2D(I,J)=CMATMUL(0,HD,DUM2D1,I,J)
     $                  +CMATMUL(0,YD,DUM2D2,I,J)
            BETA2E(I,J)=CMATMUL(0,HE,DUM2E1,I,J)
     $                  +CMATMUL(0,YE,DUM2E2,I,J)
            B2HMU(I,J)=BETA2U(I,J)
            B2HMD(I,J)=BETA2D(I,J)
            B2HME(I,J)=BETA2E(I,J)
          END IF
!
!Calculate the differentials DH
!
          DH(1,I,J)=1.D0/16.D0/PI**2*BETA1U(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2U(I,J)
          DH(2,I,J)=1.D0/16.D0/PI**2*BETA1D(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2D(I,J)
          DH(3,I,J)=1.D0/16.D0/PI**2*BETA1E(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2E(I,J)
        END DO
      END DO
!
!Convert into form readable by RKSTP (converting to BT notation...)
!Need new do loop due to the (J,I) ordering
!
      DO I=1,3
        DO J=1,3
          F(323+(I-1)*3+J)=-DH(1,J,I)
          F(332+(I-1)*3+J)=-DH(2,J,I)
          F(341+(I-1)*3+J)=-DH(3,J,I)
        END DO
      END DO
!
!Now for the (mass)^2 parameters. DUMSP are dummy matrices for the
!calculation of S'
!
      IF(SW2LP.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            DUMSP1(I,J)=3.D0*G(351)*ID(I,J)+MQM(I,J)
            DUMSP2(I,J)=3.D0*G(352)*ID(I,J)-MQM(I,J)
            DUMSP3(I,J)=G(352)*ID(I,J)+MLM(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            DUMSP1(I,J)=(0.D0,0.D0)
            DUMSP2(I,J)=(0.D0,0.D0)
            DUMSP3(I,J)=(0.D0,0.D0)
          END DO
        END DO
      END IF
!
!Dummy matrices for S, S' and SIGMAi
!
      DO I=1,3
        DO J=1,3
          DUMS(I,J)=MQM(I,J)-MLM(I,J)-2.D0*MUPM(I,J)+MDM(I,J)+MEM(I,J)
          IF(SW2LP.EQ.1)THEN
            DUMSP(I,J)=-CMATMUL(0,DUMSP1,YUDYU,I,J)+4.D0*YUDMUPMYU(I,J)
     $                 +CMATMUL(0,DUMSP2,YDDYD,I,J)-2.D0*YDDMDMYD(I,J)
     $                 +CMATMUL(0,DUMSP3,YEDYE,I,J)-2.D0*YEDMEMYE(I,J)
            DUMSIG1(I,J)=MQM(I,J)+3.D0*MLM(I,J)+8.D0*MUPM(I,J)
     $                   +2.D0*MDM(I,J)+6.D0*MEM(I,J)
            DUMSIG2(I,J)=3.D0*MQM(I,J)+MLM(I,J)
            DUMSIG3(I,J)=2.D0*MQM(I,J)+MUPM(I,J)+MDM(I,J)
          ELSE
            DUMSP(I,J)=(0.D0,0.D0)
            DUMSIG1(I,J)=(0.D0,0.D0)
            DUMSIG2(I,J)=(0.D0,0.D0)
            DUMSIG3(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!These definitions make the equations shorter
!
      S=G(351)-G(352)+CTRACE(DUMS)
      IF(SW2LP.EQ.1)THEN
        SP=CTRACE(DUMSP)+(3.D0/2.D0*G(292)**2+3.D0/10.D0*G(291)**2)
     $     *(G(351)-G(352)-CTRACE(MLM))+(8.D0/3.D0*G(293)**2+3.D0/2.D0
     $     *G(292)**2+1.D0/30.D0*G(291)**2)*CTRACE(MQM)-(16.D0/3.D0
     $     *G(293)**2+16.D0/15.D0*G(291)**2)*CTRACE(MUPM)+(8.D0/3.D0
     $     *G(293)**2+2.D0/15.D0*G(291)**2)*CTRACE(MDM)+6.D0/5.D0
     $     *G(291)**2*CTRACE(MEM)
        SIG1=1.D0/5.D0*G(291)**2*(3.D0*(G(351)+G(352))
     $       +CTRACE(DUMSIG1))
        SIG2=G(292)**2*(G(351)+G(352)+CTRACE(DUMSIG2))
        SIG3=G(293)**2*CTRACE(DUMSIG3)
      ELSE
        SP=(0.D0,0.D0)
        SIG1=(0.D0,0.D0)
        SIG2=(0.D0,0.D0)
        SIG3=(0.D0,0.D0)
      END IF
!
!Now the dummy matrices for the higgs betas
!
      IF(SW2LP.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            DUM2HU11(I,J)=G(351)*ID(I,J)+MQM(I,J)
            DUM2HU12(I,J)=(G(351)+G(352))*ID(I,J)+MQM(I,J)
            DUM2HU21(I,J)=DUM2HU11(I,J)
            DUM2HD11(I,J)=G(352)*ID(I,J)+MQM(I,J)
            DUM2HD12(I,J)=DUM2HU12(I,J)
            DUM2HD13(I,J)=G(352)*ID(I,J)+MLM(I,J)
            DUM2HD21(I,J)=DUM2HD11(I,J)
            DUM2HD31(I,J)=DUM2HD13(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            DUM2HU11(I,J)=(0.D0,0.D0)
            DUM2HU12(I,J)=(0.D0,0.D0)
            DUM2HU21(I,J)=(0.D0,0.D0)
            DUM2HD11(I,J)=(0.D0,0.D0)
            DUM2HD12(I,J)=(0.D0,0.D0)
            DUM2HD13(I,J)=(0.D0,0.D0)
            DUM2HD21(I,J)=(0.D0,0.D0)
            DUM2HD31(I,J)=(0.D0,0.D0)
          END DO
        END DO
      END IF
      DO I=1,3
        DO J=1,3
          DUM1HU1(I,J)=G(351)*ID(I,J)+MQM(I,J)
          DUM1HD1(I,J)=G(352)*ID(I,J)+MQM(I,J)
          DUM1HD2(I,J)=G(352)*ID(I,J)+MLM(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2HU1(I,J)=6.D0*CMATMUL(0,DUM2HU11,YUDYUYUDYU,I,J)+6.D0
     $                   *YUDMUPMYUYUDYU(I,J)
     $                   +CMATMUL(0,DUM2HU12,YUDYUYDDYD,I,J)
     $                   +YUDMUPMYUYDDYD(I,J)+YUDYUMQMYDDYD(I,J)
     $                   +YUDYUYDDMDMYD(I,J)+6.D0*HUDHUYUDYU(I,J)+6.D0
     $                   *HUDYUYUDHU(I,J)+HDDHDYUDYU(I,J)
     $                   +YDDYDHUDHU(I,J)+HDDYDYUDHU(I,J)
     $                   +YDDHDHUDYU(I,J)
            DUM2HU2(I,J)=CMATMUL(0,DUM2HU21,YUDYU,I,J)+YUDMUPMYU(I,J)
     $                   +HUDHU(I,J)
            DUM2HD1(I,J)=6.D0*CMATMUL(0,DUM2HD11,YDDYDYDDYD,I,J)
     $                   +6.D0*YDDMDMYDYDDYD(I,J)
     $                   +CMATMUL(0,DUM2HD12,YUDYUYDDYD,I,J)
     $                   +YUDMUPMYUYDDYD(I,J)+YUDYUMQMYDDYD(I,J)
     $                   +YUDYUYDDMDMYD(I,J)+2.D0
     $                   *CMATMUL(0,DUM2HD13,YEDYEYEDYE,I,J)+2.D0
     $                   *YEDMEMYEYEDYE(I,J)+6.D0*HDDHDYDDYD(I,J)+6.D0
     $                   *HDDYDYDDHD(I,J)+HUDHUYDDYD(I,J)
     $                   +YUDYUHDDHD(I,J)+HUDYUYDDHD(I,J)
     $                   +YUDHUHDDYD(I,J)+2.D0*HEDHEYEDYE(I,J)+2.D0
     $                   *HEDYEYEDHE(I,J)
            DUM2HD2(I,J)=CMATMUL(0,DUM2HD21,YDDYD,I,J)+YDDMDMYD(I,J)
     $                   +HDDHD(I,J)
            DUM2HD3(I,J)=CMATMUL(0,DUM2HD31,YEDYE,I,J)+YEDMEMYE(I,J)
     $                   +HEDHE(I,J)
          ELSE
            DUM2HU1(I,J)=(0.D0,0.D0)
            DUM2HU2(I,J)=(0.D0,0.D0)
            DUM2HD1(I,J)=(0.D0,0.D0)
            DUM2HD2(I,J)=(0.D0,0.D0)
            DUM2HD3(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          DUM1HU(I,J)=CMATMUL(0,DUM1HU1,YUDYU,I,J)+YUDMUPMYU(I,J)
     $                +HUDHU(I,J)
          DUM1HD(I,J)=6.D0*CMATMUL(0,DUM1HD1,YDDYD,I,J)+6.D0
     $                *YDDMDMYD(I,J)+2.D0*CMATMUL(0,DUM1HD2,YEDYE,I,J)
     $                +2.D0*YEDMEMYE(I,J)+6.D0*HDDHD(I,J)
     $                +2.D0*HEDHE(I,J)
        END DO
      END DO
!
!The higgs mass^2 beta parameters are not matrices:
!
      BETA1HU=6.D0*CTRACE(DUM1HU)-6.D0*G(292)**2*CMODSQ(G(322))
     $        -6.D0/5.D0*G(291)**2*CMODSQ(G(321))+3.D0/5.D0*G(291)**2*S
      BETA1HD=CTRACE(DUM1HD)-6.D0*G(292)**2*CMODSQ(G(322))-6.D0/5.D0
     $        *G(291)**2*CMODSQ(G(321))-3.D0/5.D0*G(291)**2*S
      IF(SW2LP.EQ.1)THEN
        BETA2HU=-6.D0*CTRACE(DUM2HU1)+(32.D0*G(293)**2+8.D0/5.D0
     $          *G(291)**2)*CTRACE(DUM2HU2)+32.D0*G(293)**2*(2.D0
     $          *CMODSQ(G(323))*TYUDYU-CCON(G(323))*TYUDHU-G(323)
     $          *THUDYU)+8.D0/5.D0*G(291)**2*(2.D0*CMODSQ(G(321))
     $          *TYUDYU-CCON(G(321))*TYUDHU-G(321)*THUDYU)+6.D0/5.D0
     $          *G(291)**2*SP+33.D0*G(292)**4*CMODSQ(G(322))
     $          +18.D0/5.D0*G(292)**2*G(291)**2*(CMODSQ(G(322))
     $          +CMODSQ(G(321))+CRE(2,G(321),G(322)))+621.D0/25.D0
     $          *G(291)**4*CMODSQ(G(321))+3.D0*G(292)**2*SIG2
     $          +3.D0/5.D0*G(291)**2*SIG1
        BETA2HD=-6.D0*CTRACE(DUM2HD1)+(32.D0*G(293)**2-4.D0/5.D0
     $          *G(291)**2)*CTRACE(DUM2HD2)+32.D0*G(293)**2*(2.D0
     $          *CMODSQ(G(323))*TYDDYD-CCON(G(323))*TYDDHD-G(323)
     $          *THDDYD)-4.D0/5.D0*G(291)**2*(2.D0*CMODSQ(G(321))
     $          *TYDDYD-CCON(G(321))*TYDDHD-G(321)*THDDYD)+12.D0/5.D0
     $          *G(291)**2*(CTRACE(DUM2HD3)+2.D0*CMODSQ(G(321))*TYEDYE
     $          -G(321)*THEDYE-CCON(G(321))*TYEDHE)-6.D0/5.D0
     $          *G(291)**2*SP+33.D0*G(292)**4*CMODSQ(G(322))
     $          +18.D0/5.D0*G(292)**2*G(291)**2*(CMODSQ(G(322))
     $          +CMODSQ(G(321))+CRE(2,G(321),G(322)))+621.D0/25.D0
     $          *G(291)**4*CMODSQ(G(321))+3.D0*G(292)**2*SIG2
     $          +3.D0/5.D0*G(291)**2*SIG1
         B2HUM=BETA2HU
         B2HDM=BETA2HD
      END IF
!
!Calculate the differentials F(351) and F(352)
!
      F(351)=1.D0/16.D0/PI**2*BETA1HU+1.D0/(16.D0*PI**2)**2*BETA2HU
      F(352)=1.D0/16.D0/PI**2*BETA1HD+1.D0/(16.D0*PI**2)**2*BETA2HD
!
!Now for the dummy matrices for the rest of the mass^2
!
      IF(SW2LP.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            DUM2Q31(I,J)=MQM(I,J)+4.D0*G(351)*ID(I,J)
            DUM2Q41(I,J)=MQM(I,J)+4.D0*G(352)*ID(I,J)
            DUM2L21(I,J)=MLM(I,J)+4.D0*G(352)*ID(I,J)
            DUM2U31(I,J)=MUPM(I,J)+4.D0*G(351)*ID(I,J)
            DUM2D31(I,J)=MDM(I,J)+4.D0*G(352)*ID(I,J)
            DUM2E21(I,J)=MEM(I,J)+4.D0*G(352)*ID(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            DUM2Q31(I,J)=(0.D0,0.D0)
            DUM2Q41(I,J)=(0.D0,0.D0)
            DUM2L21(I,J)=(0.D0,0.D0)
            DUM2U31(I,J)=(0.D0,0.D0)
            DUM2D31(I,J)=(0.D0,0.D0)
            DUM2E21(I,J)=(0.D0,0.D0)
          END DO
        END DO
      END IF
!
!This is the second round of dummy matrices
!
      DO I=1,3
        DO J=1,3
          DUM1Q1(I,J)=MQM(I,J)+2.D0*G(351)*ID(I,J)
          DUM1Q2(I,J)=MQM(I,J)+2.D0*G(352)*ID(I,J)
          DUM1Q3(I,J)=YUDYU(I,J)+YDDYD(I,J)
          DUM1L1(I,J)=MLM(I,J)+2.D0*G(352)*ID(I,J)
          DUM1U1(I,J)=2.D0*MUPM(I,J)+4.D0*G(351)*ID(I,J)
          DUM1D1(I,J)=2.D0*MDM(I,J)+4.D0*G(352)*ID(I,J)
          DUM1E1(I,J)=2.D0*MEM(I,J)+4.D0*G(352)*ID(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2Q1(I,J)=2.D0*DUM2Q31(I,J)
            DUM2Q2(I,J)=2.D0*DUM2Q41(I,J)
            DUM2Q3(I,J)=CMATMUL(0,DUM2Q31,YUDYU,I,J)+2.D0
     $                  *YUDMUPMYU(I,J)+YUDYUMQM(I,J)
            DUM2Q4(I,J)=CMATMUL(0,DUM2Q41,YDDYD,I,J)+2.D0
     $                  *YDDMDMYD(I,J)+YDDYDMQM(I,J)
            DUM2Q5(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2Q6(I,J)=MQMYUDYU(I,J)+YUDMUPMYU(I,J)
            DUM2Q7(I,J)=6.D0*MQMYDDYD(I,J)+6.D0*YDDMDMYD(I,J)+2.D0
     $                  *MLMYEDYE(I,J)+2.D0*YEDMEMYE(I,J)
            DUM2Q8(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2Q9(I,J)=6.D0*HDDHD(I,J)+2.D0*HEDHE(I,J)
            DUM2QA(I,J)=6.D0*YDDHD(I,J)+2.D0*YEDHE(I,J)
            DUM2QB(I,J)=6.D0*HDDYD(I,J)+2.D0*HEDYE(I,J)
            DUM2QC(I,J)=2.D0*DUM1Q1(I,J)
            DUM2QD(I,J)=MQM(I,J)+2.D0*G(352)*ID(I,J)
            DUM2L1(I,J)=2.D0*DUM2L21(I,J)
            DUM2L2(I,J)=CMATMUL(0,DUM2L21,YEDYE,I,J)+2.D0
     $                  *YEDMEMYE(I,J)+YEDYEMLM(I,J)
            DUM2L3(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2L4(I,J)=6.D0*MQMYDDYD(I,J)+6.D0*YDDMDMYD(I,J)+2.D0
     $                  *MLMYEDYE(I,J)+2.D0*YEDMEMYE(I,J)
            DUM2L5(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2L6(I,J)=6.D0*HDDHD(I,J)+2.D0*HEDHE(I,J)
            DUM2L7(I,J)=6.D0*YDDHD(I,J)+2.D0*YEDHE(I,J)
            DUM2L8(I,J)=6.D0*HDDYD(I,J)+2.D0*HEDYE(I,J)
            DUM2L9(I,J)=DUM1L1(I,J)
            DUM2U1(I,J)=2.D0*DUM2U31(I,J)
            DUM2U2(I,J)=2.D0*MUPM(I,J)+4.D0*(G(351)+G(352))*ID(I,J)
            DUM2U3(I,J)=CMATMUL(0,DUM2U31,YUYUD,I,J)+2.D0*YUMQMYUD(I,J)
     $                  +YUYUDMUPM(I,J)
            DUM2U4(I,J)=MQMYUDYU(I,J)+YUDMUPMYU(I,J)
            DUM2U5(I,J)=MUPM(I,J)+2.D0*G(351)*ID(I,J)
            DUM2D1(I,J)=2.D0*DUM2D31(I,J)
            DUM2D2(I,J)=2.D0*MDM(I,J)+4.D0*(G(351)+G(352))*ID(I,J)
            DUM2D3(I,J)=CMATMUL(0,DUM2D31,YDYDD,I,J)+2.D0*YDMQMYDD(I,J)
     $                  +YDYDDMDM(I,J)
            DUM2D4(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2D5(I,J)=3.D0*MQMYDDYD(I,J)+3.D0*YDDMDMYD(I,J)
     $                  +MLMYEDYE(I,J)+YEDMEMYE(I,J)
            DUM2D6(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2D7(I,J)=3.D0*HDDHD(I,J)+HEDHE(I,J)
            DUM2D8(I,J)=3.D0*HDDYD(I,J)+HEDYE(I,J)
            DUM2D9(I,J)=3.D0*YDDHD(I,J)+YEDHE(I,J)
            DUM2DA(I,J)=MDM(I,J)+2.D0*G(352)*ID(I,J)
            DUM2E1(I,J)=2.D0*DUM2E21(I,J)
            DUM2E2(I,J)=CMATMUL(0,DUM2E21,YEYED,I,J)+2.D0*YEMLMYED(I,J)
     $                  +YEYEDMEM(I,J)
            DUM2E3(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2E4(I,J)=3.D0*MQMYDDYD(I,J)+3.D0*YDDMDMYD(I,J)
     $                  +MLMYEDYE(I,J)+YEDMEMYE(I,J)
            DUM2E5(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2E6(I,J)=3.D0*HDDHD(I,J)+HEDHE(I,J)
            DUM2E7(I,J)=3.D0*HDDYD(I,J)+HEDYE(I,J)
            DUM2E8(I,J)=3.D0*YDDHD(I,J)+YEDHE(I,J)
            DUM2E9(I,J)=MEM(I,J)+2.D0*G(352)*ID(I,J)
          ELSE
            DUM2Q1(I,J)=(0.D0,0.D0)
            DUM2Q2(I,J)=(0.D0,0.D0)
            DUM2Q3(I,J)=(0.D0,0.D0)
            DUM2Q4(I,J)=(0.D0,0.D0)
            DUM2Q5(I,J)=(0.D0,0.D0)
            DUM2Q6(I,J)=(0.D0,0.D0)
            DUM2Q7(I,J)=(0.D0,0.D0)
            DUM2Q8(I,J)=(0.D0,0.D0)
            DUM2Q9(I,J)=(0.D0,0.D0)
            DUM2QA(I,J)=(0.D0,0.D0)
            DUM2QB(I,J)=(0.D0,0.D0)
            DUM2QC(I,J)=(0.D0,0.D0)
            DUM2QD(I,J)=(0.D0,0.D0)
            DUM2L1(I,J)=(0.D0,0.D0)
            DUM2L2(I,J)=(0.D0,0.D0)
            DUM2L3(I,J)=(0.D0,0.D0)
            DUM2L4(I,J)=(0.D0,0.D0)
            DUM2L5(I,J)=(0.D0,0.D0)
            DUM2L6(I,J)=(0.D0,0.D0)
            DUM2L7(I,J)=(0.D0,0.D0)
            DUM2L8(I,J)=(0.D0,0.D0)
            DUM2L9(I,J)=(0.D0,0.D0)
            DUM2U1(I,J)=(0.D0,0.D0)
            DUM2U2(I,J)=(0.D0,0.D0)
            DUM2U3(I,J)=(0.D0,0.D0)
            DUM2U4(I,J)=(0.D0,0.D0)
            DUM2U5(I,J)=(0.D0,0.D0)
            DUM2D1(I,J)=(0.D0,0.D0)
            DUM2D2(I,J)=(0.D0,0.D0)
            DUM2D3(I,J)=(0.D0,0.D0)
            DUM2D4(I,J)=(0.D0,0.D0)
            DUM2D5(I,J)=(0.D0,0.D0)
            DUM2D6(I,J)=(0.D0,0.D0)
            DUM2D7(I,J)=(0.D0,0.D0)
            DUM2D8(I,J)=(0.D0,0.D0)
            DUM2D9(I,J)=(0.D0,0.D0)
            DUM2DA(I,J)=(0.D0,0.D0)
            DUM2E1(I,J)=(0.D0,0.D0)
            DUM2E2(I,J)=(0.D0,0.D0)
            DUM2E3(I,J)=(0.D0,0.D0)
            DUM2E4(I,J)=(0.D0,0.D0)
            DUM2E5(I,J)=(0.D0,0.D0)
            DUM2E6(I,J)=(0.D0,0.D0)
            DUM2E7(I,J)=(0.D0,0.D0)
            DUM2E8(I,J)=(0.D0,0.D0)
            DUM2E9(I,J)=(0.D0,0.D0)
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the rest of the mass^2
!
      DO I=1,3
        DO J=1,3
          BETA1MQ(I,J)=CMATMUL(0,DUM1Q1,YUDYU,I,J)
     $                 +CMATMUL(0,DUM1Q2,YDDYD,I,J)
     $                 +CMATMUL(0,DUM1Q3,MQM,I,J)+2.D0
     $                 *YUDMUPMYU(I,J)+2.D0*YDDMDMYD(I,J)+2.D0
     $                 *HUDHU(I,J)+2.D0*HDDHD(I,J)+(-32.D0/3.D0
     $                 *G(293)**2*CMODSQ(G(323))-6.D0*G(292)**2
     $                 *CMODSQ(G(322))-2.D0/15.D0*G(291)**2
     $                 *CMODSQ(G(321))+1.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1ML(I,J)=CMATMUL(0,DUM1L1,YEDYE,I,J)+2.D0*YEDMEMYE(I,J)
     $                 +YEDYEMLM(I,J)+2.D0*HEDHE(I,J)+(-6.D0*G(292)**2
     $                 *CMODSQ(G(322))-6.D0/5.D0*G(291)**2
     $                 *CMODSQ(G(321))-3.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1MU(I,J)=CMATMUL(0,DUM1U1,YUYUD,I,J)+4.D0*YUMQMYUD(I,J)
     $                 +2.D0*YUYUDMUPM(I,J)+4.D0*HUHUD(I,J)+(-32.D0/3.D0
     $                 *G(293)**2*CMODSQ(G(323))-32.D0/15.D0*G(291)**2
     $                 *CMODSQ(G(321))-4.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1MD(I,J)=CMATMUL(0,DUM1D1,YDYDD,I,J)+4.D0*YDMQMYDD(I,J)
     $                 +2.D0*YDYDDMDM(I,J)+4.D0*HDHDD(I,J)+(-32.D0/3.D0
     $                 *G(293)**2*CMODSQ(G(323))-8.D0/15.D0*G(291)**2
     $                 *CMODSQ(G(321))+2.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1ME(I,J)=CMATMUL(0,DUM1E1,YEYED,I,J)+4.D0*YEMLMYED(I,J)
     $                 +2.D0*YEYEDMEM(I,J)+4.D0*HEHED(I,J)+(-24.D0/5.D0
     $                 *G(291)**2*CMODSQ(G(321))+6.D0/5.D0*G(291)**2*S)
     $                 *ID(I,J)
         IF(SW2LP.EQ.1)THEN
            BETA2MQ(I,J)=-CMATMUL(0,DUM2Q1,YUDYUYUDYU,I,J)-4.D0
     $                   *YUDMUPMYUYUDYU(I,J)-4.D0*YUDYUMQMYUDYU(I,J)
     $                   -4.D0*YUDYUYUDMUPMYU(I,J)
     $                   -2.D0*YUDYUYUDYUMQM(I,J)
     $                   -CMATMUL(0,DUM2Q2,YDDYDYDDYD,I,J)
     $                   -4.D0*YDDMDMYDYDDYD(I,J)
     $                   -4.D0*YDDYDMQMYDDYD(I,J)
     $                   -4.D0*YDDYDYDDMDMYD(I,J)
     $                   -2.D0*YDDYDYDDYDMQM(I,J)
     $                   -DUM2Q3(I,J)*3.D0*TYUDYU-DUM2Q4(I,J)
     $                   *CTRACE(DUM2Q5)-6.D0*YUDYU(I,J)*CTRACE(DUM2Q6)
     $                   -YDDYD(I,J)*CTRACE(DUM2Q7)-4.D0
     $                   *(YUDYUHUDHU(I,J)+HUDHUYUDYU(I,J)
     $                   +YUDHUHUDYU(I,J)+HUDYUYUDHU(I,J))-4.D0
     $                   *(YDDYDHDDHD(I,J)+HDDHDYDDYD(I,J)
     $                   +YDDHDHDDYD(I,J)+HDDYDYDDHD(I,J))-HUDHU(I,J)
     $                   *6.D0*TYUDYU-YUDYU(I,J)*6.D0*THUDHU-HUDYU(I,J)
     $                   *6.D0*TYUDHU-YUDHU(I,J)*6.D0*THUDYU-HDDHD(I,J)
     $                   *CTRACE(DUM2Q8)-YDDYD(I,J)*CTRACE(DUM2Q9)
     $                   -HDDYD(I,J)*CTRACE(DUM2QA)-YDDHD(I,J)
     $                   *CTRACE(DUM2QB)+2.D0/5.D0*G(291)**2
     $                   *(CMATMUL(0,DUM2QC,YUDYU,I,J)+4.D0
     $                   *YUDMUPMYU(I,J)+2.D0*YUDYUMQM(I,J)+4.D0
     $                   *HUDHU(I,J)-4.D0*G(321)*HUDYU(I,J)-4.D0
     $                   *CCON(G(321))*YUDHU(I,J)+8.D0*CMODSQ(G(321))
     $                   *YUDYU(I,J)+CMATMUL(0,DUM2QD,YDDYD,I,J)+2.D0
     $                   *YDDMDMYD(I,J)+YDDYDMQM(I,J)+2.D0*HDDHD(I,J)
     $                   -2.D0*G(321)*HDDYD(I,J)-2.D0*CCON(G(321))
     $                   *YDDHD(I,J)+4.D0*CMODSQ(G(321))*YDDYD(I,J))
     $                   +(2.D0/5.D0*G(291)**2*SP-128.D0/3.D0*G(293)**4
     $                   *CMODSQ(G(323))+32.D0*G(293)**2*G(292)**2
     $                   *(CMODSQ(G(323))+CMODSQ(G(322))
     $                   +CRE(2,G(322),G(323)))+32.D0/45.D0*G(293)**2
     $                   *G(291)**2*(CMODSQ(G(323))+CMODSQ(G(321))
     $                   +CRE(2,G(321),G(323)))+33.D0*G(292)**4
     $                   *CMODSQ(G(322))+2.D0/5.D0*G(292)**2
     $                   *G(291)**2*(CMODSQ(G(322))+CMODSQ(G(321))
     $                   +CRE(2,G(321),G(322)))+199.D0/75.D0*G(291)**4
     $                   *CMODSQ(G(321))+16.D0/3.D0*G(293)**2*SIG3+3.D0
     $                   *G(292)**2*SIG2+1.D0/15.D0*G(291)**2*SIG1)
     $                   *ID(I,J)
            BETA2ML(I,J)=-CMATMUL(0,DUM2L1,YEDYEYEDYE,I,J)-4.D0
     $                   *YEDMEMYEYEDYE(I,J)-4.D0*YEDYEMLMYEDYE(I,J)
     $                   -4.D0*YEDYEYEDMEMYE(I,J)
     $                   -2.D0*YEDYEYEDYEMLM(I,J)
     $                   -DUM2L2(I,J)*CTRACE(DUM2L3)-YEDYE(I,J)
     $                   *CTRACE(DUM2L4)-4.D0*(YEDYEHEDHE(I,J)
     $                   +HEDHEYEDYE(I,J)+YEDHEHEDYE(I,J)
     $                   +HEDYEYEDHE(I,J))-HEDHE(I,J)*CTRACE(DUM2L5)
     $                   -YEDYE(I,J)*CTRACE(DUM2L6)-HEDYE(I,J)
     $                   *CTRACE(DUM2L7)-YEDHE(I,J)*CTRACE(DUM2L8)
     $                   +6.D0/5.D0*G(291)**2
     $                   *(CMATMUL(0,DUM2L9,YEDYE,I,J)
     $                   +2.D0*YEDMEMYE(I,J)+YEDYEMLM(I,J)
     $                   +2.D0*HEDHE(I,J)-2.D0*G(321)*HEDYE(I,J)
     $                   -2.D0*CCON(G(321))*YEDHE(I,J)+4.D0
     $                   *CMODSQ(G(321))*YEDYE(I,J))+(-6.D0/5.D0
     $                   *G(291)**2*SP+33.D0*G(292)**4*CMODSQ(G(322))
     $                   +18.D0/5.D0*G(292)**2*G(291)**2
     $                   *(CMODSQ(G(322))+CMODSQ(G(321))
     $                   +CRE(2,G(321),G(322)))+621.D0/25.D0
     $                   *G(291)**4*CMODSQ(G(321))+3.D0*G(292)**2*SIG2
     $                   +3.D0/5.D0*G(291)**2*SIG1)*ID(I,J)
            BETA2MU(I,J)=-CMATMUL(0,DUM2U1,YUYUDYUYUD,I,J)-4.D0
     $                   *YUMQMYUDYUYUD(I,J)-4.D0*YUYUDMUPMYUYUD(I,J)
     $                   -4.D0*YUYUDYUMQMYUD(I,J)
     $                   -2.D0*YUYUDYUYUDMUPM(I,J)
     $                   -CMATMUL(0,DUM2U2,YUYDDYDYUD,I,J)-4.D0
     $                   *YUMQMYDDYDYUD(I,J)-4.D0*YUYDDMDMYDYUD(I,J)
     $                   -4.D0*YUYDDYDMQMYUD(I,J)
     $                   -2.D0*YUYDDYDYUDMUPM(I,J)
     $                   -DUM2U3(I,J)*6.D0*TYUDYU-12.D0*YUYUD(I,J)
     $                   *CTRACE(DUM2U4)-4.D0*(HUHUDYUYUD(I,J)
     $                   +YUYUDHUHUD(I,J)+HUYUDYUHUD(I,J)
     $                   +YUHUDHUYUD(I,J))-4.D0*(HUHDDYDYUD(I,J)
     $                   +YUYDDHDHUD(I,J)+HUYDDYDHUD(I,J)
     $                   +YUHDDHDYUD(I,J))-12.D0*(HUHUD(I,J)*TYUDYU
     $                   +YUYUD(I,J)*THUDHU+HUYUD(I,J)*THUDYU+YUHUD(I,J)
     $                   *TYUDHU)+(6.D0*G(292)**2-2.D0/5.D0*G(291)**2)
     $                   *(CMATMUL(0,DUM2U5,YUYUD,I,J)
     $                   +2.D0*YUMQMYUD(I,J)+YUYUDMUPM(I,J)
     $                   +2.D0*HUHUD(I,J))+12.D0*G(292)**2
     $                   *(2.D0*CMODSQ(G(322))*YUYUD(I,J)-CCON(G(322))
     $                   *HUYUD(I,J)-G(322)*YUHUD(I,J))-4.D0/5.D0
     $                   *G(291)**2*(2.D0*CMODSQ(G(321))*YUYUD(I,J)
     $                   -CCON(G(321))*HUYUD(I,J)-G(321)*YUHUD(I,J))
     $                   +(-8.D0/5.D0*G(291)**2*SP-128.D0/3.D0*G(293)**4
     $                   *CMODSQ(G(323))+512.D0/45.D0*G(293)**2
     $                   *G(291)**2*(CMODSQ(G(323))+CMODSQ(G(321))
     $                   +CRE(2,G(321),G(323)))+3424.D0/75.D0
     $                   *G(291)**4*CMODSQ(G(321))+16.D0/3.D0*G(293)**2
     $                   *SIG3+16.D0/15.D0*G(291)**2*SIG1)*ID(I,J)
            BETA2MD(I,J)=-CMATMUL(0,DUM2D1,YDYDDYDYDD,I,J)-4.D0
     $                   *YDMQMYDDYDYDD(I,J)-4.D0*YDYDDMDMYDYDD(I,J)
     $                   -4.D0*YDYDDYDMQMYDD(I,J)
     $                   -2.D0*YDYDDYDYDDMDM(I,J)
     $                   -CMATMUL(0,DUM2D2,YDYUDYUYDD,I,J)-4.D0
     $                   *YDMQMYUDYUYDD(I,J)-4.D0*YDYUDMUPMYUYDD(I,J)
     $                   -4.D0*YDYUDYUMQMYDD(I,J)
     $                   -2.D0*YDYUDYUYDDMDM(I,J)
     $                   -DUM2D3(I,J)*CTRACE(DUM2D4)-4.D0*YDYDD(I,J)
     $                   *CTRACE(DUM2D5)-4.D0*(HDHDDYDYDD(I,J)
     $                   +YDYDDHDHDD(I,J)+HDYDDYDHDD(I,J)
     $                   +YDHDDHDYDD(I,J))-4.D0*(HDHUDYUYDD(I,J)
     $                   +YDYUDHUHDD(I,J)+HDYUDYUHDD(I,J)
     $                   +YDHUDHUYDD(I,J))-4.D0*HDHDD(I,J)
     $                   *CTRACE(DUM2D6)-4.D0*YDYDD(I,J)*CTRACE(DUM2D7)
     $                   -4.D0*HDYDD(I,J)*CTRACE(DUM2D8)-4.D0*YDHDD(I,J)
     $                   *CTRACE(DUM2D9)+(6.D0*G(292)**2+2.D0/5.D0
     $                   *G(291)**2)*(CMATMUL(0,DUM2DA,YDYDD,I,J)+2.D0
     $                   *YDMQMYDD(I,J)+YDYDDMDM(I,J)+2.D0*HDHDD(I,J))
     $                   +12.D0*G(292)**2*(2.D0*CMODSQ(G(322))
     $                   *YDYDD(I,J)-CCON(G(322))*HDYDD(I,J)-G(322)
     $                   *YDHDD(I,J))+4.D0/5.D0*G(291)**2*(2.D0
     $                   *CMODSQ(G(321))*YDYDD(I,J)-CCON(G(321))
     $                   *HDYDD(I,J)-G(321)*YDHDD(I,J))+(4.D0/5.D0
     $                   *G(291)**2*SP-128.D0/3.D0*G(293)**4
     $                   *CMODSQ(G(323))+128.D0/45.D0*G(293)**2
     $                   *G(291)**2*(CMODSQ(G(323))+CMODSQ(G(321))
     $                   +CRE(2,G(321),G(323)))+808.D0/75.D0*G(291)**4
     $                   *CMODSQ(G(321))+16.D0/3.D0*G(293)**2*SIG3
     $                   +4.D0/15.D0*G(291)**2*SIG1)*ID(I,J)
            BETA2ME(I,J)=-CMATMUL(0,DUM2E1,YEYEDYEYED,I,J)-4.D0
     $                   *YEMLMYEDYEYED(I,J)-4.D0*YEYEDMEMYEYED(I,J)
     $                   -4.D0*YEYEDYEMLMYED(I,J)
     $                   -2.D0*YEYEDYEYEDMEM(I,J)
     $                   -DUM2E2(I,J)*CTRACE(DUM2E3)-4.D0*YEYED(I,J)
     $                   *CTRACE(DUM2E4)-4.D0*(HEHEDYEYED(I,J)
     $                   +YEYEDHEHED(I,J)+HEYEDYEHED(I,J)
     $                   +YEHEDHEYED(I,J))-4.D0*HEHED(I,J)
     $                   *CTRACE(DUM2E5)-4.D0*YEYED(I,J)*CTRACE(DUM2E6)
     $                   -4.D0*HEYED(I,J)*CTRACE(DUM2E7)-4.D0*YEHED(I,J)
     $                   *CTRACE(DUM2E8)+(6.D0*G(292)**2-6.D0/5.D0
     $                   *G(291)**2)*(CMATMUL(0,DUM2E9,YEYED,I,J)+2.D0
     $                   *YEMLMYED(I,J)+YEYEDMEM(I,J)+2.D0*HEHED(I,J))
     $                   +12.D0*G(292)**2*(2.D0*CMODSQ(G(322))
     $                   *YEYED(I,J)-CCON(G(322))*HEYED(I,J)-G(322)
     $                   *YEHED(I,J))-12.D0/5.D0*G(291)**2*(2.D0
     $                   *CMODSQ(G(321))*YEYED(I,J)-CCON(G(321))
     $                   *HEYED(I,J)-G(321)*YEHED(I,J))+(12.D0/5.D0
     $                   *G(291)**2*SP+2808.D0/25.D0*G(291)**4
     $                   *CMODSQ(G(321))+12.D0/5.D0*G(291)**2*SIG1)
     $                   *ID(I,J)
            B2MQM(I,J)=BETA2MQ(I,J)
            B2MLM(I,J)=BETA2ML(I,J)
            B2MUM(I,J)=BETA2MU(I,J)
            B2MDM(I,J)=BETA2MD(I,J)
            B2MEM(I,J)=BETA2ME(I,J)
          END  IF
!
!Calculate the differentials DM
!
          DM(1,I,J)=1.D0/16.D0/PI**2*BETA1MQ(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2MQ(I,J)
          DM(2,I,J)=1.D0/16.D0/PI**2*BETA1ML(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2ML(I,J)
          DM(3,I,J)=1.D0/16.D0/PI**2*BETA1MU(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2MU(I,J)
          DM(4,I,J)=1.D0/16.D0/PI**2*BETA1MD(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2MD(I,J)
          DM(5,I,J)=1.D0/16.D0/PI**2*BETA1ME(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2ME(I,J)
!
!Convert into form readable by RKSTP
!
          F(352+(I-1)*3+J)=DM(1,I,J)
          F(361+(I-1)*3+J)=DM(2,I,J)
          F(370+(I-1)*3+J)=DM(3,I,J)
          F(379+(I-1)*3+J)=DM(4,I,J)
          F(388+(I-1)*3+J)=DM(5,I,J)
        END DO
      END DO
!
!***********************************************************************
!END OF MSSM SECTION
!
!Now I can find the threshold gauge running with a change at m_H to SM
!
      DO I=1,3
        SUM=(0.D0,0.D0)
        DO J=1,3
          IF(THHH.EQ.0)THEN
            SUM=SUM+B2LPSM(I,J)*G(J)**2
          ELSE
            SUM=SUM+B2LPM(I,J)*G(J)**2
          END IF
        END DO
        IF(THHH.EQ.0)THEN
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CSM(I,1)*TLYUDLYU+CSM(I,2)*TLYDDLYD
     $         +CSM(I,3)*TLYEDLYE)))
        ELSE
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CM(I,1)*TYUDYU+CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        END IF
      END DO
!
!Next the full Yukawas
!
      DO I=1,3
        DO J=1,3
          DUMU1(I,J)=THSH*SQFTUQDFTUQ(I,J)
     $               +4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $               +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
          DUMU2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMD1(I,J)=THSH*SQFTDQDFTDQ(I,J)
     $               +1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $               +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
          DUMD2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUME1(I,J)=THSH*SLFTELDFTEL(I,J)+THSB*SEGTPESGTPET(I,J)
          DUME2(I,J)=2.D0*THSH*SEFTEEFTEED(I,J)
     $               +3.D0*THSW*SLGTLTGTLS(I,J)+THSB*SLGTPLTGTPLS(I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
!
          B1U(I,J)=1.D0/2.D0*(3.D0*FUFUDFU(I,J)+FDFDDFU(I,J))
     $             +CMATMUL(0,FU,DUMU1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUMU2,FU,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(204))
     $              *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $              *CONJG(G(184))*SQGTPQTFTUQ(I,J))
     $             -4.D0/3.D0*THSB*THSH*CONJG(G(184))
     $              *SUFTUUGTPUT(I,J)
     $             +FU(I,J)*3.D0*TFUDFU
     $             +1.D0/2.D0*THSH*FU(I,J)*(3.D0*THSW*MGTHUSQ
     $              +THSB*MGTPHUSQ)
     $             -FU(I,J)*(17.D0/20.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1D(I,J)=1.D0/2.D0*(3.D0*FDFDDFD(I,J)+FUFUDFD(I,J))
     $             +CMATMUL(0,FD,DUMD1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUMD2,FD,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(205))
     $              *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $              *CONJG(G(185))*SQGTPQTFTDQ(I,J))
     $             -2.D0/3.D0*THSB*THSH*CONJG(G(185))
     $              *SDFTDDGTPDT(I,J)
     $             +FD(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FD(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FD(I,J)*(3.D0/12.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1E(I,J)=3.D0/2.D0*FEFEDFE(I,J)
     $             +CMATMUL(0,FE,DUME1,I,J)
     $             +1.D0/4.D0*CMATMUL(0,DUME2,FE,I,J)
     $             +THSH*(-3.D0*THSW*CONJG(G(205))
     $              *SLGTLTFTEL(I,J)+THSB
     $              *CONJG(G(185))*SLGTPLTFTEL(I,J))
     $             -2.D0*THSB*THSH*CONJG(G(185))
     $              *SEFTEEGTPET(I,J)
     $             +FE(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FE(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FE(I,J)*(9.D0/4.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Convert into form readable by RKSTP. The transpose in BETA2 takes
!account of the differences in notation.
!
          F(3+(I-1)*3+J)=(1.D0/16.D0/PI**2*B1U(I,J)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I))*THHH
          F(12+(I-1)*3+J)=(1.D0/16.D0/PI**2*B1D(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YMD(J,I))*THHH
          F(21+(I-1)*3+J)=(1.D0/16.D0/PI**2*B1E(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YME(J,I))*THHH
        END DO
      END DO
!
!The lambdas use the same dummy matrices. I will reuse the betas
!Only find the lambdas if we are below m_H
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            B1U(I,J)=(0.D0,0.D0)
            B1D(I,J)=(0.D0,0.D0)
            B1E(I,J)=(0.D0,0.D0)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            IF(SW2LP.EQ.1)THEN
              DUM2U(I,J)=3.D0/2.D0*LYUDLYU2(I,J)
     $                  -LYUDLYULYDDLYD(I,J)
     $                  -1.D0/4.D0*LYDDLYDLYUDLYU(I,J)
     $                  +11.D0/4.D0*LYDDLYD2(I,J)
     $                  +Y2*(5.D0/4.D0*LYDDLYD(I,J)-9.D0/4.D0
     $                  *LYUDLYU(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(429)**2*ID(I,J)-2.D0*G(429)*(3.D0
     $                  *LYUDLYU(I,J)+LYDDLYD(I,J))+(221.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYUDLYU(I,J)-(17.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYDDLYD(I,J)+Y4*ID(I,J)+((7.D0/150.D0
     $                  +2.D0/3.D0*NG)*G(1)**4-9.D0/20.D0*G(1)**2
     $                  *G(2)**2+19.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2D(I,J)=3.D0/2.D0*LYDDLYD2(I,J)
     $                  -LYDDLYDLYUDLYU(I,J)
     $                  -1.D0/4.D0*LYUDLYULYDDLYD(I,J)
     $                  +11.D0/4.D0*LYUDLYU2(I,J)
     $                  +Y2*(5.D0/4.D0*LYUDLYU(I,J)-9.D0/4.D0
     $                  *LYDDLYD(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(429)**2*ID(I,J)-2.D0*G(429)*(3.D0
     $                  *LYDDLYD(I,J)+LYUDLYU(I,J))+(161.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYDDLYD(I,J)-(77.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYUDLYU(I,J)+Y4*ID(I,J)+(-(37.D0/300.D0
     $                  -4.D0/15.D0*NG)*G(1)**4-27.D0/20.D0*G(1)**2
     $                  *G(2)**2+31.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2E(I,J)=3.D0/2.D0*LYEDLYE2(I,J)
     $                  -Y2*9.D0/4.D0*LYEDLYE(I,J)-CHI4*ID(I,J)
     $                  +3.D0/2.D0*G(429)**2*ID(I,J)-6.D0*G(429)
     $                  *LYEDLYE(I,J)+(441.D0/80.D0*G(1)**2
     $                  +117.D0/16.D0*G(2)**2)*LYEDLYE(I,J)
     $                  +Y4*ID(I,J)+((21.D0/100.D0+8.D0/5.D0*NG)
     $                  *G(1)**4+27.D0/20.D0*G(1)**2*G(2)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4)*ID(I,J)
            ELSE
              DUM2U(I,J)=(0.D0,0.D0)
              DUM2D(I,J)=(0.D0,0.D0)
              DUM2E(I,J)=(0.D0,0.D0)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            B1U(I,J)=1.D0/2.D0*(3.D0*LULUDLU(I,J)+LDLDDLU(I,J)
     $                -4.D0*LDLDDLU(I,J))
     $               +CMATMUL(0,LU,DUMU1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUMU2,LU,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(289))
     $                *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $                *CONJG(G(287))*SQGTPQTFTUQ(I,J))
     $               -4.D0/3.D0*THSB*THSH*CONJG(G(287))
     $                *SUFTUUGTPUT(I,J)
     $               +LU(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LU(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LU(I,J)*(17.D0/20.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1D(I,J)=1.D0/2.D0*(3.D0*LDLDDLD(I,J)+LULUDLD(I,J)
     $                -4.D0*LULUDLD(I,J))
     $               +CMATMUL(0,LD,DUMD1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUMD2,LD,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(290))
     $                *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $                *CONJG(G(288))*SQGTPQTFTDQ(I,J))
     $               -2.D0/3.D0*THSB*THSH*CONJG(G(288))
     $                *SDFTDDGTPDT(I,J)
     $               +LD(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LD(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LD(I,J)*(3.D0/12.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1E(I,J)=3.D0/2.D0*LELEDLE(I,J)+CMATMUL(0,LE,DUME1,I,J)
     $               +1.D0/4.D0*CMATMUL(0,DUME2,LE,I,J)
     $               +THSH*(-3.D0*THSW*CONJG(G(290))
     $                *SLGTLTFTEL(I,J)+THSB
     $                *CONJG(G(288))*SLGTPLTFTEL(I,J))
     $               -2.D0*THSB*THSH*CONJG(G(288))
     $                *SEFTEEGTPET(I,J)
     $               +LE(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LE(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LE(I,J)*(9.D0/4.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2)
!
            IF(SW2LP.EQ.1)THEN
              BETA2U(I,J)=CMATMUL(0,LYU,DUM2U,I,J)
              BETA2D(I,J)=CMATMUL(0,LYD,DUM2D,I,J)
              BETA2E(I,J)=CMATMUL(0,LYE,DUM2E,I,J)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
!
!Convert into form readable by RKSTP.
!
            F(111+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2U(J,I)
            F(120+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2D(J,I)
            F(129+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2E(J,I)
          END DO
        END DO
      END IF
!
!Next I am going to work out the gaugino terms, \mu and M_{1,2,3}
!and the running of B in MV notation
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*SUFTUUFTUUD(I,J)+3.D0*SDFTDDFTDDD(I,J)
     $                   +SEFTEEFTEED(I,J)+3.D0*SQFTUQDFTUQ(I,J)
     $                   +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
!
          DUMM1(I,J)=1.D0/3.D0*SQGTPQDGTPQ(I,J)+SLGTPLDGTPL(I,J)
     $               +8.D0/3.D0*SUGTPUGTPUD(I,J)
     $               +2.D0/3.D0*SDGTPDGTPDD(I,J)+2.D0*SEGTPEGTPED(I,J)
!
          DUMM2(I,J)=3.D0*SQGTQDGTQ(I,J)+SLGTLDGTL(I,J)
!
          DUMM3(I,J)=2.D0*SQGTSQDGTSQ(I,J)+SUGTSUGTSUD(I,J)
     $               +SDGTSDGTSDD(I,J)
!     
        END DO
      END DO
      TDUMGRKMU=CTRACE(DUMGRKMU1)
      TDUMM1=CTRACE(DUMM1)
      TDUMM2=CTRACE(DUMM2)
      TDUMM3=CTRACE(DUMM3)
      SQTSUMTSFUSAUT=TCSFMUL(THSQ,SUMTSFUSAUT)
      SQTSDMTSFDSADT=TCSFMUL(THSQ,SDMTSFDSADT)
      SLTSEMTSFESAET=TCSFMUL(THSL,SEMTSFESAET)
!
      IF(THHH.EQ.0)THEN
        BM(1)=G(31)*(TDUMM1+THSH*(MSGTPHUSQ+MCGTPHDSQ))
     $        -2.D0*THSH*(G(287)*CONJG(G(108))*G(288)
     $                    +CONJG(G(287))*G(108)*CONJG(G(288)))
!
        BM(2)=G(32)*(TDUMM2+THSH*(MSGTHUSQ+MCGTHDSQ))
     $        -2.D0*THSH*(G(289)*CONJG(G(108))*G(290)
     $                    +CONJG(G(289))*G(108)*CONJG(G(290)))
     $       -12.D0*G(32)*G(2)**2
!
        BM(3)=G(33)*TDUMM3-18.D0*G(33)*G(3)**2
!
        BGRKMU=1.D0/2.D0*G(108)*THSH*TDUMGRKMU
     $       +1.D0/4.D0*G(108)*THSH*((3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $        +(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ))
     $       -(3.D0*THSW*G(289)*(G(32)+(0.D0,1.D0)*G(600))*G(290)
     $        +THSB*G(287)*(G(31)+(0.D0,1.D0)*G(599))*G(288))
     $       -G(108)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
!
        BETA1B=(0.D0,0.D0)
      ELSE
        BM(1)=G(31)*(TDUMM1+THSH*(MGTPHUSQ+MGTPHDSQ))
!
        BM(2)=G(32)*(TDUMM2+THSH*(MGTHUSQ+MGTHDSQ))
     $        -12.D0*G(32)*G(2)**2
!
        BM(3)=G(33)*TDUMM3-18.D0*G(33)*G(3)**2
!
        BGRKMU=1.D0/2.D0*G(108)*THSH*TDUMGRKMU
     $         +1.D0/4.D0*G(108)*THSH*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ
     $          +3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $         -G(108)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
        BETA1B=-1.D0/2.D0*(G(457)**2+3.D0*G(458)**2)*G(109)
     $         +6.D0*SQTSUMTSFUSAUT
     $         +6.D0*SQTSDMTSFDSADT
     $         +2.D0*SLTSEMTSFESAET
     $         -(2.D0*THSH*THSB*G(108)*CONJG(G(184))*CONJG(G(185))
     $          *(G(31)-(0.D0,1.D0)*G(599))+6.D0*THSH*THSW*G(108)
     $          *CONJG(G(204))*CONJG(G(205))*(G(32)-(0.D0,1.D0)*G(600)))
     $         -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(109)
     $         +1.D0/2.D0*(6.D0*TFUDFU+6.D0*TFDDFD+2.D0*TFEDFE
     $          +THSB*THSH*(MGTPHUSQ+MGTPHDSQ)+3.D0*THSW*THSH
     $          *(MGTHUSQ+MGTHDSQ))*G(109)
      END IF
!
!The RKSTP compatible derivatives are: (with a minus sign
!to convert the two loop mu running to BT notation)
!NB: There is some work here converting the MSSM MV notation
!    M_{1,2,3} 2-lp beta functions to BT notation.
!
      F(31)=THSB*(1.D0/16.D0/PI**2*BM(1)
     $      +1.D0/(16.D0*PI**2)**2*1.D0/2.D0*(B2M(1)+CONJG(B2M(1))))
      F(32)=THSW*(1.D0/16.D0/PI**2*BM(2)
     $      +1.D0/(16.D0*PI**2)**2*1.D0/2.D0*(B2M(2)+CONJG(B2M(2))))
      F(33)=THGL*(1.D0/16.D0/PI**2*BM(3)
     $      +1.D0/(16.D0*PI**2)**2*1.D0/2.D0*(B2M(3)+CONJG(B2M(3))))
!
      F(108)=(1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2*B2GRKMUM)
     $                                                            *THSH
      IF(THHH.EQ.1)THEN
        F(109)=1.D0/16.D0/PI**2*BETA1B-1.D0/(16.D0*PI**2)**2*BETA2BM
      END IF
!
!V_U and V_D - from PRD 49,4882 (1994)
!
      BETA1VU=3.D0/4.D0*(1.D0/5.D0*G(291)**2+G(292)**2)-3.D0*TYUDYU
      BETA1VD=3.D0/4.D0*(1.D0/5.D0*G(291)**2+G(292)**2)-3.D0*TYDDYD
     $        -TYEDYE
      IF(SW2LP.EQ.1)THEN
        BETA2VU=3.D0/4.D0*(3.D0*TYUDYUYUDYU+3.D0*TYUDYUYDDYD)-(19.D0
     $          /10.D0*G(291)**2+9.D0/2.D0*G(292)**2+20.D0*G(293)**2)
     $          *TYUDYU-(279.D0/800.D0+1803.D0/1600.D0*3.D0)*G(291)**4
     $          -(207.D0/32.D0+357.D0/64.D0*3.D0)*G(292)**4-(27.D0/80.D0
     $          +9.D0/80.D0*3.D0)*G(291)**2*G(292)**2
        BETA2VD=3.D0/4.D0*(3.D0*TYDDYDYDDYD+3.D0*TYDDYDYUDYU
     $          +TYEDYEYEDYE)-(2.D0/5.D0*G(291)**2+9.D0/2.D0*G(292)**2
     $          +20.D0*G(293)**2)*TYDDYD-(9.D0/5.D0*G(291)**2+3.D0/2.D0
     $          *G(292)**2)*TYEDYE-(279.D0/800.D0+1803.D0/1600.D0*3.D0)
     $          *G(291)**4-(207.D0/32.D0+357.D0/64.D0*3.D0)*G(292)**4
     $          -(27.D0/80.D0+9.D0/80.D0*3.D0)*G(291)**2*G(292)**2
      END IF
!
      F(110)=THHH*(G(110)*(1.D0/16.D0/PI**2*BETA1VU
     $             +1.D0/(16.D0*PI**2)**2*BETA2VU))
      F(111)=THHH*(G(111)*(1.D0/16.D0/PI**2*BETA1VD
     $             +1.D0/(16.D0*PI**2)**2*BETA2VD))
!
!Now the trilinear couplings - it should be ok to use the same
!dummy matrices, but the betas must be reset
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=(0.D0,0.D0)
          BETA2U(I,J)=(0.D0,0.D0)
          BETA1D(I,J)=(0.D0,0.D0)
          BETA2D(I,J)=(0.D0,0.D0)
          BETA1E(I,J)=(0.D0,0.D0)
          BETA2E(I,J)=(0.D0,0.D0)
        END DO
      END DO
!
!Some common terms
!
      DO I=1,3
        DO J=1,3
          DUM1U1(I,J)=8.D0/9.D0*THSB*GTPUTGTPUS(I,J)
     $               +8.D0/3.D0*THGL*GTSUTGTSUS(I,J)
     $               +2.D0*THSH*FTUUDFTUU(I,J)
          DUM1U2(I,J)=THSH*FTUQFTUQD(I,J)+THSH*FTDQFTDQD(I,J)
     $               +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $               +3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $               +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
!
          DUM1D1(I,J)=DUM1U2(I,J)
          DUM1D2(I,J)=2.D0/9.D0*THSB*GTPDTGTPDS(I,J)
     $               +8.D0/3.D0*THGL*GTSDTGTSDS(I,J)
     $               +2.D0*THSH*FTDDDFTDD(I,J)
!
          DUM1E1(I,J)=THSH*FTELFTELD(I,J)
     $               +1.D0/2.D0*THSB*GTPLSGTPLT(I,J)
     $               +3.D0/2.D0*THSW*GTLSGTLT(I,J)
          DUM1E2(I,J)=2.D0*THSB*GTPETGTPES(I,J)
     $               +2.D0*THSH*FTEEDFTEE(I,J)
        END DO
      END DO
!
!Running above m_H first
!
      IF(THHH.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            SUAUID(I,J)=CSFMUL(THSU,AU,ID,I,J)
            SUAUFUQDFUQ(I,J)=CSFMUL(THSU,AU,FUQDFUQ,I,J)
            SQIDAU(I,J)=CSFMUL(THSQ,ID,AU,I,J)
            SQFUHUDAU(I,J)=CSFMUL(THSQ,FUHUD,AU,I,J)
            SQFUURFUURDAU(I,J)=CSFMUL(THSQ,FUURFUURD,AU,I,J)
            SDADFDQDFUQ(I,J)=CSFMUL(THSD,AD,FDQDFUQ,I,J)
!
            SQIDAD(I,J)=CSFMUL(THSQ,ID,AD,I,J)
            SQFDDRFDDRDAD(I,J)=CSFMUL(THSQ,FDDRFDDRD,AD,I,J)
            SQFDHDDAD(I,J)=CSFMUL(THSQ,FDHDD,AD,I,J)
            SLFEHDDAE(I,J)=CSFMUL(THSL,FEHDD,AE,I,J)
            SDADID(I,J)=CSFMUL(THSD,AD,ID,I,J)
            SDADFDQDFDQ(I,J)=CSFMUL(THSD,AD,FDQDFDQ,I,J)
            SUAUFUQDFDQ(I,J)=CSFMUL(THSU,AU,FUQDFDQ,I,J)
!
            SLIDAE(I,J)=CSFMUL(THSL,ID,AE,I,J)
            SLFEERFEERDAE(I,J)=CSFMUL(THSL,FEERFEERD,AE,I,J)
            SEAEID(I,J)=CSFMUL(THSE,AE,ID,I,J)
            SEAEFELDFEL(I,J)=CSFMUL(THSE,AE,FELDFEL,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SUSQIDAUID(I,J)=CSFMUL(THSU,SQIDAU,ID,I,J)
            SDSQIDADID(I,J)=CSFMUL(THSD,SQIDAD,ID,I,J)
            SESLIDAEID(I,J)=CSFMUL(THSE,SLIDAE,ID,I,J)
          END DO
        END DO
        SUTSQFUHUDAU=TCSFMUL(THSU,SQFUHUDAU)
        SDTSQFDHDDAD=TCSFMUL(THSD,SQFDHDDAD)
        SETSLFEHDDAE=TCSFMUL(THSE,SLFEHDDAE)
        DO I=1,3
          DO J=1,3
            BETA1U(I,J)=-2.D0/3.D0*G(457)**2*SUAUID(I,J)
     $                 +2.D0*SUAUFUQDFUQ(I,J)
     $                 -2.D0*(G(457)**2/9.D0+4.D0*G(459)**2/3.D0)
     $                  *SUSQIDAUID(I,J)
     $                 +6.D0*FUHU(I,J)*SUTSQFUHUDAU
     $                 +(G(457)**2/6.D0-3.D0*G(458)**2/2.D0)*SQIDAU(I,J)
     $                 +4.D0*SQFUURFUURDAU(I,J)
     $                 +2.D0*SDADFDQDFUQ(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*CONJG(G(184))*GTPQSFTUU(I,J)
     $                  -4.D0/3.D0*GTPQSFUGTPUS(I,J)
     $                  -4.D0*THSH*FTUQGTPUS(I,J)*CONJG(G(184)))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSFUGTSUS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *CONJG(G(204))*GTQSFTUU(I,J)
     $                 +CSFMUL(THSU,AU,DUM1U1,I,J)
     $                 +(3.D0*TFUDFU+1.D0/2.D0*THSH*THSB*MGTPHUSQ
     $                  +3.D0/2.D0*THSH*THSW*MGTHUSQ)*AU(I,J)
     $                 +CSFMUL(THSQ,DUM1U2,AU,I,J)
     $                 -(13.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*AU(I,J)
!
            BETA1D(I,J)=-(G(457)**2/6.D0+3.D0*G(458)**2/2.D0)
     $                  *SQIDAD(I,J)
     $                 +4.D0*SQFDDRFDDRDAD(I,J)
     $                 +2.D0*(G(457)**2/18.D0-4.D0*G(459)**2/3.D0)
     $                  *SDSQIDADID(I,J)
     $                 +6.D0*FDHD(I,J)*SDTSQFDHDDAD
     $                 +2.D0*FDHD(I,J)*SETSLFEHDDAE
     $                 -1.D0/3.D0*G(457)**2*SDADID(I,J)
     $                 +2.D0*SDADFDQDFDQ(I,J)
     $                 +2.D0*SUAUFUQDFDQ(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(-THSH*CONJG(G(185))*GTPQSFTDD(I,J)
     $                  +2.D0/3.D0*GTPQSFDGTPDS(I,J)
     $                  -2.D0*THSH*FTDQGTPDS(I,J)*CONJG(G(185)))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSFDGTSDS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *CONJG(G(205))*GTQSFTDD(I,J)
     $                 +CSFMUL(THSQ,DUM1D1,AD,I,J)
     $                 +(3.D0*TFDDFD+TFEDFE
     $                  +1.D0/2.D0*THSH*THSB*MGTPHDSQ
     $                  +3.D0/2.D0*THSH*THSW*MGTHDSQ)*AD(I,J)
     $                 +CSFMUL(THSD,AD,DUM1D2,I,J)
     $                 -(7.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*AD(I,J)
!
            BETA1E(I,J)=(G(457)**2/2.D0-3.D0*G(458)**2/2.D0)*SLIDAE(I,J)
     $                 +4.D0*SLFEERFEERDAE(I,J)
     $                 -G(457)**2*SESLIDAEID(I,J)
     $                 +2.D0*FEHD(I,J)*SETSLFEHDDAE
     $                 +6.D0*FEHD(I,J)*SDTSQFDHDDAD
     $                 -G(457)**2*SEAEID(I,J)
     $                 +2.D0*SEAEFELDFEL(I,J)
     $                 +2.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*CONJG(G(185))*GTPLSFTEE(I,J)
     $                  -2.D0*GTPLSFEGTPES(I,J)
     $                  -2.D0*THSH*FTELGTPES(I,J)*CONJG(G(185)))
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *CONJG(G(205))*GTLSFTEE(I,J)
     $                 +CSFMUL(THSL,DUM1E1,AE,I,J)
     $                 +(3.D0*TFDDFD+TFEDFE
     $                  +1.D0/2.D0*THSH*THSB*MGTPHDSQ
     $                  +3.D0/2.D0*THSH*THSW*MGTHDSQ)*AE(I,J)
     $                 +CSFMUL(THSE,AE,DUM1E2,I,J)
     $                 -(27.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $                  *AE(I,J)
          END DO
        END DO
!
!Convert into form readable by RKSTP (converting 2-lp to BT notation...)
!
        DO I=1,3
          DO J=1,3
            F(33+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BETA1U(I,J)
     $                      -1.D0/(16.D0*PI**2)**2*B2HMU(J,I))*THSU(J)
            F(42+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BETA1D(I,J)
     $                      -1.D0/(16.D0*PI**2)**2*B2HMD(J,I))*THSD(J)
            F(51+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BETA1E(I,J)
     $                      -1.D0/(16.D0*PI**2)**2*B2HME(J,I))*THSE(J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            SUTRIUID(I,J)=CSFMUL(THSU,TRIU,ID,I,J)
            SUTRIUSFUQDSFUQ(I,J)=CSFMUL(THSU,TRIU,SFUQDSFUQ,I,J)
            SQIDTRIU(I,J)=CSFMUL(THSQ,ID,TRIU,I,J)
            SQFUHUDTRIU(I,J)=CSFMUL(THSQ,FUHUD,TRIU,I,J)
            SQSFUURSFUURDTRIU(I,J)=CSFMUL(THSQ,SFUURSFUURD,TRIU,I,J)
            SQCFDDRCFDDRDTRIU(I,J)=CSFMUL(THSQ,CFDDRCFDDRD,TRIU,I,J)
!
            SQIDTRID(I,J)=CSFMUL(THSQ,ID,TRID,I,J)
            SQSFUURSFUURDTRID(I,J)=CSFMUL(THSQ,SFUURSFUURD,TRID,I,J)
            SQCFDDRCFDDRDTRID(I,J)=CSFMUL(THSQ,CFDDRCFDDRD,TRID,I,J)
            SQFDHDDTRID(I,J)=CSFMUL(THSQ,FDHDD,TRID,I,J)
            SLFEHDDTRIE(I,J)=CSFMUL(THSL,FEHDD,TRIE,I,J)
            SDTRIDID(I,J)=CSFMUL(THSD,TRID,ID,I,J)
            SDTRIDCFDQDCFDQ(I,J)=CSFMUL(THSD,TRID,CFDQDCFDQ,I,J)
!
            SLIDTRIE(I,J)=CSFMUL(THSL,ID,TRIE,I,J)
            SLCFEERCFEERDTRIE(I,J)=CSFMUL(THSL,CFEERCFEERD,TRIE,I,J)
            SETRIEID(I,J)=CSFMUL(THSE,TRIE,ID,I,J)
            SETRIECFELDCFEL(I,J)=CSFMUL(THSE,TRIE,CFELDCFEL,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SUSQIDTRIUID(I,J)=CSFMUL(THSU,SQIDTRIU,ID,I,J)
            SDSQIDTRIDID(I,J)=CSFMUL(THSD,SQIDTRID,ID,I,J)
            SESLIDTRIEID(I,J)=CSFMUL(THSE,SLIDTRIE,ID,I,J)
          END DO
        END DO
        SUTSQFUHUDTRIU=TCSFMUL(THSU,SQFUHUDTRIU)
        SDTSQFDHDDTRID=TCSFMUL(THSD,SQFDHDDTRID)
        SETSLFEHDDTRIE=TCSFMUL(THSE,SLFEHDDTRIE)
        DO I=1,3
          DO J=1,3
            B1TRIU(I,J)=2.D0*G(541)**2/3.D0*SUTRIUID(I,J)
     $                 +2.D0*SUTRIUSFUQDSFUQ(I,J)
     $                 -2.D0*(G(457)**2/9.D0+4.D0*G(459)**2/3.D0)
     $                  *SUSQIDTRIUID(I,J)
     $                 +6.D0*FUHU(I,J)*SUTSQFUHUDTRIU
     $                 -2.D0*(G(541)**2/12.D0-3.D0*G(542)**2/4.D0)
     $                  *SQIDTRIU(I,J)
     $                 +4.D0*SQSFUURSFUURDTRIU(I,J)
     $                 -2.D0*SQCFDDRCFDDRDTRIU(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*CONJG(G(287))*GTPQSFTUU(I,J)
     $                  -4.D0/3.D0*GTPQSLUGTPUS(I,J)
     $                  -4.D0*THSH*FTUQGTPUS(I,J)*CONJG(G(287)))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSLUGTSUS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *CONJG(G(289))*GTQSFTUU(I,J)
     $                 +2.D0/3.D0*THSB*THSH*CONJG(G(108))*G(288)
     $                  *(4.D0*FTUQGTPUS(I,J)-GTPQSFTUU(I,J))
     $                 +6.D0*THSH*THSW*CONJG(G(108))*G(290)
     $                  *GTQSFTUU(I,J)
     $                 -4.D0*THSH*CONJG(G(108))*FTDQLDDFTUU(I,J)
     $                 +CSFMUL(THSU,TRIU,DUM1U1,I,J)
     $                 +(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)*TRIU(I,J)
     $                 +1.D0/2.D0*(THSB*MCGTPHDSQ+3.D0*THSW*MCGTHDSQ
     $                  +THSB*MSGTPHUSQ+3.D0*MSGTHUSQ)*THSH*TRIU(I,J)
     $                 +CSFMUL(THSQ,DUM1U2,TRIU,I,J)
     $                 -(13.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*TRIU(I,J)
!
            B1TRID(I,J)=-2.D0*(G(541)**2/12.D0+3.D0*G(542)**2/4.D0)
     $                  *SQIDTRID(I,J)
     $                 -2.D0*SQSFUURSFUURDTRID(I,J)
     $                 +4.D0*SQCFDDRCFDDRDTRID(I,J)
     $                 +2.D0*(G(457)**2/18.D0-4.D0*G(459)**2/3.D0)
     $                  *SDSQIDTRIDID(I,J)
     $                 +6.D0*FDHD(I,J)*SDTSQFDHDDTRID
     $                 +2.D0*FDHD(I,J)*SETSLFEHDDTRIE
     $                 -G(541)**2/3.D0*SDTRIDID(I,J)
     $                 +2.D0*SDTRIDCFDQDCFDQ(I,J)
     $                 -4.D0*THSH*CONJG(G(108))*FTUQLUDFTDD(I,J)
     $                 +2.D0/3.D0*THSB*THSH*CONJG(G(108))*G(287)
     $                  *(2.D0*FTDQGTPDS(I,J)+GTPQSFTDD(I,J))
     $                 +6.D0*THSH*THSW*CONJG(G(108))*G(289)
     $                  *GTQSFTDD(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(-THSH*CONJG(G(288))*GTPQSFTDD(I,J)
     $                  +2.D0/3.D0*GTPQSLDGTPDS(I,J)
     $                  -2.D0*THSH*FTDQGTPDS(I,J)*CONJG(G(288)))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSLDGTSDS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *CONJG(G(290))*GTQSFTDD(I,J)
     $                 +CSFMUL(THSQ,DUM1D1,TRID,I,J)
     $                 +(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)*TRID(I,J)
     $                 +1.D0/2.D0*(THSB*MCGTPHDSQ+3.D0*THSW*MCGTHDSQ
     $                  +THSB*MSGTPHUSQ+3.D0*MSGTHUSQ)*THSH*TRID(I,J)
     $                 +CSFMUL(THSD,TRID,DUM1D2,I,J)
     $                 -(7.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*TRID(I,J)
!
            B1TRIE(I,J)=2.D0*(G(541)**2/4.D0-3.D0*G(542)**2/4.D0)
     $                  *SLIDTRIE(I,J)
     $                 +4.D0*SLCFEERCFEERDTRIE(I,J)
     $                 -G(457)**2*SESLIDTRIEID(I,J)
     $                 +2.D0*FEHD(I,J)*SETSLFEHDDTRIE
     $                 +6.D0*FEHD(I,J)*SDTSQFDHDDTRID
     $                 -G(541)**2*SETRIEID(I,J)
     $                 +2.D0*SETRIECFELDCFEL(I,J)
     $                 +2.D0*THSB*THSH*CONJG(G(108))*G(287)
     $                  *(2.D0*FTELGTPES(I,J)-GTPLSFTEE(I,J))
     $                 +6.D0*THSH*THSW*CONJG(G(108))*G(289)
     $                  *GTLSFTEE(I,J)
     $                 +2.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*CONJG(G(288))*GTPLSFTEE(I,J)
     $                  -2.D0*GTPLSLEGTPES(I,J)
     $                  -2.D0*THSH*FTELGTPES(I,J)*CONJG(G(288)))
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *CONJG(G(290))*GTLSFTEE(I,J)
     $                 +CSFMUL(THSL,DUM1E1,TRIE,I,J)
     $                 +(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)*TRIE(I,J)
     $                 +1.D0/2.D0*(THSB*MCGTPHDSQ+3.D0*THSW*MCGTHDSQ
     $                  +THSB*MSGTPHUSQ+3.D0*MSGTHUSQ)*THSH*TRIE(I,J)
     $                 +CSFMUL(THSE,TRIE,DUM1E2,I,J)
     $                 -(27.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $                  *TRIE(I,J)
          END DO
        END DO
!
!Convert into form readable by RKSTP (converting 2-lp to BT notation...)
!
        DO I=1,3
          DO J=1,3
            F(399+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*B1TRIU(I,J)
     $                       +1.D0/(16.D0*PI**2)**2*(-SINB
     $                       *B2HMU(J,I)-COSB*(CONJG(G(398))
     $                       *B2YMU(J,I)+YU(J,I)*CONJG(B2GRKMUM))))
     $                       *THSU(J)
            F(408+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*B1TRID(I,J)
     $                       +1.D0/(16.D0*PI**2)**2*(-COSB
     $                       *B2HMD(J,I)-SINB*(CONJG(G(398))
     $                       *B2YMD(J,I)+YD(J,I)*CONJG(B2GRKMUM))))
     $                       *THSD(J)
            F(417+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*B1TRIE(I,J)
     $                       +1.D0/(16.D0*PI**2)**2*(-COSB
     $                       *B2HME(J,I)-SINB*(CONJG(G(398))
     $                       *B2YME(J,I)+YE(J,I)*CONJG(B2GRKMUM))))
     $                       *THSE(J)
          END DO
        END DO
      END IF
!
!Now for the (mass)^2 parameters. DUMSP are dummy matrices for the
!calculation of S'
!
      DO I=1,3
        DO J=1,3
          SUIDMU(I,J)=CSFMUL(THSU,ID,MUP,I,J)
          SQIDMQ(I,J)=CSFMUL(THSQ,ID,MQ,I,J)
          SDIDMD(I,J)=CSFMUL(THSD,ID,MD,I,J)
          SLIDML(I,J)=CSFMUL(THSL,ID,ML,I,J)
          SEIDME(I,J)=CSFMUL(THSE,ID,ME,I,J)
          SUFUHUSMU(I,J)=CSFMUL(THSU,FUHUS,MUP,I,J)
          SDFDHDSMD(I,J)=CSFMUL(THSD,FDHDS,MD,I,J)
          SQFUHUTMQ(I,J)=CSFMUL(THSQ,FUHUT,MQ,I,J)
          SQFDHDTMQ(I,J)=CSFMUL(THSQ,FDHDT,MQ,I,J)
          SEFEHDSME(I,J)=CSFMUL(THSE,FEHDS,ME,I,J)
          SLFEHDTML(I,J)=CSFMUL(THSL,FEHDT,ML,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          SQSQIDMQID(I,J)=CSFMUL(THSQ,SQIDMQ,ID,I,J)
          SUSUFUHUSMUFUHUT(I,J)=CSFMUL(THSU,SUFUHUSMU,FUHUT,I,J)
          SDSDFDHDSMDFDHDT(I,J)=CSFMUL(THSD,SDFDHDSMD,FDHDT,I,J)
          SUAUSAUT(I,J)=CSFMUL(THSU,AUS,AUT,I,J)
          DUM1Q1(I,J)=1.D0/18.D0*THSB*GTPQGTPQD(I,J)
     $               +3.D0/2.D0*THSW*GTQGTQD(I,J)
     $               +8.D0/3.D0*THGL*GTSQGTSQD(I,J)
     $               +THSH*FTUQSFTUQT(I,J)
     $               +THSH*FTDQSFTDQT(I,J)
!
          SUSUIDMUID(I,J)=CSFMUL(THSU,SUIDMU,ID,I,J)
          SQSQFUHUTMQFUHUS(I,J)=CSFMUL(THSQ,SQFUHUTMQ,FUHUS,I,J)
          SQAUTAUS(I,J)=CSFMUL(THSQ,AUT,AUS,I,J)
          DUM1U1(I,J)=8.D0/9.D0*THSB*GTPUDGTPU(I,J)
     $               +8.D0/3.D0*THGL*GTSUDGTSU(I,J)
     $               +2.D0*THSH*FTUUTFTUUS(I,J)
!
          SDSDIDMDID(I,J)=CSFMUL(THSD,SDIDMD,ID,I,J)
          SQSQFDHDTMQFDHDS(I,J)=CSFMUL(THSQ,SQFDHDTMQ,FDHDS,I,J)
          SQADTADS(I,J)=CSFMUL(THSQ,ADT,ADS,I,J)
          DUM1D1(I,J)=2.D0/9.D0*THSB*GTPDDGTPD(I,J)
     $               +8.D0/3.D0*THGL*GTSDDGTSD(I,J)
     $               +2.D0*THSH*FTDDTFTDDS(I,J)
!
          SLSLIDMLID(I,J)=CSFMUL(THSL,SLIDML,ID,I,J)
          SESEFEHDSMEFEHDT(I,J)=CSFMUL(THSE,SEFEHDSME,FEHDT,I,J)
          SDADSADT(I,J)=CSFMUL(THSD,ADS,ADT,I,J)
          SEAESAET(I,J)=CSFMUL(THSE,AES,AET,I,J)
          DUM1L1(I,J)=1.D0/2.D0*THSB*GTPLGTPLD(I,J)
     $               +3.D0/2.D0*THSW*GTLGTLD(I,J)
     $               +THSH*FTELSFTELT(I,J)
!
          SESEIDMEID(I,J)=CSFMUL(THSE,SEIDME,ID,I,J)
          SLSLFEHDTMLFEHDS(I,J)=CSFMUL(THSL,SLFEHDTML,FEHDS,I,J)
          SLAETAES(I,J)=CSFMUL(THSL,AET,AES,I,J)
          DUM1E1(I,J)=2.D0*THSB*GTPEDGTPE(I,J)
     $               +2.D0*THSH*FTEETFTEES(I,J)
        END DO
      END DO
      SUTSUIDMU=TCSFMUL(THSU,SUIDMU)
      SQTSQIDMQ=TCSFMUL(THSQ,SQIDMQ)
      SDTSDIDMD=TCSFMUL(THSD,SDIDMD)
      SLTSLIDML=TCSFMUL(THSL,SLIDML)
      SETSEIDME=TCSFMUL(THSE,SEIDME)
      SQTMQ=TCSFMUL(THSQ,MQ)
      SUTMU=TCSFMUL(THSU,MUP)
      SDTMD=TCSFMUL(THSD,MD)
      SLTML=TCSFMUL(THSL,ML)
      SETME=TCSFMUL(THSE,ME)
!
      IF(THHH.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            SUFUQTFUQSMU(I,J)=CSFMUL(THSU,FUQTFUQS,MUP,I,J)
            SQFUURSFUURTMQ(I,J)=CSFMUL(THSQ,FUURSFUURT,MQ,I,J)
            SDMTSFDSMTSFDT(I,J)=CSFMUL(THSD,MTSFDS,MTSFDT,I,J)
            SEMTSFESMTSFET(I,J)=CSFMUL(THSE,MTSFES,MTSFET,I,J)
!
            SQFDDRSFDDRTMQ(I,J)=CSFMUL(THSQ,FDDRSFDDRT,MQ,I,J)
            SDFDQTFDQSMD(I,J)=CSFMUL(THSD,FDQTFDQS,MD,I,J)
            SLFEERSFEERTML(I,J)=CSFMUL(THSL,FEERSFEERT,ML,I,J)
            SEFELTFELSME(I,J)=CSFMUL(THSE,FELTFELS,ME,I,J)
            SUMTSFUSMTSFUT(I,J)=CSFMUL(THSU,MTSFUS,MTSFUT,I,J)
          END DO
        END DO
        SUTSUFUQTFUQSMU=TCSFMUL(THSU,SUFUQTFUQSMU)
        SQTSQFUURSFUURTMQ=TCSFMUL(THSQ,SQFUURSFUURTMQ)
        SQTSUAUSAUT=TCSFMUL(THSQ,SUAUSAUT)
        SQTSDMTSFDSMTSFDT=TCSFMUL(THSQ,SDMTSFDSMTSFDT)
        SLTSEMTSFESMTSFET=TCSFMUL(THSL,SEMTSFESMTSFET)
!
        SQTSQFDDRSFDDRTMQ=TCSFMUL(THSQ,SQFDDRSFDDRTMQ)
        SDTSDFDQTFDQSMD=TCSFMUL(THSD,SDFDQTFDQSMD)
        SLTSLFEERSFEERTML=TCSFMUL(THSL,SLFEERSFEERTML)
        SETSEFELTFELSME=TCSFMUL(THSE,SEFELTFELSME)
        SQTSUMTSFUSMTSFUT=TCSFMUL(THSQ,SUMTSFUSMTSFUT)
        SQTSDADSADT=TCSFMUL(THSQ,SDADSADT)
        SLTSEAESAET=TCSFMUL(THSL,SEAESAET)
!
        BMHUPMT=3.D0/2.D0*(G(457)**2+G(458)**2)*G(61)
     $         -G(457)**2*G(62)
     $         -2.D0*G(457)**2*SUTSUIDMU
     $         +6.D0*SUTSUFUQTFUQSMU
     $         +G(457)**2*SQTSQIDMQ
     $         +6.D0*SQTSQFUURSFUURTMQ
     $         +G(457)**2*SDTSDIDMD
     $         -G(457)**2*SLTSLIDML
     $         +G(457)**2*SETSEIDME
     $         +6.D0*SQTSUAUSAUT
     $         +6.D0*SQTSDMTSFDSMTSFDT
     $         +2.D0*SLTSEMTSFESMTSFET
     $         -2.D0*THSH*MMUSQ*(THSB*MGTPHUSQ+3.D0*THSW*MGTHUSQ)
     $         -2.D0*THSH*(THSB*M1PM1PSQ*MGTPHUSQ+3.D0*THSW
     $          *M2PM2PSQ*MGTHUSQ)
     $         -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(61)
     $         +(6.D0*TFUDFU+THSB*THSH*MGTPHUSQ+3.D0*THSW*THSH
     $          *MGTHUSQ)*G(61)
!
        BMHDPMT=-G(457)**2*G(61)
     $         +3.D0/2.D0*(G(457)**2+G(458)**2)*G(62)
     $         +2.D0*G(457)**2*SUTSUIDMU
     $         -G(457)**2*SQTSQIDMQ
     $         +6.D0*SQTSQFDDRSFDDRTMQ
     $         -G(457)**2*SDTSDIDMD
     $         +6.D0*SDTSDFDQTFDQSMD
     $         +G(457)**2*SLTSLIDML
     $         +2.D0*SLTSLFEERSFEERTML
     $         -G(457)**2*SETSEIDME
     $         +2.D0*SETSEFELTFELSME
     $         +6.D0*SQTSUMTSFUSMTSFUT
     $         +6.D0*SQTSDADSADT
     $         +2.D0*SLTSEAESAET
     $         -2.D0*THSH*MMUSQ*(THSB*MGTPHDSQ+3.D0*THSW*MGTHDSQ)
     $         -2.D0*THSH*(THSB*M1PM1PSQ*MGTPHDSQ+3.D0*THSW
     $          *M2PM2PSQ*MGTHDSQ)
     $         -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(62)
     $         +(6.D0*TFDDFD+2.D0*TFEDFE+THSB*THSH*MGTPHDSQ
     $          +3.D0*THSW*THSH*MGTHDSQ)*G(62)
!
        F(61)=1.D0/16.D0/PI**2*BMHUPMT+1.D0/(16.D0*PI**2)**2
     $       *(B2HUM+CONJG(G(398))*B2GRKMUM+G(398)*CONJG(B2GRKMUM))
        F(62)=1.D0/16.D0/PI**2*BMHDPMT+1.D0/(16.D0*PI**2)**2
     $       *(B2HDM+CONJG(G(398))*B2GRKMUM+G(398)*CONJG(B2GRKMUM))
!
!Now the soft mass matrices
!
        DO I=1,3
          DO J=1,3
            SQMTSFUTMTSFUS(I,J)=CSFMUL(THSQ,MTSFUT,MTSFUS,I,J)
            SQMTSFDTMTSFDS(I,J)=CSFMUL(THSQ,MTSFDT,MTSFDS,I,J)
            SLMTSFETMTSFES(I,J)=CSFMUL(THSL,MTSFET,MTSFES,I,J)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            BETA1MQ(I,J)=(1.D0/3.D0*G(457)**2*ID(I,J)+2.D0
     $                   *FUURSFUURT(I,J))*G(61)
     $                  +(-1.D0/3.D0*G(457)**2*ID(I,J)+2.D0
     $                   *FDDRSFDDRT(I,J))*G(62)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +(G(457)**2/18.D0+3.D0*G(458)**2/2.D0
     $                   +8.D0*G(459)**2/3.D0)*SQSQIDMQID(I,J)
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -1.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SUSUFUHUSMUFUHUT(I,J)
     $                  +2.D0*SDSDFDHDSMDFDHDT(I,J)
     $                  +2.D0*SUAUSAUT(I,J)
     $                  +2.D0*SUMTSFUSMTSFUT(I,J)
     $                  +2.D0*SDADSADT(I,J)
     $                  +2.D0*SDMTSFDSMTSFDT(I,J)
     $                  -2.D0/9.D0*THSB*M1PM1PSQ*GTPQGTPQD(I,J)
     $                  -6.D0*THSW*M2PM2PSQ*GTQGTQD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSQGTSQD(I,J)
     $                  -4.D0*THSH*MMUSQ*(FTUQSFTUQT(I,J)
     $                   +FTDQSFTDQT(I,J))
     $                  -(G(1)**2/10.D0+9.D0/2.D0*G(2)**2
     $                   +8.D0*G(3)**2)*MQ(I,J)
     $                  +CSFMUL(THSQ,DUM1Q1,MQ,I,J)
     $                  +CSFMUL(THSQ,MQ,DUM1Q1,I,J)
!
            BETA1MU(I,J)=(-4.D0/3.D0*G(457)**2*ID(I,J)+4.D0
     $                   *FUQTFUQS(I,J))*G(61)
     $                  +4.D0/3.D0*G(457)**2*ID(I,J)*G(62)
     $                  +8.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +8.D0/3.D0*(G(457)**2/3.D0+G(459)**2)
     $                   *SUSUIDMUID(I,J)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +4.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFUHUTMQFUHUS(I,J)
     $                  +4.D0*SQAUTAUS(I,J)
     $                  +4.D0*SQMTSFUTMTSFUS(I,J)
     $                  -32.D0/9.D0*THSB*M1PM1PSQ*GTPUDGTPU(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSUDGTSU(I,J)
     $                  -8.D0*THSH*MMUSQ*FTUUTFTUUS(I,J)
     $                  -(8.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MUP(I,J)
     $                  +CSFMUL(THSU,DUM1U1,MUP,I,J)
     $                  +CSFMUL(THSU,MUP,DUM1U1,I,J)
!
            BETA1MD(I,J)=2.D0/3.D0*G(457)**2*ID(I,J)*G(61)
     $                  +(-2.D0/3.D0*G(457)**2*ID(I,J)+4.D0
     $                   *FDQTFDQS(I,J))*G(62)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +2.D0/3.D0*(G(457)**2/3.D0+4.D0*G(459)**2)
     $                   *SDSDIDMDID(I,J)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFDHDTMQFDHDS(I,J)
     $                  +4.D0*SQADTADS(I,J)
     $                  +4.D0*SQMTSFDTMTSFDS(I,J)
     $                  -8.D0/9.D0*THSB*M1PM1PSQ*GTPDDGTPD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSDDGTSD(I,J)
     $                  -8.D0*THSH*MMUSQ*FTDDTFTDDS(I,J)
     $                  -(2.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MD(I,J)
     $                  +CSFMUL(THSD,DUM1D1,MD,I,J)
     $                  +CSFMUL(THSD,MD,DUM1D1,I,J)
!
            BETA1ML(I,J)=-G(457)**2*ID(I,J)*G(61)
     $                  +(G(457)**2*ID(I,J)+2.D0*FEERSFEERT(I,J))*G(62)
     $                  +2.D0*G(457)**2*ID(I,J)*SUTMU
     $                  -G(457)**2*ID(I,J)*SQTMQ
     $                  -G(457)**2*ID(I,J)*SDTMD
     $                  +G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/2.D0*(G(457)**2+3.D0*G(458)**2)
     $                   *SLSLIDMLID(I,J)
     $                  -G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SESEFEHDSMEFEHDT(I,J)
     $                  +2.D0*SEAESAET(I,J)
     $                  +2.D0*SEMTSFESMTSFET(I,J)
     $                  -2.D0*THSB*M1PM1PSQ*GTPLGTPLD(I,J)
     $                  -6.D0*THGL*M2PM2PSQ*GTLGTLD(I,J)
     $                  -4.D0*THSH*MMUSQ*FTELSFTELT(I,J)
     $                  -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*ML(I,J)
     $                  +CSFMUL(THSL,DUM1L1,ML,I,J)
     $                  +CSFMUL(THSL,ML,DUM1L1,I,J)
!
            BETA1ME(I,J)=2.D0*G(457)**2*ID(I,J)*G(61)
     $                  +(-2.D0*G(457)**2*ID(I,J)+4.D0*FELTFELS(I,J))
     $                   *G(62)
     $                  -4.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -2.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*G(457)**2*SESEIDMEID(I,J)
     $                  +4.D0*SLSLFEHDTMLFEHDS(I,J)
     $                  +4.D0*SLAETAES(I,J)
     $                  +4.D0*SLMTSFETMTSFES(I,J)
     $                  -8.D0*THSB*M1PM1PSQ*GTPEDGTPE(I,J)
     $                  -8.D0*THSH*MMUSQ*FTEETFTEES(I,J)
     $                  -18.D0*G(1)**2/5.D0*ME(I,J)
     $                  +CSFMUL(THSE,DUM1E1,ME,I,J)
     $                  +CSFMUL(THSE,ME,DUM1E1,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            F(62+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MQ(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MQM(I,J)
            F(71+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ML(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MLM(I,J)
            F(80+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MU(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MUM(I,J)
            F(89+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MD(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MDM(I,J)
            F(98+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ME(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MEM(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            SUSFUQTSFUQSMU(I,J)=CSFMUL(THSU,SFUQTSFUQS,MUP,I,J)
            SQSFUURSSFUURTMQ(I,J)=CSFMUL(THSQ,SFUURSSFUURT,MQ,I,J)
            SQCFDDRSCFDDRTMQ(I,J)=CSFMUL(THSQ,CFDDRSCFDDRT,MQ,I,J)
            SDCFDQTCFDQSMD(I,J)=CSFMUL(THSD,CFDQTCFDQS,MD,I,J)
            SLCFEERSCFEERTML(I,J)=CSFMUL(THSL,CFEERSCFEERT,ML,I,J)
            SECFELTCFELSME(I,J)=CSFMUL(THSE,CFELTCFELS,ME,I,J)
            SUTRIUTRIUD(I,J)=CSFMUL(THSU,TRIU,TRIUD,I,J)
            SDTRIDTRIDD(I,J)=CSFMUL(THSD,TRID,TRIDD,I,J)
            SETRIETRIED(I,J)=CSFMUL(THSE,TRIE,TRIED,I,J)
          END DO
        END DO
        SUTSUSFUQTSFUQSMU=TCSFMUL(THSU,SUSFUQTSFUQSMU)
        SQTSQSFUURSSFUURTMQ=TCSFMUL(THSQ,SQSFUURSSFUURTMQ)
        SQTSQCFDDRSCFDDRTMQ=TCSFMUL(THSQ,SQCFDDRSCFDDRTMQ)
        SDTSDCFDQTCFDQSMD=TCSFMUL(THSD,SDCFDQTCFDQSMD)
        SLTSLCFEERSCFEERTML=TCSFMUL(THSL,SLCFEERSCFEERTML)
        SETSECFELTCFELSME=TCSFMUL(THSE,SECFELTCFELSME)
        SQTSUTRIUTRIUD=TCSFMUL(THSQ,SUTRIUTRIUD)
        SQTSDTRIDTRIDD=TCSFMUL(THSQ,SDTRIDTRIDD)
        SLTSETRIETRIED=TCSFMUL(THSL,SETRIETRIED)
!      
        BMHUD=3.D0/2.D0*(G(543)**2+G(544)**2)*G(427)
     $       +2.D0*G(541)**2*SUTSUIDMU
     $       +6.D0*SUTSUSFUQTSFUQSMU
     $       -G(541)**2*SQTSQIDMQ
     $       +6.D0*SQTSQSFUURSSFUURTMQ
     $       +6.D0*SQTSQCFDDRSCFDDRTMQ
     $       -G(541)**2*SDTSDIDMD
     $       +6.D0*SDTSDCFDQTCFDQSMD
     $       +G(541)**2*SLTSLIDML
     $       +2.D0*SLTSLCFEERSCFEERTML
     $       -G(541)**2*SETSEIDME
     $       +2.D0*SETSECFELTCFELSME
     $       +6.D0*SQTSUTRIUTRIUD+6.D0*SQTSDTRIDTRIDD
     $        +2.D0*SLTSETRIETRIED
     $       -2.D0*THSH*MMUSQ*(THSB*(MSGTPHUSQ+MCGTPHDSQ)
     $        +3.D0*THSW*(MSGTHUSQ+MCGTHDSQ))
     $       -2.D0*THSH*(THSB*M1PM1PSQ*(MSGTPHUSQ+MCGTPHDSQ)
     $        +3.D0*THSW*M2PM2PSQ*(MSGTHUSQ+MCGTHDSQ))
     $       -1.D0/2.D0*(-4.D0*THSH*THSB*CONJG(G(108))*G(287)*G(288)
     $        *(G(31)+(0.D0,1.D0)*G(599))-12.D0*THSH*THSW*CONJG(G(108))
     $        *G(289)*G(290)*(G(32)+(0.D0,1.D0)*G(600)))
     $       -1.D0/2.D0*(-4.D0*THSH*THSB*G(108)*CONJG(G(287))
     $        *CONJG(G(288))*(G(31)-(0.D0,1.D0)*G(599))-12.D0*THSH*THSW
     $        *G(108)*CONJG(G(289))*CONJG(G(290))*(G(32)-(0.D0,1.D0)
     $        *G(600)))
     $       -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(427)
     $       +(6.D0*TLUDLU+THSB*THSH*MSGTPHUSQ+3.D0*THSW*THSH*MSGTHUSQ
     $        +6.D0*TLDDLD+2.D0*TLEDLE+THSB*THSH*MCGTPHDSQ
     $        +3.D0*THSW*THSH*MCGTHDSQ)*G(427)
!
        F(427)=1.D0/16.D0/PI**2*BMHUD+1.D0/(16.D0*PI**2)**2
     $        *(SINB**2*B2HUM+COSB**2*B2HDM+CONJG(G(398))*B2GRKMUM
     $         +G(398)*CONJG(B2GRKMUM)+SINB*COSB*(BETA2B+CONJG(BETA2B)))
!
!Now the soft mass matrices
!
        DO I=1,3
          DO J=1,3
            SUTRIUSTRIUT(I,J)=CSFMUL(THSU,TRIUS,TRIUT,I,J)
            SDTRIDSTRIDT(I,J)=CSFMUL(THSD,TRIDS,TRIDT,I,J)
            SQTRIUTTRIUS(I,J)=CSFMUL(THSQ,TRIUT,TRIUS,I,J)
            SQTRIDTTRIDS(I,J)=CSFMUL(THSQ,TRIDT,TRIDS,I,J)
            SETRIESTRIET(I,J)=CSFMUL(THSE,TRIES,TRIET,I,J)
            SLTRIETTRIES(I,J)=CSFMUL(THSL,TRIET,TRIES,I,J)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            BETA1MQ(I,J)=(-1.D0/3.D0*G(541)**2*ID(I,J)+2.D0
     $                   *SFUURSSFUURT(I,J)+2.D0
     $                   *CFDDRSCFDDRT(I,J))*G(427)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +(G(457)**2/18.D0+3.D0*G(458)**2/2.D0
     $                   +8.D0*G(459)**2/3.D0)*SQSQIDMQID(I,J)
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -1.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SUSUFUHUSMUFUHUT(I,J)
     $                  +2.D0*SDSDFDHDSMDFDHDT(I,J)
     $                  +2.D0*SUTRIUSTRIUT(I,J)
     $                  +2.D0*SDTRIDSTRIDT(I,J)
     $                  -2.D0/9.D0*THSB*M1PM1PSQ*GTPQGTPQD(I,J)
     $                  -6.D0*THSW*M2PM2PSQ*GTQGTQD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSQGTSQD(I,J)
     $                  -4.D0*THSH*MMUSQ*(FTUQSFTUQT(I,J)
     $                   +FTDQSFTDQT(I,J))
     $                  -(G(1)**2/10.D0+9.D0/2.D0*G(2)**2
     $                   +8.D0*G(3)**2)*MQ(I,J)
     $                  +CSFMUL(THSQ,DUM1Q1,MQ,I,J)
     $                  +CSFMUL(THSQ,MQ,DUM1Q1,I,J)
!
            BETA1MU(I,J)=(4.D0/3.D0*G(541)**2*ID(I,J)+4.D0
     $                   *SFUQTSFUQS(I,J))*G(427)
     $                  +8.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +8.D0/3.D0*(G(457)**2/3.D0+G(459)**2)
     $                   *SUSUIDMUID(I,J)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +4.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFUHUTMQFUHUS(I,J)
     $                  +4.D0*SQTRIUTTRIUS(I,J)
     $                  -32.D0/9.D0*THSB*M1PM1PSQ*GTPUDGTPU(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSUDGTSU(I,J)
     $                  -8.D0*THSH*MMUSQ*FTUUTFTUUS(I,J)
     $                  -(8.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MUP(I,J)
     $                  +CSFMUL(THSU,DUM1U1,MUP,I,J)
     $                  +CSFMUL(THSU,MUP,DUM1U1,I,J)
!
            BETA1MD(I,J)=(-2.D0/3.D0*G(541)**2*ID(I,J)+4.D0
     $                   *CFDQTCFDQS(I,J))*G(427)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +2.D0/3.D0*(G(457)**2/3.D0+4.D0*G(459)**2)
     $                   *SDSDIDMDID(I,J)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFDHDTMQFDHDS(I,J)
     $                  +4.D0*SQTRIDTTRIDS(I,J)
     $                  -8.D0/9.D0*THSB*M1PM1PSQ*GTPDDGTPD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSDDGTSD(I,J)
     $                  -8.D0*THSH*MMUSQ*FTDDTFTDDS(I,J)
     $                  -(2.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MD(I,J)
     $                  +CSFMUL(THSD,DUM1D1,MD,I,J)
     $                  +CSFMUL(THSD,MD,DUM1D1,I,J)
!
            BETA1ML(I,J)=(G(541)**2*ID(I,J)+2.D0*CFEERSCFEERT(I,J))
     $                   *G(427)
     $                  +2.D0*G(457)**2*ID(I,J)*SUTMU
     $                  -G(457)**2*ID(I,J)*SQTMQ
     $                  -G(457)**2*ID(I,J)*SDTMD
     $                  +G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/2.D0*(G(457)**2+3.D0*G(458)**2)
     $                   *SLSLIDMLID(I,J)
     $                  -G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SESEFEHDSMEFEHDT(I,J)
     $                  +2.D0*SETRIESTRIET(I,J)
     $                  -2.D0*THSB*M1PM1PSQ*GTPLGTPLD(I,J)
     $                  -6.D0*THGL*M2PM2PSQ*GTLGTLD(I,J)
     $                  -4.D0*THSH*MMUSQ*FTELSFTELT(I,J)
     $                  -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*ML(I,J)
     $                  +CSFMUL(THSL,DUM1L1,ML,I,J)
     $                  +CSFMUL(THSL,ML,DUM1L1,I,J)
!
            BETA1ME(I,J)=(-2.D0*G(541)**2*ID(I,J)
     $                   +4.D0*CFELTCFELS(I,J))*G(427)
     $                  -4.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -2.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*G(457)**2*SESEIDMEID(I,J)
     $                  +4.D0*SLSLFEHDTMLFEHDS(I,J)
     $                  +4.D0*SLTRIETTRIES(I,J)
     $                  -8.D0*THSB*M1PM1PSQ*GTPEDGTPE(I,J)
     $                  -8.D0*THSH*MMUSQ*FTEETFTEES(I,J)
     $                  -18.D0*G(1)**2/5.D0*ME(I,J)
     $                  +CSFMUL(THSE,DUM1E1,ME,I,J)
     $                  +CSFMUL(THSE,ME,DUM1E1,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            F(62+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MQ(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MQM(I,J)
            F(71+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ML(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MLM(I,J)
            F(80+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MU(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MUM(I,J)
            F(89+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MD(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MDM(I,J)
            F(98+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ME(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MEM(I,J)
          END DO
        END DO
      END IF
!
      IF(THHH.EQ.0)THEN
!
!Finally we have the running of the Higgs Quartic Coupling and SM VEV.
!Programmed here is the MS-bar running. It therefore needs the MS-bar
!gauge and Yukawas.
!The gauge couplings and Yukawas needed to be converted to MS-bar using
!the Martin and Vaughn conversion in hep-ph/9308222.
!The following is after the conversion, so all Yukawas and Gauge
!couplings are still in the DR-bar scheme.
!
        BETALAM1=12*G(429)**2-(9.D0/5.D0*G(1)**2+9.D0*G(2)**2)
     $           *G(429)+9.D0/4.D0*(3.D0/25.D0*G(1)**4+2.D0/5.D0
     $           *G(1)**2*G(2)**2+G(2)**4)+4.D0*Y2*G(429)-4*H
        IF(SW2LP.EQ.1)THEN
          BETALAM2=-78.D0*G(429)**3+18.D0*(3.D0/5.D0*G(1)**2
     $            +3.D0*G(2)**2)*G(429)**2-((265.D0/8.D0-10*NG)
     $            *G(2)**4-117.D0/20.D0*G(1)**2*G(2)**2
     $            -9.D0/25.D0*(229.D0/24.D0+50.D0/9.D0*NG)*G(1)**4)
     $            *G(429)+(473.D0/8.D0-8.D0*NG)*G(2)**6-3.D0/5.D0
     $            *(121.D0/24.D0+8.D0/3.D0*NG)*G(1)**2*G(2)**4
     $            -9.D0/25.D0*(239.D0/24.D0+40.D0/9.D0*NG)
     $            *G(1)**4*G(2)**2-27.D0/125.D0*(59.D0/24.D0
     $            +40.D0/9.D0*NG)*G(1)**6+(-14.D0/5.D0*G(1)**2
     $            +18.D0*G(2)**2-128.D0*G(3)**2)*TLYUDLYU2
     $            +(34.D0/5.D0*G(1)**2+18.D0*G(2)**2-128.D0
     $            *G(3)**2)*TLYDDLYD2+(-42.D0/5.D0*G(1)**2
     $            +6.D0*G(2)**2)*TLYEDLYE2-3.D0/2.D0*G(2)**4
     $            *Y2+G(429)*((83.D0/10.D0*G(1)**2+27.D0/2.D0
     $            *G(2)**2+112.D0*G(3)**2)*TLYUDLYU+(-1.D0/10.D0
     $            *G(1)**2+27.D0/2.D0*G(2)**2+112.D0*G(3)**2)
     $            *TLYDDLYD+(93.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $            *TLYEDLYE)+3.D0/5.D0*G(1)**2*((-57.D0/10.D0
     $            *G(1)**2+21.D0*G(2)**2)*TLYUDLYU+(3.D0/2.D0
     $            *G(1)**2+9.D0*G(2)**2)*TLYDDLYD+(-15.D0/2.D0
     $            *G(1)**2+11.D0*G(2)**2)*TLYEDLYE)-24.D0
     $            *G(429)**2*Y2-G(429)*H+6.D0*G(429)
     $            *TLYUDLYULYDDLYD+20.D0*(3.D0*TLYUDLYU3+3.D0*TLYDDLYD3
     $            +TLYEDLYE3)-12.D0*TDUMLUD
        END IF
!
        F(429)=1.D0/(16.D0*PI**2)*BETALAM1+1.D0/(16.D0*PI**2)**2
     $         *BETALAM2
!
!Calculate the betas for the standard model vev.
!As with lambda this is the MS-bar running with DR-bar inputs except
!v and lambda
!
      BETAVEV1=9.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-Y2
      IF(SW2LP.EQ.1)THEN
        BETAVEV2=-3.D0/2.D0*G(429)**2-(83.D0/40.D0*G(1)**2+27.D0/8.D0
     $           *G(2)**2+28.D0*G(3)**2)*TYUDYU-(-1.D0/40.D0*G(1)**2
     $           +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TYDDYD
     $           -(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TYEDYE+CHI4
     $           -27.D0/80.D0*G(1)**2*G(2)**2-(93.D0/800.D0+1.D0/2.D0
     $           *NG)*G(1)**4+(463.D0/32.D0-5.D0/2.D0*NG)*G(2)**4
      END IF
!
      F(428)=G(428)*(1.D0/(16.D0*PI**2)*BETAVEV1+1.D0/(16.D0*PI**2)**2
     $      *BETAVEV2)
!
      END IF
!
!Now the tilde terms are calculated.
!
!
      DO I=1,3
        DO J=1,3
          DUMTUQ1(I,J)=4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $                 +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
     $                 +THSH*SQFTUQDFTUQ(I,J)
          DUMTUQ2(I,J)=SUFTUUFTUUD(I,J)+SQFTUQDFTUQ(I,J)
          IF(THHH.EQ.0)THEN
            DUMTUQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(289)
     $                   +1.D0/3.D0*THSB*GTPQS(I,J)*G(287)
          ELSE
            DUMTUQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(204)
     $                   +1.D0/3.D0*THSB*GTPQS(I,J)*G(184)
          END IF
          DUMTUQ4(I,J)=3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $                 +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $                 +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
          DUMTUQ5(I,J)=FTUQFTUQD(I,J)+FTDQFTDQD(I,J)
!
          DUMTDQ1(I,J)=1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $                 +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
     $                 +THSH*SQFTDQDFTDQ(I,J)
          DUMTDQ2(I,J)=3.D0*SDFTDDFTDDD(I,J)+SEFTEEFTEED(I,J)
     $                 +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
          IF(THHH.EQ.0)THEN
            DUMTDQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(290)
     $                   -1.D0/3.D0*THSB*GTPQS(I,J)*G(288)
          ELSE
            DUMTDQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(205)
     $                   -1.D0/3.D0*THSB*GTPQS(I,J)*G(185)
          END IF
          DUMTDQ4(I,J)=3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $                 +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $                 +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
          DUMTDQ5(I,J)=FTUQFTUQD(I,J)+FTDQFTDQD(I,J)
!
          DUMTEL1(I,J)=THSB*SEGTPESGTPET(I,J)
     $                 +THSH*SLFTELDFTEL(I,J)
          IF(THHH.EQ.0)THEN
            DUMTEL3(I,J)=-3.D0*THSW*GTLS(I,J)*G(290)
     $                   +THSB*GTPLS(I,J)*G(288)
          ELSE
            DUMTEL3(I,J)=-3.D0*THSW*GTLS(I,J)*G(205)
     $                   +THSB*GTPLS(I,J)*G(185)
          END IF
          DUMTEL4(I,J)=3.D0*THSW*GTLSGTLT(I,J)
     $                 +THSB*GTPLSGTPLT(I,J)
          DUMTEL5(I,J)=FTELFTELD(I,J)
!
          IF(THHH.EQ.0)THEN
            DUMTUU2(I,J)=LULUD(I,J)+LDLDD(I,J)
          ELSE
            DUMTUU2(I,J)=FUFUD(I,J)+FDFDD(I,J)
          END IF
          DUMTUU3(I,J)=THSH*(SUFTUUFTUUD(I,J)+SDFTDDFTDDD(I,J))
     $                 +3.D0/2.D0*THSW*SQGTQTGTQS(I,J)
     $                 +1.D0/18.D0*THSB*SQGTPQTGTPQS(I,J)
     $                 +8.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMTUU4(I,J)=8.D0/9.D0*THSB*GTPUTGTPUS(I,J)
     $                 +8.D0/3.D0*THGL*GTSUTGTSUS(I,J)
!
          DUMTDD2(I,J)=DUMTUU2(I,J)
          DUMTDD3(I,J)=DUMTUU3(I,J)
          DUMTDD4(I,J)=2.D0/9.D0*THSB*GTPDTGTPDS(I,J)
     $                 +8.D0/3.D0*THGL*GTSDTGTSDS(I,J)
!
          DUMTEE3(I,J)=THSH*SEFTEEFTEED(I,J)
     $                 +3.D0/2.D0*THSW*SLGTLTGTLS(I,J)
     $                 +1.D0/2.D0*THSB*SLGTPLTGTPLS(I,J)
        END DO
      END DO
      TDUMTUQ2=CTRACE(DUMTUQ2)
      TDUMTDQ2=CTRACE(DUMTDQ2)
      TDUMTEL2=TDUMTDQ2
      TDUMTUU1=TDUMTUQ2
      TDUMTDD1=TDUMTDQ2
      TDUMTEE1=TDUMTEL2
      DO I=1,3
        DO J=1,3
          IF(THHH.EQ.0)THEN
            BFTUQ(I,J)=FTUQLUDLU(I,J)
     $                 +CMATMUL(0,FTUQ,DUMTUQ1,I,J)
     $                 +3.D0/2.D0*THSH*FTUQ(I,J)*TDUMTUQ2
     $                 +1.D0/4.D0*THSH*FTUQ(I,J)
     $                  *(3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $                 +CMATMUL(0,DUMTUQ3,LU,I,J)
     $                 -4.D0/9.D0*THSB
     $                  *CMATMUL(0,GTPQS,SUFTUUGTPUT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *CMATMUL(0,GTSQS,SUFTUUGTSUT,I,J)
     $                 +CSFMUL(THSQ,DUMTUQ4,FTUQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTUQ5,FTUQ,I,J)
     $                 -FTUQ(I,J)*(5.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDQ(I,J)=FTDQLDDLD(I,J)
     $                 +CMATMUL(0,FTDQ,DUMTDQ1,I,J)
     $                 +1.D0/2.D0*THSH*FTDQ(I,J)*TDUMTDQ2
     $                 +1.D0/4.D0*THSH*FTDQ(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +CMATMUL(0,DUMTDQ3,LD,I,J)
     $                 +2.D0/9.D0*THSB
     $                  *CMATMUL(0,GTPQS,SDFTDDGTPDT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *CMATMUL(0,GTSQS,SDFTDDGTSDT,I,J)
     $                 +CSFMUL(THSQ,DUMTDQ4,FTDQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTDQ5,FTDQ,I,J)
     $                 -FTDQ(I,J)*(13.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEL(I,J)=FTELLEDLE(I,J)
     $                 +CMATMUL(0,FTEL,DUMTEL1,I,J)
     $                 +1.D0/2.D0*THSH*FTEL(I,J)*TDUMTEL2
     $                 +1.D0/4.D0*THSH*FTEL(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +CMATMUL(0,DUMTEL3,LE,I,J)
     $                 -2.D0*THSB*CMATMUL(0,GTPLS,SEFTEEGTPET,I,J)
     $                 +1.D0/2.D0*CSFMUL(THSL,DUMTEL4,FTEL,I,J)
     $                 +THSH*CSFMUL(THSL,DUMTEL5,FTEL,I,J)
     $                 -FTEL(I,J)*(9.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BFTUU(I,J)=3.D0/2.D0*THSH*FTUU(I,J)*TDUMTUU1
     $                 +1.D0/4.D0*THSH*FTUU(I,J)
     $                  *(3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTUU2,FTUU,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTUU3,FTUU,I,J)
     $                 -4.D0/9.D0*THSB*CMATMUL(0,SQGTPQTFTUQ,GTPUS,I,J)
     $                 -16.D0/3.D0*THGL*CMATMUL(0,SQGTSQTFTUQ,GTSUS,I,J)
     $                 -4.D0/3.D0*THSB*G(287)*LUGTPUS(I,J)
     $                 +2.D0*THSH*CMATMUL(0,SUFTUUFTUUD,FTUU,I,J)
     $                 +CSFMUL(THSU,FTUU,DUMTUU4,I,J)
     $                 -FTUU(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDD(I,J)=1.D0/2.D0*THSH*FTDD(I,J)*TDUMTDD1
     $                 +1.D0/4.D0*THSH*FTDD(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTDD2,FTDD,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTDD3,FTDD,I,J)
     $                 +2.D0/9.D0*THSB*CMATMUL(0,SQGTPQTFTDQ,GTPDS,I,J)
     $                 -16.D0/3.D0*THGL*CMATMUL(0,SQGTSQTFTDQ,GTSDS,I,J)
     $                 -2.D0/3.D0*THSB*G(288)*LDGTPDS(I,J)
     $                 +2.D0*THSH*CMATMUL(0,SDFTDDFTDDD,FTDD,I,J)
     $                 +CSFMUL(THSD,FTDD,DUMTDD4,I,J)
     $                 -FTDD(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEE(I,J)=1.D0/2.D0*THSH*FTEE(I,J)*TDUMTEE1
     $                 +1.D0/4.D0*THSH*FTEE(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,LELED,FTEE,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTEE3,FTEE,I,J)
     $                 -2.D0*THSB*CMATMUL(0,SLGTPLTFTEL,GTPES,I,J)
     $                 -2.D0*THSB*G(288)*LEGTPES(I,J)
     $                 +2.D0*THSH*CMATMUL(0,SEFTEEFTEED,FTEE,I,J)
     $                 +2.D0*THSB*CMATMUL(0,SEFTEEGTPET,GTPES,I,J)
     $                 -FTEE(I,J)*(9.D0/10.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2)
          ELSE
            BFTUQ(I,J)=FTUQFUDFU(I,J)
     $                 +CMATMUL(0,FTUQ,DUMTUQ1,I,J)
     $                 +3.D0/2.D0*THSH*FTUQ(I,J)*TDUMTUQ2
     $                 +1.D0/4.D0*THSH*FTUQ(I,J)
     $                  *(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $                 +CMATMUL(0,DUMTUQ3,FU,I,J)
     $                 -4.D0/9.D0*THSB
     $                  *CMATMUL(0,GTPQS,SUFTUUGTPUT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *CMATMUL(0,GTSQS,SUFTUUGTSUT,I,J)
     $                 +CSFMUL(THSQ,DUMTUQ4,FTUQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTUQ5,FTUQ,I,J)
     $                 -FTUQ(I,J)*(5.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDQ(I,J)=FTDQFDDFD(I,J)
     $                 +CMATMUL(0,FTDQ,DUMTDQ1,I,J)
     $                 +1.D0/2.D0*THSH*FTDQ(I,J)*TDUMTDQ2
     $                 +1.D0/4.D0*THSH*FTDQ(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +CMATMUL(0,DUMTDQ3,FD,I,J)
     $                 +2.D0/9.D0*THSB
     $                  *CMATMUL(0,GTPQS,SDFTDDGTPDT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *CMATMUL(0,GTSQS,SDFTDDGTSDT,I,J)
     $                 +CSFMUL(THSQ,DUMTDQ4,FTDQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTDQ5,FTDQ,I,J)
     $                 -FTDQ(I,J)*(13.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEL(I,J)=FTELFEDFE(I,J)
     $                 +CMATMUL(0,FTEL,DUMTEL1,I,J)
     $                 +1.D0/2.D0*THSH*FTEL(I,J)*TDUMTEL2
     $                 +1.D0/4.D0*THSH*FTEL(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +CMATMUL(0,DUMTEL3,FE,I,J)
     $                 -2.D0*THSB*CMATMUL(0,GTPLS,SEFTEEGTPET,I,J)
     $                 +1.D0/2.D0*CSFMUL(THSL,DUMTEL4,FTEL,I,J)
     $                 +THSH*CSFMUL(THSL,DUMTEL5,FTEL,I,J)
     $                 -FTEL(I,J)*(9.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BFTUU(I,J)=3.D0/2.D0*THSH*FTUU(I,J)*TDUMTUU1
     $                 +1.D0/4.D0*THSH*FTUU(I,J)
     $                  *(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTUU2,FTUU,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTUU3,FTUU,I,J)
     $                 -4.D0/9.D0*THSB*CMATMUL(0,SQGTPQTFTUQ,GTPUS,I,J)
     $                 -16.D0/3.D0*THGL*CMATMUL(0,SQGTSQTFTUQ,GTSUS,I,J)
     $                 -4.D0/3.D0*THSB*G(184)*FUGTPUS(I,J)
     $                 +2.D0*THSH*CMATMUL(0,SUFTUUFTUUD,FTUU,I,J)
     $                 +CSFMUL(THSU,FTUU,DUMTUU4,I,J)
     $                 -FTUU(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDD(I,J)=1.D0/2.D0*THSH*FTDD(I,J)*TDUMTDD1
     $                 +1.D0/4.D0*THSH*FTDD(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTDD2,FTDD,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTDD3,FTDD,I,J)
     $                 +2.D0/9.D0*THSB*CMATMUL(0,SQGTPQTFTDQ,GTPDS,I,J)
     $                 -16.D0/3.D0*THGL*CMATMUL(0,SQGTSQTFTDQ,GTSDS,I,J)
     $                 -2.D0/3.D0*THSB*G(185)*FDGTPDS(I,J)
     $                 +2.D0*THSH*CMATMUL(0,SDFTDDFTDDD,FTDD,I,J)
     $                 +CSFMUL(THSD,FTDD,DUMTDD4,I,J)
     $                 -FTDD(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEE(I,J)=1.D0/2.D0*THSH*FTEE(I,J)*TDUMTEE1
     $                 +1.D0/4.D0*THSH*FTEE(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,FEFED,FTEE,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,DUMTEE3,FTEE,I,J)
     $                 -2.D0*THSB*CMATMUL(0,SLGTPLTFTEL,GTPES,I,J)
     $                 -2.D0*THSB*G(185)*FEGTPES(I,J)
     $                 +2.D0*THSH*CMATMUL(0,SEFTEEFTEED,FTEE,I,J)
     $                 +2.D0*THSB*CMATMUL(0,SEFTEEGTPET,GTPES,I,J)
     $                 -FTEE(I,J)*(9.D0/10.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2)
          END IF
!
!Convert into form readable by RKSTP. Two loop running is based
!on MV notation RGEs so I must take the transpose
!
          F(232+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BFTUQ(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMU(J,I))
          F(241+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BFTDQ(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMD(J,I))
          F(250+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BFTEL(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YME(J,I))
          F(259+(I-1)*3+J)=(1.D0/16.D0/PI**2*BFTUU(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMU(J,I))*THSU(J)
          F(268+(I-1)*3+J)=(1.D0/16.D0/PI**2*BFTDD(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMD(J,I))*THSD(J)
          F(277+(I-1)*3+J)=(1.D0/16.D0/PI**2*BFTEE(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YME(J,I))*THSE(J)
        END DO
      END DO
!
!That's the yukawas done, now for the gauge couplings
!
      DO I=1,3
        DO J=1,3
          DUMTQ1(I,J)=3.D0*SQGTQDGTQ(I,J)+SLGTLDGTL(I,J)
          IF(THHH.EQ.0)THEN
            DUMTQ2(I,J)=LUSLUT(I,J)+LDSLDT(I,J)
          ELSE
            DUMTQ2(I,J)=FUSFUT(I,J)+FDSFDT(I,J)
          END IF
          DUMTQ3(I,J)=THSH*(SUFTUUSFTUUT(I,J)+SDFTDDSFTDDT(I,J))
     $                +3.D0/2.D0*THSW*SQGTQDGTQ(I,J)
     $                +1.D0/18.D0*THSB*SQGTPQDGTPQ(I,J)
     $                +8.D0/3.D0*THGL*SQGTSQDGTSQ(I,J)
          DUMTQ4(I,J)=3.D0*THSW*GTQGTQD(I,J)+1.D0/9.D0*THSB
     $                 *GTPQGTPQD(I,J)+16.D0/3.D0*THGL*GTSQGTSQD(I,J) 
          DUMTQ5(I,J)=FTUQSFTUQT(I,J)+FTDQSFTDQT(I,J)
!
          DUMTL2(I,J)=THSH*SEFTEESFTEET(I,J)
     $                +3.D0/2.D0*THSW*SLGTLDGTL(I,J)
     $                +1.D0/2.D0*THSB*SLGTPLDGTPL(I,J)
          DUMTL3(I,J)=3.D0*THSW*GTLGTLD(I,J)+THSB*GTPLGTPLD(I,J)
!
          DUMTPQ1(I,J)=1.D0/3.D0*SQGTPQDGTPQ(I,J)+SLGTPLDGTPL(I,J)
     $                 +8.D0/3.D0*SUGTPUGTPUD(I,J)
     $                 +2.D0/3.D0*SDGTPDGTPDD(I,J)
     $                 +2.D0*SEGTPEGTPED(I,J)
          DUMTPQ2(I,J)=DUMTQ2(I,J)
          DUMTPQ3(I,J)=DUMTQ3(I,J)
          DUMTPQ4(I,J)=DUMTQ4(I,J)
          DUMTPQ5(I,J)=DUMTQ5(I,J)
!
          DUMTPL2(I,J)=DUMTL2(I,J)
          DUMTPL3(I,J)=DUMTL3(I,J)
!
          DUMTPU1(I,J)=4.D0/9.D0*THSB*SUGTPUGTPUD(I,J)
     $                 +4.D0/3.D0*THGL*SUGTSUGTSUD(I,J)
     $                 +THSH*SQFTUQTFTUQS(I,J)
          DUMTPU3(I,J)=8.D0/9.D0*THSB*GTPUDGTPU(I,J)
     $                 +8.D0/3.D0*THGL*GTSUDGTSU(I,J)
!
          DUMTPD1(I,J)=1.D0/9.D0*THSB*SDGTPDGTPDD(I,J)
     $                 +4.D0/3.D0*THGL*SDGTSDGTSDD(I,J)
     $                 +THSH*SQFTDQTFTDQS(I,J)
          DUMTPD3(I,J)=2.D0/9.D0*THSB*GTPDDGTPD(I,J)
     $                 +8.D0/3.D0*THGL*GTSDDGTSD(I,J)
!
          DUMTPE1(I,J)=THSB*SEGTPEGTPED(I,J)+THSH*SLFTELTFTELS(I,J)
        END DO
      END DO
      TDUMTQ1=CTRACE(DUMTQ1)
      TDUMTL1=TDUMTQ1
      TDUMTPQ1=CTRACE(DUMTPQ1)
      TDUMTPL1=TDUMTPQ1
      TDUMTPU2=TDUMTPQ1
      TDUMTPD2=TDUMTPU2
      TDUMTPE2=TDUMTPD2
      DO I=1,3
        DO J=1,3
          IF(THHH.EQ.0)THEN
            BGTQ(I,J)=1.D0/2.D0*THSW*GTQ(I,J)*TDUMTQ1
     $                +1.D0/2.D0*THSH*THSW*GTQ(I,J)
     $                 *(MSGTHUSQ+MCGTHDSQ)
     $                +1.D0/2.D0*CMATMUL(0,GTQ,DUMTQ2,I,J)
     $                +1.D0/2.D0*CMATMUL(0,GTQ,DUMTQ3,I,J)
     $                -2.D0*THSH*(G(289)*FTUQSLUT(I,J)
     $                 +G(290)*FTDQSLDT(I,J))
     $                +1.D0/2.D0*CSFMUL(THSQ,DUMTQ4,GTQ,I,J)
     $                +THSH*CSFMUL(THSQ,DUMTQ5,GTQ,I,J)
     $                -GTQ(I,J)*(1.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTL(I,J)=1.D0/2.D0*THSW*GTL(I,J)*TDUMTL1
     $                +1.D0/2.D0*THSH*THSW*GTL(I,J)
     $                 *(MSGTHUSQ+MCGTHDSQ)
     $                +1.D0/2.D0*CMATMUL(0,GTL,LESLET,I,J)
     $                +1.D0/2.D0*CMATMUL(0,GTL,DUMTL2,I,J)
     $                -2.D0*THSH*G(290)*FTELSLET(I,J)
     $                +1.D0/2.D0*CSFMUL(THSL,DUMTL3,GTL,I,J)
     $                +THSH*CSFMUL(THSL,FTELSFTELT,GTL,I,J)
     $                -GTL(I,J)*(9.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2)
!
            BGTPQ(I,J)=1.D0/2.D0*THSB*GTPQ(I,J)*TDUMTPQ1
     $                 +1.D0/2.D0*THSH*THSB*GTPQ(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,GTPQ,DUMTPQ2,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTPQ,DUMTPQ3,I,J)
     $                 +4.D0*THSH*(-2.D0
     $                  *CMATMUL(0,FTUQS,SUGTPUFTUUT,I,J)
     $                  +CMATMUL(0,FTDQS,SDGTPDFTDDT,I,J))
     $                 +6.D0*THSH*(G(287)*FTUQSLUT(I,J)
     $                  -G(288)*FTDQSLDT(I,J))
     $                 +1.D0/2.D0*CSFMUL(THSQ,DUMTPQ4,GTPQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTPQ5,GTPQ,I,J)
     $                 -GTPQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTPL(I,J)=1.D0/2.D0*THSB*GTPL(I,J)*TDUMTPL1
     $                 +1.D0/2.D0*THSH*THSB*GTPL(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,GTPL,LESLET,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTPL,DUMTPL2,I,J)
     $                 -4.D0*THSH*CMATMUL(0,FTELS,SEGTPEFTEET,I,J)
     $                 +2.D0*THSH*G(288)*FTELSLET(I,J)
     $                 +1.D0/2.D0*CSFMUL(THSL,DUMTPL3,GTPL,I,J)
     $                 +THSH*CSFMUL(THSL,FTELSFTELT,GTPL,I,J)
     $                 -GTPL(I,J)*(9.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BGTPU(I,J)=LUTLUSGTPU(I,J)
     $                 +CMATMUL(0,DUMTPU1,GTPU,I,J)
     $                 +1.D0/2.D0*THSB*GTPU(I,J)*TDUMTPU2
     $                 +1.D0/2.D0*THSB*THSH*GTPU(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 -3.D0*THSB*THSH*LUTFTUUS(I,J)*G(287)
     $                 -THSH*CMATMUL(0,SQFTUQTGTPQ,FTUUS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SUGTPUFTUUT,FTUUS,I,J)
     $                 +CSFMUL(THSU,GTPU,DUMTPU3,I,J)
     $                 -GTPU(I,J)*(4.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPD(I,J)=LDTLDSGTPD(I,J)
     $                 +CMATMUL(0,DUMTPD1,GTPD,I,J)
     $                 +1.D0/2.D0*THSB*GTPD(I,J)*TDUMTPD2
     $                 +1.D0/2.D0*THSB*THSH*GTPD(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 -6.D0*THSB*THSH*LDTFTDDS(I,J)*G(288)
     $                 +2.D0*THSH*CMATMUL(0,SQFTDQTGTPQ,FTDDS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SDGTPDFTDDT,FTDDS,I,J)
     $                 +CSFMUL(THSD,GTPD,DUMTPD3,I,J)
     $                 -GTPD(I,J)*(1.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPE(I,J)=LETLESGTPE(I,J)
     $                 +CMATMUL(0,DUMTPE1,GTPE,I,J)
     $                 +1.D0/2.D0*THSB*GTPE(I,J)*TDUMTPE2
     $                 +1.D0/2.D0*THSB*THSH*GTPE(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 -2.D0*THSB*THSH*LETFTEES(I,J)*G(288)
     $                 -2.D0*THSH*CMATMUL(0,SLFTELTGTPL,FTEES,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SEGTPEFTEET,FTEES,I,J)
     $                 +2.D0*THSB*CMATMUL(0,SEGTPEGTPED,GTPE,I,J)
     $                 -GTPE(I,J)*9.D0/5.D0*G(1)**2
          ELSE 
            BGTQ(I,J)=1.D0/2.D0*THSW*GTQ(I,J)*TDUMTQ1
     $                +1.D0/2.D0*THSH*THSW*GTQ(I,J)
     $                 *(MGTHUSQ+MGTHDSQ)
     $                +1.D0/2.D0*CMATMUL(0,GTQ,DUMTQ2,I,J)
     $                +1.D0/2.D0*CMATMUL(0,GTQ,DUMTQ3,I,J)
     $                -2.D0*THSH*(G(204)*FTUQSFUT(I,J)
     $                 +G(205)*FTDQSFDT(I,J))
     $                +1.D0/2.D0*CSFMUL(THSQ,DUMTQ4,GTQ,I,J)
     $                +THSH*CSFMUL(THSQ,DUMTQ5,GTQ,I,J)
     $                -GTQ(I,J)*(1.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTL(I,J)=1.D0/2.D0*THSW*GTL(I,J)*TDUMTL1
     $                +1.D0/2.D0*THSH*THSW*GTL(I,J)
     $                 *(MGTHUSQ+MGTHDSQ)
     $                +1.D0/2.D0*CMATMUL(0,GTL,FESFET,I,J)
     $                +1.D0/2.D0*CMATMUL(0,GTL,DUMTL2,I,J)
     $                -2.D0*THSH*G(205)*FTELSFET(I,J)
     $                +1.D0/2.D0*CSFMUL(THSL,DUMTL3,GTL,I,J)
     $                +THSH*CSFMUL(THSL,FTELSFTELT,GTL,I,J)
     $                -GTL(I,J)*(9.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2)
!
            BGTPQ(I,J)=1.D0/2.D0*THSB*GTPQ(I,J)*TDUMTPQ1
     $                 +1.D0/2.D0*THSH*THSB*GTPQ(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,GTPQ,DUMTPQ2,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTPQ,DUMTPQ3,I,J)
     $                 +4.D0*THSH*(-2.D0
     $                  *CMATMUL(0,FTUQS,SUGTPUFTUUT,I,J)
     $                  +CMATMUL(0,FTDQS,SDGTPDFTDDT,I,J))
     $                 +6.D0*THSH*(G(184)*FTUQSFUT(I,J)
     $                  -G(185)*FTDQSFDT(I,J))
     $                 +1.D0/2.D0*CSFMUL(THSQ,DUMTPQ4,GTPQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTPQ5,GTPQ,I,J)
     $                 -GTPQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTPL(I,J)=1.D0/2.D0*THSB*GTPL(I,J)*TDUMTPL1
     $                 +1.D0/2.D0*THSH*THSB*GTPL(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 +1.D0/2.D0*CMATMUL(0,GTPL,FESFET,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTPL,DUMTPL2,I,J)
     $                 -4.D0*THSH*CMATMUL(0,FTELS,SEGTPEFTEET,I,J)
     $                 +2.D0*THSH*G(185)*FTELSFET(I,J)
     $                 +1.D0/2.D0*CSFMUL(THSL,DUMTPL3,GTPL,I,J)
     $                 +THSH*CSFMUL(THSL,FTELSFTELT,GTPL,I,J)
     $                 -GTPL(I,J)*(9.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BGTPU(I,J)=FUTFUSGTPU(I,J)
     $                 +CMATMUL(0,DUMTPU1,GTPU,I,J)
     $                 +1.D0/2.D0*THSB*GTPU(I,J)*TDUMTPU2
     $                 +1.D0/2.D0*THSB*THSH*GTPU(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 -3.D0*THSB*THSH*FUTFTUUS(I,J)*G(184)
     $                 -THSH*CMATMUL(0,SQFTUQTGTPQ,FTUUS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SUGTPUFTUUT,FTUUS,I,J)
     $                 +CSFMUL(THSU,GTPU,DUMTPU3,I,J)
     $                 -GTPU(I,J)*(4.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPD(I,J)=FDTFDSGTPD(I,J)
     $                 +CMATMUL(0,DUMTPD1,GTPD,I,J)
     $                 +1.D0/2.D0*THSB*GTPD(I,J)*TDUMTPD2
     $                 +1.D0/2.D0*THSB*THSH*GTPD(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 -6.D0*THSB*THSH*FDTFTDDS(I,J)*G(185)
     $                 +2.D0*THSH*CMATMUL(0,SQFTDQTGTPQ,FTDDS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SDGTPDFTDDT,FTDDS,I,J)
     $                 +CSFMUL(THSD,GTPD,DUMTPD3,I,J)
     $                 -GTPD(I,J)*(1.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPE(I,J)=FETFESGTPE(I,J)
     $                 +CMATMUL(0,DUMTPE1,GTPE,I,J)
     $                 +1.D0/2.D0*THSB*GTPE(I,J)*TDUMTPE2
     $                 +1.D0/2.D0*THSB*THSH*GTPE(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 -2.D0*THSB*THSH*FETFTEES(I,J)*G(185)
     $                 -2.D0*THSH*CMATMUL(0,SLFTELTGTPL,FTEES,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SEGTPEFTEET,FTEES,I,J)
     $                 +2.D0*THSB*CMATMUL(0,SEGTPEGTPED,GTPE,I,J)
     $                 -GTPE(I,J)*9.D0/5.D0*G(1)**2
          END IF
!
!Convert into form readable by RKSTP
!
          F(185+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BGTQ(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(2)*ID(I,J))
          F(194+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BGTL(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(2)*ID(I,J))
          F(138+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BGTPQ(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
          F(147+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BGTPL(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
          F(156+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTPU(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
     $       *THSU(J)
          F(165+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTPD(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
     $       *THSD(J)
          F(174+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTPE(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
     $       *THSE(J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          DUMTSQ1(I,J)=2.D0*SQGTSQDGTSQ(I,J)+SUGTSUGTSUD(I,J)
     $                 +SDGTSDGTSDD(I,J)
          IF(THHH.EQ.0)THEN
            DUMTSQ2(I,J)=LUSLUT(I,J)+LDSLDT(I,J)
          ELSE
            DUMTSQ2(I,J)=FUSFUT(I,J)+FDSFDT(I,J)
          END IF
          DUMTSQ3(I,J)=DUMTQ3(I,J)
          DUMTSQ4(I,J)=DUMTQ4(I,J)
          DUMTSQ5(I,J)=DUMTQ5(I,J)
!
          DUMTSU1(I,J)=DUMTPU1(I,J)
          DUMTSU3(I,J)=SUGTSUGTSUD(I,J)+SDGTSDGTSDD(I,J)
          DUMTSU4(I,J)=DUMTPU3(I,J)
!
          DUMTSD1(I,J)=DUMTPD1(I,J)
          DUMTSD4(I,J)=DUMTPD3(I,J)
        END DO
      END DO
      TDUMTSQ1=CTRACE(DUMTSQ1)
      TDUMTSU2=CTRACE(SQGTSQDGTSQ)
      TDUMTSU3=CTRACE(DUMTSU3)
      TDUMTSD2=TDUMTSU2
      TDUMTSD3=TDUMTSU3
      DO I=1,3
        DO J=1,3
          IF(THHH.EQ.0)THEN
            BGTSQ(I,J)=1.D0/2.D0*THGL*TDUMTSQ1*GTSQ(I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTSQ,DUMTSQ2,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTSQ,DUMTSQ3,I,J)
     $                 -2.D0*THSH*(CMATMUL(0,FTUQS,SUGTSUFTUUT,I,J)
     $                  +CMATMUL(0,FTDQS,SDGTSDFTDDT,I,J))
     $                 +1.D0/2.D0*CSFMUL(THSQ,DUMTSQ4,GTSQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTSQ5,GTSQ,I,J)
     $                 -GTSQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+13.D0*G(3)**2)
!
            BGTSU(I,J)=LUTLUSGTSU(I,J)
     $                 +CMATMUL(0,DUMTSU1,GTSU,I,J)
     $                 +THGL*GTSU(I,J)*TDUMTSU2
     $                 +1.D0/2.D0*THGL*GTSU(I,J)*TDUMTSU3
     $                 -4.D0*THSH*CMATMUL(0,SQFTUQTGTSQ,FTUUS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SUGTSUFTUUT,FTUUS,I,J)
     $                 +CSFMUL(THSU,GTSU,DUMTSU4,I,J)
     $                 -GTSU(I,J)*(4.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)
!
            BGTSD(I,J)=LDTLDSGTSD(I,J)
     $                 +CMATMUL(0,DUMTSD1,GTSD,I,J)
     $                 +THGL*GTSD(I,J)*TDUMTSD2
     $                 +1.D0/2.D0*THGL*GTSD(I,J)*TDUMTSD3
     $                 -4.D0*THSH*CMATMUL(0,SQFTDQTGTSQ,FTDDS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SDGTSDFTDDT,FTDDS,I,J)
     $                 +CSFMUL(THSD,GTSD,DUMTSD4,I,J)
     $                 -GTSD(I,J)*(1.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)
          ELSE
            BGTSQ(I,J)=1.D0/2.D0*THGL*TDUMTSQ1*GTSQ(I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTSQ,DUMTSQ2,I,J)
     $                 +1.D0/2.D0*CMATMUL(0,GTSQ,DUMTSQ3,I,J)
     $                 -2.D0*THSH*(CMATMUL(0,FTUQS,SUGTSUFTUUT,I,J)
     $                  +CMATMUL(0,FTDQS,SDGTSDFTDDT,I,J))
     $                 +1.D0/2.D0*CSFMUL(THSQ,DUMTSQ4,GTSQ,I,J)
     $                 +THSH*CSFMUL(THSQ,DUMTSQ5,GTSQ,I,J)
     $                 -GTSQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+13.D0*G(3)**2)
!
            BGTSU(I,J)=FUTFUSGTSU(I,J)
     $                 +CMATMUL(0,DUMTSU1,GTSU,I,J)
     $                 +THGL*GTSU(I,J)*TDUMTSU2
     $                 +1.D0/2.D0*THGL*GTSU(I,J)*TDUMTSU3
     $                 -4.D0*THSH*CMATMUL(0,SQFTUQTGTSQ,FTUUS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SUGTSUFTUUT,FTUUS,I,J)
     $                 +CSFMUL(THSU,GTSU,DUMTSU4,I,J)
     $                 -GTSU(I,J)*(4.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)
!
            BGTSD(I,J)=FDTFDSGTSD(I,J)
     $                 +CMATMUL(0,DUMTSD1,GTSD,I,J)
     $                 +THGL*GTSD(I,J)*TDUMTSD2
     $                 +1.D0/2.D0*THGL*GTSD(I,J)*TDUMTSD3
     $                 -4.D0*THSH*CMATMUL(0,SQFTDQTGTSQ,FTDDS,I,J)
     $                 +2.D0*THSH*CMATMUL(0,SDGTSDFTDDT,FTDDS,I,J)
     $                 +CSFMUL(THSD,GTSD,DUMTSD4,I,J)
     $                 -GTSD(I,J)*(1.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)

          END IF
!
!Convert into form readable by RKSTP
!
          F(205+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BGTSQ(I,J)
     $           +1.D0/(16.D0*PI**2)**2*B2GM(3)*ID(I,J))
          F(214+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTSU(I,J)
     $           +1.D0/(16.D0*PI**2)**2*B2GM(3)*ID(I,J))*THSU(J)
          F(223+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTSD(I,J)
     $           +1.D0/(16.D0*PI**2)**2*B2GM(3)*ID(I,J))*THSD(J)
        END DO
      END DO
!
!The Higgs tilde gauge couplings are calculated twice,
!once for below m_H when s and c are absorbed into the definition.
!Some terms are valid in both regimes and do not need to be
!calculated twice.
!
      DO I=1,3
        DO J=1,3
          DUMGTHU1(I,J)=CMATMUL(0,SQFTUQTGTQ,FUS,I,J)
          DUMGTHD1(I,J)=3.D0*CMATMUL(0,SQFTDQTGTQ,FDS,I,J)
     $                  +CMATMUL(0,SLFTELTGTL,FES,I,J)
          DUMGTPHU1(I,J)=CMATMUL(0,SQFTUQTGTPQ,FUS,I,J)
          DUMGTPHU2(I,J)=CMATMUL(0,SUGTPUFTUUT,FUS,I,J)
          DUMGTPHD1(I,J)=-CMATMUL(0,SQFTDQTGTPQ,FDS,I,J)
     $                   +CMATMUL(0,SLFTELTGTPL,FES,I,J)
          DUMGTPHD2(I,J)=CMATMUL(0,SDGTPDFTDDT,FDS,I,J)
     $                   +CMATMUL(0,SEGTPEFTEET,FES,I,J)
        END DO
      END DO
      TSQGTQDGTQ=CTRACE(SQGTQDGTQ)
      TSLGTLDGTL=CTRACE(SLGTLDGTL)
      TSUFTUUFTUUD=CTRACE(SUFTUUFTUUD)
      TSQFTUQDFTUQ=CTRACE(SQFTUQDFTUQ)
      TDUMGTHU1=CTRACE(DUMGTHU1)
      TFUSFUT=CTRACE(FUSFUT)
      TFDSFDT=CTRACE(FDSFDT)
      TFESFET=CTRACE(FESFET)
      TSDFTDDFTDDD=CTRACE(SDFTDDFTDDD)
      TSEFTEEFTEED=CTRACE(SEFTEEFTEED)
      TSQFTDQDFTDQ=CTRACE(SQFTDQDFTDQ)
      TSLFTELDFTEL=CTRACE(SLFTELDFTEL)
      TDUMGTHD1=CTRACE(DUMGTHD1)
      TSQGTPQDGTPQ=CTRACE(SQGTPQDGTPQ)
      TSLGTPLDGTPL=CTRACE(SLGTPLDGTPL)
      TSUGTPUGTPUD=CTRACE(SUGTPUGTPUD)
      TSDGTPDGTPDD=CTRACE(SDGTPDGTPDD)
      TSEGTPEGTPED=CTRACE(SEGTPEGTPED)
      TDUMGTPHU1=CTRACE(DUMGTPHU1)
      TDUMGTPHU2=CTRACE(DUMGTPHU2)
      TDUMGTPHD1=CTRACE(DUMGTPHD1)
      TDUMGTPHD2=CTRACE(DUMGTPHD2)
!
      BGTHU=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(204)
     $      +1.D0/2.D0*THSH*THSW*(MGTHUSQ+MGTHDSQ)*G(204)
     $      +1.D0/2.D0*THSH*(3.D0*G(204)*TSUFTUUFTUUD
     $       +3.D0*G(204)*TSQFTUQDFTUQ)
     $      +1.D0/4.D0*THSH*G(204)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $      -6.D0*TDUMGTHU1
     $      +G(204)*3.D0*TFUSFUT
     $      +1.D0/2.D0*THSH*G(204)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $      -G(204)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
      BGTHD=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(205)
     $      +1.D0/2.D0*THSH*THSW*(MGTHUSQ+MGTHDSQ)*G(205)
     $      +1.D0/2.D0*THSH*(3.D0*G(205)*TSDFTDDFTDDD
     $       +G(205)*TSEFTEEFTEED+3.D0*G(205)*TSQFTDQDFTDQ
     $       +G(205)*TSLFTELDFTEL)
     $      +1.D0/4.D0*THSH*G(205)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $      -2.D0*TDUMGTHD1
     $      +G(205)*(3.D0*TFDSFDT+TFESFET)
     $      +1.D0/2.D0*THSH*G(205)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $      -G(205)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
      BGTPHU=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $        +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $        +2.D0*TSEGTPEGTPED)*G(184)
     $       +1.D0/2.D0*THSH*THSB*(MGTPHUSQ+MGTPHDSQ)*G(184)
     $       +1.D0/2.D0*THSH*(3.D0*G(184)*TSUFTUUFTUUD
     $        +3.D0*G(184)*TSQFTUQDFTUQ)
     $       +1.D0/4.D0*THSH*G(184)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $       +2.D0*TDUMGTPHU1-8.D0*TDUMGTPHU2
     $       +G(184)*3.D0*TFUSFUT
     $       +1.D0/2.D0*THSH*G(184)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $       -G(184)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
      BGTPHD=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $        +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $        +2.D0*TSEGTPEGTPED)*G(185)
     $       +1.D0/2.D0*THSH*THSB*(MGTPHUSQ+MGTPHDSQ)*G(185)
     $       +1.D0/2.D0*THSH*(3.D0*G(185)*TSDFTDDFTDDD
     $        +G(185)*TSEFTEEFTEED+3.D0*G(185)*TSQFTDQDFTDQ
     $        +G(185)*TSLFTELDFTEL)
     $       +1.D0/4.D0*THSH*G(185)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       +2.D0*TDUMGTPHD1-4.D0*TDUMGTPHD2
     $       +G(185)*(3.D0*TFDSFDT+TFESFET)
     $       +1.D0/2.D0*THSH*G(185)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       -G(185)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
!Convert into form readable by RKSTP
!
      F(184)=THHH*(1.D0/16.D0/PI**2*BGTPHU
     $           +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0))
      F(185)=THHH*(1.D0/16.D0/PI**2*BGTPHD
     $           +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0))
      F(204)=THHH*(1.D0/16.D0/PI**2*BGTHU
     $           +1.D0/(16.D0*PI**2)**2*B2GM(2))
      F(205)=THHH*(1.D0/16.D0/PI**2*BGTHD
     $           +1.D0/(16.D0*PI**2)**2*B2GM(2))
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            DUMGTHU1(I,J)=CMATMUL(0,SQFTUQTGTQ,LUS,I,J)
            DUMGTHD1(I,J)=3.D0*CMATMUL(0,SQFTDQTGTQ,LDS,I,J)
     $                    +CMATMUL(0,SLFTELTGTL,LES,I,J)
            DUMGTPHU1(I,J)=CMATMUL(0,SQFTUQTGTPQ,LUS,I,J)
            DUMGTPHU2(I,J)=CMATMUL(0,SUGTPUFTUUT,LUS,I,J)
            DUMGTPHD1(I,J)=-CMATMUL(0,SQFTDQTGTPQ,LDS,I,J)
     $                     +CMATMUL(0,SLFTELTGTPL,LES,I,J)
            DUMGTPHD2(I,J)=CMATMUL(0,SDGTPDFTDDT,LDS,I,J)
     $                     +CMATMUL(0,SEGTPEFTEET,LES,I,J)
          END DO
        END DO
        TDUMGTHU1=CTRACE(DUMGTHU1)
        TLUSLUT=CTRACE(LUSLUT)
        TLDSLDT=CTRACE(LDSLDT)
        TLESLET=CTRACE(LESLET)
        TDUMGTHD1=CTRACE(DUMGTHD1)
        TDUMGTPHU1=CTRACE(DUMGTPHU1)
        TDUMGTPHU2=CTRACE(DUMGTPHU2)
        TDUMGTPHD1=CTRACE(DUMGTPHD1)
        TDUMGTPHD2=CTRACE(DUMGTPHD2)
!
        BSGTHU=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(289)
     $        +1.D0/2.D0*THSH*THSW*(MSGTHUSQ+MCGTHDSQ)*G(289)
     $        +1.D0/2.D0*THSH*(3.D0*G(289)*TSUFTUUFTUUD
     $         +3.D0*G(289)*TSQFTUQDFTUQ)
     $        +1.D0/4.D0*THSH*G(289)
     $         *(3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $        -6.D0*TDUMGTHU1
****
     $        +THSB*THSH*G(290)*CONJG(G(288))*G(287)
     $        -THSW*THSH*MCGTHDSQ*G(289)
****
     $        +G(289)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $        +1.D0/2.D0*THSH*G(289)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $         +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $        -G(289)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
        BCGTHD=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(290)
     $        +1.D0/2.D0*THSH*THSW*(MSGTHUSQ+MCGTHDSQ)*G(290)
     $        +1.D0/2.D0*THSH*(3.D0*G(290)*TSDFTDDFTDDD
     $         +G(290)*TSEFTEEFTEED+3.D0*G(290)*TSQFTDQDFTDQ
     $         +G(290)*TSLFTELDFTEL)
     $        +1.D0/4.D0*THSH*G(290)*(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $        -2.D0*TDUMGTHD1
****
     $        +THSB*THSH*G(289)*CONJG(G(287))*G(288)
     $        -THSW*THSH*MSGTHUSQ*G(290)
****
     $        +G(290)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $        +1.D0/2.D0*THSH*G(290)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $         +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $        -G(290)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
        BSGTPHU=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $          +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $          +2.D0*TSEGTPEGTPED)*G(287)
     $         +1.D0/2.D0*THSH*THSB*(MSGTPHUSQ+MCGTPHDSQ)*G(287)
     $         +1.D0/2.D0*THSH*(3.D0*G(287)*TSUFTUUFTUUD
     $          +3.D0*G(287)*TSQFTUQDFTUQ)
     $         +1.D0/4.D0*THSH*G(287)*(3.D0*THSW*MSGTHUSQ
     $          +THSB*MSGTPHUSQ)
     $         +2.D0*TDUMGTPHU1-8.D0*TDUMGTPHU2
****
     $         +THSB*THSH*MCGTPHDSQ*G(287)
     $         +3.D0*THSW*THSH*G(288)*CONJG(G(290))*G(289)
****
     $         +G(287)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $         +1.D0/2.D0*THSH*G(287)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $          +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $         -G(287)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
        BCGTPHD=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $          +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $          +2.D0*TSEGTPEGTPED)*G(288)
     $         +1.D0/2.D0*THSH*THSB*(MSGTPHUSQ+MCGTPHDSQ)*G(288)
     $         +1.D0/2.D0*THSH*(3.D0*G(288)*TSDFTDDFTDDD
     $          +G(288)*TSEFTEEFTEED+3.D0*G(288)*TSQFTDQDFTDQ
     $          +G(288)*TSLFTELDFTEL)
     $         +1.D0/4.D0*THSH*G(288)*(3.D0*THSW*MCGTHDSQ
     $          +THSB*MCGTPHDSQ)
     $         +2.D0*TDUMGTPHD1-4.D0*TDUMGTPHD2
****
     $         +THSB*THSH*MSGTPHUSQ*G(288)
     $         +3.D0*THSW*THSH*G(287)*CONJG(G(289))*G(290)
****
     $         +G(288)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $         +1.D0/2.D0*THSH*G(288)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $          +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $         -G(288)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
!The two loop part of the running will be the same as the
!ordinary g running multiplied by an s (or corresponding c)
!
!
!Convert into form readable by RKSTP
!
        F(287)=1.D0/16.D0/PI**2*BSGTPHU
     $             +1.D0/(16.D0*PI**2)**2*SINB*B2GM(1)*DSQRT(3.D0/5.D0)
        F(288)=1.D0/16.D0/PI**2*BCGTPHD
     $             +1.D0/(16.D0*PI**2)**2*COSB*B2GM(1)*DSQRT(3.D0/5.D0)
        F(289)=1.D0/16.D0/PI**2*BSGTHU
     $             +1.D0/(16.D0*PI**2)**2*SINB*B2GM(2)
        F(290)=1.D0/16.D0/PI**2*BCGTHD
     $             +1.D0/(16.D0*PI**2)**2*COSB*B2GM(2)
      END IF
!
!Now I work out the \tilde{\mu}^* f^{h_u}_u terms. These only exist
!above m_H
!
      IF(THHH.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            SQIDMTSFU(I,J)=CSFMUL(THSQ,ID,MTSFU,I,J)
            SQIDMTSFD(I,J)=CSFMUL(THSQ,ID,MTSFD,I,J)
            SQFUHUDMTSFU(I,J)=CSFMUL(THSQ,FUHUD,MTSFU,I,J)
            SQFDHDDMTSFD(I,J)=CSFMUL(THSQ,FDHDD,MTSFD,I,J)
            SLIDMTSFE(I,J)=CSFMUL(THSL,ID,MTSFE,I,J)
            SLFEHDDMTSFE(I,J)=CSFMUL(THSL,FEHDD,MTSFE,I,J)
            SUMTSFUID(I,J)=CSFMUL(THSU,MTSFU,ID,I,J)
            SUMTSFUFUQDFDQ(I,J)=CSFMUL(THSU,MTSFU,FUQDFDQ,I,J)
            SDMTSFDID(I,J)=CSFMUL(THSD,MTSFD,ID,I,J)
            SDMTSFDFDQDFUQ(I,J)=CSFMUL(THSD,MTSFD,FDQDFUQ,I,J)
            SEMTSFEID(I,J)=CSFMUL(THSE,MTSFE,ID,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SUSQIDMTSFUID(I,J)=CSFMUL(THSU,SQIDMTSFU,ID,I,J)
            SDSQIDMTSFDID(I,J)=CSFMUL(THSD,SQIDMTSFD,ID,I,J)
            SESLIDMTSFEID(I,J)=CSFMUL(THSE,SLIDMTSFE,ID,I,J)
          END DO
        END DO
        SUTSQFUHUDMTSFU=TCSFMUL(THSU,SQFUHUDMTSFU)
        SDTSQFDHDDMTSFD=TCSFMUL(THSD,SQFDHDDMTSFD)
        SETSLFEHDDMTSFE=TCSFMUL(THSE,SLFEHDDMTSFE)
        DO I=1,3
          DO J=1,3
            DUMMTSFU1(I,J)=(G(457)**2/6.D0-3.D0*G(458)**2/2.D0)*ID(I,J)
     $                    +2.D0*FDDRFDDRD(I,J)
            DUMMTSFU2(I,J)=8.D0/9.D0*THSB*GTPUTGTPUS(I,J)
     $                    +8.D0/3.D0*THGL*GTSUTGTSUS(I,J)
     $                    +2.D0*THSH*FTUUDFTUU(I,J)
            DUMMTSFU3(I,J)=THSH*FTUQFTUQD(I,J)+THSH*FTDQFTDQD(I,J)
     $                    +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $                    +3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $                    +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
!
            DUMMTSFD1(I,J)=(G(457)**2/6.D0+3.D0*G(458)**2/2.D0)*ID(I,J)
     $                    -2.D0*FUURFUURD(I,J)
            DUMMTSFD2(I,J)=2.D0/9.D0*THSB*GTPDTGTPDS(I,J)
     $                    +8.D0/3.D0*THGL*GTSDTGTSDS(I,J)
     $                    +2.D0*THSH*FTDDDFTDD(I,J)
            DUMMTSFD3(I,J)=DUMMTSFU3(I,J)
!
            DUMMTSFE1(I,J)=(G(457)**2/2.D0-3.D0*G(458)**2/2.D0)*ID(I,J)
            DUMMTSFE2(I,J)=2.D0*THSB*GTPETGTPES(I,J)
     $                    +2.D0*THSH*FTEEDFTEE(I,J)
            DUMMTSFE3(I,J)=THSH*FTELFTELD(I,J)
     $                    +1.D0/2.D0*THSB*GTPLSGTPLT(I,J)
     $                    +3.D0/2.D0*THSW*GTLSGTLT(I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            BMTSFU(I,J)=2.D0/3.D0*G(457)**2*SUMTSFUID(I,J)
     $                 -2.D0*(G(457)**2/9.D0+4.D0*G(459)**2/3.D0)
     $                  *SUSQIDMTSFUID(I,J)
     $                 +6.D0*FUHU(I,J)*SUTSQFUHUDMTSFU
     $                 -CSFMUL(THSQ,DUMMTSFU1,MTSFU,I,J)
     $                 -2.D0*SDMTSFDFDQDFUQ(I,J)
     $                 -2.D0/3.D0*THSB*THSH*CONJG(G(108))*G(185)
     $                  *(4.D0*FTUQGTPUS(I,J)-GTPQSFTUU(I,J))
     $                 -6.D0*THSH*THSW*CONJG(G(108))*G(205)
     $                  *GTQSFTUU(I,J)
     $                 +4.D0*THSH*CONJG(G(108))*FTDQFDDFTUU(I,J)
     $                 +CSFMUL(THSU,MTSFU,DUMMTSFU2,I,J)
     $                 +(3.D0*TFDDFD+TFEDFE+1.D0/2.D0*THSB*THSH*MGTPHDSQ
     $                  +3.D0/2.D0*THSW*THSH*MGTHDSQ)*MTSFU(I,J)
     $                 +CSFMUL(THSQ,DUMMTSFU3,MTSFU,I,J)
     $                 -(13.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*MTSFU(I,J)
!
            BMTSFD(I,J)=1.D0/3.D0*G(457)**2*SDMTSFDID(I,J)
     $                 +2.D0*(G(457)**2/18.D0-4.D0*G(459)**2/3.D0)
     $                  *SDSQIDMTSFDID(I,J)
     $                 +6.D0*FDHD(I,J)*SDTSQFDHDDMTSFD
     $                 +2.D0*FDHD(I,J)*SETSLFEHDDMTSFE
     $                 +CSFMUL(THSQ,DUMMTSFD1,MTSFD,I,J)
     $                 -2.D0*SUMTSFUFUQDFDQ(I,J)
     $                 -2.D0/3.D0*THSB*THSH*CONJG(G(108))*G(184)
     $                  *(2.D0*FTDQGTPDS(I,J)+GTPQSFTDD(I,J))
     $                 -6.D0*THSH*THSW*CONJG(G(108))*G(204)
     $                  *GTQSFTDD(I,J)
     $                 +4.D0*THSH*CONJG(G(108))*FTUQFUDFTDD(I,J)
     $                 +CSFMUL(THSD,MTSFD,DUMMTSFD2,I,J)
     $                 +(3.D0*TFUDFU+1.D0/2.D0*THSB*THSH*MGTPHUSQ
     $                  +3.D0/2.D0*THSW*THSH*MGTHUSQ)*MTSFD(I,J)
     $                 +CSFMUL(THSQ,DUMMTSFD3,MTSFD,I,J)
     $                 -(7.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*MTSFD(I,J)
!
            BMTSFE(I,J)=G(457)**2*SEMTSFEID(I,J)
     $                 -G(457)**2*SESLIDMTSFEID(I,J)
     $                 +2.D0*FEHD(I,J)*SETSLFEHDDMTSFE
     $                 +6.D0*FEHD(I,J)*SDTSQFDHDDMTSFD
     $                 -CSFMUL(THSL,DUMMTSFE1,MTSFE,I,J)
     $                 -2.D0*THSB*THSH*CONJG(G(108))*G(184)
     $                  *(2.D0*FTELGTPES(I,J)-GTPLSFTEE(I,J))
     $                 -6.D0*THSH*THSW*CONJG(G(108))*G(204)
     $                  *GTLSFTEE(I,J)
     $                 +CSFMUL(THSE,MTSFE,DUMMTSFE2,I,J)
     $                 +(3.D0*TFUDFU+1.D0/2.D0*THSB*THSH*MGTPHUSQ
     $                  +3.D0/2.D0*THSW*THSH*MGTHUSQ)*MTSFE(I,J)
     $                 +CSFMUL(THSL,DUMMTSFE3,MTSFE,I,J)
     $                 -(27.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $                  *MTSFE(I,J)
!
!Convert into form readable by RKSTP
!
            F(429+(I-1)*3+J)=1.D0/16.D0/PI**2*BMTSFU(I,J)
     $           +1.D0/(16.D0*PI**2)**2*(CONJG(G(398))*B2YMU(J,I)
     $                                    +YU(J,I)*CONJG(B2GRKMUM))
            F(438+(I-1)*3+J)=1.D0/16.D0/PI**2*BMTSFD(I,J)
     $           +1.D0/(16.D0*PI**2)**2*(CONJG(G(398))*B2YMD(J,I)
     $                                    +YD(J,I)*CONJG(B2GRKMUM))
            F(447+(I-1)*3+J)=1.D0/16.D0/PI**2*BMTSFE(I,J)
     $           +1.D0/(16.D0*PI**2)**2*(CONJG(G(398))*B2YME(J,I)
     $                                    +YE(J,I)*CONJG(B2GRKMUM))

          END DO
        END DO
      END IF
!
!
!Next I am going to work out the primed gaugino terms, M'_{1,2,3}
!
!
      IF(THHH.EQ.0)THEN
        BMP(1)=G(599)*(TDUMM1+THSH*(MSGTPHUSQ+MCGTPHDSQ))
     $        -2.D0*(0.D0,1.D0)*THSH*(G(287)*CONJG(G(108))*G(288)
     $                            -CONJG(G(287))*G(108)*CONJG(G(288)))
!
        BMP(2)=G(600)*(TDUMM2+THSH*(MSGTHUSQ+MCGTHDSQ))
     $        -2.D0*(0.D0,1.D0)*THSH*(G(289)*CONJG(G(108))*G(290)
     $                            -CONJG(G(289))*G(108)*CONJG(G(290)))
     $        -12.D0*G(600)*G(2)**2
!
        BMP(3)=G(601)*TDUMM3-18.D0*G(601)*G(3)**2
      ELSE
        BMP(1)=G(599)*(TDUMM1+THSH*(MGTPHUSQ+MGTPHDSQ))
!
        BMP(2)=G(600)*(TDUMM2+THSH*(MGTHUSQ+MGTHDSQ))
     $        -12.D0*G(600)*G(2)**2
!
        BMP(3)=G(601)*TDUMM3-18.D0*G(601)*G(3)**2
      END IF
!
!NB: There is some work here converting the MSSM MV notation
!    M_{1,2,3} 2-lp beta functions to BT notation.
!The RKSTP compatible derivatives are:
!
      F(599)=THSB*(1.D0/16.D0/PI**2*BMP(1)
     $         +1.D0/(16.D0*PI**2)**2*(0.D0,1.D0)/2.D0
     $                                *(B2M(1)-CONJG(B2M(1))))
      F(600)=THSW*(1.D0/16.D0/PI**2*BMP(2)
     $         +1.D0/(16.D0*PI**2)**2*(0.D0,1.D0)/2.D0
     $                                *(B2M(2)-CONJG(B2M(2))))
      F(601)=THGL*(1.D0/16.D0/PI**2*BMP(3)
     $         +1.D0/(16.D0*PI**2)**2*(0.D0,1.D0)/2.D0
     $                                *(B2M(3)-CONJG(B2M(3))))
!
!Rotate F to original basis.
!
      DO I=1,601
        FCURR(I)=F(I)
      END DO
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15)THEN
!
!If necessary, rotate the Fs back into the old current basis.
!
        CALL ROTSQBACK(F,FCURR)
      END IF
!
!Set the OLDN for the next step
!
      OLDNSQ=NSQ
      OLDNSU=NSU
      OLDNSD=NSD
      OLDNSL=NSL
      OLDNSE=NSE
!
      RETURN
      END
!
      SUBROUTINE CRKSTP(N,H,X,Y,SUB,W)
C-----------------------------------------------------------------------
C     A double complex copy of the subroutine rkstp 
C     From CERN Program Library, routine D209, with error message for
C     N.LT.1 replaced by STOP 99 to eliminate Kernlib error routine.
C
!      IMPLICIT NONE
!
!      dimension y(n),w(n,3)
      DOUBLE COMPLEX Y(N),W(N,3)
      DOUBLE PRECISION H,HLOCAL,H2,X,H6,XH,XH2
      INTEGER N,NLOCAL,J
      LOGICAL MFLAG,RFLAG
      EXTERNAL SUB
C
C     ******************************************************************
C
C     THIS SUBROUTINE REPLACES X BY X+H AND ADVANCES THE SOLUTION OF THE
C     SYSTEM OF DIFFERENTIAL EQUATIONS DY/DX=F(X,Y) FROM Y(X) TO Y(X+H)
C     USING A FIFTH-ORDER RUNGE-KUTTA METHOD.
C
C     SUB IS THE NAME OF A SUBROUTINE SUB(X,Y,F) WHICH SETS THE VECTOR F
C     TO THE DERIVATIVE AT X OF THE VECTOR Y.
C
C     W IS A WORKING-SPACE ARRAY, TREATED AS CONSISTING OF THREE CONSEC-
C     UTIVE WORKING VECTORS OF LENGTH N.
C
C     ******************************************************************
C
C  START.
      IF (N.LT.1) STOP 99
      NLOCAL=N
      HLOCAL=H
      H2=0.5d0*HLOCAL
      H6=HLOCAL/6.d0
      XH=X+HLOCAL
      XH2=X+H2
      CALL SUB(X,Y,W(1,1))
      DO 1 J=1,NLOCAL
         W(J,2)=Y(J)+H2*W(J,1)
    1 CONTINUE
      CALL SUB(XH2,W(1,2),W(1,3))
      DO 2 J=1,NLOCAL
         W(J,1)=W(J,1)+2.d0*W(J,3)
         W(J,2)=Y(J)+H2*W(J,3)
    2 CONTINUE
      CALL SUB(XH2,W(1,2),W(1,3))
      DO 3 J=1,NLOCAL
         W(J,1)=W(J,1)+2.d0*W(J,3)
         W(J,2)=Y(J)+HLOCAL*W(J,3)
    3 CONTINUE
      CALL SUB(XH,W(1,2),W(1,3))
      DO 4 J=1,NLOCAL
         Y(J)=Y(J)+H6*(W(J,1)+W(J,3))
    4 CONTINUE
      X=XH
      RETURN
      END
!
      FUNCTION CSFMUL(TH,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B but taking
!         into account restrictions on the multiple from
!         thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A, B input as (3x3) complex
!         I, J integer
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),B(3,3),CSFMUL
      INTEGER LOOP,TH(3),I,J
!
      CSFMUL=(0.D0,0.D0)
      DO LOOP=1,3
        CSFMUL=CSFMUL+A(I,LOOP)*B(LOOP,J)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
!
      FUNCTION CTRACE(A)
!
!Purpose: To compute the trace of matrix A
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),CTRACE
      INTEGER I
!
      CTRACE=(0.D0,0.D0)
      DO I=1,3
        CTRACE=CTRACE+A(I,I)
      END DO
!
      RETURN
      END
!
      SUBROUTINE DAGGER(A,DAGA)
!
!Purpose: To compute the dagger of matrix A
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION DAGA(3,3)
      INTEGER I,J
!
      DO I=1,3
        DO J=1,3
          DAGA(I,J)=A(J,I)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE DECOUTCOUP(Q,G)
!
!Purpose: To output m^2_U in the mass basis as it runs down.
!
      IMPLICIT NONE
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE COMPLEX GTMP(601),G(601)
      DOUBLE PRECISION Q,COSB,SINB
      INTEGER I,SWABS,SWROT
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!These switches can be used to print all complex entries or
!turn off the rotation.
!
      SWABS=1
      SWROT=1
!
!Remember the old gs
!
      DO I=1,601
        GTMP(I)=G(I)
      END DO
!
!Perform the rotation to the mass basis if nec.
!
      IF(SWROT.EQ.1)CALL STROTBACK(GTMP,G,0)
!
!Make some adjustments to facilitate comparisons
!
      IF(SWABS.EQ.1)THEN
        DO I=1,601
          IF(ABS(G(I)).LT.1.D-30)THEN
            G(I)=(1.D-30,1.D-30)
          END IF
        END DO
      END IF
!
!Print out the terms - absolute squares first
!
      IF(SWABS.EQ.1)THEN
        I=81
        WRITE(63,52)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
!
!Now twice the size with real and imaginary parts
!
      ELSE
        I=81
        WRITE(77,54)Q,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5)
     $    ,G(I+6),G(I+7),G(I+8)
      END IF
!
!Return the gs to their original values
!
      DO I=1,601
        G(I)=GTMP(I)
      END DO
!
   52 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
   54 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10)
!
      RETURN
      END
!
      SUBROUTINE DECRUN(GIN,GOUT,QIN,SVLQ)
!
!Purpose: To run the inputs to from m_H to the required scale.
!
!         Uses knowledge of the decoupling from RGEFLAV such as
!         decoupling points and frozen values of couplings, masses
!         and squark mixings.
!
      IMPLICIT NONE
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
      DOUBLE COMPLEX GIN(601),GRUN(601),GOUT(601),GFIX(601),W(1803)
      DOUBLE PRECISION DGRUN(601),DW(1803),QIN,TQIN
      DOUBLE PRECISION DT,T,TMH,TTH(20),TMDEC,Q,MDEC,PI
      INTEGER I,J,II,NSTEP,STTH,FINTH,NSTEPTH,BELOWMS,SVLQ
      EXTERNAL CRGE601,DRGE601
!
!Reinstate this (and certain lines below) to write out the running
!of m^2_U according to DECRUN using DECOUTCOUP.
!
!      OPEN(63,FILE='out/stdown.dat',STATUS='UNKNOWN')
      IF(QIN.LT.RGEMS)THEN
        BELOWMS=1
      ELSE
        BELOWMS=0
      END IF
!
      MDEC=QTHQL(1)
      IF(QTHUR(1).LT.QTHQL(1).AND.SVLQ.EQ.1)MDEC=QTHUR(1)
      IF(QTHDR(1).LT.QTHQL(1).AND.SVLQ.EQ.0)MDEC=QTHDR(1)
!
!Prevent running from going below mt
!
      IF(MDEC.LT.MT)THEN
        WRITE(33,*)'ATTEMPTED TO RUN BELOW MT IN DECRUN'
        MDEC=MT
      END IF
!
      PI=4.D0*DATAN(1.D0)
!
      CALL STROTATE(GIN,GRUN,0)
!
!The number of steps is fixed and large, but could be allowed
!to float or be equal to the numbers in RGEFLAV.
!
      NSTEP=1000
      NSTEPTH=500
      TMH=0.D0
      TQIN=LOG(QIN/QNH)
      TMDEC=LOG(MDEC/QNH)
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/QNH)
      END DO
!
!Find start and end points. FINTH is threshold below MDEC
!
      STTH=LOCMH
      FINTH=20
      DO I=20,1,-1
        IF(MDEC.LT.QTHSORT(I))FINTH=FINTH-1
      END DO
!
      IF(MDEC.LT.QIN)THEN
!
!Running down past the relevant thresholds
!
        IF(STTH.NE.FINTH+1)THEN
          DO I=STTH,FINTH+2,-1
            IF(QTHSORT(I).EQ.QTHSORT(I-1))GOTO 41
            DT=(TTH(I-1)-TTH(I))/FLOAT(NSTEPTH)
            DO II=1,NSTEPTH
              T=TTH(I)+(TTH(I-1)-TTH(I))*FLOAT(II-1)/FLOAT(NSTEPTH)
              SSQSTEP=QNH*EXP(T)
!
!Finite Pierce correction
!
              IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
                BELOWMS=1
                CALL STROTBACK(GRUN,GFIX,1)
                GFIX(12)=GFIX(12)*(1.D0-DBLE(RTISA))
                GFIX(21)=GFIX(21)*(1.D0-DBLE(RBISA))
                GFIX(30)=GFIX(30)*(1.D0-DBLE(RLISA))
                IF(RGEMS.LT.QNH)THEN
                  GFIX(120)=GFIX(120)*(1.D0-DBLE(RTISA))
                  GFIX(129)=GFIX(129)*(1.D0-DBLE(RBISA))
                  GFIX(138)=GFIX(138)*(1.D0-DBLE(RLISA))
                END IF
                CALL STROTATE(GFIX,GRUN,1)
              END IF
!
              IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
              CALL SQDIAG(GRUN)
              CALL CHDEC(SSQSTEP,GRUN,1)
              IF(COMP.EQ.0)THEN
                DO J=1,601
                  DGRUN(J)=DBLE(GRUN(J))
                END DO
                CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
                DO J=1,601
                  GRUN(J)=DCMPLX(DGRUN(J))
                END DO
              ELSE
                CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
              END IF
              EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!              CALL DECOUTCOUP(QNH*EXP(T),GRUN)
            END DO
  41      END DO
        END IF
!
!Now the final run to MDEC
!
        DT=(TMDEC-TTH(FINTH+1))/FLOAT(NSTEP)
        DO II=1,NSTEP
          T=TTH(FINTH+1)+(TMDEC-TTH(FINTH+1))*FLOAT(II-1)/FLOAT(NSTEP)
          SSQSTEP=QNH*EXP(T)
          Q=QNH*EXP(T)
!
          IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
            BELOWMS=1
            CALL STROTBACK(GRUN,GFIX,1)
            GFIX(12)=GFIX(12)*(1.D0-DBLE(RTISA))
            GFIX(21)=GFIX(21)*(1.D0-DBLE(RBISA))
            GFIX(30)=GFIX(30)*(1.D0-DBLE(RLISA))
            IF(RGEMS.LT.QNH)THEN
              GFIX(120)=GFIX(120)*(1.D0-DBLE(RTISA))
              GFIX(129)=GFIX(129)*(1.D0-DBLE(RBISA))
              GFIX(138)=GFIX(138)*(1.D0-DBLE(RLISA))
            END IF
            CALL STROTATE(GFIX,GRUN,1)
          END IF
!
          IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          CALL SQDIAG(GRUN)
          CALL CHDEC(SSQSTEP,GRUN,1)
          IF(COMP.EQ.0)THEN
            DO J=1,601
              DGRUN(J)=DBLE(GRUN(J))
            END DO
            CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
            DO J=1,601
              GRUN(J)=DCMPLX(DGRUN(J))
            END DO
          ELSE
            CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
          END IF
          EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!          CALL DECOUTCOUP(QNH*EXP(T),GRUN)
        END DO
      ELSE IF(MDEC.GT.QIN)THEN
!
!Running up past the relevant thresholds
!
        IF(STTH.NE.FINTH)THEN
          DO I=STTH,FINTH-1
            IF(QTHSORT(I).EQ.QTHSORT(I+1))GOTO 42
            DT=(TTH(I+1)-TTH(I))/FLOAT(NSTEPTH)
            DO II=1,NSTEPTH
              T=TTH(I)+(TTH(I+1)-TTH(I))*FLOAT(II-1)/FLOAT(NSTEPTH)
              SSQSTEP=QNH*EXP(T)
!
              IF(BELOWMS.EQ.1.AND.QNH*EXP(T+DT).GT.RGEMS)THEN
                BELOWMS=0
                CALL STROTBACK(GRUN,GFIX,1)
                GFIX(12)=GFIX(12)/(1.D0-DBLE(RTISA))
                GFIX(21)=GFIX(21)/(1.D0-DBLE(RBISA))
                GFIX(30)=GFIX(30)/(1.D0-DBLE(RLISA))
                IF(RGEMS.LT.QNH)THEN
                  GFIX(120)=GFIX(120)/(1.D0-DBLE(RTISA))
                  GFIX(129)=GFIX(129)/(1.D0-DBLE(RBISA))
                  GFIX(138)=GFIX(138)/(1.D0-DBLE(RLISA))
                END IF
                CALL STROTATE(GFIX,GRUN,1)
              END IF
!
              IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
              CALL CHINT(SSQSTEP,GRUN)
              IF(COMP.EQ.0)THEN
                DO J=1,601
                  DGRUN(J)=DBLE(GRUN(J))
                END DO
                CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
                DO J=1,601
                  GRUN(J)=DCMPLX(DGRUN(J))
                END DO
              ELSE
                CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
              END IF
              EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            END DO
  42      END DO
        END IF
!
!Now the final run to MDEC
!
        DT=(TMDEC-TTH(FINTH))/FLOAT(NSTEP)
        DO II=1,NSTEP
          T=TTH(FINTH)+(TMDEC-TTH(FINTH))*FLOAT(II-1)/FLOAT(NSTEP)
          SSQSTEP=QNH*EXP(T)
          Q=QNH*EXP(T)
!
          IF(BELOWMS.EQ.1.AND.QNH*EXP(T+DT).GT.RGEMS)THEN
            BELOWMS=0
            CALL STROTBACK(GRUN,GFIX,1)
            GFIX(12)=GFIX(12)/(1.D0-DBLE(RTISA))
            GFIX(21)=GFIX(21)/(1.D0-DBLE(RBISA))
            GFIX(30)=GFIX(30)/(1.D0-DBLE(RLISA))
            IF(RGEMS.LT.QNH)THEN
              GFIX(120)=GFIX(120)/(1.D0-DBLE(RTISA))
              GFIX(129)=GFIX(129)/(1.D0-DBLE(RBISA))
              GFIX(138)=GFIX(138)/(1.D0-DBLE(RLISA))
            END IF
            CALL STROTATE(GFIX,GRUN,1)
          END IF
!
          IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          CALL CHINT(SSQSTEP,GRUN)
          IF(COMP.EQ.0)THEN
            DO J=1,601
              DGRUN(J)=DBLE(GRUN(J))
            END DO
            CALL DRKSTP(601,DT,T,DGRUN,DRGE601,DW)
            DO J=1,601
              GRUN(J)=DCMPLX(DGRUN(J))
            END DO
          ELSE
            CALL CRKSTP(601,DT,T,GRUN,CRGE601,W)
          END IF
          EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
      ELSE
        T=TQIN
      END IF
!
      QIN=QNH*EXP(T)
!
!Work out the final squark mass matrices in the quark mass basis
!
      CALL SQDIAG(GRUN)
      CALL MASSSQM(GRUN)
!
      CALL STROTBACK(GRUN,GOUT,0)
!
!Save the new values of GHIK if m_H<MDEC
!
      IF(QNH.LT.QIN)THEN
        DO I=1,601
          GHIK(I)=GOUT(I)
        END DO
      END IF
!
      CLOSE(63)
!
      RETURN
      END
!
      SUBROUTINE DIAGSQM(QMSQ,UMSQ,DMSQ,LMSQ,EMSQ)
!
!Purpose: To use the results of MASSSQM for the squark mass
!         squared matrices in the quark mass basis to find the
!         eigenvalues correspoding to u,c,t ; d,s,b ; e,nu,tau
!         in that order. Used by EWSB routine.
!
      IMPLICIT NONE
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      DOUBLE COMPLEX QVA(3),UVA(3),DVA(3),LVA(3),EVA(3)
      DOUBLE COMPLEX QVE(3,3),UVE(3,3),DVE(3,3),LVE(3,3),EVE(3,3)
      DOUBLE COMPLEX QIN(3,3),UIN(3,3),DIN(3,3),EIN(3,3),LIN(3,3)
      DOUBLE COMPLEX EVERTMP(3,3),CWORK1(99),CWORK2(6)
      DOUBLE PRECISION VALUEQ,VALUEU,VALUED,VALUEL,VALUEE
      DOUBLE PRECISION QMSQ(3),UMSQ(3),DMSQ(3),LMSQ(3),EMSQ(3)
      INTEGER I,J
      INTEGER ENTRYQ,ENTRYU,ENTRYD,ENTRYL,ENTRYE
      INTEGER CIERR1,CIERR2,CIERR3,CIERR4,CIERR5
!
      DO I=1,3
        DO J=1,3
          QIN(I,J)=MQQMASS(I,J)
          UIN(I,J)=MUQMASS(I,J)
          DIN(I,J)=MDQMASS(I,J)
          LIN(I,J)=MLQMASS(I,J)
          EIN(I,J)=MEQMASS(I,J)
        END DO
      END DO
!
!Diagonalise and recombine into complex quantities
!
      CALL ZGEEV('V','V',3,QIN,3,QVA,QVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR1)
      CALL ZGEEV('V','V',3,UIN,3,UVA,UVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR2)
      CALL ZGEEV('V','V',3,DIN,3,DVA,DVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR3)
      CALL ZGEEV('V','V',3,LIN,3,LVA,LVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR4)
      CALL ZGEEV('V','V',3,EIN,3,EVA,EVE,3,EVERTMP,3,CWORK1,99
     $           ,CWORK2,CIERR5)
!
!Now work out the correct order of the eigenvalues
!depending on the eigenvector.
!
      DO J=1,3
        ENTRYQ=1
        VALUEQ=ABS(QVE(1,J))
        ENTRYU=1
        VALUEU=ABS(UVE(1,J))
        ENTRYD=1
        VALUED=ABS(DVE(1,J))
        ENTRYL=1
        VALUEL=ABS(LVE(1,J))
        ENTRYE=1
        VALUEE=ABS(EVE(1,J))
        DO I=2,3
          IF(ABS(QVE(I,J)).GT.VALUEQ)THEN
            VALUEQ=ABS(QVE(I,J))
            ENTRYQ=I
          END IF
          IF(ABS(UVE(I,J)).GT.VALUEU)THEN
            VALUEU=ABS(UVE(I,J))
            ENTRYU=I
          END IF
          IF(ABS(DVE(I,J)).GT.VALUED)THEN
            VALUED=ABS(DVE(I,J))
            ENTRYD=I
          END IF
          IF(ABS(LVE(I,J)).GT.VALUEL)THEN
            VALUEL=ABS(LVE(I,J))
            ENTRYL=I
          END IF
          IF(ABS(EVE(I,J)).GT.VALUEE)THEN
            VALUEE=ABS(EVE(I,J))
            ENTRYE=I
          END IF
        END DO
        QMSQ(ENTRYQ)=DSQRT(ABS(QVA(J))**2)
        UMSQ(ENTRYU)=DSQRT(ABS(UVA(J))**2)
        DMSQ(ENTRYD)=DSQRT(ABS(DVA(J))**2)
        LMSQ(ENTRYL)=DSQRT(ABS(LVA(J))**2)
        EMSQ(ENTRYE)=DSQRT(ABS(EVA(J))**2)
      END DO
!
      RETURN
      END
!
      SUBROUTINE DOWNMHCOND
!
!Purpose: Apply the matching conditions at m_H when running down.
!         If m_H is greater than m_SUSY, the values of decoupled
!         operators are frozen here for use later.
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      DOUBLE PRECISION SINB,COSB
      DOUBLE COMPLEX B,VU,VD
      INTEGER I,J
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
      B=G(109)
      VU=G(110)
      VD=G(111)
!
      G(429)=(3.D0/5.D0*(G(1))**2+(G(2))**2)/4.D0
     $                                  *((1+TANB**2)/(1-TANB**2))**2
      G(428)=DSQRT(DBLE(VU*CONJG(VU))+DBLE(VD*CONJG(VD)))
      DO I=1,3
        DO J=1,3
          G(111+(I-1)*3+J)=SINB*G(3+(I-1)*3+J)
          G(120+(I-1)*3+J)=COSB*G(12+(I-1)*3+J)
          G(129+(I-1)*3+J)=COSB*G(21+(I-1)*3+J)
          G(399+(I-1)*3+J)=SINB*G(33+(I-1)*3+J)-COSB*G(429+(I-1)*3+J)
          G(408+(I-1)*3+J)=COSB*G(42+(I-1)*3+J)-SINB*G(438+(I-1)*3+J)
          G(417+(I-1)*3+J)=COSB*G(51+(I-1)*3+J)-SINB*G(447+(I-1)*3+J)
        END DO
      END DO
      G(287)=SINB*G(184)
      G(288)=COSB*G(185)
      G(289)=SINB*G(204)
      G(290)=COSB*G(205)
      G(427)=SINB**2*G(61)+COSB**2*G(62)-2.D0*SINB*COSB*B
!
!Save the combinations necessary for the EWSB conditions if m_H > m_SUSY
!
      IF(RGEMS.LT.QNH)THEN
        CSHSH=SINB**2*G(61)+COSB**2*G(62)-2.D0*SINB*COSB*B
        CSHLH=SINB*COSB*G(61)-SINB*COSB*G(62)+(SINB**2-COSB**2)*B
        CLHLH=COSB**2*G(61)+SINB**2*G(62)+2.D0*SINB*COSB*B
        DO I=1,3
          DO J=1,3
            CULSHUR(I,J)=SINB*G(33+(I-1)*3+J)-COSB*G(429+(I-1)*3+J)
            CDLSHDR(I,J)=COSB*G(42+(I-1)*3+J)-SINB*G(438+(I-1)*3+J)
            CELSHER(I,J)=COSB*G(51+(I-1)*3+J)-SINB*G(447+(I-1)*3+J)
            CULLHUR(I,J)=COSB*G(33+(I-1)*3+J)+SINB*G(429+(I-1)*3+J)
            CDLLHDR(I,J)=-SINB*G(42+(I-1)*3+J)-COSB*G(438+(I-1)*3+J)
            CELLHER(I,J)=-SINB*G(51+(I-1)*3+J)-COSB*G(447+(I-1)*3+J)
          END DO
        END DO
      END IF
!
      RETURN
      END
!
      SUBROUTINE DOWNMHIGHMH(QEND,NSTEP)
!
!Purpose: The final run down before stop decay routine is called.
!         NOTE: If m_H is less than M_t, only run to M_t and 
!               make sure that this value is passed to SQSIX
!
      IMPLICIT NONE
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      DOUBLE PRECISION TTH(20),THIGH,TMT,T,QEND,SINB,COSB
      DOUBLE PRECISION DW(1803),THSP,PI,DT
      DOUBLE COMPLEX B,VU,VD
      DOUBLE COMPLEX W(1803),ID(3,3)
      INTEGER NSTEP,I,J,II,NHTMT,BELOWMS
      EXTERNAL CRGE601,DRGE601
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $ ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $ ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $ ,ID(3,3)/(1.D0,0.D0)/
!
      BELOWMS=0
      PI=4.D0*DATAN(1.D0)
!
      THIGH=0.D0
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/MHIGH)
      END DO
      TMT=LOG(MT/MHIGH)
!
!Check location of lowest threshold above m_t
!
      NHTMT=NLTMT+1
!
      IF(NHTMT.NE.21)THEN
        DT=(TTH(20)-THIGH)/FLOAT(NSTEP)
      ELSE
        DT=(TMT-THIGH)/FLOAT(NSTEP)
      END IF
!
      DO II=1,NSTEP
        IF(NHTMT.NE.21)THEN
          T=THIGH+(TTH(20)-THIGH)*FLOAT(II-1)/FLOAT(NSTEP)
        ELSE
          T=THIGH+(TMT-THIGH)*FLOAT(II-1)/FLOAT(NSTEP)
        END IF
        SSQSTEP=MHIGH*EXP(T)
!
        IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
          BELOWMS=1
          CALL ROTBACK(1)
          G(12)=G(12)*(1.D0-DBLE(RTISA))
          G(21)=G(21)*(1.D0-DBLE(RBISA))
          G(30)=G(30)*(1.D0-DBLE(RLISA))
          IF(RGEMS.LE.QNH)THEN
            G(120)=G(120)*(1.D0-DBLE(RTISA))
            G(129)=G(129)*(1.D0-DBLE(RBISA))
            G(138)=G(138)*(1.D0-DBLE(RLISA))
          END IF
          CALL ROTATE(1)
        END IF
!
        IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!
        CALL SQDIAG(G)
        CALL CHDEC(SSQSTEP,G,1)
!
        IF(COMP.EQ.0)THEN
          DO I=1,601
            DG(I)=DBLE(G(I))
          END DO
          CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
          DO I=1,601
            G(I)=DCMPLX(DG(I))
          END DO
        ELSE
          CALL CRKSTP(601,DT,T,G,CRGE601,W)
        END IF
!
        EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(NHTMT.NE.21)THEN
        DO I=19,1,-1
          IF(NSTEPTHRESH(I).EQ.0)GOTO 40
          IF(I.LT.LOCMH)GOTO 50
          DT=(TTH(I)-TTH(I+1))/FLOAT(NSTEPTHRESH(I))
!
          DO II=1,NSTEPTHRESH(I)
            T=TTH(I+1)+(TTH(I)-TTH(I+1))*FLOAT(II-1)
     $                                  /FLOAT(NSTEPTHRESH(I))
            SSQSTEP=MHIGH*EXP(T)
!
            IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
              BELOWMS=1
              CALL ROTBACK(1)
              G(12)=G(12)*(1.D0-DBLE(RTISA))
              G(21)=G(21)*(1.D0-DBLE(RBISA))
              G(30)=G(30)*(1.D0-DBLE(RLISA))
              IF(RGEMS.LE.QNH)THEN
                G(120)=G(120)*(1.D0-DBLE(RTISA))
                G(129)=G(129)*(1.D0-DBLE(RBISA))
                G(138)=G(138)*(1.D0-DBLE(RLISA))
              END IF
              CALL ROTATE(1)
            END IF
!
            IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            CALL SQDIAG(G)
            CALL CHDEC(SSQSTEP,G,1)
            IF(COMP.EQ.0)THEN
              DO J=1,601
                DG(J)=DBLE(G(J))
              END DO
              CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
              DO J=1,601
                G(J)=DCMPLX(DG(J))
              END DO
            ELSE
              CALL CRKSTP(601,DT,T,G,CRGE601,W)
            END IF
            EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          END DO
!
  40    END DO
  50    CONTINUE
      END IF
      IF(NHTMT.EQ.21)THEN
        WRITE(*,*)'WARNING: M_H IS LESS THAN M_T'
      END IF
!
!SM assignments at m_H
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
      B=G(109)
      VU=G(110)
      VD=G(111)
!
      G(429)=(3.D0/5.D0*(G(1))**2+(G(2))**2)/4.D0
     $                                  *((1+TANB**2)/(1-TANB**2))**2
      G(428)=DSQRT(DBLE(VU*CONJG(VU))+DBLE(VD*CONJG(VD)))
      DO I=1,3
        DO J=1,3
          G(111+(I-1)*3+J)=SINB*G(3+(I-1)*3+J)
          G(120+(I-1)*3+J)=COSB*G(12+(I-1)*3+J)
          G(129+(I-1)*3+J)=COSB*G(21+(I-1)*3+J)
          G(399+(I-1)*3+J)=SINB*G(33+(I-1)*3+J)-COSB*G(429+(I-1)*3+J)
          G(408+(I-1)*3+J)=COSB*G(42+(I-1)*3+J)-SINB*G(438+(I-1)*3+J)
          G(417+(I-1)*3+J)=COSB*G(51+(I-1)*3+J)-SINB*G(447+(I-1)*3+J)
        END DO
      END DO
      G(287)=SINB*G(184)
      G(288)=COSB*G(185)
      G(289)=SINB*G(204)
      G(290)=COSB*G(205)
      G(427)=SINB**2*G(61)+COSB**2*G(62)-2.D0*SINB*COSB*B
!
      QEND=MHIGH*EXP(T)
!
      RETURN
      END
!
      SUBROUTINE DOWNMHIGHMZ(QEND,NSTEP,THCNG)
!
!Purpose: To run the couplings down from M_HIGH to M_H and
!         rotate back into the mass basis
!
!         THCNG is used to switch off the changing of the
!         thresholds before the final run to allow the Yukawas
!         to become diagonal at M_t.
!
      IMPLICIT NONE
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
!DEC:
!The common block for changing the decoupling points during each run.
!The order of the entries in NEWTH and BELOW are: mq(1-3), mu(1-3),
!md(1-3), ml(1-3), me(1-3), higgsino (mu), m_gluino, m_A0, M_1, M_2.
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      DOUBLE PRECISION TTH(20),THIGH,DT,T,QEND,DW(1803)
      DOUBLE PRECISION TMT,TZ,PI,DGSM(32),DWSM(96)
      DOUBLE COMPLEX W(1803),ID(3,3),GSM(32),WSM(96)
      INTEGER NSTEP,I,J,K,II,NHTMT,NSTEPMT,LOOPNSTEP,BELOWMS,THCNG
      INTEGER NSTEPSM
      EXTERNAL CRGE601,DRGE601,CSMRGEDR,DSMRGEDR
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $ ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $ ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $ ,ID(3,3)/(1.D0,0.D0)/
!
!USED BY OUTCOUP. CLOSE AT THE END OF THE PROG.
!
!      OPEN(71,FILE='out/dau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(72,FILE='out/dad.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(73,FILE='out/dae.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(74,FILE='out/dmhumu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(75,FILE='out/dmhdmu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(76,FILE='out/dmq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(77,FILE='out/dmup.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(78,FILE='out/dmd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(79,FILE='out/dml.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(80,FILE='out/dme.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(81,FILE='out/dmu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(82,FILE='out/db.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(83,FILE='out/dmtsfu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(84,FILE='out/dmtsfd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(85,FILE='out/dmtsfe.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(86,FILE='out/dtriu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(87,FILE='out/dtrid.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(88,FILE='out/dtrie.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(89,FILE='out/dmhud.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(90,FILE='out/daum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(91,FILE='out/dadm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(92,FILE='out/daem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(93,FILE='out/dmhumum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(94,FILE='out/dmhdmum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(95,FILE='out/dmqm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(96,FILE='out/dmupm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(97,FILE='out/dmdm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(98,FILE='out/dmlm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(99,FILE='out/dmem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(70,FILE='out/dmum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(69,FILE='out/dbm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(68,FILE='out/dg1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(67,FILE='out/dg2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(66,FILE='out/dg3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(65,FILE='out/dytau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(64,FILE='out/dyb.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(63,FILE='out/dyu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(62,FILE='out/dm1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(61,FILE='out/dm2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(60,FILE='out/dm3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(59,FILE='out/dgtpq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(58,FILE='out/dgtpu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(57,FILE='out/dgtq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(56,FILE='out/dftuq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(55,FILE='out/dvu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(54,FILE='out/dvd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(53,FILE='out/dmuflavcomp.dat',STATUS='UNKNOWN'
!     $ ,FORM='FORMATTED')
!      OPEN(52,FILE='out/dm12comp.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!
      BELOWMS=0
      PI=4.D0*DATAN(1.D0)
!
!Set the initial values for changing the thresholds.
!The new thresholds are set to the old ones in case they
!are not changed.
!
      DO I=1,20
        BELOW(I)=0
      END DO
      DO I=1,3
        NEWTH(I)=QTHQL(I)
        NEWTH(I+3)=QTHUR(I)
        NEWTH(I+6)=QTHDR(I)
        NEWTH(I+9)=QTHLL(I)
        NEWTH(I+12)=QTHER(I)
      END DO
      NEWTH(16)=QNSH
      NEWTH(17)=QNSG
      NEWTH(18)=QNH
      NEWTH(19)=QTHSB
      NEWTH(20)=QTHSW
!
      THIGH=0.D0
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/MHIGH)
      END DO
      TMT=LOG(MT/MHIGH)
!
!Check location of lowest threshold above m_t
!
      NHTMT=NLTMT+1
!
      IF(NHTMT.NE.21)THEN
        DT=(TTH(20)-THIGH)/FLOAT(NSTEP)
      ELSE
        DT=(TMT-THIGH)/FLOAT(NSTEP)
      END IF
!
      DO II=1,NSTEP
        IF(NHTMT.NE.21)THEN
          T=THIGH+(TTH(20)-THIGH)*FLOAT(II-1)/FLOAT(NSTEP)
        ELSE
          T=THIGH+(TMT-THIGH)*FLOAT(II-1)/FLOAT(NSTEP)
        END IF
        SSQSTEP=MHIGH*EXP(T)
!
        IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!
!Carry out the squark rotation business to make sure that we are
!in the correct basis for the squark thresholds
!
        CALL SQDIAG(G)
        CALL CHDEC(SSQSTEP,G,THCNG)
!        CALL OUTCOUP(MHIGH*EXP(T),0)
!
!Finite corrections to the Yukawas from Isajet and other
!m_SUSY conditions.
!
        IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
          BELOWMS=1
          CALL DOWNMSCOND
        END IF
!
        IF(COMP.EQ.0)THEN
          DO I=1,601
            DG(I)=DBLE(G(I))
          END DO
          CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
          DO I=1,601
            G(I)=DCMPLX(DG(I))
          END DO
        ELSE
          CALL CRKSTP(601,DT,T,G,CRGE601,W)
        END IF
        EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(NHTMT.NE.21)THEN
!
        IF(LOCMH.EQ.20)CALL DOWNMHCOND
!
        DO I=19,NHTMT,-1
          IF(NSTEPTHRESH(I).EQ.0)GOTO 40
          DT=(TTH(I)-TTH(I+1))/FLOAT(NSTEPTHRESH(I))
!
          DO II=1,NSTEPTHRESH(I)
            T=TTH(I+1)+(TTH(I)-TTH(I+1))*FLOAT(II-1)
     $                                  /FLOAT(NSTEPTHRESH(I))
            SSQSTEP=MHIGH*EXP(T)
            IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            CALL SQDIAG(G)
            CALL CHDEC(SSQSTEP,G,THCNG)
!            CALL OUTCOUP(MHIGH*EXP(T),0)
            IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
              BELOWMS=1
              CALL DOWNMSCOND
            END IF
            IF(COMP.EQ.0)THEN
              DO J=1,601
                DG(J)=DBLE(G(J))
              END DO
              CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
              DO J=1,601
                G(J)=DCMPLX(DG(J))
              END DO
            ELSE
              CALL CRKSTP(601,DT,T,G,CRGE601,W)
            END IF
            EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          END DO
!
   40     IF(LOCMH.EQ.I)CALL DOWNMHCOND
        END DO
!
!Continue running to m_t
!
        IF(NHTMT.NE.1)THEN
          DT=(TMT-TTH(NHTMT))/FLOAT(NSTEPTHRESH(NHTMT-1))
          LOOPNSTEP=NSTEPTHRESH(NHTMT-1)
        ELSE
          NSTEPMT=INT(ABS(DLOG(MT/QTHSORT(1)))*NSTEP/DLOG(MHIGH/MT)*25)
          DT=(TMT-TTH(NHTMT))/FLOAT(NSTEPMT)
          LOOPNSTEP=NSTEPMT
        END IF
!
        DO II=1,LOOPNSTEP
          T=TTH(NHTMT)+(TMT-TTH(NHTMT))*FLOAT(II-1)/FLOAT(LOOPNSTEP)
          SSQSTEP=MHIGH*EXP(T)
          IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          CALL SQDIAG(G)
          CALL CHDEC(SSQSTEP,G,THCNG)
!          CALL OUTCOUP(MHIGH*EXP(T),0)
          IF(BELOWMS.EQ.0.AND.SSQSTEP.LT.RGEMS)THEN
            BELOWMS=1
            CALL DOWNMSCOND
          END IF
          IF(COMP.EQ.0)THEN
            DO I=1,601
              DG(I)=DBLE(G(I))
            END DO
            CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
            DO I=1,601
              G(I)=DCMPLX(DG(I))
            END DO
          ELSE
            CALL CRKSTP(601,DT,T,G,CRGE601,W)
          END IF
          EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
!
      ELSE
        WRITE(*,*)'WARNING: ALL THRESHOLDS ARE LESS THAN MT'
        CALL DOWNMHCOND
      END IF
!
!Finish the running by continuing to M_Z. This is primarily done
!since on the first run up, the value of the quartic coupling in the
!SM was zero, as was lambda_t.
!
      NSTEPSM=50
      TZ=DLOG(MZ/MHIGH)
      DT=(TZ-TMT)/FLOAT(NSTEPSM)
      NU=3 !Decouple the top at M_Z
      DO I=1,3
        GSM(I)=G(I)
      END DO
      DO I=4,30
        GSM(I)=G(108+I)
      END DO
      GSM(31)=G(429)
      GSM(32)=G(428)
!
!Rotate into the mass basis - we will not run the full flavour
!SM RGEs
!
      CALL ROTBACKSM(GSM)
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          DGSM(I)=DBLE(GSM(I))
        END DO
      END IF
!
      DO II=1,NSTEPSM
        T=TMT+(TZ-TMT)*FLOAT(II-1)/FLOAT(NSTEPSM)
        SMQSTEP=MHIGH*EXP(T)
        IF(II.EQ.1)EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        IF(COMP.EQ.0)THEN
          CALL DRKSTP(32,DT,T,DGSM,DSMRGEDR,DWSM)
        ELSE
          CALL CRKSTP(32,DT,T,GSM,CSMRGEDR,WSM)
        END IF
        EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          GSM(I)=DCMPLX(DGSM(I))
        END DO
      END IF
!
      QEND=MHIGH*EXP(T)
!
!If any of the entries in BELOW are still 0, then
!inform the user.
!
!      DO I=1,20
!        IF(BELOW(I).EQ.0)WRITE(*,*)'THRESH ',I,' IS BELOW MT'
!      END DO
!
!Save the weak scale SM values of the quartic coupling and up Yukawa
!
      LAMBDAMZ=GSM(31)
      LAMTMZ=GSM(12)
!
!Now replace the old thresholds with the new ones
!
      DO I=1,3
        QTHQL(I)=NEWTH(I)
        QTHUR(I)=NEWTH(I+3)
        QTHDR(I)=NEWTH(I+6)
        QTHLL(I)=NEWTH(I+9)
        QTHER(I)=NEWTH(I+12)
      END DO
      QNSH=NEWTH(16)
      QNSG=NEWTH(17)
      QNH=NEWTH(18)
      QTHSB=NEWTH(19)
      QTHSW=NEWTH(20)
      CALL SORTTH
!
!Recheck the location of the first threshold above m_t
!
      NLTMT=0
      DO I=1,20
        IF(QTHSORT(I).LE.MT)NLTMT=I
      END DO
!
!Put out warning for all sfermion masses less than 0.D0
!
      DO I=1,3
        IF(QTHQL(I).LT.MT+1.D0.AND.DBLE(MQVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE UP-LEFT SQUARK EIGENVALUE'
        IF(QTHUR(I).LT.MT+1.D0.AND.DBLE(MUPVA(I)).LE.0.d0)
     $       WRITE(*,*)'WARNING: NEGATIVE UP-RIGHT SQUARK EIGENVALUE'
        IF(QTHDR(I).LT.MT+1.D0.AND.DBLE(MDVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE DOWN-RIGHT SQUARK EIGENVALUE'
        IF(QTHLL(I).LT.MT+1.D0.AND.DBLE(MLVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE L-SLEPTON EIGENVALUE'
        IF(QTHER(I).LT.MT+1.D0.AND.DBLE(MEVA(I)).LE.0.D0)
     $       WRITE(*,*)'WARNING: NEGATIVE E-SLEPTON EIGENVALUE'
      END DO
!
!Fix the number of steps in case the thresholds have changed.
!
      IF(NLTMT.LT.19)THEN
        IF(NLTMT.NE.0)THEN
          NSTEPTHRESH(NLTMT)=INT(ABS(DLOG(MT/QTHSORT(NLTMT+1)))
     $                                     *NSTEP/DLOG(MHIGH/MT)*25)
        END IF
        DO I=NLTMT+1,19
          NSTEPTHRESH(I)=INT(ABS(DLOG(QTHSORT(I)/QTHSORT(I+1)))
     $                                     *NSTEP/DLOG(MHIGH/MT)*25)
          IF(NSTEPTHRESH(I).EQ.0.AND.
     $                            ABS(QTHSORT(I)-QTHSORT(I+1)).GT.1D-10)
     $                                              NSTEPTHRESH(I)=10
        END DO
      END IF
!
!If the location of the thresholds is changing, it may be useful to the
!user to know their values. The following lines can be reinstated to
!write them to the screen, or changed to print to a file.
!
!      IF(THCNG.EQ.0)THEN
!        WRITE(*,*)
!        WRITE(*,*)'RGEMS',RGEMS
!        WRITE(*,*)'MU (16)',QNSH
!        WRITE(*,*)'GLUINO (17)',QNSG
!        WRITE(*,*)'HIGGS (18)',QNH
!        WRITE(*,*)'BINO (19)',QTHSB
!        WRITE(*,*)'WINO (20)',QTHSW
!        DO I=1,3
!        WRITE(*,*)'SQUARKS (Q,UR,DR) (1,4,7)',I,' '
!     $                                    ,QTHQL(I),QTHUR(I),QTHDR(I)
!        END DO
!        DO I=1,3
!          WRITE(*,*)'SLEPTONS (L,ER) (10,13)',I,' ',QTHLL(I),QTHER(I)
!        END DO
!      END IF
!
      IF(NHTMT.GT.LOCMH)CALL DOWNMHCOND
!
!      CLOSE(52)
!      CLOSE(53)
!      CLOSE(54)
!      CLOSE(55)
!      CLOSE(56)
!      CLOSE(57)
!      CLOSE(58)
!      CLOSE(59)
!      CLOSE(60)
!      CLOSE(61)
!      CLOSE(62)
!      CLOSE(63)
!      CLOSE(64)
!      CLOSE(65)
!      CLOSE(66)
!      CLOSE(67)
!      CLOSE(68)
!      CLOSE(69)
!      CLOSE(70)
!      CLOSE(71)
!      CLOSE(72)
!      CLOSE(73)
!      CLOSE(74)
!      CLOSE(75)
!      CLOSE(76)
!      CLOSE(77)
!      CLOSE(78)
!      CLOSE(79)
!      CLOSE(80)
!      CLOSE(81)
!      CLOSE(82)
!      CLOSE(83)
!      CLOSE(84)
!      CLOSE(85)
!      CLOSE(86)
!      CLOSE(87)
!      CLOSE(88)
!      CLOSE(89)
!      CLOSE(90)
!      CLOSE(91)
!      CLOSE(92)
!      CLOSE(93)
!      CLOSE(94)
!      CLOSE(95)
!      CLOSE(96)
!      CLOSE(97)
!      CLOSE(98)
!      CLOSE(99)
!
      RETURN
      END
!
      SUBROUTINE DOWNMSCOND
!
!Purpose: Apply various conditions at M_SUSY
!         The EWSB relations are imposed to derive various values.
!         Two versions depending on if m_H is less than m_SUSY
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      REAL G0(601),QIN,MUIN,SIG1,SIG2 !The variables used by SUGEFFFL
      DOUBLE COMPLEX BMHM1SQ,BMHM2SQ,BMHB,BMHAU(3,3),BMHAD(3,3)
      DOUBLE COMPLEX BMHAE(3,3),BMHMTSFU(3,3),BMHMTSFD(3,3)
      DOUBLE COMPLEX BMHMTSFE(3,3),GSAV(601),MDIFF
      DOUBLE COMPLEX BMHM1SQNEW,BMHM2SQNEW
      DOUBLE PRECISION QMSQ(3),UMSQ(3),DMSQ(3),LMSQ(3),EMSQ(3)
      DOUBLE PRECISION SINB,COSB,VU,VD,VUSQ,VDSQ,VSM,MUMSQ
      INTEGER I,J
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
      VU=DBLE(G(110))
      VD=DBLE(G(111))
      VSM=DBLE(G(428))
      VUSQ=VU**2
      VDSQ=VD**2
!
!Save the original G entries
!
      DO I=1,601
        GSAV(I)=G(I)
      END DO
!
!Solve the simultaneous equations
!
      IF(RGEMS.LT.QNH)THEN
        CSHSH=G(427)
        DO I=1,3
          DO J=1,3
            CULSHUR(I,J)=G(399+(I-1)*3+J)
            CDLSHDR(I,J)=G(408+(I-1)*3+J)
            CELSHER(I,J)=G(417+(I-1)*3+J)
          END DO
        END DO
        BMHM1SQ=COSB**2*CSHSH-2.D0*SINB*COSB*CSHLH+SINB**2*CLHLH
        BMHM2SQ=SINB**2*CSHSH+2.D0*SINB*COSB*CSHLH+COSB**2*CLHLH
        BMHB=SINB*COSB*(-CSHSH+CLHLH)+(SINB**2-COSB**2)*CSHLH
        DO I=1,3
          DO J=1,3
            BMHAU(I,J)=SINB*CULSHUR(I,J)+COSB*CULLHUR(I,J)
            BMHAD(I,J)=COSB*CDLSHDR(I,J)-SINB*CDLLHDR(I,J)
            BMHAE(I,J)=COSB*CELSHER(I,J)-SINB*CELLHER(I,J)
            BMHMTSFU(I,J)=-COSB*CULSHUR(I,J)+SINB*CULLHUR(I,J)
            BMHMTSFD(I,J)=-SINB*CDLSHDR(I,J)-COSB*CDLLHDR(I,J)
            BMHMTSFE(I,J)=-SINB*CELSHER(I,J)-COSB*CELLHER(I,J)
!
            G(33+(I-1)*3+J)=BMHAU(I,J)
            G(42+(I-1)*3+J)=BMHAD(I,J)
            G(51+(I-1)*3+J)=BMHAE(I,J)
            G(429+(I-1)*3+J)=BMHMTSFU(I,J)
            G(438+(I-1)*3+J)=BMHMTSFD(I,J)
            G(447+(I-1)*3+J)=BMHMTSFE(I,J)
          END DO
        END DO
      END IF
!
!Find the squark mass eigenvalues in the quark mass basis
!
      CALL MASSSQM(G)
      CALL DIAGSQM(QMSQ,UMSQ,DMSQ,LMSQ,EMSQ)
!
!Apply the finite shifts to the third generation Yukawas
!
      CALL ROTBACK(1)
      G(12)=G(12)*(1.D0-DBLE(RTISA))
      G(21)=G(21)*(1.D0-DBLE(RBISA))
      G(30)=G(30)*(1.D0-DBLE(RLISA))
      IF(RGEMS.LE.QNH)THEN
        G(120)=G(120)*(1.D0-DBLE(RTISA))
        G(129)=G(129)*(1.D0-DBLE(RBISA))
        G(138)=G(138)*(1.D0-DBLE(RLISA))
      END IF
!
!Set up and call the 1-loop radiative corrections to the scalar potential
!Since Isajet has differing sign conventions, we flip the signs
!of \mu and the a-parameters.
!
      QIN=REAL(RGEMS)
      MUIN=-REAL(REAL(G(108)))
      DO I=1,601
        G0(I)=REAL(REAL(G(I))) !The important quantities have Im()~0
      END DO
      DO I=1,27
        G0(33+I)=-G0(33+I)
        G0(323+I)=-G0(323+I)
        G0(399+I)=-G0(399+I)
        G0(429+I)=-G0(429+I)
      END DO
!
!Fix the values of the sfermion masses
!
      DO I=1,3
        G0(62+(I-1)*3+I)=REAL(QMSQ(I))
        G0(80+(I-1)*3+I)=REAL(UMSQ(I))
        G0(89+(I-1)*3+I)=REAL(DMSQ(I))
        G0(71+(I-1)*3+I)=REAL(LMSQ(I))
        G0(98+(I-1)*3+I)=REAL(EMSQ(I))
      END DO
!
      IF(RGEMS.LT.QNH)THEN
        G0(110)=REAL(SINB*VSM)
        G0(111)=REAL(COSB*VSM)
      END IF
!
!NOTE: Following Isajet's notation, SIG1 is the correction
!      to the down-type Higgs field.
!
      CALL SUGEFFFL(QIN,G0,MUIN,SIG1,SIG2)
!
!Now use EWSB conditions to fix b and (m^2_1+m^2_2).
!Also fix b and mu for the MSSM running couplings.
!
      IF(RGEMS.LT.QNH)THEN
        MDIFF=BMHM1SQ-BMHM2SQ
        BMHM1SQNEW=((-(MDIFF+SIG1-SIG2)*1.D0/2.D0*(1.D0+TANB**2)+1.D0/4.D0
     $          *(G(1)**2+G(2)**2)*(SINB**2-COSB**2)/COSB**2*VSM**2)
     $          /(1.D0-TANB**2))+1.D0/2.D0*(MDIFF+SIG1-SIG2)-SIG1
        BMHM2SQNEW=BMHM1SQNEW-MDIFF
!
!Set the new value of b with the old m1 and m2
!
        BMHB=(BMHM1SQ+BMHM2SQ+SIG1+SIG2)*SINB*COSB
!
!Now replace the values of m1 and m2
!
        BMHM1SQ=BMHM1SQNEW
        BMHM2SQ=BMHM2SQNEW
!
!The MSSM mu parameter is fixed next
!
        MUMSQ=((DBLE(G(351))+SIG2)*TANB**2-DBLE(G(352))-SIG1)
     $       /(1.D0-TANB**2)-1.D0/4.D0*(DBLE(G(291))**2+DBLE(G(292))**2)
     $       *VSM**2
        IF(MUMSQ.LT.0.D0)MUMSQ=MZ**2
!
!Now reset the MSSM mu and b
!
        G(399)=(G(351)+G(352)+SIG1+SIG2+2.D0*G(398)**2)*SINB*COSB
        G(398)=DSQRT(MUMSQ)*EXP((0.D0,1.D0)*PHASEMU)
!
!Reset the frozen parameters
!
        CSHSH=SINB**2*BMHM2SQ+COSB**2*BMHM1SQ-2.D0*SINB*COSB*BMHB
        CSHLH=SINB*COSB*BMHM2SQ-SINB*COSB*BMHM1SQ+(SINB**2-COSB**2)*BMHB
        CLHLH=COSB**2*BMHM2SQ+SINB**2*BMHM1SQ+2.D0*SINB*COSB*BMHB
!
!Reset the running parameter
!
        G(427)=CSHSH
      ELSE
!
!Carry out the greater than m_H conditions in a similar manner
!
        MDIFF=G(62)-G(61)
        BMHM1SQNEW=((-(MDIFF+SIG1-SIG2)*1.D0/2.D0*(1.D0+VUSQ/VDSQ)
     $       +1.D0/4.D0*(G(1)**2+G(2)**2)*(VUSQ-VDSQ)*(VUSQ+VDSQ)/VDSQ)
     $       /(1.D0-VUSQ/VDSQ))+1.D0/2.D0*(MDIFF+SIG1-SIG2)-SIG1
        BMHM2SQNEW=BMHM1SQNEW-MDIFF
        G(109)=(G(62)+G(61)+SIG1+SIG2)*VU*VD/(VUSQ+VDSQ)
        G(62)=BMHM1SQNEW
        G(61)=BMHM2SQNEW
        MUMSQ=((DBLE(G(351))+SIG2)*VUSQ/VDSQ-DBLE(G(352))-SIG1)
     $         /(1.D0-VUSQ/VDSQ)+1.D0/4.D0*(DBLE(G(291))**2
     $         +DBLE(G(292))**2)*(VUSQ+VDSQ)
        IF(MUMSQ.LT.0.D0)MUMSQ=MZ**2
        G(399)=(G(351)+G(352)+SIG1+SIG2+2.D0*G(398)**2)
     $                                              *VU*VD/(VUSQ+VDSQ)
        G(398)=DSQRT(MUMSQ)*EXP((0.D0,1.D0)*PHASEMU)
      END IF
      CALL ROTATE(1)
!
!Replace the original G entries, excluding the fixed couplings
!
      DO I=1,3
        G(I)=GSAV(I)
      END DO
      DO I=31,60
        G(I)=GSAV(I)
      END DO
      DO I=63,108
        G(I)=GSAV(I)
      END DO
      DO I=110,111
        G(I)=GSAV(I)
      END DO
      DO I=139,397
        G(I)=GSAV(I)
      END DO
      DO I=400,426
        G(I)=GSAV(I)
      END DO
      DO I=428,601
        G(I)=GSAV(I)
      END DO
!
      RETURN
      END
!
      SUBROUTINE DOWNSQM(G,Q,OUT)
!
!Purpose: Construct the 6x6 mass matrix for down-type squarks.
!         The basis is d_l,s_l,b_l,d_r,s_r,b_r.
!
      IMPLICIT NONE
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE COMPLEX G(601),OUT(6,6)
!
      DOUBLE COMPLEX AD(3,3),MQ(3,3),MD(3,3),MTSFD(3,3)
      DOUBLE COMPLEX FDQ(3,3),FDDR(3,3)
      DOUBLE COMPLEX GLP,G2L
      DOUBLE COMPLEX ADT(3,3),MTSFDT(3,3),FDQT(3,3),FDDRT(3,3)
      DOUBLE COMPLEX ADS(3,3),MTSFDS(3,3),FDQS(3,3),FDDRS(3,3)
      DOUBLE COMPLEX TRID(3,3),CFDQ(3,3),CFDDR(3,3)
      DOUBLE COMPLEX TRIDT(3,3),CFDQT(3,3),CFDDRT(3,3)
      DOUBLE COMPLEX TRIDS(3,3),CFDQS(3,3),CFDDRS(3,3)
      DOUBLE COMPLEX ID(3,3),CMATMUL
      DOUBLE PRECISION TANB,SINB,COSB,Q,VU,VD,VSM,VUSQ,VDSQ,VSMSQ
      INTEGER I,J
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $     ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $     ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $     ,ID(3,3)/(1.D0,0.D0)/
!
      DO I=1,6
        DO J=1,6
          OUT(I,J)=(0.D0,0.D0)
        END DO
      END DO
!
!Insert the soft mass matrices in the quark mass basis
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=MQQMASS(I,J)
          MD(I,J)=MDQMASS(I,J)
        END DO
      END DO
!
!Calculate cos (beta) and sin (beta)
!
      VU=DBLE(G(110))
      VD=DBLE(G(111))
      TANB=VU/VD
      COSB=DSQRT(1.D0/(1.D0+TANB**2))
      SINB=TANB*COSB
      VUSQ=VU**2
      VDSQ=VD**2
!
      VSM=DBLE(G(428))
      VSMSQ=VSM**2
!
!Convert the other entries in G(601) to the matrices needed
!
      DO I=1,3
        DO J=1,3
          AD(I,J)=G(42+(I-1)*3+J)
          MTSFD(I,J)=G(438+(I-1)*3+J)
          FDQ(I,J)=G(12+(I-1)*3+J) !Quartics are not run independently
          FDDR(I,J)=G(12+(I-1)*3+J)
!
          TRID(I,J)=G(408+(I-1)*3+J)
          CFDQ(I,J)=G(120+(I-1)*3+J)
          CFDDR(I,J)=G(120+(I-1)*3+J)
        END DO
      END DO
!
!Quartics are not run independently at this time.
!Set other quartics to be their non-tilde counterparts.
!This can be fixed if the quartic running is introduced.
!*****NB: IF THE QUARTIC RUNNING IS INTRODUCED CARE MUST
!         BE TAKEN WITH G(541) AND G(542) SINCE THEY ARE
!         DEFINED DIFFERENTLY IN THE COMPLEX VERSION AS
!         OPPOSED TO THE REAL VERSION. SEE THE NOTE AT
!         THE BEGINNING OF drge601.f
!
      GLP=DSQRT(3.D0/5.D0)*G(1)
      G2L=G(2)
      G(541)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*DSQRT(3.D0/5.D0)*G(1)
      G(542)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*G(2)
!
      DO I=1,3
        DO J=1,3
          ADT(I,J)=AD(J,I)
          MTSFDT(I,J)=MTSFD(J,I)
          FDQT(I,J)=FDQ(J,I)
          FDDRT(I,J)=FDDR(J,I)
          ADS(I,J)=CONJG(AD(I,J))
          MTSFDS(I,J)=CONJG(MTSFD(I,J))
          FDQS(I,J)=CONJG(FDQ(I,J))
          FDDRS(I,J)=CONJG(FDDR(I,J))
!
          TRIDT(I,J)=TRID(J,I)
          CFDQT(I,J)=CFDQ(J,I)
          CFDDRT(I,J)=CFDDR(J,I)
          TRIDS(I,J)=CONJG(TRID(I,J))
          CFDQS(I,J)=CONJG(CFDQ(I,J))
          CFDDRS(I,J)=CONJG(CFDDR(I,J))
        END DO
      END DO
!
!Split matrix into 3x3 blocks.
!First, the top left block
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J)=VDSQ*CMATMUL(0,FDDRS,FDDRT,I,J)+MQ(I,J)
     $               +(GLP**2/12.D0+G2L**2/4.D0)*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I,J)=VSMSQ*CMATMUL(0,CFDDRS,CFDDRT,I,J)+MQ(I,J)
     $               -VSMSQ*(G(541)**2/12.D0+G(542)**2/4.D0)*ID(I,J)
          END IF
        END DO
      END DO
!
!Next the bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J)=VU*MTSFDT(I,J)-VD*ADT(I,J)
          ELSE
            OUT(I+3,J)=-(VSM*TRIDT(I,J))
          END IF
        END DO
      END DO
!
!Top right is the dagger of bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J+3)=VU*MTSFDS(I,J)-VD*ADS(I,J)
          ELSE
            OUT(I,J+3)=-(VSM*TRIDS(I,J))
          END IF
        END DO
      END DO
!
!Finally bottom right
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J+3)=VDSQ*CMATMUL(0,FDQT,FDQS,I,J)+MD(I,J)
     $                   +GLP**2/6.D0*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I+3,J+3)=VSMSQ*CMATMUL(0,CFDQT,CFDQS,I,J)+MD(I,J)
     $                   -VSMSQ*G(541)**2/6.D0*ID(I,J)
          END IF
        END DO
      END DO
!
      RETURN
      END
!
      FUNCTION DRE(STAR,A,B)
!
!Purpose: To find the real part of A*B with:
!         If STAR=1 compute conjugate of A
!         If STAR=2 compute conjugate of B
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A,B
      DOUBLE PRECISION CON
      DOUBLE PRECISION DRE
      INTEGER STAR
!
      DRE=0.D0
      IF(STAR.EQ.1)THEN
        DRE=DBLE(CON(A)*B)
      ELSE IF(STAR.EQ.2)THEN
        DRE=DBLE(A*CON(B))
      ELSE IF(STAR.EQ.0)THEN
        DRE=DBLE(A*B)
      ELSE
        WRITE(*,*)'WRONG STAR'
      END IF
!
      RETURN
      END
!
      SUBROUTINE DRGE215(T,G,F)
!
!Contains threshold RGEs for gauge and yukawas.
!
!DRGE215 IS FOR RUNNNG THRESHOLDS, NO TILDES, FOR THE FIRST TIME
!The sfermion thresholds are not all distinct. All left-squarks
!are at the same point, all up-right-squarks are at another point,
!all down-right-squarks at another, all left-sleptons another and
!all right-sleptons another.
!
!     G(  1) = g_1         G(  2) = g_2         G(  3) = g_3
!     G(  4) = FU(1,1)     G(  5) = FU(1,2)     G( 12) = FU(3,3)
!     G( 13) = FD(1,1)     G( 22) = FE(1,1)     G( 30) = FE(3,3)
!     G( 31) = mu          G( 32) = V_U         G( 33) = V_D
!
!     G( 34) = LU(1,1)     G( 43) = LD(1,1)     G( 52) = LE(1,1)
!
!     G( 61) = VEV_SM      G( 62) = LAM_SM
!
!     G( 63) = GTPQ(1,1)   G( 72) = GTPL(1,1)   G( 81) = GTPU(1,1)
!     G( 90) = GTPD(1,1)   G( 99) = GTPE(1,1)   G(108) = GTP_Hu
!     G(109) = GTP_Hd      G(110) = GTQ(1,1)    G(119) = GTL(1,1)
!     G(128) = GT_Hu       G(129) = GT_Hd       G(130) = GTSQ(1,1)
!     G(139) = GTSU(1,1)   G(148) = GTSD(1,1)   G(157) = FTUQ(1,1)
!     G(166) = FTDQ(1,1)   G(175) = FTEL(1,1)   G(184) = FTUU(1,1)
!     G(193) = FTDD(1,1)   G(202) = FTEE(1,1)   G(211) = sGTP_Hu
!     G(212) = cGTP_Hd     G(213) = sGT_Hu      G(214) = cGT_Hd
!
!     G(215) = mu(M)
!
!This is the BT version which receives G in book notation
!
      IMPLICIT NONE
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE PRECISION T
      DOUBLE PRECISION G(215)
      DOUBLE PRECISION F(215)
!
      DOUBLE PRECISION FU(3,3),FD(3,3),FE(3,3)
      DOUBLE PRECISION YU(3,3),YD(3,3),YE(3,3)
      DOUBLE PRECISION AU(3,3),AD(3,3),AE(3,3)
      DOUBLE PRECISION HU(3,3),HD(3,3),HE(3,3),DH(3,3,3)
      DOUBLE PRECISION MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE PRECISION DM(5,3,3)
      DOUBLE PRECISION LU(3,3),LD(3,3),LE(3,3),LYU(3,3),LYD(3,3)
      DOUBLE PRECISION LYE(3,3)
      DOUBLE PRECISION GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3),GTPE(3,3)
      DOUBLE PRECISION GTQ(3,3),GTL(3,3),GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE PRECISION FTUQ(3,3),FTDQ(3,3),FTEL(3,3)
      DOUBLE PRECISION FTUU(3,3),FTDD(3,3),FTEE(3,3)
!
      DOUBLE PRECISION DMATMUL,SFMUL,DTRACE,SUM
!
!These are used in the calculation of the RGEs which contain
!thresholds
!
      DOUBLE PRECISION DUMU1(3,3),DUMU2(3,3),DUMD1(3,3),DUMD2(3,3)
      DOUBLE PRECISION DUME1(3,3),DUME2(3,3),DUMGRKMU1(3,3),TDUMGRKMU
      DOUBLE PRECISION DUMLUD1(3,3),DUMLUD2(3,3),DUMLUD(3,3),TDUMLUD
!
      DOUBLE PRECISION FUS(3,3),FDS(3,3),FES(3,3)
      DOUBLE PRECISION FTUQD(3,3),FTUUD(3,3),FTDDD(3,3),FTDQD(3,3)
      DOUBLE PRECISION FTELD(3,3),FTEED(3,3)
      DOUBLE PRECISION GTPUS(3,3),GTSUS(3,3),GTQS(3,3),GTLS(3,3)
      DOUBLE PRECISION GTPQS(3,3),GTPLS(3,3),GTSQS(3,3),GTPDS(3,3)
      DOUBLE PRECISION GTPES(3,3),GTSDS(3,3),GTPUT(3,3),GTQT(3,3)
      DOUBLE PRECISION GTLT(3,3),GTPQT(3,3),GTPLT(3,3),GTPDT(3,3)
      DOUBLE PRECISION GTPET(3,3),GTSUT(3,3),GTSQT(3,3),GTSDT(3,3)
!
      DOUBLE PRECISION FUFUD(3,3),FDFDD(3,3),FEFED(3,3)
      DOUBLE PRECISION YUDYU(3,3),YDDYD(3,3),YEDYE(3,3)
      DOUBLE PRECISION LULUD(3,3),LDLDD(3,3),LELED(3,3)
      DOUBLE PRECISION LYUDLYU(3,3),LYDDLYD(3,3),LYEDLYE(3,3)
      DOUBLE PRECISION LYUDLYU2(3,3),LYDDLYD2(3,3),LYEDLYE2(3,3)
      DOUBLE PRECISION LYUDLYU3(3,3),LYDDLYD3(3,3),LYEDLYE3(3,3)
      DOUBLE PRECISION TLYUDLYU,TLYDDLYD,TLYEDLYE
      DOUBLE PRECISION TLYUDLYU2,TLYDDLYD2,TLYEDLYE2
      DOUBLE PRECISION TLYUDLYU3,TLYDDLYD3,TLYEDLYE3
      DOUBLE PRECISION LYUDLYULYDDLYD(3,3),LYDDLYDLYUDLYU(3,3)
      DOUBLE PRECISION TLYUDLYULYDDLYD
      DOUBLE PRECISION FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
      DOUBLE PRECISION YUYUD(3,3),YDYDD(3,3),YEYED(3,3)
      DOUBLE PRECISION YUYDD(3,3),YDYUD(3,3)
      DOUBLE PRECISION TFUDFU,TFDDFD,TFEDFE
      DOUBLE PRECISION TYUDYU,TYDDYD,TYEDYE
      DOUBLE PRECISION TLUDLU,TLDDLD,TLEDLE
      DOUBLE PRECISION FUFUDFU(3,3),FUFUDFD(3,3),FDFDDFU(3,3)
      DOUBLE PRECISION FDFDDFD(3,3),FEFEDFE(3,3)
      DOUBLE PRECISION LULUDLU(3,3),LULUDLD(3,3),LDLDDLU(3,3)
      DOUBLE PRECISION LDLDDLD(3,3),LELEDLE(3,3)
      DOUBLE PRECISION SQFTUQDFTUQ(3,3),SQFTDQDFTDQ(3,3),SQGTQTGTQS(3,3)
      DOUBLE PRECISION SQGTQTFTUQ(3,3),SQGTQTFTDQ(3,3),SQGTPQTGTPQS(3,3)
      DOUBLE PRECISION SQGTPQTFTUQ(3,3),SQGTPQTFTDQ(3,3)
      DOUBLE PRECISION SQGTSQTGTSQS(3,3),SUGTPUSGTPUT(3,3)
      DOUBLE PRECISION SUGTSUSGTSUT(3,3),SUFTUUFTUUD(3,3)
      DOUBLE PRECISION SUFTUUGTPUT(3,3),SDFTDDFTDDD(3,3)
      DOUBLE PRECISION SDFTDDGTPDT(3,3),SDGTPDSGTPDT(3,3)
      DOUBLE PRECISION SDGTSDSGTSDT(3,3),SLFTELDFTEL(3,3)
      DOUBLE PRECISION SLGTLTGTLS(3,3),SLGTLTFTEL(3,3)
      DOUBLE PRECISION SLGTPLTGTPLS(3,3),SLGTPLTFTEL(3,3)
      DOUBLE PRECISION SEGTPESGTPET(3,3),SEFTEEFTEED(3,3)
      DOUBLE PRECISION SEFTEEGTPET(3,3)
      DOUBLE PRECISION MGTPHUSQ,MGTPHDSQ,MGTHUSQ,MGTHDSQ
      DOUBLE PRECISION MSGTPHUSQ,MCGTPHDSQ,MSGTHUSQ,MCGTHDSQ
      DOUBLE PRECISION Y2,H,Y4,CHI4
!
      DOUBLE PRECISION B1U(3,3),B1D(3,3),B1E(3,3)
      DOUBLE PRECISION BGRKMU,BETALAM1,BETALAM2

!
!The following are for two loop and soft RGEs
!
      DOUBLE PRECISION DUM2U1(3,3),DUM2U2(3,3)
      DOUBLE PRECISION DUM2D1(3,3),DUM2D2(3,3)
      DOUBLE PRECISION DUM2E1(3,3),DUM2E2(3,3)
      DOUBLE PRECISION DUM2U(3,3),DUM2D(3,3),DUM2E(3,3)
      DOUBLE PRECISION DUM1GRKMU(3,3),DUM2GRKMU(3,3)
!
      DOUBLE PRECISION YUYUDYUYUD(3,3),YDYDDYDYDD(3,3)
      DOUBLE PRECISION YEYEDYEYED(3,3),YUDYUYUDYU(3,3)
      DOUBLE PRECISION YDDYDYDDYD(3,3),YEDYEYEDYE(3,3)
      DOUBLE PRECISION YDDYDYUDYU(3,3),YUDYUYDDYD(3,3)
      DOUBLE PRECISION YUYDDYDYUD(3,3)
!
      DOUBLE PRECISION TYUDYUYUDYU,TYDDYDYDDYD,TYEDYEYEDYE
      DOUBLE PRECISION TYUDYUYDDYD,TYDDYDYUDYU
!
      DOUBLE PRECISION BETA1U(3,3),BETA2U(3,3),BETA1D(3,3),BETA2D(3,3)
      DOUBLE PRECISION BETA1E(3,3),BETA2E(3,3)
      DOUBLE PRECISION BETA2GRKMU,BETA1VU,BETA1VD,BETA2VU,BETA2VD
      DOUBLE PRECISION B2YMU(3,3),B2YMD(3,3),B2YME(3,3)
      DOUBLE PRECISION BETAVEV1,BETAVEV2
!
      DOUBLE PRECISION ID(3,3),MVMU,MVMUM
      DOUBLE PRECISION B1LP(3),B1LPM(3),B2LPSM(3,3),B2LPM(3,3)
      DOUBLE PRECISION CM(3,3),CSM(3,3)
      DOUBLE PRECISION PI,Q
!
      INTEGER I,J,NG,ND,NE,NNU,NU,NSQ,NSU,NSD,NSL,NSE,NSH,NH,NSW,NSG
      INTEGER THLH,THHH,THSH,THSQ(3),THSU(3),THSD(3),THSL(3)
      INTEGER THSE(3),THSB,THSW,THGL
!
      DATA ID(1,1)/1.D0/,ID(1,2)/0.D0/,ID(1,3)/0.D0/
      DATA ID(2,1)/0.D0/,ID(2,2)/1.D0/,ID(2,3)/0.D0/
      DATA ID(3,1)/0.D0/,ID(3,2)/0.D0/,ID(3,3)/1.D0/
      DATA B1LPM(1)/6.6D0/,B1LPM(2)/1.D0/,B1LPM(3)/-3.D0/
      DATA B2LPM(1,1)/7.96D0/,B2LPM(1,2)/5.4D0/,B2LPM(1,3)/17.6D0/
      DATA B2LPM(2,1)/1.8D0/,B2LPM(2,2)/25.D0/,B2LPM(2,3)/24.D0/
      DATA B2LPM(3,1)/2.2D0/,B2LPM(3,2)/9.D0/,B2LPM(3,3)/14.D0/
      DATA CM(1,1)/5.2D0/,CM(1,2)/2.8D0/,CM(1,3)/3.6D0/
      DATA CM(2,1)/6.D0/,CM(2,2)/6.D0/,CM(2,3)/2.D0/
      DATA CM(3,1)/4.D0/,CM(3,2)/4.D0/,CM(3,3)/0.D0/
!
!Set all F's and betas to zero
!
      DO I=1,215
        F(I)=0.D0
      END DO
      DO I=1,3
        DO J=1,3
          B1U(I,J)=0.D0
          B1D(I,J)=0.D0
          B1E(I,J)=0.D0
          BETA1U(I,J)=0.D0
          BETA2U(I,J)=0.D0
          BETA1D(I,J)=0.D0
          BETA2D(I,J)=0.D0
          BETA1E(I,J)=0.D0
          BETA2E(I,J)=0.D0
          B2YMU(I,J)=0.D0
          B2YMD(I,J)=0.D0
          B2YME(I,J)=0.D0
        END DO
      END DO
      BGRKMU=0.D0
      BETA2GRKMU=0.D0
      BETA1VU=0.D0
      BETA1VD=0.D0
      BETA2VU=0.D0
      BETA2VD=0.D0
      BETALAM1=0.D0
      BETALAM2=0.D0
      BETAVEV1=0.D0
      BETAVEV2=0.D0
!
      PI=4.D0*ATAN(1.D0)
      Q=SSQSTEP
      IF(Q.LT.1.D0)THEN
        WRITE(*,*)'ERROR IN Q: ',Q
        STOP
      END IF
      NG=3.D0
      NU=3
      ND=3
      NE=3
      NNU=3
!
      NSQ=3
      IF((Q-QTHQL(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHQL(1)).LT.ABS(EPS).AND.EPS.LT.0))NSQ=0
      NSU=3
      IF((Q-QTHUR(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHUR(1)).LT.ABS(EPS).AND.EPS.LT.0))NSU=0
      NSD=3
      IF((Q-QTHDR(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHDR(1)).LT.ABS(EPS).AND.EPS.LT.0))NSD=0
      NSL=3
      IF((Q-QTHLL(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHLL(1)).LT.ABS(EPS).AND.EPS.LT.0))NSL=0
      NSE=3
      IF((Q-QTHER(1)).LT.-EPS.OR.
     $       (ABS(Q-QTHER(1)).LT.ABS(EPS).AND.EPS.LT.0))NSE=0
!
      IF ((Q-QNSH).GT.EPS.OR.
     $         (ABS(Q-QNSH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSH=2
      ELSE
        NSH=0
      END IF
      IF ((Q-QNSG).GT.EPS.OR.
     $         (ABS(Q-QNSG).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF ((Q-QNH).GT.EPS.OR.
     $         (ABS(Q-QNH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF ((Q-QTHSB).GT.EPS.OR.
     $         (ABS(Q-QTHSB).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        THSB=1
      ELSE
        THSB=0
      END IF
      IF ((Q-QTHSW).GT.EPS.OR.
     $         (ABS(Q-QTHSW).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSW=1
        THSW=1
      ELSE
        NSW=0
        THSW=0
      END IF
      THLH=1
      THHH=NH/2 !This works so long as THHH is an integer variable
      IF(THHH.NE.1)THEN !Perform check
        IF(THHH.NE.0)WRITE(*,*)'ERROR IN THHH'
      END IF
      THSH=NSH/2
      DO I=1,3
        THSQ(I)=1
        IF((Q-QTHQL(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHQL(1)).LT.ABS(EPS).AND.EPS.LT.0))THSQ(I)=0
        THSU(I)=1
        IF((Q-QTHUR(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHUR(1)).LT.ABS(EPS).AND.EPS.LT.0))THSU(I)=0
        THSD(I)=1
        IF((Q-QTHDR(1)).LT.-EPS.OR.
     $         (ABS(Q-QTHDR(1)).LT.ABS(EPS).AND.EPS.LT.0))THSD(I)=0
        THSL(I)=NSL/3
        THSE(I)=NSE/3
      END DO
      THGL=NSG
!
!Convert input into 3x3 matrices.
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
!
!MV notation Yukawas
!
          YU(I,J)=G(3+(J-1)*3+I)
          YD(I,J)=G(12+(J-1)*3+I)
          YE(I,J)=G(21+(J-1)*3+I)
!
          LU(I,J)=G(33+(I-1)*3+J)
          LD(I,J)=G(42+(I-1)*3+J)
          LE(I,J)=G(51+(I-1)*3+J)
          GTPQ(I,J)=G(62+(I-1)*3+J)
          GTPL(I,J)=G(71+(I-1)*3+J)
          GTPU(I,J)=G(80+(I-1)*3+J)
          GTPD(I,J)=G(89+(I-1)*3+J)
          GTPE(I,J)=G(98+(I-1)*3+J)
          GTQ(I,J)=G(109+(I-1)*3+J)
          GTL(I,J)=G(118+(I-1)*3+J)
          GTSQ(I,J)=G(129+(I-1)*3+J)
          GTSU(I,J)=G(138+(I-1)*3+J)
          GTSD(I,J)=G(147+(I-1)*3+J)
          FTUQ(I,J)=G(156+(I-1)*3+J)
          FTDQ(I,J)=G(165+(I-1)*3+J)
          FTEL(I,J)=G(174+(I-1)*3+J)
          FTUU(I,J)=G(183+(I-1)*3+J)
          FTDD(I,J)=G(192+(I-1)*3+J)
          FTEE(I,J)=G(201+(I-1)*3+J)
        END DO
      END DO
!
      MVMU=G(31)
      MVMUM=G(215)
!
!SET THE TILDE TERMS TO THEIR SM COUNTERPARTS.
!
      DO I=1,3
        DO J=1,3
          GTPQ(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(62+(I-1)*3+J)=GTPQ(I,J)
          GTPL(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(71+(I-1)*3+J)=GTPL(I,J)
          GTPU(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(80+(I-1)*3+J)=GTPU(I,J)
          GTPD(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(89+(I-1)*3+J)=GTPD(I,J)
          GTPE(I,J)=G(1)*SQRT(3.D0/5.D0)*ID(I,J)
          G(98+(I-1)*3+J)=GTPE(I,J)
          GTQ(I,J)=G(2)*ID(I,J)
          G(109+(I-1)*3+J)=GTQ(I,J)
          GTL(I,J)=G(2)*ID(I,J)
          G(118+(I-1)*3+J)=GTL(I,J)
          GTSQ(I,J)=G(3)*ID(I,J)
          G(129+(I-1)*3+J)=GTSQ(I,J)
          GTSU(I,J)=G(3)*ID(I,J)
          G(138+(I-1)*3+J)=GTSU(I,J)
          GTSD(I,J)=G(3)*ID(I,J)
          G(147+(I-1)*3+J)=GTSD(I,J)
          FTUQ(I,J)=YU(J,I)
          G(156+(I-1)*3+J)=FTUQ(I,J)
          FTDQ(I,J)=YD(J,I)
          G(165+(I-1)*3+J)=FTDQ(I,J)
          FTEL(I,J)=YE(J,I)
          G(174+(I-1)*3+J)=FTEL(I,J)
          FTUU(I,J)=YU(J,I)
          G(183+(I-1)*3+J)=FTUU(I,J)
          FTDD(I,J)=YD(J,I)
          G(192+(I-1)*3+J)=FTDD(I,J)
          FTEE(I,J)=YE(J,I)
          G(201+(I-1)*3+J)=FTEE(I,J)
        END DO
      END DO
      G(108)=SQRT(3.D0/5.D0)*G(1)
      G(109)=SQRT(3.D0/5.D0)*G(1)
      G(128)=G(2)
      G(129)=G(2)
      G(211)=0.D0
      G(212)=0.D0
      G(213)=0.D0
      G(214)=0.D0
!
!The separated out contributions are from PRD 49 4882 (1194),
!Castano,Piard,Ramond
!
      B1LP(1)=2.D0/5.D0*(17.D0/12.D0*DBLE(NU)+5.D0/12.D0*DBLE(ND)
     $        +5.D0/4.D0*DBLE(NE)+1.D0/4.D0*DBLE(NNU))+1.D0/30.D0
     $        *DBLE(NSQ)+4.D0/15.D0*DBLE(NSU)+1.D0/15.D0*DBLE(NSD)
     $        +1.D0/10.D0*DBLE(NSL)+1.D0/5.D0*DBLE(NSE)
     $        +1.D0/5.D0*DBLE(NSH)+1.D0/10.D0*DBLE(NH)
      B1LP(2)=-22.D0/3.D0+1.D0/2.D0*(DBLE(NU)+DBLE(ND))+1.D0/6.D0
     $        *(DBLE(NE)+DBLE(NNU))+1.D0/2.D0*DBLE(NSQ)
     $        +1.D0/6.D0*DBLE(NSL)+1.D0/3.D0*DBLE(NSH)+1.D0/6.D0
     $        *DBLE(NH)+4.D0/3.D0*DBLE(NSW)
      B1LP(3)=-11.D0+2.D0/3.D0*(DBLE(NU)+DBLE(ND))+1.D0/3.D0
     $        *DBLE(NSQ)+1.D0/6.D0*DBLE(NSU)+1.D0/6.D0*DBLE(NSD)
     $        +2.D0*DBLE(NSG)
      IF(THHH.EQ.0)THEN
        B2LPSM(1,1)=-(-NG*19.D0/15.D0-9.D0/50.D0)
        B2LPSM(1,2)=-(-NG*3.D0/5.D0-9.D0/10.D0)
        B2LPSM(1,3)=-(-NG*44.D0/15.D0)
        B2LPSM(2,1)=-(-NG*1.D0/5.D0-3.D0/10.D0)
        B2LPSM(2,2)=-(136.D0/3.D0-NG*49.D0/3.D0-13.D0/6.D0)
        B2LPSM(2,3)=-(-NG*4.D0)
        B2LPSM(3,1)=-(-NG*11.D0/30.D0)
        B2LPSM(3,2)=-(-NG*3.D0/2.D0)
        B2LPSM(3,3)=-(102.D0-NG*76.D0/3.D0)
        CSM(1,1)=1.7D0
        CSM(1,2)=.5D0
        CSM(1,3)=1.5D0
        CSM(2,1)=1.5D0
        CSM(2,2)=1.5D0
        CSM(2,3)=.5D0
        CSM(3,1)=2.D0
        CSM(3,2)=2.D0
        CSM(3,3)=0.D0
      END IF
!
!I need many variations on the 3x3 matrices.
!
      CALL DAGGER(FTUQ,FTUQD)
      CALL DAGGER(FTUU,FTUUD)
      CALL DAGGER(FTDD,FTDDD)
      CALL DAGGER(FTDQ,FTDQD)
      CALL DAGGER(FTEL,FTELD)
      CALL DAGGER(FTEE,FTEED)
      DO I=1,3
        DO J=1,3
          GTPUS(I,J)=GTPU(I,J)
          GTSUS(I,J)=GTSU(I,J)
          GTQS(I,J)=GTQ(I,J)
          GTLS(I,J)=GTL(I,J)
          GTPQS(I,J)=GTPQ(I,J)
          GTPLS(I,J)=GTPL(I,J)
          GTSQS(I,J)=GTSQ(I,J)
          GTPDS(I,J)=GTPD(I,J)
          GTPES(I,J)=GTPE(I,J)
          GTSDS(I,J)=GTSD(I,J)
          GTPUT(I,J)=GTPU(J,I)
          GTQT(I,J)=GTQ(J,I)
          GTLT(I,J)=GTL(J,I)
          GTPQT(I,J)=GTPQ(J,I)
          GTPLT(I,J)=GTPL(J,I)
          GTPDT(I,J)=GTPD(J,I)
          GTPET(I,J)=GTPE(J,I)
          GTSUT(I,J)=GTSU(J,I)
          GTSQT(I,J)=GTSQ(J,I)
          GTSDT(I,J)=GTSD(J,I)
        END DO
      END DO
!
!Now all the matrix multiples
!
      DO I=1,3
        DO J=1,3
          YUDYU(I,J)=DMATMUL(1,YU,YU,I,J)
          YDDYD(I,J)=DMATMUL(1,YD,YD,I,J)
          YEDYE(I,J)=DMATMUL(1,YE,YE,I,J)
          YUYUD(I,J)=DMATMUL(2,YU,YU,I,J)
          YDYDD(I,J)=DMATMUL(2,YD,YD,I,J)
          YEYED(I,J)=DMATMUL(2,YE,YE,I,J)
          YUYDD(I,J)=DMATMUL(2,YU,YD,I,J)
          YDYUD(I,J)=DMATMUL(2,YD,YU,I,J)
!
          FUFUD(I,J)=DMATMUL(2,FU,FU,I,J)
          FDFDD(I,J)=DMATMUL(2,FD,FD,I,J)
          FEFED(I,J)=DMATMUL(2,FE,FE,I,J)
          FUDFU(I,J)=DMATMUL(1,FU,FU,I,J)
          FDDFD(I,J)=DMATMUL(1,FD,FD,I,J)
          FEDFE(I,J)=DMATMUL(1,FE,FE,I,J)
          LULUD(I,J)=DMATMUL(2,LU,LU,I,J)
          LDLDD(I,J)=DMATMUL(2,LD,LD,I,J)
          LELED(I,J)=DMATMUL(2,LE,LE,I,J)
        END DO
      END DO
!
      IF(THHH.EQ.0)THEN
        MSGTPHUSQ=G(211)**2
        MCGTPHDSQ=G(212)**2
        MSGTHUSQ=G(213)**2
        MCGTHDSQ=G(214)**2
      END IF
      MGTPHUSQ=G(108)**2
      MGTPHDSQ=G(109)**2
      MGTHUSQ=G(128)**2
      MGTHDSQ=G(129)**2
!
      TYUDYU=DTRACE(YUDYU)
      TYDDYD=DTRACE(YDDYD)
      TYEDYE=DTRACE(YEDYE)
!
      TFUDFU=DTRACE(FUDFU)
      TFDDFD=DTRACE(FDDFD)
      TFEDFE=DTRACE(FEDFE)
      TLUDLU=DTRACE(LULUD)
      TLDDLD=DTRACE(LDLDD)
      TLEDLE=DTRACE(LELED)
!
      DO I=1,3
        DO J=1,3
          FUFUDFU(I,J)=DMATMUL(0,FUFUD,FU,I,J)
          FUFUDFD(I,J)=DMATMUL(0,FUFUD,FD,I,J)
          FDFDDFU(I,J)=DMATMUL(0,FDFDD,FU,I,J)
          FDFDDFD(I,J)=DMATMUL(0,FDFDD,FD,I,J)
          FEFEDFE(I,J)=DMATMUL(0,FEFED,FE,I,J)
          LULUDLU(I,J)=DMATMUL(0,LULUD,LU,I,J)
          LULUDLD(I,J)=DMATMUL(0,LULUD,LD,I,J)
          LDLDDLU(I,J)=DMATMUL(0,LDLDD,LU,I,J)
          LDLDDLD(I,J)=DMATMUL(0,LDLDD,LD,I,J)
          LELEDLE(I,J)=DMATMUL(0,LELED,LE,I,J)
          SQFTUQDFTUQ(I,J)=SFMUL(THSQ,FTUQD,FTUQ,I,J)
          SQFTDQDFTDQ(I,J)=SFMUL(THSQ,FTDQD,FTDQ,I,J)
          SQGTQTGTQS(I,J)=SFMUL(THSQ,GTQT,GTQS,I,J)
          SQGTQTFTUQ(I,J)=SFMUL(THSQ,GTQT,FTUQ,I,J)
          SQGTQTFTDQ(I,J)=SFMUL(THSQ,GTQT,FTDQ,I,J)
          SQGTPQTGTPQS(I,J)=SFMUL(THSQ,GTPQT,GTPQS,I,J)
          SQGTPQTFTUQ(I,J)=SFMUL(THSQ,GTPQT,FTUQ,I,J)
          SQGTPQTFTDQ(I,J)=SFMUL(THSQ,GTPQT,FTDQ,I,J)
          SQGTSQTGTSQS(I,J)=SFMUL(THSQ,GTSQT,GTSQS,I,J)
          SUFTUUFTUUD(I,J)=SFMUL(THSU,FTUU,FTUUD,I,J)
          SUFTUUGTPUT(I,J)=SFMUL(THSU,FTUU,GTPUT,I,J)
          SUGTPUSGTPUT(I,J)=SFMUL(THSU,GTPUS,GTPUT,I,J)
          SUGTSUSGTSUT(I,J)=SFMUL(THSU,GTSUS,GTSUT,I,J)
          SDFTDDFTDDD(I,J)=SFMUL(THSD,FTDD,FTDDD,I,J)
          SDFTDDGTPDT(I,J)=SFMUL(THSD,FTDD,GTPDT,I,J)
          SDGTPDSGTPDT(I,J)=SFMUL(THSD,GTPDS,GTPDT,I,J)
          SDGTSDSGTSDT(I,J)=SFMUL(THSD,GTSDS,GTSDT,I,J)
          SLFTELDFTEL(I,J)=SFMUL(THSL,FTELD,FTEL,I,J)
          SLGTLTGTLS(I,J)=SFMUL(THSL,GTLT,GTLS,I,J)
          SLGTLTFTEL(I,J)=SFMUL(THSL,GTLT,FTEL,I,J)
          SLGTPLTGTPLS(I,J)=SFMUL(THSL,GTPLT,GTPLS,I,J)
          SLGTPLTFTEL(I,J)=SFMUL(THSL,GTPLT,FTEL,I,J)
          SEFTEEFTEED(I,J)=SFMUL(THSE,FTEE,FTEED,I,J)
          SEFTEEGTPET(I,J)=SFMUL(THSE,FTEE,GTPET,I,J)
          SEGTPESGTPET(I,J)=SFMUL(THSE,GTPES,GTPET,I,J)
        END DO
      END DO
!
!These are the two loop terms. All in MV notation
!
      IF(SW2LP.EQ.1)THEN
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=DMATMUL(0,YUYUD,YUYUD,I,J)
            YDYDDYDYDD(I,J)=DMATMUL(0,YDYDD,YDYDD,I,J)
            YEYEDYEYED(I,J)=DMATMUL(0,YEYED,YEYED,I,J)
            YUDYUYUDYU(I,J)=DMATMUL(0,YUDYU,YUDYU,I,J)
            YDDYDYDDYD(I,J)=DMATMUL(0,YDDYD,YDDYD,I,J)
            YEDYEYEDYE(I,J)=DMATMUL(0,YEDYE,YEDYE,I,J)
            YDDYDYUDYU(I,J)=DMATMUL(0,YDDYD,YUDYU,I,J)
            YUDYUYDDYD(I,J)=DMATMUL(0,YUDYU,YDDYD,I,J)
            YUYDDYDYUD(I,J)=DMATMUL(0,YUYDD,YDYUD,I,J)
          END DO
        END DO
!
        TYUDYUYUDYU=DTRACE(YUDYUYUDYU)
        TYUDYUYDDYD=DTRACE(YUDYUYDDYD)
        TYDDYDYUDYU=DTRACE(YDDYDYUDYU)
        TYDDYDYDDYD=DTRACE(YDDYDYDDYD)
        TYEDYEYEDYE=DTRACE(YEDYEYEDYE)
!
!These are SM terms for the two loop running below m_H
!I am going to use LYU for the MV notation SM Yukawa
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=DMATMUL(1,LYU,LYU,I,J)
              LYDDLYD(I,J)=DMATMUL(1,LYD,LYD,I,J)
              LYEDLYE(I,J)=DMATMUL(1,LYE,LYE,I,J)
            END DO
          END DO
          TLYUDLYU=DTRACE(LYUDLYU)
          TLYDDLYD=DTRACE(LYDDLYD)
          TLYEDLYE=DTRACE(LYEDLYE)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=DMATMUL(0,LYUDLYU,LYUDLYU,I,J)
              LYDDLYD2(I,J)=DMATMUL(0,LYDDLYD,LYDDLYD,I,J)
              LYEDLYE2(I,J)=DMATMUL(0,LYEDLYE,LYEDLYE,I,J)
              LYUDLYULYDDLYD(I,J)=DMATMUL(0,LYUDLYU,LYDDLYD,I,J)
              LYDDLYDLYUDLYU(I,J)=DMATMUL(0,LYDDLYD,LYUDLYU,I,J)
              DUMLUD1(I,J)=LYUDLYU(I,J)+LYDDLYD(I,J)
            END DO
          END DO
          TLYUDLYU2=DTRACE(LYUDLYU2)
          TLYDDLYD2=DTRACE(LYDDLYD2)
          TLYEDLYE2=DTRACE(LYEDLYE2)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=DMATMUL(0,DUMLUD1,LYDDLYD,I,J)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=DMATMUL(0,LYUDLYU2,LYUDLYU,I,J)
              LYDDLYD3(I,J)=DMATMUL(0,LYDDLYD2,LYDDLYD,I,J)
              LYEDLYE3(I,J)=DMATMUL(0,LYEDLYE2,LYEDLYE,I,J)
              DUMLUD(I,J)=DMATMUL(0,LYUDLYU,DUMLUD2,I,J)
            END DO
          END DO
          TLYUDLYU3=DTRACE(LYUDLYU3)
          TLYDDLYD3=DTRACE(LYDDLYD3)
          TLYEDLYE3=DTRACE(LYEDLYE3)
          TLYUDLYULYDDLYD=DTRACE(LYUDLYULYDDLYD)
          TDUMLUD=DTRACE(DUMLUD)
!
          Y2=3.D0*TLYUDLYU+3.D0*TLYDDLYD+TLYEDLYE
          H=3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
          Y4=(83.D0/40.D0*G(1)**2+27.D0/8.D0*G(2)**2
     $       +28.D0*G(3)**2)*TLYUDLYU+(-1.D0/40.D0*G(1)**2
     $       +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TLYDDLYD
     $       +(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TLYEDLYE
          CHI4=9.D0/4.D0*(3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
     $         -2.D0/3.D0*TLYUDLYULYDDLYD)
        END IF
      ELSE
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=0.D0
            YDYDDYDYDD(I,J)=0.D0
            YEYEDYEYED(I,J)=0.D0
            YUDYUYUDYU(I,J)=0.D0
            YDDYDYDDYD(I,J)=0.D0
            YEDYEYEDYE(I,J)=0.D0
            YDDYDYUDYU(I,J)=0.D0
            YUDYUYDDYD(I,J)=0.D0
            YUYDDYDYUD(I,J)=0.D0
          END DO
        END DO
        TYUDYUYUDYU=0.D0
        TYUDYUYDDYD=0.D0
        TYDDYDYUDYU=0.D0
        TYDDYDYDDYD=0.D0
        TYEDYEYEDYE=0.D0
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=0.D0
              LYDDLYD(I,J)=0.D0
              LYEDLYE(I,J)=0.D0
            END DO
          END DO
          TLYUDLYU=0.D0
          TLYDDLYD=0.D0
          TLYEDLYE=0.D0
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=0.D0
              LYDDLYD2(I,J)=0.D0
              LYEDLYE2(I,J)=0.D0
              LYUDLYULYDDLYD(I,J)=0.D0
              LYDDLYDLYUDLYU(I,J)=0.D0
              DUMLUD1(I,J)=0.D0
            END DO
          END DO
          TLYUDLYU2=0.D0
          TLYDDLYD2=0.D0
          TLYEDLYE2=0.D0
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=0.D0
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=0.D0
              LYDDLYD3(I,J)=0.D0
              LYEDLYE3(I,J)=0.D0
              DUMLUD(I,J)=0.D0
            END DO
          END DO
          TLYUDLYU3=0.D0
          TLYDDLYD3=0.D0
          TLYEDLYE3=0.D0
          TLYUDLYULYDDLYD=0.D0
          TDUMLUD=0.D0
!
          Y2=0.D0
          H=0.D0
          Y4=0.D0
          CHI4=0.D0
        END IF
      END IF
!
!The threshold gauge running with a change at m_H to SM
!
      DO I=1,3
        SUM=0.D0
        DO J=1,3
          IF(THHH.EQ.0)THEN
            SUM=SUM+B2LPSM(I,J)*G(J)**2
          ELSE
            SUM=SUM+B2LPM(I,J)*G(J)**2
          END IF
        END DO
        IF(THHH.EQ.0)THEN
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CSM(I,1)*TLYUDLYU+CSM(I,2)*TLYDDLYD
     $         +CSM(I,3)*TLYEDLYE)))
        ELSE
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CM(I,1)*TYUDYU+CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        END IF
      END DO
!
!Next the full Yukawas
!
      DO I=1,3
        DO J=1,3
          DUMU1(I,J)=THSH*SQFTUQDFTUQ(I,J)
     $               +4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $               +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
          DUMU2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMD1(I,J)=THSH*SQFTDQDFTDQ(I,J)
     $               +1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $               +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
          DUMD2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUME1(I,J)=THSH*SLFTELDFTEL(I,J)+THSB*SEGTPESGTPET(I,J)
          DUME2(I,J)=2.D0*THSH*SEFTEEFTEED(I,J)
     $               +3.D0*THSW*SLGTLTGTLS(I,J)+THSB*SLGTPLTGTPLS(I,J)
!
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U2(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2D1(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)+
     $                  YEYEDYEYED(I,J)
            DUM2D2(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E1(I,J)=DUM2D1(I,J)
            DUM2E2(I,J)=DUM2U2(I,J)
          ELSE
            DUM2U1(I,J)=0.D0
            DUM2U2(I,J)=0.D0
            DUM2D1(I,J)=0.D0
            DUM2D2(I,J)=0.D0
            DUM2E1(I,J)=0.D0
            DUM2E2(I,J)=0.D0
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Here are the two loop terms
!
          IF(SW2LP.EQ.1)THEN
            DUM2U(I,J)=-3.D0*DTRACE(DUM2U1)*ID(I,J)-YDDYD(I,J)
     $                 *DTRACE(DUM2U2)-9.D0*YUDYU(I,J)*TYUDYU
     $                 -4.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)-2.D0
     $                 *YDDYDYUDYU(I,J)+(16.D0*G(3)**2+4.D0/5.D0
     $                 *G(1)**2)*TYUDYU*ID(I,J)+(6.D0*G(2)**2
     $                 +2.D0/5.D0*G(1)**2)*YUDYU(I,J)+2.D0/5.D0
     $                 *G(1)**2*YDDYD(I,J)+(-16.D0/9.D0*G(3)**4
     $                 +8.D0*G(3)**2*G(2)**2+136.D0/45.D0
     $                 *G(3)**2*G(1)**2+15.D0/2.D0*G(2)**4
     $                 +G(2)**2*G(1)**2+2743.D0/450.D0
     $                 *G(1)**4)*ID(I,J)
            DUM2D(I,J)=-3.D0*DTRACE(DUM2D1)*ID(I,J)-3.D0*YUDYU(I,J)
     $                 *TYUDYU-3.D0*YDDYD(I,J)*DTRACE(DUM2D2)
     $                 -4.D0*YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)
     $                 -2.D0*YUDYUYDDYD(I,J)+(16.D0*G(3)**2
     $                 -2.D0/5.D0*G(1)**2)*TYDDYD*ID(I,J)+6.D0/5.D0
     $                 *G(1)**2*TYEDYE*ID(I,J)+4.D0/5.D0*G(1)**2
     $                 *YUDYU(I,J)+(6.D0*G(2)**2+4.D0/5.D0
     $                 *G(1)**2)*YDDYD(I,J)+(-16.D0/9.D0*G(3)**4
     $                 +8.D0*G(3)**2*G(2)**2+8.D0/9.D0
     $                 *G(3)**2*G(1)**2+15.D0/2.D0*G(2)**4
     $                 +G(2)**2*G(1)**2+287.D0/90.D0
     $                 *G(1)**4)*ID(I,J)
            DUM2E(I,J)=-3.D0*DTRACE(DUM2E1)*ID(I,J)-3.D0*YEDYE(I,J)
     $                 *DTRACE(DUM2E2)-4.D0*YEDYEYEDYE(I,J)+(16.D0
     $                 *G(3)**2-2.D0/5.D0*G(1)**2)*TYDDYD
     $                 *ID(I,J)+6.D0/5.D0*G(1)**2*TYEDYE*ID(I,J)
     $                 +6.D0*G(2)**2*YEDYE(I,J)+(15.D0/2.D0
     $                 *G(2)**4+9.D0/5.D0*G(2)**2*G(1)**2
     $                 +27.D0/2.D0*G(1)**4)*ID(I,J)
          ELSE
            DUM2U(I,J)=0.D0
            DUM2D(I,J)=0.D0
            DUM2E(I,J)=0.D0
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the Yukawas
!
      DO I=1,3
        DO J=1,3
!
          B1U(I,J)=1.D0/2.D0*(3.D0*FUFUDFU(I,J)+FDFDDFU(I,J))
     $             +DMATMUL(0,FU,DUMU1,I,J)
     $             +1.D0/4.D0*DMATMUL(0,DUMU2,FU,I,J)
     $             +THSH*(-3.D0*THSW*G(128)
     $              *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $              *G(108)*SQGTPQTFTUQ(I,J))
     $             -4.D0/3.D0*THSB*THSH*G(108)
     $              *SUFTUUGTPUT(I,J)
     $             +FU(I,J)*3.D0*TFUDFU
     $             +1.D0/2.D0*THSH*FU(I,J)*(3.D0*THSW*MGTHUSQ
     $              +THSB*MGTPHUSQ)
     $             -FU(I,J)*(17.D0/20.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1D(I,J)=1.D0/2.D0*(3.D0*FDFDDFD(I,J)+FUFUDFD(I,J))
     $             +DMATMUL(0,FD,DUMD1,I,J)
     $             +1.D0/4.D0*DMATMUL(0,DUMD2,FD,I,J)
     $             +THSH*(-3.D0*THSW*G(129)
     $              *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $              *G(109)*SQGTPQTFTDQ(I,J))
     $             -2.D0/3.D0*THSB*THSH*G(109)
     $              *SDFTDDGTPDT(I,J)
     $             +FD(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FD(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FD(I,J)*(3.D0/12.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1E(I,J)=3.D0/2.D0*FEFEDFE(I,J)
     $             +DMATMUL(0,FE,DUME1,I,J)
     $             +1.D0/4.D0*DMATMUL(0,DUME2,FE,I,J)
     $             +THSH*(-3.D0*THSW*G(129)
     $              *SLGTLTFTEL(I,J)+THSB
     $              *G(109)*SLGTPLTFTEL(I,J))
     $             -2.D0*THSB*THSH*G(109)
     $              *SEFTEEGTPET(I,J)
     $             +FE(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FE(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FE(I,J)*(9.D0/4.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2)
!
          IF(SW2LP.EQ.1)THEN
            B2YMU(I,J)=DMATMUL(0,YU,DUM2U,I,J)
            B2YMD(I,J)=DMATMUL(0,YD,DUM2D,I,J)
            B2YME(I,J)=DMATMUL(0,YE,DUM2E,I,J)
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Convert into form readable by RKSTP. The transpose in BETA2 takes
!account of the differences in notation.
!
          F(3+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I)
          F(12+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YMD(J,I)
          F(21+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YME(J,I)
        END DO
      END DO
!
!The lambdas use the same dummy matrices. I will reuse the betas
!Only find the lambdas if we are below m_H
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            B1U(I,J)=0.D0
            B1D(I,J)=0.D0
            B1E(I,J)=0.D0
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            IF(SW2LP.EQ.1)THEN
              DUM2U(I,J)=3.D0/2.D0*LYUDLYU2(I,J)
     $                  -LYUDLYULYDDLYD(I,J)
     $                  -1.D0/4.D0*LYDDLYDLYUDLYU(I,J)
     $                  +11.D0/4.D0*LYDDLYD2(I,J)
     $                  +Y2*(5.D0/4.D0*LYDDLYD(I,J)-9.D0/4.D0
     $                  *LYUDLYU(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(62)**2*ID(I,J)-2.D0*G(62)*(3.D0
     $                  *LYUDLYU(I,J)+LYDDLYD(I,J))+(221.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYUDLYU(I,J)-(17.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYDDLYD(I,J)+Y4*ID(I,J)+((7.D0/150.D0
     $                  +2.D0/3.D0*NG)*G(1)**4-9.D0/20.D0*G(1)**2
     $                  *G(2)**2+19.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2D(I,J)=3.D0/2.D0*LYDDLYD2(I,J)
     $                  -LYDDLYDLYUDLYU(I,J)
     $                  -1.D0/4.D0*LYUDLYULYDDLYD(I,J)
     $                  +11.D0/4.D0*LYUDLYU2(I,J)
     $                  +Y2*(5.D0/4.D0*LYUDLYU(I,J)-9.D0/4.D0
     $                  *LYDDLYD(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(62)**2*ID(I,J)-2.D0*G(62)*(3.D0
     $                  *LYDDLYD(I,J)+LYUDLYU(I,J))+(161.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYDDLYD(I,J)-(77.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYUDLYU(I,J)+Y4*ID(I,J)+(-(37.D0/300.D0
     $                  -4.D0/15.D0*NG)*G(1)**4-27.D0/20.D0*G(1)**2
     $                  *G(2)**2+31.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2E(I,J)=3.D0/2.D0*LYEDLYE2(I,J)
     $                  -Y2*9.D0/4.D0*LYEDLYE(I,J)-CHI4*ID(I,J)
     $                  +3.D0/2.D0*G(62)**2*ID(I,J)-6.D0*G(62)
     $                  *LYEDLYE(I,J)+(441.D0/80.D0*G(1)**2
     $                  +117.D0/16.D0*G(2)**2)*LYEDLYE(I,J)
     $                  +Y4*ID(I,J)+((21.D0/100.D0+8.D0/5.D0*NG)
     $                  *G(1)**4+27.D0/20.D0*G(1)**2*G(2)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4)*ID(I,J)
            ELSE
              DUM2U(I,J)=0.D0
              DUM2D(I,J)=0.D0
              DUM2E(I,J)=0.D0
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            B1U(I,J)=1.D0/2.D0*(3.D0*LULUDLU(I,J)+LDLDDLU(I,J)
     $                -4.D0*LDLDDLU(I,J))
     $               +DMATMUL(0,LU,DUMU1,I,J)
     $               +1.D0/4.D0*DMATMUL(0,DUMU2,LU,I,J)
     $               +THSH*(-3.D0*THSW*G(213)
     $                *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $                *G(211)*SQGTPQTFTUQ(I,J))
     $               -4.D0/3.D0*THSB*THSH*G(211)
     $                *SUFTUUGTPUT(I,J)
     $               +LU(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LU(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LU(I,J)*(17.D0/20.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1D(I,J)=1.D0/2.D0*(3.D0*LDLDDLD(I,J)+LULUDLD(I,J)
     $                -4.D0*LULUDLD(I,J))
     $               +DMATMUL(0,LD,DUMD1,I,J)
     $               +1.D0/4.D0*DMATMUL(0,DUMD2,LD,I,J)
     $               +THSH*(-3.D0*THSW*G(214)
     $                *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $                *G(212)*SQGTPQTFTDQ(I,J))
     $               -2.D0/3.D0*THSB*THSH*G(212)
     $                *SDFTDDGTPDT(I,J)
     $               +LD(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LD(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LD(I,J)*(3.D0/12.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1E(I,J)=3.D0/2.D0*LELEDLE(I,J)+DMATMUL(0,LE,DUME1,I,J)
     $               +1.D0/4.D0*DMATMUL(0,DUME2,LE,I,J)
     $               +THSH*(-3.D0*THSW*G(214)
     $                *SLGTLTFTEL(I,J)+THSB
     $                *G(212)*SLGTPLTFTEL(I,J))
     $               -2.D0*THSB*THSH*G(212)
     $                *SEFTEEGTPET(I,J)
     $               +LE(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LE(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LE(I,J)*(9.D0/4.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2)
!
            IF(SW2LP.EQ.1)THEN
              BETA2U(I,J)=DMATMUL(0,LYU,DUM2U,I,J)
              BETA2D(I,J)=DMATMUL(0,LYD,DUM2D,I,J)
              BETA2E(I,J)=DMATMUL(0,LYE,DUM2E,I,J)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
!
!Convert into form readable by RKSTP.
!
            F(33+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2U(J,I)
            F(42+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2D(J,I)
            F(51+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2E(J,I)
          END DO
        END DO
      END IF
!
!Next I am going to work out the gaugino terms, \mu and M_{1,2,3}
!and the running of B in MV notation
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*SUFTUUFTUUD(I,J)+3.D0*SDFTDDFTDDD(I,J)
     $                   +SEFTEEFTEED(I,J)+3.D0*SQFTUQDFTUQ(I,J)
     $                   +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
!
          IF(SW2LP.EQ.1)THEN
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2GRKMU(I,J)=0.D0
          END IF
        END DO
      END DO
      TDUMGRKMU=DTRACE(DUMGRKMU1)
!
      BGRKMU=1.D0/2.D0*G(31)*THSH*TDUMGRKMU
     $       +1.D0/4.D0*G(31)*THSH*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ
     $        +3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       -G(31)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMU*(-3.D0*DTRACE(DUM2GRKMU)+(16.D0*G(3)**2
     $             +4.D0/5.D0*G(1)**2)*TYUDYU+(16.D0*G(3)**2
     $             -2.D0/5.D0*G(1)**2)*TYDDYD+6.D0/5.D0*G(1)**2
     $             *TYEDYE+15.D0/2.D0*G(2)**4+9.D0/5.D0*G(1)**2
     $             *G(2)**2+207.D0/50.D0*G(1)**4)
      END IF
!
!The RKSTP compatible derivative is:
!
      F(31)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $      *BETA2GRKMU
!
!V_U and V_D - from PRD 49,4882 (1994)
!
      BETA1VU=3.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-3.D0*TYUDYU
      BETA1VD=3.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-3.D0*TYDDYD
     $        -TYEDYE
      IF(SW2LP.EQ.1)THEN
        BETA2VU=3.D0/4.D0*(3.D0*TYUDYUYUDYU+3.D0*TYUDYUYDDYD)-(19.D0
     $          /10.D0*G(1)**2+9.D0/2.D0*G(2)**2+20.D0*G(3)**2)
     $          *TYUDYU-(279.D0/800.D0+1803.D0/1600.D0*3.D0)*G(1)**4
     $          -(207.D0/32.D0+357.D0/64.D0*3.D0)*G(2)**4-(27.D0/80.D0
     $          +9.D0/80.D0*3.D0)*G(1)**2*G(2)**2
        BETA2VD=3.D0/4.D0*(3.D0*TYDDYDYDDYD+3.D0*TYDDYDYUDYU
     $          +TYEDYEYEDYE)-(2.D0/5.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $          +20.D0*G(3)**2)*TYDDYD-(9.D0/5.D0*G(1)**2+3.D0/2.D0
     $          *G(2)**2)*TYEDYE-(279.D0/800.D0+1803.D0/1600.D0*3.D0)
     $          *G(1)**4-(207.D0/32.D0+357.D0/64.D0*3.D0)*G(2)**4
     $          -(27.D0/80.D0+9.D0/80.D0*3.D0)*G(1)**2*G(2)**2
      END IF
!
      F(32)=G(32)*(1.D0/16.D0/PI**2*BETA1VU+1.D0/(16.D0*PI**2)**2
     $       *BETA2VU)
      F(33)=G(33)*(1.D0/16.D0/PI**2*BETA1VD+1.D0/(16.D0*PI**2)**2
     $       *BETA2VD)
!
      IF(THHH.EQ.0)THEN
!
!Finally we have the running of the Higgs Quartic Coupling and SM VEV.
!Programmed here is the MS-bar running. It therefore needs the MS-bar
!gauge and Yukawas.
!The gauge couplings and Yukawas needed to be converted to MS-bar using
!the Martin and Vaughn conversion in hep-ph/9308222.
!The following is after the conversion, so all Yukawas and Gauge
!couplings are still in the DR-bar scheme.
!
        BETALAM1=12*G(62)**2-(9.D0/5.D0*G(1)**2+9.D0*G(2)**2)
     $           *G(62)+9.D0/4.D0*(3.D0/25.D0*G(1)**4+2.D0/5.D0
     $           *G(1)**2*G(2)**2+G(2)**4)+4.D0*Y2*G(62)-4*H
        IF(SW2LP.EQ.1)THEN
          BETALAM2=-78.D0*G(62)**3+18.D0*(3.D0/5.D0*G(1)**2
     $            +3.D0*G(2)**2)*G(62)**2-((265.D0/8.D0-10*NG)
     $            *G(2)**4-117.D0/20.D0*G(1)**2*G(2)**2
     $            -9.D0/25.D0*(229.D0/24.D0+50.D0/9.D0*NG)*G(1)**4)
     $            *G(62)+(473.D0/8.D0-8.D0*NG)*G(2)**6-3.D0/5.D0
     $            *(121.D0/24.D0+8.D0/3.D0*NG)*G(1)**2*G(2)**4
     $            -9.D0/25.D0*(239.D0/24.D0+40.D0/9.D0*NG)
     $            *G(1)**4*G(2)**2-27.D0/125.D0*(59.D0/24.D0
     $            +40.D0/9.D0*NG)*G(1)**6+(-14.D0/5.D0*G(1)**2
     $            +18.D0*G(2)**2-128.D0*G(3)**2)*TLYUDLYU2
     $            +(34.D0/5.D0*G(1)**2+18.D0*G(2)**2-128.D0
     $            *G(3)**2)*TLYDDLYD2+(-42.D0/5.D0*G(1)**2
     $            +6.D0*G(2)**2)*TLYEDLYE2-3.D0/2.D0*G(2)**4
     $            *Y2+G(62)*((83.D0/10.D0*G(1)**2+27.D0/2.D0
     $            *G(2)**2+112.D0*G(3)**2)*TLYUDLYU+(-1.D0/10.D0
     $            *G(1)**2+27.D0/2.D0*G(2)**2+112.D0*G(3)**2)
     $            *TLYDDLYD+(93.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $            *TLYEDLYE)+3.D0/5.D0*G(1)**2*((-57.D0/10.D0
     $            *G(1)**2+21.D0*G(2)**2)*TLYUDLYU+(3.D0/2.D0
     $            *G(1)**2+9.D0*G(2)**2)*TLYDDLYD+(-15.D0/2.D0
     $            *G(1)**2+11.D0*G(2)**2)*TLYEDLYE)-24.D0
     $            *G(62)**2*Y2-G(62)*H+6.D0*G(62)
     $            *TLYUDLYULYDDLYD+20.D0*(3.D0*TLYUDLYU3+3.D0*TLYDDLYD3
     $            +TLYEDLYE3)-12.D0*TDUMLUD
        END IF
!
        F(62)=1.D0/(16.D0*PI**2)*BETALAM1+1.D0/(16.D0*PI**2)**2
     $         *BETALAM2
!
!Calculate the betas for the standard model vev.
!As with lambda this is the MS-bar running with DR-bar inputs except
!v and lambda
!
      BETAVEV1=9.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-Y2
      IF(SW2LP.EQ.1)THEN
        BETAVEV2=-3.D0/2.D0*G(62)**2-(83.D0/40.D0*G(1)**2+27.D0/8.D0
     $           *G(2)**2+28.D0*G(3)**2)*TYUDYU-(-1.D0/40.D0*G(1)**2
     $           +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TYDDYD
     $           -(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TYEDYE+CHI4
     $           -27.D0/80.D0*G(1)**2*G(2)**2-(93.D0/800.D0+1.D0/2.D0
     $           *NG)*G(1)**4+(463.D0/32.D0-5.D0/2.D0*NG)*G(2)**4
      END IF
!
      F(61)=G(61)*(1.D0/(16.D0*PI**2)*BETAVEV1+1.D0/(16.D0*PI**2)**2
     $      *BETAVEV2)
!
      END IF
!
!Finally, the MSSM mu parameter
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2GRKMU(I,J)=0.D0
          END IF
        END DO
      END DO
!
      BGRKMU=MVMUM*(DTRACE(DUMGRKMU1)-3.D0*G(2)**2-3.D0/5.D0
     $           *G(1)**2)
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMUM*(-3.D0*DTRACE(DUM2GRKMU)+(16.D0*G(3)**2+4.D0
     $             /5.D0*G(1)**2)*TYUDYU+(16.D0*G(3)**2-2.D0/5.D0
     $             *G(1)**2)*TYDDYD+6.D0/5.D0*G(1)**2*TYEDYE
     $             +15.D0/2.D0*G(2)**4+9.D0/5.D0*G(1)**2*G(2)**2
     $             +207.D0/50.D0*G(1)**4)
      END IF
!
!The RKSTP compatible derivative is.
!
      F(215)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $       *BETA2GRKMU
!
      RETURN
      END
!
      SUBROUTINE DRGE601(T,GCURR,FCURR)
!
!Contains the full threshold RGEs.
!
!DRGE601 IS FOR THRESHOLDS WITH THETAS.
!
!     G(  1) = g_1         G(  2) = g_2         G(  3) = g_3
!     G(  4) = FU(1,1)     G(  5) = FU(1,2)     G( 12) = FU(3,3)
!     G( 13) = FD(1,1)     G( 22) = FE(1,1)     G( 30) = FE(3,3)
!     G( 31) = M_1         G( 32) = M_2         G( 33) = M_3
!     G( 34) = AU(1,1)     G( 43) = AD(1,1)     G( 52) = AE(1,1)
!     G( 61) = M_Hu^2+MT^2 G( 62) = M_Hd^2+MT^2 G( 63) = MQ(1,1)
!     G( 72) = ML(1,1)     G( 81) = MUP(1,1)    G( 90) = MD(1,1)
!     G( 99) = ME(1,1)     G(108) = mu          G(109) = b
!     G(110) = V_U         G(111) = V_D
!
!     G(112) = LU(1,1)     G(121) = LD(1,1)     G(130) = LE(1,1)
!
!     G(139) = GTPQ(1,1)   G(148) = GTPL(1,1)   G(157) = GTPU(1,1)
!     G(166) = GTPD(1,1)   G(175) = GTPE(1,1)   G(184) = GTP_Hu
!     G(185) = GTP_Hd      G(186) = GTQ(1,1)    G(195) = GTL(1,1)
!     G(204) = GT_Hu       G(205) = GT_Hd       G(206) = GTSQ(1,1)
!     G(215) = GTSU(1,1)   G(224) = GTSD(1,1)   G(233) = FTUQ(1,1)
!     G(242) = FTDQ(1,1)   G(251) = FTEL(1,1)   G(260) = FTUU(1,1)
!     G(269) = FTDD(1,1)   G(278) = FTEE(1,1)   G(287) = sGTP_Hu
!     G(288) = cGTP_Hd     G(289) = sGT_Hu      G(290) = cGT_Hd
!
!     G(291) = g(M)_1      G(292) = g(M)_2      G(293) = g(M)_3
!     G(294) = f(M)_U(1,1) G(303) = f(M)_D(1,1) G(312) = f(M)_E(1,1)
!     G(321) = M^(MV)_1(M) G(322) = M^(MV)_2(M) G(323) = M^(MV)_3(M)
!     G(324) = AU(M)(1,1)  G(333) = AD(M)(1,1)  G(342) = AE(M)(1,1)
!     G(351) = M(M)_Hu^2   G(352) = M(M)_Hd^2   G(353) = MQ(M)(1,1)
!     G(362) = ML(M)(1,1)  G(371) = MUP(M)(1,1) G(380) = MD(M)(1,1)
!     G(389) = ME(M)(1,1)  G(398) = mu(M)       G(399) = b(M)
!
!     G(400) = TRI_U(1,1)  G(409) = TRI_D(1,1)  G(418) = TRI_E(1,1)
!     G(427) = M_HUD       G(428) = VEV_SM      G(429) = LAM_SM
!
!     G(430) = MTSF_U(1,1) G(439) = MTSF_D(1,1) G(448) = MTSF_E(1,1)
!     G(457) = GLP         G(458) = G2L         G(459) = G3L
!     G(460) = FUHU(1,1)   G(469) = FDHD(1,1)   G(478) = FEHD(1,1)
!     G(487) = FUQ(1,1)    G(496) = FUUR(1,1)   G(505) = FDQ(1,1)
!     G(514) = FDDR(1.1)   G(523) = FEL(1,1)    G(532) = FEER(1,1)
!     G(541) = GLP^2*(C^2-S^2)      G(542) = G2L^2*(C^2-S^2)
!     G(543) = GLP*(C^2-S^2)        G(544) = G2L*(C^2-S^2)
!     G(545) = sFUQ(1,1)  G(554) = sFUUR(1,1)   G(563) = cFDQ(1,1)
!     G(572) = cFDDR(1.1) G(581) = cFEL(1,1)    G(590) = cFEER(1,1)
!     G(599) = M'_1       G(600) = M'_2         G(601) = M'_3
!
!NOTE: MQ, MU, ETC... ARE USED TO DENOTE *SQUARED* SOFT MASSES
!NOTE: BOTH G(31-33) AND G(599-601) SHOULD BE REAL. THEY ARE DEFINED
!      AS COMPLEX FOR EASE OF USE. ALSO, THE MSSM GAUGINO MASS TERMS
!      ARE MV NOTATION, ie M-iM'
!NOTE: FOR THE QUARTIC MATRICES, THE VALUES IN G ARE NOT READ FROM
!      OR WRITTEN TO BY THIS PROGRAMME. THE FINAL NUMBERS FOR THESE
!      ENTRIES IN G SHOULD BE IGNORED. IN ADDITION, IN THE CASE THAT
!      THIS CODE IS EXTENDED TO INCLUDE QUARTIC RUNNING SINCE IT IS
!      EXPECTED THAT QUARTIC RGES WILL ONLY BE DERIVABLE FOR THE
!      SQUARE OF THE TERMS ABOVE. EG FUQ*FUQ.
!NOTE: G(541) AND G(542) ARE DEFINED DIFFERENTLY IN THE DOUBLE PRECISION
!      VERSION DUE TO ISSUES WITH KEEPING MINUS SIGNS IN THE SQUARE ROOT
!
!This is the BT version which receives G in book notation
!
      IMPLICIT NONE
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE PRECISION T
      DOUBLE PRECISION G(601),GCURR(601)
      DOUBLE PRECISION F(601),FCURR(601)
!
      DOUBLE PRECISION FU(3,3),FD(3,3),FE(3,3)
      DOUBLE PRECISION YU(3,3),YD(3,3),YE(3,3)
      DOUBLE PRECISION AU(3,3),AD(3,3),AE(3,3)
      DOUBLE PRECISION HU(3,3),HD(3,3),HE(3,3),DH(3,3,3)
      DOUBLE PRECISION MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE PRECISION MQM(3,3),MLM(3,3),MUPM(3,3),MDM(3,3),MEM(3,3)
      DOUBLE PRECISION DM(5,3,3)
      DOUBLE PRECISION LU(3,3),LD(3,3),LE(3,3)
      DOUBLE PRECISION LYU(3,3),LYD(3,3),LYE(3,3)
      DOUBLE PRECISION GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3),GTPE(3,3)
      DOUBLE PRECISION GTQ(3,3),GTL(3,3)
      DOUBLE PRECISION GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE PRECISION FTUQ(3,3),FTDQ(3,3),FTEL(3,3)
      DOUBLE PRECISION FTUU(3,3),FTDD(3,3),FTEE(3,3)
!
      DOUBLE PRECISION TRIU(3,3),TRID(3,3),TRIE(3,3)
      DOUBLE PRECISION MTSFU(3,3),MTSFD(3,3),MTSFE(3,3)
      DOUBLE PRECISION FUHU(3,3),FDHD(3,3),FEHD(3,3)
      DOUBLE PRECISION FUQ(3,3),FUUR(3,3)
      DOUBLE PRECISION FDQ(3,3),FDDR(3,3)
      DOUBLE PRECISION FEL(3,3),FEER(3,3)
      DOUBLE PRECISION SFUQ(3,3),CFDQ(3,3),SFUUR(3,3)
      DOUBLE PRECISION CFDDR(3,3),CFEL(3,3),CFEER(3,3)
!
      DOUBLE PRECISION DMATMUL,SFMUL,DTRACE,TSFMUL,SUM
      DOUBLE PRECISION MODSQ,CON,DRE
!
!These are used in the calculation of the RGEs which contain
!thresholds
!
      DOUBLE PRECISION DUMU1(3,3),DUMU2(3,3),DUMD1(3,3),DUMD2(3,3)
      DOUBLE PRECISION DUME1(3,3),DUME2(3,3)
      DOUBLE PRECISION DUMM1(3,3),TDUMM1,DUMM2(3,3),TDUMM2,DUMM3(3,3)
      DOUBLE PRECISION TDUMM3,DUMGRKMU1(3,3),TDUMGRKMU
      DOUBLE PRECISION DUMTUQ1(3,3),DUMTUQ2(3,3),DUMTUQ3(3,3),TDUMTUQ2
      DOUBLE PRECISION DUMTUQ4(3,3),DUMTUQ5(3,3),DUMTDQ1(3,3)
      DOUBLE PRECISION DUMTDQ2(3,3)
      DOUBLE PRECISION TDUMTDQ2,DUMTDQ3(3,3),DUMTDQ4(3,3),DUMTDQ5(3,3)
      DOUBLE PRECISION DUMTEL1(3,3),TDUMTEL2,DUMTEL3(3,3),DUMTEL4(3,3)
      DOUBLE PRECISION DUMTEL5(3,3),TDUMTUU1,DUMTUU2(3,3),DUMTUU3(3,3)
      DOUBLE PRECISION DUMTUU4(3,3),TDUMTDD1,DUMTDD2(3,3),DUMTDD3(3,3)
      DOUBLE PRECISION DUMTDD4(3,3),TDUMTEE1,DUMTEE3(3,3)
      DOUBLE PRECISION DUMTQ1(3,3),TDUMTQ1,DUMTQ2(3,3),DUMTQ3(3,3)
      DOUBLE PRECISION DUMTQ4(3,3),DUMTQ5(3,3),TDUMTL1,DUMTL2(3,3)
      DOUBLE PRECISION DUMTL3(3,3),DUMTPQ1(3,3),TDUMTPQ1,DUMTPQ2(3,3)
      DOUBLE PRECISION DUMTPQ3(3,3),DUMTPQ4(3,3),DUMTPQ5(3,3),TDUMTPL1
      DOUBLE PRECISION DUMTPL2(3,3),DUMTPL3(3,3),DUMTPU1(3,3),TDUMTPU2
      DOUBLE PRECISION DUMTPU3(3,3),DUMTPD1(3,3),TDUMTPD2,DUMTPD3(3,3)
      DOUBLE PRECISION DUMTPE1(3,3),TDUMTPE2
      DOUBLE PRECISION DUMTSQ1(3,3),TDUMTSQ1,DUMTSQ2(3,3),DUMTSQ3(3,3)
      DOUBLE PRECISION DUMTSQ4(3,3),DUMTSQ5(3,3),DUMTSU1(3,3),TDUMTSU2
      DOUBLE PRECISION DUMTSU3(3,3),TDUMTSU3,DUMTSU4(3,3),DUMTSD1(3,3)
      DOUBLE PRECISION TDUMTSD2,TDUMTSD3,DUMTSD4(3,3)
      DOUBLE PRECISION DUMGTHU1(3,3),TDUMGTHU1,DUMGTHD1(3,3),TDUMGTHD1
      DOUBLE PRECISION DUMGTPHU1(3,3),TDUMGTPHU1,DUMGTPHU2(3,3)
      DOUBLE PRECISION TDUMGTPHU2,DUMGTPHD1(3,3),TDUMGTPHD1
      DOUBLE PRECISION DUMGTPHD2(3,3),TDUMGTPHD2
      DOUBLE PRECISION DUMLUD1(3,3),DUMLUD2(3,3),DUMLUD(3,3),TDUMLUD
      DOUBLE PRECISION DUMMTSFU1(3,3),DUMMTSFU2(3,3),DUMMTSFU3(3,3)
      DOUBLE PRECISION DUMMTSFD1(3,3),DUMMTSFD2(3,3),DUMMTSFD3(3,3)
      DOUBLE PRECISION DUMMTSFE1(3,3),DUMMTSFE2(3,3),DUMMTSFE3(3,3)
!
      DOUBLE PRECISION FUS(3,3),FDS(3,3),FES(3,3)
      DOUBLE PRECISION FUT(3,3),FDT(3,3),FET(3,3)
      DOUBLE PRECISION LUS(3,3),LDS(3,3),LES(3,3)
      DOUBLE PRECISION LUT(3,3),LDT(3,3),LET(3,3)
      DOUBLE PRECISION AUS(3,3),ADS(3,3),AES(3,3)
      DOUBLE PRECISION AUT(3,3),ADT(3,3),AET(3,3)
      DOUBLE PRECISION FTUQS(3,3),FTUQT(3,3),FTUQD(3,3),FTUUD(3,3)
      DOUBLE PRECISION FTUUS(3,3),FTUUT(3,3),FTDQS(3,3),FTDQT(3,3)
      DOUBLE PRECISION FTDDD(3,3),FTDDS(3,3),FTDDT(3,3),FTDQD(3,3)
      DOUBLE PRECISION FTELD(3,3),FTELS(3,3),FTELT(3,3),FTEED(3,3)
      DOUBLE PRECISION FTEES(3,3),FTEET(3,3),GTQD(3,3),GTLD(3,3)
      DOUBLE PRECISION GTPQD(3,3),GTSQD(3,3),GTPUD(3,3),GTPDD(3,3)
      DOUBLE PRECISION GTPED(3,3),GTPLD(3,3),GTSUD(3,3),GTSDD(3,3)
      DOUBLE PRECISION GTPUS(3,3),GTSUS(3,3),GTQS(3,3),GTLS(3,3)
      DOUBLE PRECISION GTPQS(3,3),GTPLS(3,3),GTSQS(3,3),GTPDS(3,3)
      DOUBLE PRECISION GTPES(3,3),GTSDS(3,3),GTPUT(3,3),GTQT(3,3)
      DOUBLE PRECISION GTLT(3,3),GTPQT(3,3),GTPLT(3,3),GTPDT(3,3)
      DOUBLE PRECISION GTPET(3,3),GTSUT(3,3),GTSQT(3,3),GTSDT(3,3)
      DOUBLE PRECISION FUHUD(3,3),FDHDD(3,3),FEHDD(3,3)
      DOUBLE PRECISION MTSFUS(3,3),MTSFDS(3,3),MTSFES(3,3)
      DOUBLE PRECISION MTSFUT(3,3),MTSFDT(3,3),MTSFET(3,3)
      DOUBLE PRECISION FUHUS(3,3),FDHDS(3,3),FEHDS(3,3)
      DOUBLE PRECISION TRIUS(3,3),TRIDS(3,3),TRIES(3,3)
      DOUBLE PRECISION FUHUT(3,3),FDHDT(3,3),FEHDT(3,3)
      DOUBLE PRECISION TRIUT(3,3),TRIDT(3,3),TRIET(3,3)
!
      DOUBLE PRECISION FUFUD(3,3),FDFDD(3,3),FEFED(3,3)
      DOUBLE PRECISION YUDYU(3,3),YDDYD(3,3),YEDYE(3,3)
      DOUBLE PRECISION LULUD(3,3),LDLDD(3,3),LELED(3,3)
      DOUBLE PRECISION LUDLU(3,3),LDDLD(3,3),LEDLE(3,3)
      DOUBLE PRECISION LYUDLYU(3,3),LYDDLYD(3,3),LYEDLYE(3,3)
      DOUBLE PRECISION LYUDLYU2(3,3),LYDDLYD2(3,3),LYEDLYE2(3,3)
      DOUBLE PRECISION LYUDLYU3(3,3),LYDDLYD3(3,3),LYEDLYE3(3,3)
      DOUBLE PRECISION TLYUDLYU,TLYDDLYD,TLYEDLYE
      DOUBLE PRECISION TLYUDLYU2,TLYDDLYD2,TLYEDLYE2
      DOUBLE PRECISION TLYUDLYU3,TLYDDLYD3,TLYEDLYE3
      DOUBLE PRECISION LYUDLYULYDDLYD(3,3),LYDDLYDLYUDLYU(3,3)
      DOUBLE PRECISION TLYUDLYULYDDLYD
      DOUBLE PRECISION FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
      DOUBLE PRECISION YUYUD(3,3),YDYDD(3,3),YEYED(3,3)
      DOUBLE PRECISION YUYDD(3,3),YDYUD(3,3)
      DOUBLE PRECISION FUGTPUS(3,3),FDGTPDS(3,3),FEGTPES(3,3)
      DOUBLE PRECISION FUGTSUS(3,3),FDGTSDS(3,3)
      DOUBLE PRECISION LUGTPUS(3,3),LDGTPDS(3,3),LEGTPES(3,3)
      DOUBLE PRECISION LUGTSUS(3,3),LDGTSDS(3,3)
      DOUBLE PRECISION FTUQFTUQD(3,3),FTUQSFTUQT(3,3),FTUQSFUT(3,3)
      DOUBLE PRECISION FTUQSLUT(3,3)
      DOUBLE PRECISION FTDQFTDQD(3,3),FTDQSFTDQT(3,3),FTDQSFDT(3,3)
      DOUBLE PRECISION FTDQSLDT(3,3)
      DOUBLE PRECISION FTELFTELD(3,3),FTELSFTELT(3,3),FTELSFET(3,3)
      DOUBLE PRECISION FTELSLET(3,3)
      DOUBLE PRECISION GTQSGTQT(3,3),GTQGTQD(3,3),GTPQSGTPQT(3,3)
      DOUBLE PRECISION GTPQGTPQD(3,3),GTPUTGTPUS(3,3),GTPUDGTPU(3,3)
      DOUBLE PRECISION GTPDDGTPD(3,3),GTPDTGTPDS(3,3),GTSQSGTSQT(3,3)
      DOUBLE PRECISION GTSQGTSQD(3,3),GTSUTGTSUS(3,3),GTSUDGTSU(3,3)
      DOUBLE PRECISION GTSDDGTSD(3,3),GTSDTGTSDS(3,3),GTLSGTLT(3,3)
      DOUBLE PRECISION GTLGTLD(3,3),GTPLSGTPLT(3,3),GTPLGTPLD(3,3)
      DOUBLE PRECISION GTPETGTPES(3,3),GTPEDGTPE(3,3)
      DOUBLE PRECISION FDDRFDDRD(3,3),FDQDFUQ(3,3),FTUQGTPUS(3,3)
      DOUBLE PRECISION FTDQGTPDS(3,3),FTELGTPES(3,3),FTUQFUDFTDD(3,3)
      DOUBLE PRECISION FTDQFDDFTUU(3,3),FTDQLDDFTUU(3,3)
      DOUBLE PRECISION FTUQLUDFTDD(3,3),FUURFUURD(3,3),FUQDFDQ(3,3)
      DOUBLE PRECISION FUQDFUQ(3,3),FDQDFDQ(3,3)
      DOUBLE PRECISION FUQTFUQS(3,3),FUURSFUURT(3,3)
      DOUBLE PRECISION FDQTFDQS(3,3),FDDRSFDDRT(3,3),FELTFELS(3,3)
      DOUBLE PRECISION FEERSFEERT(3,3)
      DOUBLE PRECISION GTPQSFTUU(3,3),GTQSFTUU(3,3),FDDFTUU(3,3)
      DOUBLE PRECISION LUDFTDD(3,3),LDDFTUU(3,3)
      DOUBLE PRECISION FUDFTDD(3,3),FEERFEERD(3,3),FELDFEL(3,3)
      DOUBLE PRECISION FTUUDFTUU(3,3),FTDDDFTDD(3,3),FTEEDFTEE(3,3)
      DOUBLE PRECISION GTPQSFTDD(3,3),GTQSFTDD(3,3)
      DOUBLE PRECISION GTPLSFTEE(3,3),GTLSFTEE(3,3)
      DOUBLE PRECISION FUSFUT(3,3),FDSFDT(3,3),FESFET(3,3)
      DOUBLE PRECISION LUSLUT(3,3),LDSLDT(3,3),LESLET(3,3)
      DOUBLE PRECISION FUTFUS(3,3),FDTFDS(3,3),FETFES(3,3)
      DOUBLE PRECISION LUTLUS(3,3),LDTLDS(3,3),LETLES(3,3)
      DOUBLE PRECISION FUTFTUUS(3,3),FDTFTDDS(3,3),FETFTEES(3,3)
      DOUBLE PRECISION LUTFTUUS(3,3),LDTFTDDS(3,3),LETFTEES(3,3)
      DOUBLE PRECISION TFUDFU,TFDDFD,TFEDFE
      DOUBLE PRECISION TYUDYU,TYDDYD,TYEDYE
      DOUBLE PRECISION TLUDLU,TLDDLD,TLEDLE
      DOUBLE PRECISION FUFUDFU(3,3),FUFUDFD(3,3),FUTFUSGTPU(3,3)
      DOUBLE PRECISION LUTLUSGTPU(3,3)
      DOUBLE PRECISION FUTFUSGTSU(3,3),LUTLUSGTSU(3,3)
      DOUBLE PRECISION FDTFDSGTSD(3,3),LDTLDSGTSD(3,3)
      DOUBLE PRECISION FDFDDFU(3,3)
      DOUBLE PRECISION FDTFDSGTPD(3,3),LDTLDSGTPD(3,3)
      DOUBLE PRECISION FDFDDFD(3,3),FEFEDFE(3,3)
      DOUBLE PRECISION FETFESGTPE(3,3),LETLESGTPE(3,3)
      DOUBLE PRECISION LULUDLU(3,3),LULUDLD(3,3),LDLDDLU(3,3)
      DOUBLE PRECISION LDLDDLD(3,3),LELEDLE(3,3)
      DOUBLE PRECISION FTUQFUDFU(3,3),FTDQFDDFD(3,3),FTELFEDFE(3,3)
      DOUBLE PRECISION FTUQLUDLU(3,3),FTDQLDDLD(3,3),FTELLEDLE(3,3)
      DOUBLE PRECISION SQFTUQDFTUQ(3,3),SQFTUQTFTUQS(3,3)
      DOUBLE PRECISION SQFTUQTGTQ(3,3),SQFTUQTGTPQ(3,3),SQFTUQTGTSQ(3,3)
      DOUBLE PRECISION SQFTDQDFTDQ(3,3),SQFTDQTFTDQS(3,3)
      DOUBLE PRECISION SQFTDQTGTPQ(3,3),SQFTDQTGTQ(3,3)
      DOUBLE PRECISION SQFTDQTGTSQ(3,3),SQGTQTGTQS(3,3)
      DOUBLE PRECISION SQGTQTFTUQ(3,3),SQGTQTFTDQ(3,3)
      DOUBLE PRECISION SQGTQDGTQ(3,3),SQGTPQTGTPQS(3,3)
      DOUBLE PRECISION SQGTPQTFTUQ(3,3),SQGTPQTFTDQ(3,3)
      DOUBLE PRECISION SQGTPQDGTPQ(3,3),SQGTSQTGTSQS(3,3)
      DOUBLE PRECISION SQGTSQTFTUQ(3,3),SQGTSQTFTDQ(3,3)
      DOUBLE PRECISION SQGTSQDGTSQ(3,3),SQIDMTSFU(3,3),SQIDMTSFD(3,3)
      DOUBLE PRECISION SQFUHUDMTSFU(3,3),SQFDHDDMTSFD(3,3)
      DOUBLE PRECISION SUGTPUFTUUT(3,3),SUGTPUGTPUD(3,3)
      DOUBLE PRECISION SUGTPUSGTPUT(3,3),SUGTSUGTSUD(3,3)
      DOUBLE PRECISION SUGTSUSGTSUT(3,3),SUGTSUFTUUT(3,3)
      DOUBLE PRECISION SUFTUUFTUUD(3,3),SUFTUUGTPUT(3,3)
      DOUBLE PRECISION SUFTUUSFTUUT(3,3),SUFTUUGTSUT(3,3)
      DOUBLE PRECISION SUMTSFUID(3,3),SUSQIDMTSFUID(3,3)
      DOUBLE PRECISION SUMTSFUFUQDFDQ(3,3)
      DOUBLE PRECISION SDFTDDFTDDD(3,3),SDFTDDSFTDDT(3,3)
      DOUBLE PRECISION SDFTDDGTPDT(3,3),SDFTDDGTSDT(3,3)
      DOUBLE PRECISION SDGTPDFTDDT(3,3),SDGTPDSGTPDT(3,3)
      DOUBLE PRECISION SDGTPDGTPDD(3,3),SDGTSDGTSDD(3,3)
      DOUBLE PRECISION SDGTSDSGTSDT(3,3)
      DOUBLE PRECISION SDGTSDFTDDT(3,3),SDMTSFDFDQDFUQ(3,3)
      DOUBLE PRECISION SDMTSFDID(3,3),SDSQIDMTSFDID(3,3)
      DOUBLE PRECISION SLFTELDFTEL(3,3),SLFTELTFTELS(3,3)
      DOUBLE PRECISION SLFTELTGTL(3,3),SLFTELTGTPL(3,3),SLGTLTGTLS(3,3)
      DOUBLE PRECISION SLGTLTFTEL(3,3),SLGTLDGTL(3,3),SLGTPLTGTPLS(3,3)
      DOUBLE PRECISION SLGTPLDGTPL(3,3),SLGTPLTFTEL(3,3)
      DOUBLE PRECISION SLFEHDDMTSFE(3,3),SLIDMTSFE(3,3)
      DOUBLE PRECISION SEGTPESGTPET(3,3)
      DOUBLE PRECISION SEGTPEGTPED(3,3),SEGTPEFTEET(3,3)
      DOUBLE PRECISION SEFTEEFTEED(3,3)
      DOUBLE PRECISION SEFTEESFTEET(3,3),SEFTEEGTPET(3,3),SEMTSFEID(3,3)
      DOUBLE PRECISION SESLIDMTSFEID(3,3)
      DOUBLE PRECISION TFUSFUT,TFDSFDT,TFESFET,TLUSLUT,TLDSLDT,TLESLET
      DOUBLE PRECISION TSQGTQDGTQ,TSLGTLDGTL,TSUFTUUFTUUD,TSDFTDDFTDDD
      DOUBLE PRECISION TSEFTEEFTEED,TSQFTUQDFTUQ,TSQFTDQDFTDQ
      DOUBLE PRECISION TSLFTELDFTEL,TSQGTPQDGTPQ,TSLGTPLDGTPL
      DOUBLE PRECISION TSDGTPDGTPDD,TSUGTPUGTPUD,TSEGTPEGTPED
      DOUBLE PRECISION SUTSQFUHUDMTSFU,SDTSQFDHDDMTSFD
      DOUBLE PRECISION SETSLFEHDDMTSFE
      DOUBLE PRECISION SQTSUMTSFUSAUT,SUMTSFUSAUT(3,3)
      DOUBLE PRECISION SQTSDMTSFDSADT,SDMTSFDSADT(3,3)
      DOUBLE PRECISION SLTSEMTSFESAET,SEMTSFESAET(3,3)
      DOUBLE PRECISION MGTPHUSQ,MGTPHDSQ,MGTHUSQ,MGTHDSQ
      DOUBLE PRECISION MSGTPHUSQ,MCGTPHDSQ,MSGTHUSQ,MCGTHDSQ
      DOUBLE PRECISION MMUSQ,M1PM1PSQ,M2PM2PSQ,M3PM3PSQ
      DOUBLE PRECISION Y2,H,Y4,CHI4
!
      DOUBLE PRECISION SUAUID(3,3),SUAUFUQDFUQ(3,3),SUSQIDAUID(3,3)
      DOUBLE PRECISION SQIDAU(3,3),SUTSQFUHUDAU,SQFUHUDAU(3,3)
      DOUBLE PRECISION SQFUURFUURDAU(3,3),GTPQSFUGTPUS(3,3)
      DOUBLE PRECISION SDADFDQDFUQ(3,3),GTSQSFUGTSUS(3,3)
      DOUBLE PRECISION SQIDAD(3,3),SQFDDRFDDRDAD(3,3),SDSQIDADID(3,3)
      DOUBLE PRECISION SDTSQFDHDDAD,SQFDHDDAD(3,3)
      DOUBLE PRECISION SETSLFEHDDAE,SLFEHDDAE(3,3),SDADID(3,3)
      DOUBLE PRECISION SDADFDQDFDQ(3,3),GTPQSFDGTPDS(3,3)
      DOUBLE PRECISION GTSQSFDGTSDS(3,3),SUAUFUQDFDQ(3,3)
      DOUBLE PRECISION SLIDAE(3,3),SLFEERFEERDAE(3,3),SEAEID(3,3)
      DOUBLE PRECISION SESLIDAEID(3,3),SEAEFELDFEL(3,3)
      DOUBLE PRECISION GTPLSFEGTPES(3,3)
!
      DOUBLE PRECISION SFUQDSFUQ(3,3),SFUURSFUURD(3,3)
      DOUBLE PRECISION CFDDRCFDDRD(3,3),CFDQDSFUQ(3,3)
      DOUBLE PRECISION CFDQDCFDQ(3,3)
      DOUBLE PRECISION SFUQDCFDQ(3,3)
      DOUBLE PRECISION CFEERCFEERD(3,3)
      DOUBLE PRECISION CFELDCFEL(3,3),SUTRIUID(3,3)
      DOUBLE PRECISION SUTRIUSFUQDSFUQ(3,3),SQIDTRIU(3,3)
      DOUBLE PRECISION SUSQIDTRIUID(3,3),SUTSQFUHUDTRIU
      DOUBLE PRECISION SQFUHUDTRIU(3,3),SQSFUURSFUURDTRIU(3,3)
      DOUBLE PRECISION SQCFDDRCFDDRDTRIU(3,3)
      DOUBLE PRECISION GTPQSLUGTPUS(3,3)
      DOUBLE PRECISION GTSQSLUGTSUS(3,3),SQIDTRID(3,3)
      DOUBLE PRECISION SQSFUURSFUURDTRID(3,3)
      DOUBLE PRECISION SQCFDDRCFDDRDTRID(3,3),SDSQIDTRIDID(3,3)
      DOUBLE PRECISION SDTSQFDHDDTRID,SQFDHDDTRID(3,3)
      DOUBLE PRECISION SETSLFEHDDTRIE,SLFEHDDTRIE(3,3),SDTRIDID(3,3)
      DOUBLE PRECISION SDTRIDCFDQDCFDQ(3,3)
      DOUBLE PRECISION GTPQSLDGTPDS(3,3),GTSQSLDGTSDS(3,3)
      DOUBLE PRECISION SLIDTRIE(3,3),SLCFEERCFEERDTRIE(3,3)
      DOUBLE PRECISION SESLIDTRIEID(3,3),SETRIEID(3,3)
      DOUBLE PRECISION SETRIECFELDCFEL(3,3),GTPLSLEGTPES(3,3)
!
      DOUBLE PRECISION SUTSUIDMU,SUIDMU(3,3),SUTSUFUQTFUQSMU
      DOUBLE PRECISION SUFUQTFUQSMU(3,3),SQTSQIDMQ,SQIDMQ(3,3)
      DOUBLE PRECISION SQTSQFUURSFUURTMQ,SQFUURSFUURTMQ(3,3)
      DOUBLE PRECISION SDTSDIDMD,SDIDMD(3,3),SLTSLIDML,SLIDML(3,3)
      DOUBLE PRECISION SETSEIDME,SEIDME(3,3),SQTSUAUSAUT,SUAUSAUT(3,3)
      DOUBLE PRECISION SQTSDMTSFDSMTSFDT,SDMTSFDSMTSFDT(3,3)
      DOUBLE PRECISION SLTSEMTSFESMTSFET,SEMTSFESMTSFET(3,3)
      DOUBLE PRECISION SQTSQFDDRSFDDRTMQ,SQFDDRSFDDRTMQ(3,3)
      DOUBLE PRECISION SDTSDFDQTFDQSMD,SDFDQTFDQSMD(3,3)
      DOUBLE PRECISION SLTSLFEERSFEERTML,SLFEERSFEERTML(3,3)
      DOUBLE PRECISION SETSEFELTFELSME,SEFELTFELSME(3,3)
      DOUBLE PRECISION SQTSUMTSFUSMTSFUT,SUMTSFUSMTSFUT(3,3)
      DOUBLE PRECISION SQTSDADSADT,SDADSADT(3,3),SLTSEAESAET
      DOUBLE PRECISION SEAESAET(3,3)
!
      DOUBLE PRECISION SUTSUSFUQTSFUQSMU,SUSFUQTSFUQSMU(3,3)
      DOUBLE PRECISION SFUQTSFUQS(3,3),SQTSQSFUURSSFUURTMQ
      DOUBLE PRECISION SQSFUURSSFUURTMQ(3,3),SFUURSSFUURT(3,3)
      DOUBLE PRECISION SQTSQCFDDRSCFDDRTMQ,SQCFDDRSCFDDRTMQ(3,3)
      DOUBLE PRECISION CFDDRSCFDDRT(3,3),SDTSDCFDQTCFDQSMD
      DOUBLE PRECISION SDCFDQTCFDQSMD(3,3),CFDQTCFDQS(3,3)
      DOUBLE PRECISION SLTSLCFEERSCFEERTML,SLCFEERSCFEERTML(3,3)
      DOUBLE PRECISION CFEERSCFEERT(3,3),SETSECFELTCFELSME
      DOUBLE PRECISION SECFELTCFELSME(3,3),CFELTCFELS(3,3)
      DOUBLE PRECISION SQTSUTRIUTRIUD,SUTRIUTRIUD(3,3),TRIUD(3,3)
      DOUBLE PRECISION SQTSDTRIDTRIDD,SDTRIDTRIDD(3,3),TRIDD(3,3)
      DOUBLE PRECISION SLTSETRIETRIED,SETRIETRIED(3,3),TRIED(3,3)
!
      DOUBLE PRECISION SUTMU,SQTMQ,SDTMD,SLTML,SETME
      DOUBLE PRECISION SQSQIDMQID(3,3)
      DOUBLE PRECISION SUSUFUHUSMUFUHUT(3,3),SUFUHUSMU(3,3)
      DOUBLE PRECISION SDSDFDHDSMDFDHDT(3,3),SDFDHDSMD(3,3)
      DOUBLE PRECISION SUTRIUSTRIUT(3,3),SDTRIDSTRIDT(3,3)
      DOUBLE PRECISION SUSUIDMUID(3,3),SQSQFUHUTMQFUHUS(3,3)
      DOUBLE PRECISION SQFUHUTMQ(3,3),SQAUTAUS(3,3),SQMTSFUTMTSFUS(3,3)
      DOUBLE PRECISION SQTRIUTTRIUS(3,3),FTUUTFTUUS(3,3)
      DOUBLE PRECISION SDSDIDMDID(3,3)
      DOUBLE PRECISION SQSQFDHDTMQFDHDS(3,3),SQFDHDTMQ(3,3)
      DOUBLE PRECISION SQADTADS(3,3)
      DOUBLE PRECISION SQMTSFDTMTSFDS(3,3),FTDDTFTDDS(3,3)
      DOUBLE PRECISION SQTRIDTTRIDS(3,3)
      DOUBLE PRECISION SLSLIDMLID(3,3),SEFEHDSME(3,3)
      DOUBLE PRECISION SESEFEHDSMEFEHDT(3,3)
      DOUBLE PRECISION SETRIESTRIET(3,3)
      DOUBLE PRECISION SESEIDMEID(3,3),SLSLFEHDTMLFEHDS(3,3)
      DOUBLE PRECISION SLFEHDTML(3,3)
      DOUBLE PRECISION SLAETAES(3,3),SLMTSFETMTSFES(3,3)
      DOUBLE PRECISION FTEETFTEES(3,3),SLTRIETTRIES(3,3)
!
      DOUBLE PRECISION B1U(3,3),B1D(3,3),B1E(3,3)
      DOUBLE PRECISION BFTUQ(3,3),BFTDQ(3,3),BFTEL(3,3)
      DOUBLE PRECISION BFTUU(3,3),BFTDD(3,3),BFTEE(3,3)
      DOUBLE PRECISION BGTQ(3,3),BGTL(3,3),BGTPQ(3,3),BGTPL(3,3)
      DOUBLE PRECISION BGTPU(3,3),BGTPD(3,3),BGTPE(3,3)
      DOUBLE PRECISION BGTSQ(3,3),BGTSU(3,3),BGTSD(3,3)
      DOUBLE PRECISION BGRKMU,BM(3),BMP(3),BGTHU,BGTHD,BGTPHU,BGTPHD
      DOUBLE PRECISION BSGTHU,BCGTHD,BSGTPHU,BCGTPHD,BETALAM1,BETALAM2
      DOUBLE PRECISION BMTSFU(3,3),BMTSFD(3,3),BMTSFE(3,3)
      DOUBLE PRECISION BMHUPMT,BMHDPMT,BMHUD
!
!The following are for two loop and soft RGEs
!
      DOUBLE PRECISION S,SP,SIG1,SIG2,SIG3
      DOUBLE PRECISION DUM2U11(3,3),DUM2U12(3,3),DUM2U21(3,3)
      DOUBLE PRECISION DUM2U22(3,3),DUM2U23(3,3),DUM1D11(3,3)
      DOUBLE PRECISION DUM1D12(3,3),DUM2D11(3,3),DUM2D12(3,3)
      DOUBLE PRECISION DUM2D21(3,3),DUM2D22(3,3),DUM2D23(3,3)
      DOUBLE PRECISION DUM1E11(3,3),DUM1E21(3,3),DUM2E11(3,3)
      DOUBLE PRECISION DUM2E12(3,3),DUM2E21(3,3),DUM2E22(3,3)
      DOUBLE PRECISION DUM2E23(3,3)
      DOUBLE PRECISION DUM1U1(3,3),DUM1U2(3,3),DUM2U1(3,3),DUM2U2(3,3)
      DOUBLE PRECISION DUM1D1(3,3),DUM1D2(3,3),DUM2D1(3,3),DUM2D2(3,3)
      DOUBLE PRECISION DUM1E1(3,3),DUM1E2(3,3),DUM2E1(3,3),DUM2E2(3,3)
      DOUBLE PRECISION DUM1U(3,3),DUM2U(3,3),DUM1D(3,3),DUM2D(3,3)
      DOUBLE PRECISION DUM1E(3,3),DUM2E(3,3)
      DOUBLE PRECISION DUMS(3,3),DUMSP1(3,3),DUMSP2(3,3),DUMSP3(3,3)
      DOUBLE PRECISION DUMSIG1(3,3),DUMSIG2(3,3),DUMSIG3(3,3)
      DOUBLE PRECISION DUM1HU1(3,3)
      DOUBLE PRECISION DUMSP(3,3)
      DOUBLE PRECISION DUM2HU11(3,3),DUM2HU12(3,3),DUM2HU21(3,3)
      DOUBLE PRECISION DUM2HD11(3,3)
      DOUBLE PRECISION DUM2HD12(3,3),DUM2HD13(3,3),DUM2HD21(3,3)
      DOUBLE PRECISION DUM2HD31(3,3)
      DOUBLE PRECISION DUM2HU1(3,3),DUM2HU2(3,3),DUM1HD1(3,3)
      DOUBLE PRECISION DUM1HD2(3,3)
      DOUBLE PRECISION DUM2HD1(3,3),DUM2HD2(3,3),DUM2HD3(3,3)
      DOUBLE PRECISION DUM1HU(3,3),DUM1HD(3,3)
      DOUBLE PRECISION DUM2Q31(3,3),DUM2Q41(3,3),DUM2L21(3,3)
      DOUBLE PRECISION DUM2U31(3,3),DUM2D31(3,3)
      DOUBLE PRECISION DUM1Q1(3,3),DUM1Q2(3,3),DUM1Q3(3,3),DUM2Q1(3,3)
      DOUBLE PRECISION DUM2Q2(3,3),DUM2Q3(3,3),DUM2Q4(3,3),DUM2Q5(3,3)
      DOUBLE PRECISION DUM2Q6(3,3),DUM2Q7(3,3),DUM2Q8(3,3),DUM2Q9(3,3)
      DOUBLE PRECISION DUM2QA(3,3),DUM2QB(3,3),DUM2QC(3,3),DUM2QD(3,3)
      DOUBLE PRECISION DUM1L1(3,3),DUM2L1(3,3),DUM2L2(3,3),DUM2L3(3,3)
      DOUBLE PRECISION DUM2L4(3,3),DUM2L5(3,3),DUM2L6(3,3),DUM2L7(3,3)
      DOUBLE PRECISION DUM2L8(3,3),DUM2L9(3,3),DUM2U3(3,3),DUM2U4(3,3)
      DOUBLE PRECISION DUM2U5(3,3),DUM2D3(3,3),DUM2D4(3,3),DUM2D5(3,3)
      DOUBLE PRECISION DUM2D6(3,3),DUM2D7(3,3),DUM2D8(3,3),DUM2D9(3,3)
      DOUBLE PRECISION DUM2DA(3,3),DUM2E3(3,3),DUM2E4(3,3),DUM2E5(3,3)
      DOUBLE PRECISION DUM2E6(3,3),DUM2E7(3,3),DUM2E8(3,3),DUM2E9(3,3)
      DOUBLE PRECISION DUM1GRKMU(3,3),DUM2GRKMU(3,3),DUM1B1(3,3)
      DOUBLE PRECISION DUM1B2(3,3)
      DOUBLE PRECISION DUM2B1(3,3),DUM2B2(3,3)
!
      DOUBLE PRECISION TYUDHU,TYDDHD,TYEDHE,THUYUD
      DOUBLE PRECISION THDYDD,THEYED,THUDYU,THDDYD,THEDYE,THUDHU
      DOUBLE PRECISION YUYUDYUYUD(3,3),YDYDDYDYDD(3,3)
      DOUBLE PRECISION YEYEDYEYED(3,3),YUDYUYUDYU(3,3)
      DOUBLE PRECISION YDDYDYDDYD(3,3),YEDYEYEDYE(3,3)
      DOUBLE PRECISION YDDYDYUDYU(3,3),YUDYUYDDYD(3,3)
      DOUBLE PRECISION YUYDDYDYUD(3,3),YDYUDYUYDD(3,3)
      DOUBLE PRECISION HUYUD(3,3),HDYDD(3,3),HEYED(3,3)
      DOUBLE PRECISION YUDHU(3,3),YDDHD(3,3),YEDHE(3,3)
      DOUBLE PRECISION HUDHU(3,3),HDDHD(3,3),HEDHE(3,3)
      DOUBLE PRECISION HUDYU(3,3),HDDYD(3,3),HEDYE(3,3)
      DOUBLE PRECISION HUHUD(3,3),HDHDD(3,3),HEHED(3,3)
      DOUBLE PRECISION YUHUD(3,3),YDHDD(3,3),YEHED(3,3)
      DOUBLE PRECISION HUYDD(3,3),HDYUD(3,3),HUHDD(3,3)
      DOUBLE PRECISION HDHUD(3,3),YDHUD(3,3),YUHDD(3,3)
      DOUBLE PRECISION MUPMYU(3,3),MDMYD(3,3),MEMYE(3,3)
      DOUBLE PRECISION YUMQM(3,3),YDMQM(3,3),YEMLM(3,3)
      DOUBLE PRECISION YUDMUPMYU(3,3),YDDMDMYD(3,3),YEDMEMYE(3,3)
      DOUBLE PRECISION YEDYEMLM(3,3),YUMQMYUD(3,3),YUYUDMUPM(3,3)
      DOUBLE PRECISION YDMQMYDD(3,3),YDYDDMDM(3,3),YEMLMYED(3,3)
      DOUBLE PRECISION YEYEDMEM(3,3),YUDYUMQM(3,3),YDDYDMQM(3,3)
      DOUBLE PRECISION MQMYUDYU(3,3),MQMYDDYD(3,3),MLMYEDYE(3,3)
      DOUBLE PRECISION YUDYUYDD(3,3),YUMQMYDD(3,3),YUYDDMDM(3,3)
      DOUBLE PRECISION YDMQMYUD(3,3),YDYUDMUPM(3,3)
      DOUBLE PRECISION HUYUDYUYUD(3,3),HUYDDYDYUD(3,3),HDYUDYUYDD(3,3)
      DOUBLE PRECISION HDYDDYDYDD(3,3),HEYEDYEYED(3,3),YUDYUYUDHU(3,3)
      DOUBLE PRECISION YUDHUYUDYU(3,3),YDDYDYDDHD(3,3),YDDHDYDDYD(3,3)
      DOUBLE PRECISION YDDYDYUDHU(3,3),YDDHDYUDYU(3,3),YUDHUYDDYD(3,3)
      DOUBLE PRECISION YUDYUYDDHD(3,3),YEDYEYEDHE(3,3),YEDHEYEDYE(3,3)
      DOUBLE PRECISION HUDHUYUDYU(3,3),HUDYUYUDHU(3,3),HDDHDYDDYD(3,3)
      DOUBLE PRECISION HDDYDYDDHD(3,3),HDDHDYUDYU(3,3),YDDYDHUDHU(3,3)
      DOUBLE PRECISION HDDYDYUDHU(3,3),YDDHDHUDYU(3,3),HUDHUYDDYD(3,3)
      DOUBLE PRECISION YUDYUHDDHD(3,3),HUDYUYDDHD(3,3),YUDHUHDDYD(3,3)
      DOUBLE PRECISION HEDHEYEDYE(3,3),HEDYEYEDHE(3,3),YUDYUHUDHU(3,3)
      DOUBLE PRECISION YDDYDHDDHD(3,3),YUDHUHUDYU(3,3),YDDHDHDDYD(3,3)
      DOUBLE PRECISION YEDYEHEDHE(3,3),YEDHEHEDYE(3,3),HUHUDYUYUD(3,3)
      DOUBLE PRECISION YUYUDHUHUD(3,3),HUYUDYUHUD(3,3),YUHUDHUYUD(3,3)
      DOUBLE PRECISION HUHDDYDYUD(3,3),YUYDDHDHUD(3,3),HUYDDYDHUD(3,3)
      DOUBLE PRECISION YUHDDHDYUD(3,3),HDHDDYDYDD(3,3),YDYDDHDHDD(3,3)
      DOUBLE PRECISION HDYDDYDHDD(3,3),YDHDDHDYDD(3,3),HDHUDYUYDD(3,3)
      DOUBLE PRECISION YDYUDHUHDD(3,3),HDYUDYUHDD(3,3),YDHUDHUYDD(3,3)
      DOUBLE PRECISION HEHEDYEYED(3,3),YEYEDHEHED(3,3),HEYEDYEHED(3,3)
      DOUBLE PRECISION YEHEDHEYED(3,3),YUDMUPMYUYUDYU(3,3)
      DOUBLE PRECISION YUDMUPMYUYDDYD(3,3),YUDYUMQMYDDYD(3,3)
      DOUBLE PRECISION YUDYUYDDMDMYD(3,3),YDDMDMYDYDDYD(3,3)
      DOUBLE PRECISION YEDMEMYEYEDYE(3,3),YUDYUMQMYUDYU(3,3)
      DOUBLE PRECISION YUDYUYUDMUPMYU(3,3),YUDYUYUDYUMQM(3,3)
      DOUBLE PRECISION YDDYDMQMYDDYD(3,3),YDDYDYDDMDMYD(3,3)
      DOUBLE PRECISION YDDYDYDDYDMQM(3,3),YEDYEMLMYEDYE(3,3)
      DOUBLE PRECISION YEDYEYEDMEMYE(3,3),YEDYEYEDYEMLM(3,3)
      DOUBLE PRECISION YUMQMYUDYUYUD(3,3),YUYUDMUPMYUYUD(3,3)
      DOUBLE PRECISION YUYUDYUMQMYUD(3,3),YUYUDYUYUDMUPM(3,3)
      DOUBLE PRECISION YUMQMYDDYDYUD(3,3),YUYDDMDMYDYUD(3,3)
      DOUBLE PRECISION YUYDDYDMQMYUD(3,3),YUYDDYDYUDMUPM(3,3)
      DOUBLE PRECISION YDMQMYDDYDYDD(3,3),YDYDDMDMYDYDD(3,3)
      DOUBLE PRECISION YDYDDYDMQMYDD(3,3),YDYDDYDYDDMDM(3,3)
      DOUBLE PRECISION YDMQMYUDYUYDD(3,3),YDYUDMUPMYUYDD(3,3)
      DOUBLE PRECISION YDYUDYUMQMYDD(3,3),YDYUDYUYDDMDM(3,3)
      DOUBLE PRECISION YEMLMYEDYEYED(3,3),YEYEDMEMYEYED(3,3)
      DOUBLE PRECISION YEYEDYEMLMYED(3,3),YEYEDYEYEDMEM(3,3)
      DOUBLE PRECISION TYUDYUYUDYU,TYDDYDYDDYD,TYEDYEYEDYE
      DOUBLE PRECISION TYUDYUYDDYD,TYDDYDYUDYU
!
      DOUBLE PRECISION BETA1U(3,3),BETA2U(3,3),BETA1D(3,3),BETA2D(3,3)
      DOUBLE PRECISION BETA1E(3,3),BETA2E(3,3)
      DOUBLE PRECISION BETAT2U(3,3),BETAT2D(3,3),BETAT2E(3,3)
      DOUBLE PRECISION BETA1HU,BETA2HU,BETA1HD
      DOUBLE PRECISION BETA2HD,BETA1MQ(3,3),BETA2MQ(3,3),BETA1ML(3,3)
      DOUBLE PRECISION BETA2ML(3,3),BETA1MU(3,3),BETA2MU(3,3)
      DOUBLE PRECISION BETA1MD(3,3)
      DOUBLE PRECISION BETA2MD(3,3),BETA1ME(3,3),BETA2ME(3,3)
      DOUBLE PRECISION BETA2GRKMU,BETA1B,BETA2B,BETA1VU,BETA1VD,BETA2VU
      DOUBLE PRECISION BETA2VD,B2GM(3),B1YMU(3,3),B1YMD(3,3),B1YME(3,3)
      DOUBLE PRECISION B1TRIU(3,3),B1TRID(3,3),B1TRIE(3,3)
      DOUBLE PRECISION B2YMU(3,3),B2YMD(3,3),B2YME(3,3),BETAVEV1
      DOUBLE PRECISION BETAVEV2,B2GRKMUM,B2M(3),B2HMU(3,3),B2HMD(3,3)
      DOUBLE PRECISION B2HME(3,3),B2HUM,B2HDM,B2MQM(3,3),B2MLM(3,3)
      DOUBLE PRECISION B2MUM(3,3),B2MDM(3,3),B2MEM(3,3),BETA2BM
!
      DOUBLE PRECISION ID(3,3),MVMU,MVB
      DOUBLE PRECISION B1LP(3),B1LPM(3),B2LPSM(3,3),B2LPM(3,3)
      DOUBLE PRECISION CM(3,3),CSM(3,3),SINB,COSB
      DOUBLE PRECISION PI,Q
!
      INTEGER I,J,NG,ND,NE,NNU,NU,NSQ,NSU,NSD,NSL,NSE,NSH,NH,NSW,NSG
      INTEGER THLH,THHH,THSH,THSB,THSW,THGL
!
      DATA ID(1,1)/1.D0/,ID(1,2)/0.D0/,ID(1,3)/0.D0/
      DATA ID(2,1)/0.D0/,ID(2,2)/1.D0/,ID(2,3)/0.D0/
      DATA ID(3,1)/0.D0/,ID(3,2)/0.D0/,ID(3,3)/1.D0/
      DATA B1LPM(1)/6.6D0/,B1LPM(2)/1.D0/,B1LPM(3)/-3.D0/
      DATA B2LPM(1,1)/7.96D0/,B2LPM(1,2)/5.4D0/,B2LPM(1,3)/17.6D0/
      DATA B2LPM(2,1)/1.8D0/,B2LPM(2,2)/25.D0/,B2LPM(2,3)/24.D0/
      DATA B2LPM(3,1)/2.2D0/,B2LPM(3,2)/9.D0/,B2LPM(3,3)/14.D0/
      DATA CM(1,1)/5.2D0/,CM(1,2)/2.8D0/,CM(1,3)/3.6D0/
      DATA CM(2,1)/6.D0/,CM(2,2)/6.D0/,CM(2,3)/2.D0/
      DATA CM(3,1)/4.D0/,CM(3,2)/4.D0/,CM(3,3)/0.D0/
!
      Q=SSQSTEP !This means that thresholds depend on the q being
                !calculated by the main programme instead of the
                !integration subroutine
!
      DO I=1,601
        G(I)=GCURR(I)
      END DO
!
!Set all F's and betas to zero
!
      DO I=1,601
        F(I)=0.D0
      END DO
      DO I=1,3
        DO J=1,3
          B1U(I,J)=0.D0
          B1D(I,J)=0.D0
          B1E(I,J)=0.D0
          BETA1U(I,J)=0.D0
          BETA2U(I,J)=0.D0
          BETAT2U(I,J)=0.D0
          BETA1D(I,J)=0.D0
          BETA2D(I,J)=0.D0
          BETA1E(I,J)=0.D0
          BETA2E(I,J)=0.D0
          BETA1MQ(I,J)=0.D0
          BETA2MQ(I,J)=0.D0
          BETA1ML(I,J)=0.D0
          BETA2ML(I,J)=0.D0
          BETA1MU(I,J)=0.D0
          BETA2MU(I,J)=0.D0
          BETA1MD(I,J)=0.D0
          BETA2MD(I,J)=0.D0
          BETA1ME(I,J)=0.D0
          BETA2ME(I,J)=0.D0
          BFTUQ(I,J)=0.D0
          BFTDQ(I,J)=0.D0
          BFTEL(I,J)=0.D0
          BFTUU(I,J)=0.D0
          BFTDD(I,J)=0.D0
          BFTEE(I,J)=0.D0
          BGTQ(I,J)=0.D0
          BGTL(I,J)=0.D0
          BGTPQ(I,J)=0.D0
          BGTPL(I,J)=0.D0
          BGTPU(I,J)=0.D0
          BGTPD(I,J)=0.D0
          BGTPE(I,J)=0.D0
          BGTSQ(I,J)=0.D0
          BGTSU(I,J)=0.D0
          BGTSD(I,J)=0.D0
          BMTSFU(I,J)=0.D0
          BMTSFD(I,J)=0.D0
          BMTSFE(I,J)=0.D0
!          
          B1YMU(I,J)=0.D0
          B1YMD(I,J)=0.D0
          B1YME(I,J)=0.D0
          B1TRIU(I,J)=0.D0
          B1TRID(I,J)=0.D0
          B1TRIE(I,J)=0.D0
          B2YMU(I,J)=0.D0
          B2YMD(I,J)=0.D0
          B2YME(I,J)=0.D0
          B2HMU(I,J)=0.D0
          B2HMD(I,J)=0.D0
          B2HME(I,J)=0.D0
          B2MQM(I,J)=0.D0
          B2MLM(I,J)=0.D0
          B2MUM(I,J)=0.D0
          B2MDM(I,J)=0.D0
          B2MEM(I,J)=0.D0
        END DO
        BM(I)=0.D0
        BMP(I)=0.D0
        B2GM(I)=0.D0
        B2M(I)=0.D0
      END DO
      BETA1HU=0.D0
      BETA2HU=0.D0
      BETA1HD=0.D0
      BETA2HD=0.D0
      B2HUM=0.D0
      B2HDM=0.D0
      BGRKMU=0.D0
      BETA2GRKMU=0.D0
      B2GRKMUM=0.D0
      BETA1B=0.D0
      BETA2B=0.D0
      BETA2BM=0.D0
      BETA1VU=0.D0
      BETA1VD=0.D0
      BETA2VU=0.D0
      BETA2VD=0.D0
      BGTHU=0.D0
      BGTHD=0.D0
      BGTPHU=0.D0
      BGTPHD=0.D0
      BSGTHU=0.D0
      BCGTHD=0.D0
      BSGTPHU=0.D0
      BCGTPHD=0.D0
      BETALAM1=0.D0
      BETALAM2=0.D0
      BETAVEV1=0.D0
      BETAVEV2=0.D0
!
      BMHUPMT=0.D0
      BMHDPMT=0.D0
      BMHUD=0.D0
!
      PI=4.D0*ATAN(1.D0)
      SINB=DSQRT(DBLE(TANB)**2/(1+DBLE(TANB)**2))
      COSB=SINB/TANB
      IF(Q.LT.1.D0)THEN
        WRITE(*,*)'ERROR IN Q: ',Q
        STOP
      END IF
      NG=3.D0
      NU=3
      ND=3
      NE=3
      NNU=3
      IF ((Q-QNSH).GT.ABS(EPS).OR.
     $         (ABS(Q-QNSH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSH=2
      ELSE
        NSH=0
      END IF
      IF ((Q-QNSG).GT.ABS(EPS).OR.
     $         (ABS(Q-QNSG).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSG=1
      ELSE
        NSG=0
      END IF
      IF ((Q-QNH).GT.ABS(EPS).OR.
     $         (ABS(Q-QNH).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NH=2
      ELSE
        NH=1
      END IF
      IF ((Q-QTHSB).GT.ABS(EPS).OR.
     $         (ABS(Q-QTHSB).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        THSB=1
      ELSE
        THSB=0
      END IF
      IF ((Q-QTHSW).GT.ABS(EPS).OR.
     $         (ABS(Q-QTHSW).LT.ABS(EPS).AND.EPS.GT.0)) THEN
        NSW=1
        THSW=1
      ELSE
        NSW=0
        THSW=0
      END IF
      THLH=1
      THHH=NH/2 !This works so long as THHH is an integer variable
      IF(THHH.NE.1)THEN !Perform check
        IF(THHH.NE.0)WRITE(*,*)'ERROR IN THHH'
      END IF
      THSH=NSH/2
      THGL=NSG
!
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      DO I=1,3
        NSQ=NSQ+THSQ(I)
        NSU=NSU+THSU(I)
        NSD=NSD+THSD(I)
        NSL=NSL+THSL(I)
        NSE=NSE+THSE(I)
      END DO
!
!Since I remove entries in the sfermion mass matrix
!as they decouple, I need to continue rotating even
!when all NS_ are zero.
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15)THEN
!
!Rotate to the mass basis.
!
        CALL DROTSQ(GCURR,G)
!
      END IF
!
!Convert input into 3x3 matrices. For the Yukawas and a-parameters
!I also convert notation and define Y and H as Martin and Vaughn
!parameters.
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
          AU(I,J)=G(33+(I-1)*3+J)
          AD(I,J)=G(42+(I-1)*3+J)
          AE(I,J)=G(51+(I-1)*3+J)
          MQ(I,J)=G(62+(I-1)*3+J)
          ML(I,J)=G(71+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
          ME(I,J)=G(98+(I-1)*3+J)
!
          YU(I,J)=G(293+(J-1)*3+I)
          YD(I,J)=G(302+(J-1)*3+I)
          YE(I,J)=G(311+(J-1)*3+I)
          HU(I,J)=-G(323+(J-1)*3+I)
          HD(I,J)=-G(332+(J-1)*3+I)
          HE(I,J)=-G(341+(J-1)*3+I)
          MQM(I,J)=G(352+(I-1)*3+J)
          MLM(I,J)=G(361+(I-1)*3+J)
          MUPM(I,J)=G(370+(I-1)*3+J)
          MDM(I,J)=G(379+(I-1)*3+J)
          MEM(I,J)=G(388+(I-1)*3+J)
!
          LU(I,J)=G(111+(I-1)*3+J)
          LD(I,J)=G(120+(I-1)*3+J)
          LE(I,J)=G(129+(I-1)*3+J)
!
          GTPQ(I,J)=G(138+(I-1)*3+J)
          GTPL(I,J)=G(147+(I-1)*3+J)
          GTPU(I,J)=G(156+(I-1)*3+J)
          GTPD(I,J)=G(165+(I-1)*3+J)
          GTPE(I,J)=G(174+(I-1)*3+J)
          GTQ(I,J)=G(185+(I-1)*3+J)
          GTL(I,J)=G(194+(I-1)*3+J)
          GTSQ(I,J)=G(205+(I-1)*3+J)
          GTSU(I,J)=G(214+(I-1)*3+J)
          GTSD(I,J)=G(223+(I-1)*3+J)
          FTUQ(I,J)=G(232+(I-1)*3+J)
          FTDQ(I,J)=G(241+(I-1)*3+J)
          FTEL(I,J)=G(250+(I-1)*3+J)
          FTUU(I,J)=G(259+(I-1)*3+J)
          FTDD(I,J)=G(268+(I-1)*3+J)
          FTEE(I,J)=G(277+(I-1)*3+J)
!
          TRIU(I,J)=G(399+(I-1)*3+J)
          TRID(I,J)=G(408+(I-1)*3+J)
          TRIE(I,J)=G(417+(I-1)*3+J)
          MTSFU(I,J)=G(429+(I-1)*3+J)
          MTSFD(I,J)=G(438+(I-1)*3+J)
          MTSFE(I,J)=G(447+(I-1)*3+J)
!
!For now the quartics are not run independently
!If you wish to introduce independent quartic running
!see the note at the top of this subroutine.
!Care should be taken with quartics which are
!equal by SU(2) invariance.
!
          FUHU(I,J)=FU(I,J)
          FDHD(I,J)=FD(I,J)
          FEHD(I,J)=FE(I,J)
          FUQ(I,J)=FU(I,J)
          FDQ(I,J)=FD(I,J)
          FUUR(I,J)=FU(I,J)
          FDDR(I,J)=FD(I,J)
          FEL(I,J)=FE(I,J)
          FEER(I,J)=FE(I,J)
          IF(THHH.EQ.0)THEN
            SFUQ(I,J)=LU(I,J)
            CFDQ(I,J)=LD(I,J)
            SFUUR(I,J)=LU(I,J)
            CFDDR(I,J)=LD(I,J)
            CFEL(I,J)=LE(I,J)
            CFEER(I,J)=LE(I,J)
          ELSE
            SFUQ(I,J)=0.D0
            CFDQ(I,J)=0.D0
            SFUUR(I,J)=0.D0
            CFDDR(I,J)=0.D0
            CFEL(I,J)=0.D0
            CFEER(I,J)=0.D0
          END IF
        END DO
      END DO
!
!Set other quartics to their counterparts
!
      G(457)=DSQRT(3.D0/5.D0)*G(1)
      G(458)=G(2)
      G(459)=G(3)
      IF(THHH.EQ.0)THEN
        G(541)=(COSB**2-SINB**2)*3.D0/5.D0*G(1)**2
        G(542)=(COSB**2-SINB**2)*G(2)**2
        G(543)=(COSB**2-SINB**2)*DSQRT(3.D0/5.D0)*G(1)
        G(544)=(COSB**2-SINB**2)*G(2)
      ELSE
        G(541)=0.D0
        G(542)=0.D0
        G(543)=0.D0
        G(544)=0.D0
      END IF
!
!The mu and b running is partly in MV notation, so I must convert
!NB: MVMU is the MSSM MV notation \mu.
!
      MVMU=G(398)
      MVB=-G(399)
!      
!The separated out contributions are from PRD 49 4882 (1194),
!Castano,Piard,Ramond
!
      B1LP(1)=2.D0/5.D0*(17.D0/12.D0*DBLE(NU)+5.D0/12.D0*DBLE(ND)
     $        +5.D0/4.D0*DBLE(NE)+1.D0/4.D0*DBLE(NNU))+1.D0/30.D0
     $        *DBLE(NSQ)+4.D0/15.D0*DBLE(NSU)+1.D0/15.D0*DBLE(NSD)
     $        +1.D0/10.D0*DBLE(NSL)+1.D0/5.D0*DBLE(NSE)
     $        +1.D0/5.D0*DBLE(NSH)+1.D0/10.D0*DBLE(NH)
      B1LP(2)=-22.D0/3.D0+1.D0/2.D0*(DBLE(NU)+DBLE(ND))+1.D0/6.D0
     $        *(DBLE(NE)+DBLE(NNU))+1.D0/2.D0*DBLE(NSQ)
     $        +1.D0/6.D0*DBLE(NSL)+1.D0/3.D0*DBLE(NSH)+1.D0/6.D0
     $        *DBLE(NH)+4.D0/3.D0*DBLE(NSW)
      B1LP(3)=-11.D0+2.D0/3.D0*(DBLE(NU)+DBLE(ND))+1.D0/3.D0
     $        *DBLE(NSQ)+1.D0/6.D0*DBLE(NSU)+1.D0/6.D0*DBLE(NSD)
     $        +2.D0*DBLE(NSG)
      B2LPSM(1,1)=-(-NG*19.D0/15.D0-9.D0/50.D0)
      B2LPSM(1,2)=-(-NG*3.D0/5.D0-9.D0/10.D0)
      B2LPSM(1,3)=-(-NG*44.D0/15.D0)
      B2LPSM(2,1)=-(-NG*1.D0/5.D0-3.D0/10.D0)
      B2LPSM(2,2)=-(136.D0/3.D0-NG*49.D0/3.D0-13.D0/6.D0)
      B2LPSM(2,3)=-(-NG*4.D0)
      B2LPSM(3,1)=-(-NG*11.D0/30.D0)
      B2LPSM(3,2)=-(-NG*3.D0/2.D0)
      B2LPSM(3,3)=-(102.D0-NG*76.D0/3.D0)
      CSM(1,1)=1.7D0
      CSM(1,2)=.5D0
      CSM(1,3)=1.5D0
      CSM(2,1)=1.5D0
      CSM(2,2)=1.5D0
      CSM(2,3)=.5D0
      CSM(3,1)=2.D0
      CSM(3,2)=2.D0
      CSM(3,3)=0.D0
!
!I need many variations on the 3x3 matrices.
!
      CALL DAGGER(FTUQ,FTUQD)
      CALL DAGGER(FTUU,FTUUD)
      CALL DAGGER(FTDD,FTDDD)
      CALL DAGGER(FTDQ,FTDQD)
      CALL DAGGER(FTEL,FTELD)
      CALL DAGGER(FTEE,FTEED)
      CALL DAGGER(GTQ,GTQD)
      CALL DAGGER(GTL,GTLD)
      CALL DAGGER(GTPQ,GTPQD)
      CALL DAGGER(GTPU,GTPUD)
      CALL DAGGER(GTPD,GTPDD)
      CALL DAGGER(GTPL,GTPLD)
      CALL DAGGER(GTPE,GTPED)
      CALL DAGGER(GTSQ,GTSQD)
      CALL DAGGER(GTSU,GTSUD)
      CALL DAGGER(GTSD,GTSDD)
      CALL DAGGER(FUHU,FUHUD)
      CALL DAGGER(FDHD,FDHDD)
      CALL DAGGER(FEHD,FEHDD)
      CALL DAGGER(TRIU,TRIUD)
      CALL DAGGER(TRID,TRIDD)
      CALL DAGGER(TRIE,TRIED)
      DO I=1,3
        DO J=1,3
          FUS(I,J)=FU(I,J)
          FDS(I,J)=FD(I,J)
          FES(I,J)=FE(I,J)
          LUS(I,J)=LU(I,J)
          LDS(I,J)=LD(I,J)
          LES(I,J)=LE(I,J)
          AUS(I,J)=AU(I,J)
          ADS(I,J)=AD(I,J)
          AES(I,J)=AE(I,J)
          FTUQS(I,J)=FTUQ(I,J)
          FTUUS(I,J)=FTUU(I,J)
          FTDQS(I,J)=FTDQ(I,J)
          FTDDS(I,J)=FTDD(I,J)
          FTELS(I,J)=FTEL(I,J)
          FTEES(I,J)=FTEE(I,J)
          GTPUS(I,J)=GTPU(I,J)
          GTSUS(I,J)=GTSU(I,J)
          GTQS(I,J)=GTQ(I,J)
          GTLS(I,J)=GTL(I,J)
          GTPQS(I,J)=GTPQ(I,J)
          GTPLS(I,J)=GTPL(I,J)
          GTSQS(I,J)=GTSQ(I,J)
          GTPDS(I,J)=GTPD(I,J)
          GTPES(I,J)=GTPE(I,J)
          GTSDS(I,J)=GTSD(I,J)
          MTSFUS(I,J)=MTSFU(I,J)
          MTSFDS(I,J)=MTSFD(I,J)
          MTSFES(I,J)=MTSFE(I,J)
          FUHUS(I,J)=FUHU(I,J)
          FDHDS(I,J)=FDHD(I,J)
          FEHDS(I,J)=FEHD(I,J)
          TRIUS(I,J)=TRIU(I,J)
          TRIDS(I,J)=TRID(I,J)
          TRIES(I,J)=TRIE(I,J)
          FUT(I,J)=FU(J,I)
          FDT(I,J)=FD(J,I)
          FET(I,J)=FE(J,I)
          LUT(I,J)=LU(J,I)
          LDT(I,J)=LD(J,I)
          LET(I,J)=LE(J,I)
          AUT(I,J)=AU(J,I)
          ADT(I,J)=AD(J,I)
          AET(I,J)=AE(J,I)
          FTUQT(I,J)=FTUQ(J,I)
          FTUUT(I,J)=FTUU(J,I)
          FTDQT(I,J)=FTDQ(J,I)
          FTDDT(I,J)=FTDD(J,I)
          FTELT(I,J)=FTEL(J,I)
          FTEET(I,J)=FTEE(J,I)
          GTPUT(I,J)=GTPU(J,I)
          GTQT(I,J)=GTQ(J,I)
          GTLT(I,J)=GTL(J,I)
          GTPQT(I,J)=GTPQ(J,I)
          GTPLT(I,J)=GTPL(J,I)
          GTPDT(I,J)=GTPD(J,I)
          GTPET(I,J)=GTPE(J,I)
          GTSUT(I,J)=GTSU(J,I)
          GTSQT(I,J)=GTSQ(J,I)
          GTSDT(I,J)=GTSD(J,I)
          MTSFUT(I,J)=MTSFU(J,I)
          MTSFDT(I,J)=MTSFD(J,I)
          MTSFET(I,J)=MTSFE(J,I)
          FUHUT(I,J)=FUHU(J,I)
          FDHDT(I,J)=FDHD(J,I)
          FEHDT(I,J)=FEHD(J,I)
          TRIUT(I,J)=TRIU(J,I)
          TRIDT(I,J)=TRID(J,I)
          TRIET(I,J)=TRIE(J,I)
        END DO
      END DO
!
!Now all the matrix multiples
!
      DO I=1,3
        DO J=1,3
          YUDYU(I,J)=DMATMUL(1,YU,YU,I,J)
          YDDYD(I,J)=DMATMUL(1,YD,YD,I,J)
          YEDYE(I,J)=DMATMUL(1,YE,YE,I,J)
          YUYUD(I,J)=DMATMUL(2,YU,YU,I,J)
          YDYDD(I,J)=DMATMUL(2,YD,YD,I,J)
          YEYED(I,J)=DMATMUL(2,YE,YE,I,J)
          YUYDD(I,J)=DMATMUL(2,YU,YD,I,J)
          YDYUD(I,J)=DMATMUL(2,YD,YU,I,J)
!
          HUYUD(I,J)=DMATMUL(2,HU,YU,I,J)
          HDYDD(I,J)=DMATMUL(2,HD,YD,I,J)
          HEYED(I,J)=DMATMUL(2,HE,YE,I,J)
          YUDHU(I,J)=DMATMUL(1,YU,HU,I,J)
          YDDHD(I,J)=DMATMUL(1,YD,HD,I,J)
          YEDHE(I,J)=DMATMUL(1,YE,HE,I,J)
          HUDHU(I,J)=DMATMUL(1,HU,HU,I,J)
          HDDHD(I,J)=DMATMUL(1,HD,HD,I,J)
          HEDHE(I,J)=DMATMUL(1,HE,HE,I,J)
          HUDYU(I,J)=DMATMUL(1,HU,YU,I,J)
          HDDYD(I,J)=DMATMUL(1,HD,YD,I,J)
          HEDYE(I,J)=DMATMUL(1,HE,YE,I,J)
          HUHUD(I,J)=DMATMUL(2,HU,HU,I,J)
          HDHDD(I,J)=DMATMUL(2,HD,HD,I,J)
          HEHED(I,J)=DMATMUL(2,HE,HE,I,J)
          YUHUD(I,J)=DMATMUL(2,YU,HU,I,J)
          YDHDD(I,J)=DMATMUL(2,YD,HD,I,J)
          YEHED(I,J)=DMATMUL(2,YE,HE,I,J)
          HUYDD(I,J)=DMATMUL(2,HU,YD,I,J)
          HDYUD(I,J)=DMATMUL(2,HD,YU,I,J)
          HUHDD(I,J)=DMATMUL(2,HU,HD,I,J)
          HDHUD(I,J)=DMATMUL(2,HD,HU,I,J)
          YDHUD(I,J)=DMATMUL(2,YD,HU,I,J)
          YUHDD(I,J)=DMATMUL(2,YU,HD,I,J)
!
          MUPMYU(I,J)=DMATMUL(0,MUPM,YU,I,J)
          MDMYD(I,J)=DMATMUL(0,MDM,YD,I,J)
          MEMYE(I,J)=DMATMUL(0,MEM,YE,I,J)
          YUMQM(I,J)=DMATMUL(0,YU,MQM,I,J)
          YDMQM(I,J)=DMATMUL(0,YD,MQM,I,J)
          YEMLM(I,J)=DMATMUL(0,YE,MLM,I,J)
!
          FUFUD(I,J)=DMATMUL(2,FU,FU,I,J)
          FUSFUT(I,J)=DMATMUL(0,FUS,FUT,I,J)
          LUSLUT(I,J)=DMATMUL(0,LUS,LUT,I,J)
          FUTFUS(I,J)=DMATMUL(0,FUT,FUS,I,J)
          LUTLUS(I,J)=DMATMUL(0,LUT,LUS,I,J)
          FUTFTUUS(I,J)=DMATMUL(0,FUT,FTUUS,I,J)
          LUTFTUUS(I,J)=DMATMUL(0,LUT,FTUUS,I,J)
          FUGTPUS(I,J)=DMATMUL(0,FU,GTPUS,I,J)
          FUGTSUS(I,J)=DMATMUL(0,FU,GTSUS,I,J)
          LUGTPUS(I,J)=DMATMUL(0,LU,GTPUS,I,J)
          LUGTSUS(I,J)=DMATMUL(0,LU,GTSUS,I,J)
          FDFDD(I,J)=DMATMUL(2,FD,FD,I,J)
          FDTFDS(I,J)=DMATMUL(0,FDT,FDS,I,J)
          LDTLDS(I,J)=DMATMUL(0,LDT,LDS,I,J)
          FDTFTDDS(I,J)=DMATMUL(0,FDT,FTDDS,I,J)
          LDTFTDDS(I,J)=DMATMUL(0,LDT,FTDDS,I,J)
          FDSFDT(I,J)=DMATMUL(0,FDS,FDT,I,J)
          LDSLDT(I,J)=DMATMUL(0,LDS,LDT,I,J)
          FDGTPDS(I,J)=DMATMUL(0,FD,GTPDS,I,J)
          FDGTSDS(I,J)=DMATMUL(0,FD,GTSDS,I,J)
          LDGTPDS(I,J)=DMATMUL(0,LD,GTPDS,I,J)
          LDGTSDS(I,J)=DMATMUL(0,LD,GTSDS,I,J)
          FEFED(I,J)=DMATMUL(2,FE,FE,I,J)
          FETFES(I,J)=DMATMUL(0,FET,FES,I,J)
          LETLES(I,J)=DMATMUL(0,LET,LES,I,J)
          FETFTEES(I,J)=DMATMUL(0,FET,FTEES,I,J)
          LETFTEES(I,J)=DMATMUL(0,LET,FTEES,I,J)
          FESFET(I,J)=DMATMUL(0,FES,FET,I,J)
          LESLET(I,J)=DMATMUL(0,LES,LET,I,J)
          FEGTPES(I,J)=DMATMUL(0,FE,GTPES,I,J)
          LEGTPES(I,J)=DMATMUL(0,LE,GTPES,I,J)
          FUDFU(I,J)=DMATMUL(1,FU,FU,I,J)
          FDDFD(I,J)=DMATMUL(1,FD,FD,I,J)
          FEDFE(I,J)=DMATMUL(1,FE,FE,I,J)
          LUDLU(I,J)=DMATMUL(1,LU,LU,I,J)
          LDDLD(I,J)=DMATMUL(1,LD,LD,I,J)
          LEDLE(I,J)=DMATMUL(1,LE,LE,I,J)
          LULUD(I,J)=DMATMUL(2,LU,LU,I,J)
          LDLDD(I,J)=DMATMUL(2,LD,LD,I,J)
          LELED(I,J)=DMATMUL(2,LE,LE,I,J)
!
          FTUQFTUQD(I,J)=DMATMUL(2,FTUQ,FTUQ,I,J)
          FTUQSFTUQT(I,J)=DMATMUL(0,FTUQS,FTUQT,I,J)
          FTUQSFUT(I,J)=DMATMUL(0,FTUQS,FUT,I,J)
          FTUQSLUT(I,J)=DMATMUL(0,FTUQS,LUT,I,J)
          FTDQFTDQD(I,J)=DMATMUL(2,FTDQ,FTDQ,I,J)
          FTDQSFTDQT(I,J)=DMATMUL(0,FTDQS,FTDQT,I,J)
          FTDQSFDT(I,J)=DMATMUL(0,FTDQS,FDT,I,J)
          FTDQSLDT(I,J)=DMATMUL(0,FTDQS,LDT,I,J)
          FTELFTELD(I,J)=DMATMUL(2,FTEL,FTEL,I,J)
          FTELSFTELT(I,J)=DMATMUL(0,FTELS,FTELT,I,J)
          FTELSFET(I,J)=DMATMUL(0,FTELS,FET,I,J)
          FTELSLET(I,J)=DMATMUL(0,FTELS,LET,I,J)
          GTQSGTQT(I,J)=DMATMUL(0,GTQS,GTQT,I,J)
          GTQGTQD(I,J)=DMATMUL(2,GTQ,GTQ,I,J)
          GTPQSGTPQT(I,J)=DMATMUL(0,GTPQS,GTPQT,I,J)
          GTPQGTPQD(I,J)=DMATMUL(2,GTPQ,GTPQ,I,J)
          GTPUTGTPUS(I,J)=DMATMUL(0,GTPUT,GTPUS,I,J)
          GTPUDGTPU(I,J)=DMATMUL(1,GTPU,GTPU,I,J)
          GTPDDGTPD(I,J)=DMATMUL(1,GTPD,GTPD,I,J)
          GTPDTGTPDS(I,J)=DMATMUL(0,GTPDT,GTPDS,I,J)
          GTSQSGTSQT(I,J)=DMATMUL(0,GTSQS,GTSQT,I,J)
          GTSQGTSQD(I,J)=DMATMUL(2,GTSQ,GTSQD,I,J)
          GTSUTGTSUS(I,J)=DMATMUL(0,GTSUT,GTSUS,I,J)
          GTSUDGTSU(I,J)=DMATMUL(1,GTSU,GTSU,I,J)
          GTSDDGTSD(I,J)=DMATMUL(1,GTSD,GTSD,I,J)
          GTSDTGTSDS(I,J)=DMATMUL(0,GTSDT,GTSDS,I,J)
          GTLSGTLT(I,J)=DMATMUL(0,GTLS,GTLT,I,J)
          GTLGTLD(I,J)=DMATMUL(2,GTL,GTL,I,J)
          GTPLSGTPLT(I,J)=DMATMUL(0,GTPLS,GTPLT,I,J)
          GTPLGTPLD(I,J)=DMATMUL(2,GTPL,GTPL,I,J)
          GTPETGTPES(I,J)=DMATMUL(0,GTPET,GTPES,I,J)
          GTPEDGTPE(I,J)=DMATMUL(1,GTPE,GTPE,I,J)
          FUURFUURD(I,J)=DMATMUL(2,FUUR,FUUR,I,J)
          FUQDFUQ(I,J)=DMATMUL(1,FUQ,FUQ,I,J)
          FUQDFDQ(I,J)=DMATMUL(1,FUQ,FDQ,I,J)
          FDQDFUQ(I,J)=DMATMUL(1,FDQ,FUQ,I,J)
          FDQDFDQ(I,J)=DMATMUL(1,FDQ,FDQ,I,J)
          FDDRFDDRD(I,J)=DMATMUL(2,FDDR,FDDR,I,J)
          FEERFEERD(I,J)=DMATMUL(2,FEER,FEER,I,J)
          FELDFEL(I,J)=DMATMUL(1,FEL,FEL,I,J)
          FTUQGTPUS(I,J)=DMATMUL(0,FTUQ,GTPUS,I,J)
          FTDQGTPDS(I,J)=DMATMUL(0,FTDQ,GTPDS,I,J)
          FTELGTPES(I,J)=DMATMUL(0,FTEL,GTPES,I,J)
          GTPQSFTUU(I,J)=DMATMUL(0,GTPQS,FTUU,I,J)
          GTPQSFTDD(I,J)=DMATMUL(0,GTPQS,FTDD,I,J)
          GTPLSFTEE(I,J)=DMATMUL(0,GTPLS,FTEE,I,J)
          GTQSFTUU(I,J)=DMATMUL(0,GTQS,FTUU,I,J)
          GTQSFTDD(I,J)=DMATMUL(0,GTQS,FTDD,I,J)
          GTLSFTEE(I,J)=DMATMUL(0,GTLS,FTEE,I,J)
          FDDFTUU(I,J)=DMATMUL(1,FD,FTUU,I,J)
          LDDFTUU(I,J)=DMATMUL(1,LD,FTUU,I,J)
          FUDFTDD(I,J)=DMATMUL(1,FU,FTDD,I,J)
          LUDFTDD(I,J)=DMATMUL(1,LU,FTDD,I,J)
          FTUUDFTUU(I,J)=DMATMUL(0,FTUUD,FTUU,I,J)
          FTDDDFTDD(I,J)=DMATMUL(0,FTDDD,FTDD,I,J)
          FTEEDFTEE(I,J)=DMATMUL(0,FTEED,FTEE,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          FUQTFUQS(I,J)=FUQDFUQ(J,I)
          FUURSFUURT(I,J)=FUURFUURD(J,I)
          FDQTFDQS(I,J)=FDQDFDQ(J,I)
          FDDRSFDDRT(I,J)=FDDRFDDRD(J,I)
          FEERSFEERT(I,J)=FEERFEERD(J,I)
          FELTFELS(I,J)=FELDFEL(J,I)
          FTUUTFTUUS(I,J)=FTUUDFTUU(J,I)
          FTDDTFTDDS(I,J)=FTDDDFTDD(J,I)
          FTEETFTEES(I,J)=FTEEDFTEE(J,I)
        END DO
      END DO
!
      IF(THHH.EQ.0)THEN
        MSGTPHUSQ=G(287)**2
        MCGTPHDSQ=G(288)**2
        MSGTHUSQ=G(289)**2
        MCGTHDSQ=G(290)**2
      END IF
      MGTPHUSQ=G(184)**2
      MGTPHDSQ=G(185)**2
      MGTHUSQ=G(204)**2
      MGTHDSQ=G(205)**2
      MMUSQ=G(108)**2
      M1PM1PSQ=G(31)**2+G(599)**2
      M2PM2PSQ=G(32)**2+G(600)**2
      M3PM3PSQ=G(33)**2+G(601)**2
!
      TYUDYU=DTRACE(YUDYU)
      TYDDYD=DTRACE(YDDYD)
      TYEDYE=DTRACE(YEDYE)
      TYUDHU=DTRACE(YUDHU)
      TYDDHD=DTRACE(YDDHD)
      TYEDHE=DTRACE(YEDHE)
      THUYUD=DTRACE(HUYUD)
!
      TFUDFU=DTRACE(FUDFU)
      TFDDFD=DTRACE(FDDFD)
      TFEDFE=DTRACE(FEDFE)
      TLUDLU=DTRACE(LULUD)
      TLDDLD=DTRACE(LDLDD)
      TLEDLE=DTRACE(LELED)
!
      DO I=1,3
        DO J=1,3
          YUDMUPMYU(I,J)=DMATMUL(1,YU,MUPMYU,I,J)
          YDDMDMYD(I,J)=DMATMUL(1,YD,MDMYD,I,J)
          YEDMEMYE(I,J)=DMATMUL(1,YE,MEMYE,I,J)
          YEDYEMLM(I,J)=DMATMUL(1,YE,YEMLM,I,J)
          YUMQMYUD(I,J)=DMATMUL(2,YUMQM,YU,I,J)
          YUYUDMUPM(I,J)=DMATMUL(0,YUYUD,MUPM,I,J)
          YDMQMYDD(I,J)=DMATMUL(2,YDMQM,YD,I,J)
          YDYDDMDM(I,J)=DMATMUL(0,YDYDD,MDM,I,J)
          YEMLMYED(I,J)=DMATMUL(2,YEMLM,YE,I,J)
          YEYEDMEM(I,J)=DMATMUL(0,YEYED,MEM,I,J)
!
          FUFUDFU(I,J)=DMATMUL(0,FUFUD,FU,I,J)
          FUFUDFD(I,J)=DMATMUL(0,FUFUD,FD,I,J)
          FUTFUSGTPU(I,J)=DMATMUL(0,FUTFUS,GTPU,I,J)
          LUTLUSGTPU(I,J)=DMATMUL(0,LUTLUS,GTPU,I,J)
          FUTFUSGTSU(I,J)=DMATMUL(0,FUTFUS,GTSU,I,J)
          LUTLUSGTSU(I,J)=DMATMUL(0,LUTLUS,GTSU,I,J)
          FDFDDFU(I,J)=DMATMUL(0,FDFDD,FU,I,J)
          FDFDDFD(I,J)=DMATMUL(0,FDFDD,FD,I,J)
          FDTFDSGTPD(I,J)=DMATMUL(0,FDTFDS,GTPD,I,J)
          LDTLDSGTPD(I,J)=DMATMUL(0,LDTLDS,GTPD,I,J)
          FDTFDSGTSD(I,J)=DMATMUL(0,FDTFDS,GTSD,I,J)
          LDTLDSGTSD(I,J)=DMATMUL(0,LDTLDS,GTSD,I,J)
          FEFEDFE(I,J)=DMATMUL(0,FEFED,FE,I,J)
          FETFESGTPE(I,J)=DMATMUL(0,FETFES,GTPE,I,J)
          LETLESGTPE(I,J)=DMATMUL(0,LETLES,GTPE,I,J)
          LULUDLU(I,J)=DMATMUL(0,LULUD,LU,I,J)
          LULUDLD(I,J)=DMATMUL(0,LULUD,LD,I,J)
          LDLDDLU(I,J)=DMATMUL(0,LDLDD,LU,I,J)
          LDLDDLD(I,J)=DMATMUL(0,LDLDD,LD,I,J)
          LELEDLE(I,J)=DMATMUL(0,LELED,LE,I,J)
          FTUQFUDFU(I,J)=DMATMUL(0,FTUQ,FUDFU,I,J)
          FTUQLUDLU(I,J)=DMATMUL(0,FTUQ,LUDLU,I,J)
          FTDQFDDFD(I,J)=DMATMUL(0,FTDQ,FDDFD,I,J)
          FTDQLDDLD(I,J)=DMATMUL(0,FTDQ,LDDLD,I,J)
          FTELFEDFE(I,J)=DMATMUL(0,FTEL,FEDFE,I,J)
          FTELLEDLE(I,J)=DMATMUL(0,FTEL,LEDLE,I,J)
          FTUQFUDFTDD(I,J)=DMATMUL(0,FTUQ,FUDFTDD,I,J)
          FTUQLUDFTDD(I,J)=DMATMUL(0,FTUQ,LUDFTDD,I,J)
          FTDQFDDFTUU(I,J)=DMATMUL(0,FTDQ,FDDFTUU,I,J)
          FTDQLDDFTUU(I,J)=DMATMUL(0,FTDQ,LDDFTUU,I,J)
          GTPQSFUGTPUS(I,J)=DMATMUL(0,GTPQS,FUGTPUS,I,J)
          GTPQSLUGTPUS(I,J)=DMATMUL(0,GTPQS,LUGTPUS,I,J)
          GTSQSFUGTSUS(I,J)=DMATMUL(0,GTSQS,FUGTSUS,I,J)
          GTSQSLUGTSUS(I,J)=DMATMUL(0,GTSQS,LUGTSUS,I,J)
          GTPQSFDGTPDS(I,J)=DMATMUL(0,GTPQS,FDGTPDS,I,J)
          GTPQSLDGTPDS(I,J)=DMATMUL(0,GTPQS,LDGTPDS,I,J)
          GTSQSFDGTSDS(I,J)=DMATMUL(0,GTSQS,FDGTSDS,I,J)
          GTSQSLDGTSDS(I,J)=DMATMUL(0,GTSQS,LDGTSDS,I,J)
          GTPLSFEGTPES(I,J)=DMATMUL(0,GTPLS,FEGTPES,I,J)
          GTPLSLEGTPES(I,J)=DMATMUL(0,GTPLS,LEGTPES,I,J)
          SQFTUQDFTUQ(I,J)=SFMUL(THSQ,FTUQD,FTUQ,I,J)
          SQFTUQTFTUQS(I,J)=SFMUL(THSQ,FTUQT,FTUQS,I,J)
          SQFTUQTGTQ(I,J)=SFMUL(THSQ,FTUQT,GTQ,I,J)
          SQFTUQTGTPQ(I,J)=SFMUL(THSQ,FTUQT,GTPQ,I,J)
          SQFTUQTGTSQ(I,J)=SFMUL(THSQ,FTUQT,GTSQ,I,J)
          SQFTDQDFTDQ(I,J)=SFMUL(THSQ,FTDQD,FTDQ,I,J)
          SQFTDQTFTDQS(I,J)=SFMUL(THSQ,FTDQT,FTDQS,I,J)
          SQFTDQTGTQ(I,J)=SFMUL(THSQ,FTDQT,GTQ,I,J)
          SQFTDQTGTPQ(I,J)=SFMUL(THSQ,FTDQT,GTPQ,I,J)
          SQFTDQTGTSQ(I,J)=SFMUL(THSQ,FTDQT,GTSQ,I,J)
          SQGTQTGTQS(I,J)=SFMUL(THSQ,GTQT,GTQS,I,J)
          SQGTQTFTUQ(I,J)=SFMUL(THSQ,GTQT,FTUQ,I,J)
          SQGTQTFTDQ(I,J)=SFMUL(THSQ,GTQT,FTDQ,I,J)
          SQGTQDGTQ(I,J)=SFMUL(THSQ,GTQD,GTQ,I,J)
          SQGTPQTGTPQS(I,J)=SFMUL(THSQ,GTPQT,GTPQS,I,J)
          SQGTPQTFTUQ(I,J)=SFMUL(THSQ,GTPQT,FTUQ,I,J)
          SQGTPQTFTDQ(I,J)=SFMUL(THSQ,GTPQT,FTDQ,I,J)
          SQGTPQDGTPQ(I,J)=SFMUL(THSQ,GTPQD,GTPQ,I,J)
          SQGTSQTGTSQS(I,J)=SFMUL(THSQ,GTSQT,GTSQS,I,J)
          SQGTSQTFTUQ(I,J)=SFMUL(THSQ,GTSQT,FTUQ,I,J)
          SQGTSQTFTDQ(I,J)=SFMUL(THSQ,GTSQT,FTDQ,I,J)
          SQGTSQDGTSQ(I,J)=SFMUL(THSQ,GTSQD,GTSQ,I,J)
          SDMTSFDSADT(I,J)=SFMUL(THSD,MTSFDS,ADT,I,J)
          SUFTUUFTUUD(I,J)=SFMUL(THSU,FTUU,FTUUD,I,J)
          SUFTUUSFTUUT(I,J)=SFMUL(THSU,FTUUS,FTUUT,I,J)
          SUFTUUGTPUT(I,J)=SFMUL(THSU,FTUU,GTPUT,I,J)
          SUFTUUGTSUT(I,J)=SFMUL(THSU,FTUU,GTSUT,I,J)
          SUGTPUFTUUT(I,J)=SFMUL(THSU,GTPU,FTUUT,I,J)
          SUGTPUSGTPUT(I,J)=SFMUL(THSU,GTPUS,GTPUT,I,J)
          SUGTPUGTPUD(I,J)=SFMUL(THSU,GTPU,GTPUD,I,J)
          SUGTSUGTSUD(I,J)=SFMUL(THSU,GTSU,GTSUD,I,J)
          SUGTSUSGTSUT(I,J)=SFMUL(THSU,GTSUS,GTSUT,I,J)
          SUGTSUFTUUT(I,J)=SFMUL(THSU,GTSU,FTUUT,I,J)
          SUMTSFUSAUT(I,J)=SFMUL(THSU,MTSFUS,AUT,I,J)
          SDFTDDFTDDD(I,J)=SFMUL(THSD,FTDD,FTDDD,I,J)
          SDFTDDSFTDDT(I,J)=SFMUL(THSD,FTDDS,FTDDT,I,J)
          SDFTDDGTPDT(I,J)=SFMUL(THSD,FTDD,GTPDT,I,J)
          SDFTDDGTSDT(I,J)=SFMUL(THSD,FTDD,GTSDT,I,J)
          SDGTPDFTDDT(I,J)=SFMUL(THSD,GTPD,FTDDT,I,J)
          SDGTPDSGTPDT(I,J)=SFMUL(THSD,GTPDS,GTPDT,I,J)
          SDGTPDGTPDD(I,J)=SFMUL(THSD,GTPD,GTPDD,I,J)
          SDGTSDFTDDT(I,J)=SFMUL(THSD,GTSD,FTDDT,I,J)
          SDGTSDGTSDD(I,J)=SFMUL(THSD,GTSD,GTSDD,I,J)
          SDGTSDSGTSDT(I,J)=SFMUL(THSD,GTSDS,GTSDT,I,J)
          SLFTELDFTEL(I,J)=SFMUL(THSL,FTELD,FTEL,I,J)
          SLFTELTFTELS(I,J)=SFMUL(THSL,FTELT,FTELS,I,J)
          SLFTELTGTL(I,J)=SFMUL(THSL,FTELT,GTL,I,J)
          SLFTELTGTPL(I,J)=SFMUL(THSL,FTELT,GTPL,I,J)
          SLGTLTGTLS(I,J)=SFMUL(THSL,GTLT,GTLS,I,J)
          SLGTLTFTEL(I,J)=SFMUL(THSL,GTLT,FTEL,I,J)
          SLGTLDGTL(I,J)=SFMUL(THSL,GTLD,GTL,I,J)
          SLGTPLTGTPLS(I,J)=SFMUL(THSL,GTPLT,GTPLS,I,J)
          SLGTPLDGTPL(I,J)=SFMUL(THSL,GTPLD,GTPL,I,J)
          SLGTPLTFTEL(I,J)=SFMUL(THSL,GTPLT,FTEL,I,J)
          SEMTSFESAET(I,J)=SFMUL(THSE,MTSFES,AET,I,J)
          SEFTEEFTEED(I,J)=SFMUL(THSE,FTEE,FTEED,I,J)
          SEFTEESFTEET(I,J)=SFMUL(THSE,FTEES,FTEET,I,J)
          SEFTEEGTPET(I,J)=SFMUL(THSE,FTEE,GTPET,I,J)
          SEGTPEFTEET(I,J)=SFMUL(THSE,GTPE,FTEET,I,J)
          SEGTPESGTPET(I,J)=SFMUL(THSE,GTPES,GTPET,I,J)
          SEGTPEGTPED(I,J)=SFMUL(THSE,GTPE,GTPED,I,J)
        END DO
      END DO
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            SFUQDSFUQ(I,J)=DMATMUL(1,SFUQ,SFUQ,I,J)
            SFUURSFUURD(I,J)=DMATMUL(2,SFUUR,SFUUR,I,J)
            SFUQDCFDQ(I,J)=DMATMUL(1,SFUQ,CFDQ,I,J)
            CFDDRCFDDRD(I,J)=DMATMUL(2,CFDDR,CFDDR,I,J)
            CFDQDCFDQ(I,J)=DMATMUL(1,CFDQ,CFDQ,I,J)
            CFDQDSFUQ(I,J)=DMATMUL(1,CFDQ,SFUQ,I,J)
            CFEERCFEERD(I,J)=DMATMUL(2,CFEER,CFEER,I,J)
            CFELDCFEL(I,J)=DMATMUL(1,CFEL,CFEL,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SFUQTSFUQS(I,J)=SFUQDSFUQ(J,I)
            SFUURSSFUURT(I,J)=SFUURSFUURD(J,I)
            CFDDRSCFDDRT(I,J)=CFDDRCFDDRD(J,I)
            CFDQTCFDQS(I,J)=CFDQDCFDQ(J,I)
            CFEERSCFEERT(I,J)=CFEERCFEERD(J,I)
            CFELTCFELS(I,J)=CFELDCFEL(J,I)
          END DO
        END DO
      END IF
!
!These are the two loop terms. All in MV notation
!
      IF(SW2LP.EQ.1)THEN
!
        THDYDD=DTRACE(HDYDD)
        THEYED=DTRACE(HEYED)
        THUDYU=DTRACE(HUDYU)
        THDDYD=DTRACE(HDDYD)
        THEDYE=DTRACE(HEDYE)
        THUDHU=DTRACE(HUDHU)
!
        DO I=1,3
          DO J=1,3
            YUDYUMQM(I,J)=DMATMUL(0,YUDYU,MQM,I,J)
            YDDYDMQM(I,J)=DMATMUL(0,YDDYD,MQM,I,J)
            MQMYUDYU(I,J)=DMATMUL(0,MQM,YUDYU,I,J)
            MQMYDDYD(I,J)=DMATMUL(0,MQM,YDDYD,I,J)
            MLMYEDYE(I,J)=DMATMUL(0,MLM,YEDYE,I,J)
            YUDYUYDD(I,J)=DMATMUL(2,YUDYU,YD,I,J)
            YUMQMYDD(I,J)=DMATMUL(2,YUMQM,YD,I,J)
            YUYDDMDM(I,J)=DMATMUL(0,YUYDD,MDM,I,J)
            YDMQMYUD(I,J)=DMATMUL(2,YDMQM,YU,I,J)
            YDYUDMUPM(I,J)=DMATMUL(0,YDYUD,MUPM,I,J)
            YDMQMYUD(I,J)=DMATMUL(2,YDMQM,YU,I,J)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=DMATMUL(0,YUYUD,YUYUD,I,J)
            YDYDDYDYDD(I,J)=DMATMUL(0,YDYDD,YDYDD,I,J)
            YEYEDYEYED(I,J)=DMATMUL(0,YEYED,YEYED,I,J)
            YUDYUYUDYU(I,J)=DMATMUL(0,YUDYU,YUDYU,I,J)
            YDDYDYDDYD(I,J)=DMATMUL(0,YDDYD,YDDYD,I,J)
            YEDYEYEDYE(I,J)=DMATMUL(0,YEDYE,YEDYE,I,J)
            YDDYDYUDYU(I,J)=DMATMUL(0,YDDYD,YUDYU,I,J)
            YUDYUYDDYD(I,J)=DMATMUL(0,YUDYU,YDDYD,I,J)
            YUYDDYDYUD(I,J)=DMATMUL(0,YUYDD,YDYUD,I,J)
            YDYUDYUYDD(I,J)=DMATMUL(0,YDYUD,YUYDD,I,J)
!
            HUYUDYUYUD(I,J)=DMATMUL(0,HUYUD,YUYUD,I,J)
            HUYDDYDYUD(I,J)=DMATMUL(0,HUYDD,YDYUD,I,J)
            HDYUDYUYDD(I,J)=DMATMUL(0,HDYUD,YUYDD,I,J)
            HDYDDYDYDD(I,J)=DMATMUL(0,HDYDD,YDYDD,I,J)
            HEYEDYEYED(I,J)=DMATMUL(0,HEYED,YEYED,I,J)
            YUDYUYUDHU(I,J)=DMATMUL(0,YUDYU,YUDHU,I,J)
            YUDHUYUDYU(I,J)=DMATMUL(0,YUDHU,YUDYU,I,J)
            YDDYDYDDHD(I,J)=DMATMUL(0,YDDYD,YDDHD,I,J)
            YDDHDYDDYD(I,J)=DMATMUL(0,YDDHD,YDDYD,I,J)
            YDDYDYUDHU(I,J)=DMATMUL(0,YDDYD,YUDHU,I,J)
            YDDHDYUDYU(I,J)=DMATMUL(0,YDDHD,YUDYU,I,J)
            YUDHUYDDYD(I,J)=DMATMUL(0,YUDHU,YDDYD,I,J)
            YUDYUYDDHD(I,J)=DMATMUL(0,YUDYU,YDDHD,I,J)
            YEDYEYEDHE(I,J)=DMATMUL(0,YEDYE,YEDHE,I,J)
            YEDHEYEDYE(I,J)=DMATMUL(0,YEDHE,YEDYE,I,J)
            HUDHUYUDYU(I,J)=DMATMUL(0,HUDHU,YUDYU,I,J)
            HUDYUYUDHU(I,J)=DMATMUL(0,HUDYU,YUDHU,I,J)
            HDDHDYDDYD(I,J)=DMATMUL(0,HDDHD,YDDYD,I,J)
            HDDYDYDDHD(I,J)=DMATMUL(0,HDDYD,YDDHD,I,J)
            HDDHDYUDYU(I,J)=DMATMUL(0,HDDHD,YUDYU,I,J)
            YDDYDHUDHU(I,J)=DMATMUL(0,YDDYD,HUDHU,I,J)
            HDDYDYUDHU(I,J)=DMATMUL(0,HDDYD,YUDHU,I,J)
            YDDHDHUDYU(I,J)=DMATMUL(0,YDDHD,HUDYU,I,J)
            HUDHUYDDYD(I,J)=DMATMUL(0,HUDHU,YDDYD,I,J)
            YUDYUHDDHD(I,J)=DMATMUL(0,YUDYU,HDDHD,I,J)
            HUDYUYDDHD(I,J)=DMATMUL(0,HUDYU,YDDHD,I,J)
            YUDHUHDDYD(I,J)=DMATMUL(0,YUDHU,HDDYD,I,J)
            HEDHEYEDYE(I,J)=DMATMUL(0,HEDHE,YEDYE,I,J)
            HEDYEYEDHE(I,J)=DMATMUL(0,HEDYE,YEDHE,I,J)
            YUDYUHUDHU(I,J)=DMATMUL(0,YUDYU,HUDHU,I,J)
            YDDYDHDDHD(I,J)=DMATMUL(0,YDDYD,HDDHD,I,J)
            YUDHUHUDYU(I,J)=DMATMUL(0,YUDHU,HUDYU,I,J)
            YDDHDHDDYD(I,J)=DMATMUL(0,YDDHD,HDDYD,I,J)
            YEDYEHEDHE(I,J)=DMATMUL(0,YEDYE,HEDHE,I,J)
            YEDHEHEDYE(I,J)=DMATMUL(0,YEDHE,HEDYE,I,J)
            HUHUDYUYUD(I,J)=DMATMUL(0,HUHUD,YUYUD,I,J)
            YUYUDHUHUD(I,J)=DMATMUL(0,YUYUD,HUHUD,I,J)
            HUYUDYUHUD(I,J)=DMATMUL(0,HUYUD,YUHUD,I,J)
            YUHUDHUYUD(I,J)=DMATMUL(0,YUHUD,HUYUD,I,J)
            HUHDDYDYUD(I,J)=DMATMUL(0,HUHDD,YDYUD,I,J)
            YUYDDHDHUD(I,J)=DMATMUL(0,YUYDD,HDHUD,I,J)
            HUYDDYDHUD(I,J)=DMATMUL(0,HUYDD,YDHUD,I,J)
            YUHDDHDYUD(I,J)=DMATMUL(0,YUHDD,HDYUD,I,J)
            HDHDDYDYDD(I,J)=DMATMUL(0,HDHDD,YDYDD,I,J)
            YDYDDHDHDD(I,J)=DMATMUL(0,YDYDD,HDHDD,I,J)
            HDYDDYDHDD(I,J)=DMATMUL(0,HDYDD,YDHDD,I,J)
            YDHDDHDYDD(I,J)=DMATMUL(0,YDHDD,HDYDD,I,J)
            HDHUDYUYDD(I,J)=DMATMUL(0,HDHUD,YUYDD,I,J)
            YDYUDHUHDD(I,J)=DMATMUL(0,YDYUD,HUHDD,I,J)
            HDYUDYUHDD(I,J)=DMATMUL(0,HDYUD,YUHDD,I,J)
            YDHUDHUYDD(I,J)=DMATMUL(0,YDHUD,HUYDD,I,J)
            HEHEDYEYED(I,J)=DMATMUL(0,HEHED,YEYED,I,J)
            YEYEDHEHED(I,J)=DMATMUL(0,YEYED,HEHED,I,J)
            HEYEDYEHED(I,J)=DMATMUL(0,HEYED,YEHED,I,J)
            YEHEDHEYED(I,J)=DMATMUL(0,YEHED,HEYED,I,J)
!
            YUDMUPMYUYUDYU(I,J)=DMATMUL(0,YUDMUPMYU,YUDYU,I,J)
            YUDMUPMYUYDDYD(I,J)=DMATMUL(0,YUDMUPMYU,YDDYD,I,J)
            YUDYUMQMYDDYD(I,J)=DMATMUL(0,YUDYUMQM,YDDYD,I,J)
            YUDYUYDDMDMYD(I,J)=DMATMUL(0,YUDYUYDD,MDMYD,I,J)
            YDDMDMYDYDDYD(I,J)=DMATMUL(0,YDDMDMYD,YDDYD,I,J)
            YEDMEMYEYEDYE(I,J)=DMATMUL(0,YEDMEMYE,YEDYE,I,J)
            YUDYUMQMYUDYU(I,J)=DMATMUL(0,YUDYUMQM,YUDYU,I,J)
            YUDYUYUDMUPMYU(I,J)=DMATMUL(0,YUDYU,YUDMUPMYU,I,J)
            YUDYUYUDYUMQM(I,J)=DMATMUL(0,YUDYU,YUDYUMQM,I,J)
            YDDYDMQMYDDYD(I,J)=DMATMUL(0,YDDYDMQM,YDDYD,I,J)
            YDDYDYDDMDMYD(I,J)=DMATMUL(0,YDDYD,YDDMDMYD,I,J)
            YDDYDYDDYDMQM(I,J)=DMATMUL(0,YDDYD,YDDYDMQM,I,J)
            YEDYEMLMYEDYE(I,J)=DMATMUL(0,YEDYEMLM,YEDYE,I,J)
            YEDYEYEDMEMYE(I,J)=DMATMUL(0,YEDYE,YEDMEMYE,I,J)
            YEDYEYEDYEMLM(I,J)=DMATMUL(0,YEDYE,YEDYEMLM,I,J)
            YUMQMYUDYUYUD(I,J)=DMATMUL(0,YUMQMYUD,YUYUD,I,J)
            YUYUDMUPMYUYUD(I,J)=DMATMUL(0,YUYUDMUPM,YUYUD,I,J)
            YUYUDYUMQMYUD(I,J)=DMATMUL(0,YUYUD,YUMQMYUD,I,J)
            YUYUDYUYUDMUPM(I,J)=DMATMUL(0,YUYUD,YUYUDMUPM,I,J)
            YUMQMYDDYDYUD(I,J)=DMATMUL(0,YUMQMYDD,YDYUD,I,J)
            YUYDDMDMYDYUD(I,J)=DMATMUL(0,YUYDDMDM,YDYUD,I,J)
            YUYDDYDMQMYUD(I,J)=DMATMUL(0,YUYDD,YDMQMYUD,I,J)
            YUYDDYDYUDMUPM(I,J)=DMATMUL(0,YUYDD,YDYUDMUPM,I,J)
            YDMQMYDDYDYDD(I,J)=DMATMUL(0,YDMQMYDD,YDYDD,I,J)
            YDYDDMDMYDYDD(I,J)=DMATMUL(0,YDYDDMDM,YDYDD,I,J)
            YDYDDYDMQMYDD(I,J)=DMATMUL(0,YDYDD,YDMQMYDD,I,J)
            YDYDDYDYDDMDM(I,J)=DMATMUL(0,YDYDD,YDYDDMDM,I,J)
            YDMQMYUDYUYDD(I,J)=DMATMUL(0,YDMQMYUD,YUYDD,I,J)
            YDYUDMUPMYUYDD(I,J)=DMATMUL(0,YDYUDMUPM,YUYDD,I,J)
            YDYUDYUMQMYDD(I,J)=DMATMUL(0,YDYUD,YUMQMYDD,I,J)
            YDYUDYUYDDMDM(I,J)=DMATMUL(0,YDYUD,YUYDDMDM,I,J)
            YEMLMYEDYEYED(I,J)=DMATMUL(0,YEMLMYED,YEYED,I,J)
            YEYEDMEMYEYED(I,J)=DMATMUL(0,YEYEDMEM,YEYED,I,J)
            YEYEDYEMLMYED(I,J)=DMATMUL(0,YEYED,YEMLMYED,I,J)
            YEYEDYEYEDMEM(I,J)=DMATMUL(0,YEYED,YEYEDMEM,I,J)
          END DO
        END DO
!
        TYUDYUYUDYU=DTRACE(YUDYUYUDYU)
        TYUDYUYDDYD=DTRACE(YUDYUYDDYD)
        TYDDYDYUDYU=DTRACE(YDDYDYUDYU)
        TYDDYDYDDYD=DTRACE(YDDYDYDDYD)
        TYEDYEYEDYE=DTRACE(YEDYEYEDYE)
!
!These are SM terms for the two loop running below m_H
!I am going to use LYU for the MV notation SM Yukawa
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=LU(J,I)
              LYD(I,J)=LD(J,I)
              LYE(I,J)=LE(J,I)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=DMATMUL(1,LYU,LYU,I,J)
              LYDDLYD(I,J)=DMATMUL(1,LYD,LYD,I,J)
              LYEDLYE(I,J)=DMATMUL(1,LYE,LYE,I,J)
            END DO
          END DO
          TLYUDLYU=DTRACE(LYUDLYU)
          TLYDDLYD=DTRACE(LYDDLYD)
          TLYEDLYE=DTRACE(LYEDLYE)
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=DMATMUL(0,LYUDLYU,LYUDLYU,I,J)
              LYDDLYD2(I,J)=DMATMUL(0,LYDDLYD,LYDDLYD,I,J)
              LYEDLYE2(I,J)=DMATMUL(0,LYEDLYE,LYEDLYE,I,J)
              LYUDLYULYDDLYD(I,J)=DMATMUL(0,LYUDLYU,LYDDLYD,I,J)
              LYDDLYDLYUDLYU(I,J)=DMATMUL(0,LYDDLYD,LYUDLYU,I,J)
              DUMLUD1(I,J)=LYUDLYU(I,J)+LYDDLYD(I,J)
            END DO
          END DO
          TLYUDLYU2=DTRACE(LYUDLYU2)
          TLYDDLYD2=DTRACE(LYDDLYD2)
          TLYEDLYE2=DTRACE(LYEDLYE2)
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=DMATMUL(0,DUMLUD1,LYDDLYD,I,J)
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=DMATMUL(0,LYUDLYU2,LYUDLYU,I,J)
              LYDDLYD3(I,J)=DMATMUL(0,LYDDLYD2,LYDDLYD,I,J)
              LYEDLYE3(I,J)=DMATMUL(0,LYEDLYE2,LYEDLYE,I,J)
              DUMLUD(I,J)=DMATMUL(0,LYUDLYU,DUMLUD2,I,J)
            END DO
          END DO
          TLYUDLYU3=DTRACE(LYUDLYU3)
          TLYDDLYD3=DTRACE(LYDDLYD3)
          TLYEDLYE3=DTRACE(LYEDLYE3)
          TLYUDLYULYDDLYD=DTRACE(LYUDLYULYDDLYD)
          TDUMLUD=DTRACE(DUMLUD)
!
          Y2=3.D0*TLYUDLYU+3.D0*TLYDDLYD+TLYEDLYE
          H=3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
          Y4=(83.D0/40.D0*G(1)**2+27.D0/8.D0*G(2)**2
     $       +28.D0*G(3)**2)*TLYUDLYU+(-1.D0/40.D0*G(1)**2
     $       +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TLYDDLYD
     $       +(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TLYEDLYE
          CHI4=9.D0/4.D0*(3.D0*TLYUDLYU2+3.D0*TLYDDLYD2+TLYEDLYE2
     $         -2.D0/3.D0*TLYUDLYULYDDLYD)
        END IF
      ELSE
        THDYDD=0.D0
        THEYED=0.D0
        THUDYU=0.D0
        THDDYD=0.D0
        THEDYE=0.D0
        THUDHU=0.D0
!
        DO I=1,3
          DO J=1,3
            YUDYUMQM(I,J)=0.D0
            YDDYDMQM(I,J)=0.D0
            MQMYUDYU(I,J)=0.D0
            MQMYDDYD(I,J)=0.D0
            MLMYEDYE(I,J)=0.D0
            YUDYUYDD(I,J)=0.D0
            YUMQMYDD(I,J)=0.D0
            YUYDDMDM(I,J)=0.D0
            YDMQMYUD(I,J)=0.D0
            YDYUDMUPM(I,J)=0.D0
            YDMQMYUD(I,J)=0.D0
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            YUYUDYUYUD(I,J)=0.D0
            YDYDDYDYDD(I,J)=0.D0
            YEYEDYEYED(I,J)=0.D0
            YUDYUYUDYU(I,J)=0.D0
            YDDYDYDDYD(I,J)=0.D0
            YEDYEYEDYE(I,J)=0.D0
            YDDYDYUDYU(I,J)=0.D0
            YUDYUYDDYD(I,J)=0.D0
            YUYDDYDYUD(I,J)=0.D0
            YDYUDYUYDD(I,J)=0.D0
!
            HUYUDYUYUD(I,J)=0.D0
            HUYDDYDYUD(I,J)=0.D0
            HDYUDYUYDD(I,J)=0.D0
            HDYDDYDYDD(I,J)=0.D0
            HEYEDYEYED(I,J)=0.D0
            YUDYUYUDHU(I,J)=0.D0
            YUDHUYUDYU(I,J)=0.D0
            YDDYDYDDHD(I,J)=0.D0
            YDDHDYDDYD(I,J)=0.D0
            YDDYDYUDHU(I,J)=0.D0
            YDDHDYUDYU(I,J)=0.D0
            YUDHUYDDYD(I,J)=0.D0
            YUDYUYDDHD(I,J)=0.D0
            YEDYEYEDHE(I,J)=0.D0
            YEDHEYEDYE(I,J)=0.D0
            HUDHUYUDYU(I,J)=0.D0
            HUDYUYUDHU(I,J)=0.D0
            HDDHDYDDYD(I,J)=0.D0
            HDDYDYDDHD(I,J)=0.D0
            HDDHDYUDYU(I,J)=0.D0
            YDDYDHUDHU(I,J)=0.D0
            HDDYDYUDHU(I,J)=0.D0
            YDDHDHUDYU(I,J)=0.D0
            HUDHUYDDYD(I,J)=0.D0
            YUDYUHDDHD(I,J)=0.D0
            HUDYUYDDHD(I,J)=0.D0
            YUDHUHDDYD(I,J)=0.D0
            HEDHEYEDYE(I,J)=0.D0
            HEDYEYEDHE(I,J)=0.D0
            YUDYUHUDHU(I,J)=0.D0
            YDDYDHDDHD(I,J)=0.D0
            YUDHUHUDYU(I,J)=0.D0
            YDDHDHDDYD(I,J)=0.D0
            YEDYEHEDHE(I,J)=0.D0
            YEDHEHEDYE(I,J)=0.D0
            HUHUDYUYUD(I,J)=0.D0
            YUYUDHUHUD(I,J)=0.D0
            HUYUDYUHUD(I,J)=0.D0
            YUHUDHUYUD(I,J)=0.D0
            HUHDDYDYUD(I,J)=0.D0
            YUYDDHDHUD(I,J)=0.D0
            HUYDDYDHUD(I,J)=0.D0
            YUHDDHDYUD(I,J)=0.D0
            HDHDDYDYDD(I,J)=0.D0
            YDYDDHDHDD(I,J)=0.D0
            HDYDDYDHDD(I,J)=0.D0
            YDHDDHDYDD(I,J)=0.D0
            HDHUDYUYDD(I,J)=0.D0
            YDYUDHUHDD(I,J)=0.D0
            HDYUDYUHDD(I,J)=0.D0
            YDHUDHUYDD(I,J)=0.D0
            HEHEDYEYED(I,J)=0.D0
            YEYEDHEHED(I,J)=0.D0
            HEYEDYEHED(I,J)=0.D0
            YEHEDHEYED(I,J)=0.D0
!
            YUDMUPMYUYUDYU(I,J)=0.D0
            YUDMUPMYUYDDYD(I,J)=0.D0
            YUDYUMQMYDDYD(I,J)=0.D0
            YUDYUYDDMDMYD(I,J)=0.D0
            YDDMDMYDYDDYD(I,J)=0.D0
            YEDMEMYEYEDYE(I,J)=0.D0
            YUDYUMQMYUDYU(I,J)=0.D0
            YUDYUYUDMUPMYU(I,J)=0.D0
            YUDYUYUDYUMQM(I,J)=0.D0
            YDDYDMQMYDDYD(I,J)=0.D0
            YDDYDYDDMDMYD(I,J)=0.D0
            YDDYDYDDYDMQM(I,J)=0.D0
            YEDYEMLMYEDYE(I,J)=0.D0
            YEDYEYEDMEMYE(I,J)=0.D0
            YEDYEYEDYEMLM(I,J)=0.D0
            YUMQMYUDYUYUD(I,J)=0.D0
            YUYUDMUPMYUYUD(I,J)=0.D0
            YUYUDYUMQMYUD(I,J)=0.D0
            YUYUDYUYUDMUPM(I,J)=0.D0
            YUMQMYDDYDYUD(I,J)=0.D0
            YUYDDMDMYDYUD(I,J)=0.D0
            YUYDDYDMQMYUD(I,J)=0.D0
            YUYDDYDYUDMUPM(I,J)=0.D0
            YDMQMYDDYDYDD(I,J)=0.D0
            YDYDDMDMYDYDD(I,J)=0.D0
            YDYDDYDMQMYDD(I,J)=0.D0
            YDYDDYDYDDMDM(I,J)=0.D0
            YDMQMYUDYUYDD(I,J)=0.D0
            YDYUDMUPMYUYDD(I,J)=0.D0
            YDYUDYUMQMYDD(I,J)=0.D0
            YDYUDYUYDDMDM(I,J)=0.D0
            YEMLMYEDYEYED(I,J)=0.D0
            YEYEDMEMYEYED(I,J)=0.D0
            YEYEDYEMLMYED(I,J)=0.D0
            YEYEDYEYEDMEM(I,J)=0.D0
          END DO
        END DO
!
        TYUDYUYUDYU=0.D0
        TYUDYUYDDYD=0.D0
        TYDDYDYUDYU=0.D0
        TYDDYDYDDYD=0.D0
        TYEDYEYEDYE=0.D0
!
        IF(THHH.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              LYU(I,J)=0.D0
              LYD(I,J)=0.D0
              LYE(I,J)=0.D0
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU(I,J)=0.D0
              LYDDLYD(I,J)=0.D0
              LYEDLYE(I,J)=0.D0
            END DO
          END DO
          TLYUDLYU=0.D0
          TLYDDLYD=0.D0
          TLYEDLYE=0.D0
          DO I=1,3
            DO J=1,3
              LYUDLYU2(I,J)=0.D0
              LYDDLYD2(I,J)=0.D0
              LYEDLYE2(I,J)=0.D0
              LYUDLYULYDDLYD(I,J)=0.D0
              LYDDLYDLYUDLYU(I,J)=0.D0
              DUMLUD1(I,J)=0.D0
            END DO
          END DO
          TLYUDLYU2=0.D0
          TLYDDLYD2=0.D0
          TLYEDLYE2=0.D0
          DO I=1,3
            DO J=1,3
              DUMLUD2(I,J)=0.D0
            END DO
          END DO
          DO I=1,3
            DO J=1,3
              LYUDLYU3(I,J)=0.D0
              LYDDLYD3(I,J)=0.D0
              LYEDLYE3(I,J)=0.D0
              DUMLUD(I,J)=0.D0
            END DO
          END DO
          TLYUDLYU3=0.D0
          TLYDDLYD3=0.D0
          TLYEDLYE3=0.D0
          TLYUDLYULYDDLYD=0.D0
          TDUMLUD=0.D0
!
          Y2=0.D0
          H=0.D0
          Y4=0.D0
          CHI4=0.D0
        END IF
      END IF
!
!Now I calculate the runnng of the MSSM gauge, yukawas and soft terms
!which will be kept intact for two loop running of fs, SUSY couplings
!and tildas.
!
      DO I=1,3
        SUM=0.D0
        DO J=1,3
          SUM=SUM+B2LPM(I,J)*G(290+J)**2
        END DO
        F(290+I)=G(290+I)**3/16.D0/PI**2*(B1LPM(I)
     $           +DBLE(SW2LP)/16.D0/PI**2*(SUM-(CM(I,1)*TYUDYU
     $           +CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        B2GM(I)=G(290+I)**3*DBLE(SW2LP)*(SUM-(CM(I,1)*TYUDYU
     $          +CM(I,2)*TYDDYD+CM(I,3)*TYEDYE))
      END DO
!
      DO I=1,3
        DO J=1,3
          DUM1U1(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1D1(I,J)=DUM1U1(I,J)
          DUM1E1(I,J)=DUM1U1(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U2(I,J)=DUM1D1(I,J)
            DUM2D1(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)+
     $                  YEYEDYEYED(I,J)
            DUM2D2(I,J)=DUM1D1(I,J)
            DUM2E1(I,J)=DUM2D1(I,J)
            DUM2E2(I,J)=DUM2U2(I,J)
          ELSE
            DUM2U1(I,J)=0.D0
            DUM2U2(I,J)=0.D0
            DUM2D1(I,J)=0.D0
            DUM2D2(I,J)=0.D0
            DUM2E1(I,J)=0.D0
            DUM2E2(I,J)=0.D0
          END IF
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          DUM1U(I,J)=3.D0*TYUDYU*ID(I,J)+3.D0*YUDYU(I,J)+YDDYD(I,J)
     $               +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2
     $               -13.D0/15.D0*G(291)**2)*ID(I,J)              
          DUM1D(I,J)=DTRACE(DUM1D1)*ID(I,J)+3.D0*YDDYD(I,J)+YUDYU(I,J)
     $               +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2
     $               -7.D0/15.D0*G(291)**2)*ID(I,J)
          DUM1E(I,J)=DTRACE(DUM1E1)*ID(I,J)+3.D0*YEDYE(I,J)+(-3.D0
     $               *G(292)**2-9.D0/5.D0*G(291)**2)*ID(I,J)
!
!Here are the two loop terms
!
          IF(SW2LP.EQ.1)THEN
            DUM2U(I,J)=-3.D0*DTRACE(DUM2U1)*ID(I,J)-YDDYD(I,J)
     $                 *DTRACE(DUM2U2)-9.D0*YUDYU(I,J)*TYUDYU
     $                 -4.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)-2.D0
     $                 *YDDYDYUDYU(I,J)+(16.D0*G(293)**2+4.D0/5.D0
     $                 *G(291)**2)*TYUDYU*ID(I,J)+(6.D0*G(292)**2
     $                 +2.D0/5.D0*G(291)**2)*YUDYU(I,J)+2.D0/5.D0
     $                 *G(291)**2*YDDYD(I,J)+(-16.D0/9.D0*G(293)**4
     $                 +8.D0*G(293)**2*G(292)**2+136.D0/45.D0
     $                 *G(293)**2*G(291)**2+15.D0/2.D0*G(292)**4
     $                 +G(292)**2*G(291)**2+2743.D0/450.D0
     $                 *G(291)**4)*ID(I,J)
            DUM2D(I,J)=-3.D0*DTRACE(DUM2D1)*ID(I,J)-3.D0*YUDYU(I,J)
     $                 *TYUDYU-3.D0*YDDYD(I,J)*DTRACE(DUM2D2)
     $                 -4.D0*YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)
     $                 -2.D0*YUDYUYDDYD(I,J)+(16.D0*G(293)**2
     $                 -2.D0/5.D0*G(291)**2)*TYDDYD*ID(I,J)+6.D0/5.D0
     $                 *G(291)**2*TYEDYE*ID(I,J)+4.D0/5.D0*G(291)**2
     $                 *YUDYU(I,J)+(6.D0*G(292)**2+4.D0/5.D0
     $                 *G(291)**2)*YDDYD(I,J)+(-16.D0/9.D0*G(293)**4
     $                 +8.D0*G(293)**2*G(292)**2+8.D0/9.D0
     $                 *G(293)**2*G(291)**2+15.D0/2.D0*G(292)**4
     $                 +G(292)**2*G(291)**2+287.D0/90.D0
     $                 *G(291)**4)*ID(I,J)
            DUM2E(I,J)=-3.D0*DTRACE(DUM2E1)*ID(I,J)-3.D0*YEDYE(I,J)
     $                 *DTRACE(DUM2E2)-4.D0*YEDYEYEDYE(I,J)+(16.D0
     $                 *G(293)**2-2.D0/5.D0*G(291)**2)*TYDDYD
     $                 *ID(I,J)+6.D0/5.D0*G(291)**2*TYEDYE*ID(I,J)
     $                 +6.D0*G(292)**2*YEDYE(I,J)+(15.D0/2.D0
     $                 *G(292)**4+9.D0/5.D0*G(292)**2*G(291)**2
     $                 +27.D0/2.D0*G(291)**4)*ID(I,J)
          ELSE
            DUM2U(I,J)=0.D0
            DUM2D(I,J)=0.D0
            DUM2E(I,J)=0.D0
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the Yukawas
!
      DO I=1,3
        DO J=1,3
          B1YMU(I,J)=DMATMUL(0,YU,DUM1U,I,J)
          B1YMD(I,J)=DMATMUL(0,YD,DUM1D,I,J)
          B1YME(I,J)=DMATMUL(0,YE,DUM1E,I,J)
          IF(SW2LP.EQ.1)THEN
            B2YMU(I,J)=DMATMUL(0,YU,DUM2U,I,J)
            B2YMD(I,J)=DMATMUL(0,YD,DUM2D,I,J)
            B2YME(I,J)=DMATMUL(0,YE,DUM2E,I,J)
          END IF
        END DO
      END DO
!
!Convert to differentials
!
      DO I=1,3
        DO J=1,3
        F(293+(I-1)*3+J)=1.D0/16.D0/PI**2*B1YMU(J,I)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I)
        F(302+(I-1)*3+J)=1.D0/16.D0/PI**2*B1YMD(J,I)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMD(J,I)
        F(311+(I-1)*3+J)=1.D0/16.D0/PI**2*B1YME(J,I)
     $                   +1.D0/(16.D0*PI**2)**2*B2YME(J,I)
        END DO
      END DO
!
!Similarly, the MSSM Gaugino mass parameters are:
!
      DO I=1,3
        SUM=0.D0
        DO J=1,3
          SUM=SUM+B2LPM(I,J)*G(290+J)**2*(G(320+I)+G(320+J))
        END DO
        IF(SW2LP.EQ.1)THEN
          B2M(I)=2.D0*G(290+I)**2*(SUM+CM(I,1)*(TYUDHU-G(320+I)
     $           *TYUDYU)+CM(I,2)*(TYDDHD-G(320+I)*TYDDYD)+CM(I,3)
     $           *(TYEDHE-G(320+I)*TYEDYE))
        END IF
        F(320+I)=2.D0*G(290+I)**2/16.D0/PI**2*B1LPM(I)*G(320+I)
     $          +1.D0/(16.D0*PI**2)**2*B2M(I)
      END DO
!
!And here is the MSSM \mu
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1B1(I,J)=3.D0*YUYUD(I,J)+3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1B2(I,J)=6.D0*HUYUD(I,J)+6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2B1(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                  +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
            DUM2B2(I,J)=3.D0*HUYUDYUYUD(I,J)+3.D0*HDYDDYDYDD(I,J)
     $                  +HUYDDYDYUD(I,J)+HDYUDYUYDD(I,J)
     $                  +HEYEDYEYED(I,J)
            DUM2GRKMU(I,J)=3.D0*YUYUDYUYUD(I,J)+3.D0*YDYDDYDYDD(I,J)
     $                     +2.D0*YUYDDYDYUD(I,J)+YEYEDYEYED(I,J)
          ELSE
            DUM2B1(I,J)=0.D0
            DUM2B2(I,J)=0.D0
            DUM2GRKMU(I,J)=0.D0
          END IF
        END DO
      END DO
!
        BGRKMU=MVMU*(DTRACE(DUMGRKMU1)-3.D0*G(292)**2-3.D0/5.D0
     $           *G(291)**2)
        BETA1B=MVB*(DTRACE(DUM1B1)-3.D0*G(292)**2-3.D0/5.D0
     $         *G(291)**2)+MVMU*(DTRACE(DUM1B2)+6.D0*G(292)**2
     $         *G(322)+6.D0/5.D0*G(291)**2*G(321))
!
      IF(SW2LP.EQ.1)THEN
        BETA2GRKMU=MVMU*(-3.D0*DTRACE(DUM2GRKMU)+(16.D0*G(293)**2+4.D0
     $             /5.D0*G(291)**2)*TYUDYU+(16.D0*G(293)**2-2.D0/5.D0
     $             *G(291)**2)*TYDDYD+6.D0/5.D0*G(291)**2*TYEDYE
     $             +15.D0/2.D0*G(292)**4+9.D0/5.D0*G(291)**2*G(292)**2
     $             +207.D0/50.D0*G(291)**4)
        B2GRKMUM=BETA2GRKMU
        BETA2B=MVB*(-3.D0*DTRACE(DUM2B1)+(16.D0*G(293)**2+4.D0/5.D0
     $         *G(291)**2)*TYUDYU+(16.D0*G(293)**2-2.D0/5.D0*G(291)**2)
     $         *TYDDYD+6.D0/5.D0*G(291)**2*TYEDYE+15.D0/2.D0*G(292)**4
     $         +9.D0/5.D0*G(291)**2*G(292)**2+207.D0/50.D0*G(291)**4)
     $         +MVMU*(-12.D0*DTRACE(DUM2B2)+(32.D0*G(293)**2
     $         +8.D0/5.D0*G(291)**2)*THUYUD+(32.D0*G(293)**2-4.D0/5.D0
     $         *G(291)**2)*THDYDD+12.D0/5.D0*G(291)**2*THEYED-(32.D0
     $         *G(293)**2*G(323)+8.D0/5.D0*G(291)**2*G(321))*TYUDYU
     $         -(32.D0*G(293)**2*G(323)-4.D0/5.D0*G(291)**2*G(321))
     $         *TYDDYD-12.D0/5.D0*G(291)**2*G(321)*TYEDYE-30.D0
     $         *G(292)**4*G(322)-18.D0/5.D0*G(291)**2*G(292)**2
     $         *(G(321)+G(322))-414.D0/25.D0*G(291)**4*G(321))
         BETA2BM=BETA2B
      END IF
!
!The RKSTP compatible derivative is (with a minus sign
!to convert the running to BT notation)
!
        F(398)=1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2
     $         *BETA2GRKMU
        F(399)=-1.D0/16.D0/PI**2*BETA1B-1.D0/(16.D0*PI**2)**2*BETA2B
!
!Now the MSSM trilinear couplings - it should be ok to use the same
!dummy matrices, but the betas must be reset
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=0.D0
          BETA2U(I,J)=0.D0
          BETA1D(I,J)=0.D0
          BETA2D(I,J)=0.D0
          BETA1E(I,J)=0.D0
          BETA2E(I,J)=0.D0
        END DO
      END DO
!
!DUM2U11 is a dummy matrix for the calculation of DUM2U1
!
      DO I=1,3
        DO J=1,3
          DUM1D11(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1D12(I,J)=6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
          DUM1E11(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
          DUM1E21(I,J)=6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2U11(I,J)=3.D0*YUYUDYUYUD(I,J)+YUYDDYDYUD(I,J)
            DUM2U12(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2U21(I,J)=6.D0*HUYUDYUYUD(I,J)+HUYDDYDYUD(I,J)
     $                   +HDYUDYUYDD(I,J)
            DUM2U22(I,J)=6.D0*HDYDD(I,J)+2.D0*HEYED(I,J)
            DUM2U23(I,J)=6.D0*YDYDD(I,J)+2.D0*YEYED(I,J)
            DUM2D11(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)
     $                   +YEYEDYEYED(I,J)
            DUM2D12(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2D21(I,J)=6.D0*HDYDDYDYDD(I,J)+HUYDDYDYUD(I,J)
     $                   +HDYUDYUYDD(I,J)+2.D0*HEYEDYEYED(I,J)
            DUM2D22(I,J)=3.D0*HDYDD(I,J)+HEYED(I,J)
            DUM2D23(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E11(I,J)=3.D0*YDYDDYDYDD(I,J)+YUYDDYDYUD(I,J)
     $                   +YEYEDYEYED(I,J)
            DUM2E12(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E21(I,J)=6.D0*HDYDDYDYDD(I,J)+HUYDDYDYUD(I,J)
     $                   +HDYUDYUYDD(I,J)+2.D0*HEYEDYEYED(I,J)
            DUM2E22(I,J)=3.D0*YDYDD(I,J)+YEYED(I,J)
            DUM2E23(I,J)=3.D0*HDYDD(I,J)+HEYED(I,J)
          ELSE
            DUM2U11(I,J)=0.D0
            DUM2U12(I,J)=0.D0
            DUM2U21(I,J)=0.D0
            DUM2U22(I,J)=0.D0
            DUM2U23(I,J)=0.D0
            DUM2D11(I,J)=0.D0
            DUM2D12(I,J)=0.D0
            DUM2D21(I,J)=0.D0
            DUM2D22(I,J)=0.D0
            DUM2D23(I,J)=0.D0
            DUM2E11(I,J)=0.D0
            DUM2E12(I,J)=0.D0
            DUM2E21(I,J)=0.D0
            DUM2E22(I,J)=0.D0
            DUM2E23(I,J)=0.D0
          END IF
        END DO
      END DO
!
!DUM1U1 is a dummy matrix for the calculation of BETA1U
!
      DO I=1,3
        DO J=1,3
          DUM1U1(I,J)=3.D0*TYUDYU*ID(I,J)+5.D0*YUDYU(I,J)+YDDYD(I,J)
     $                +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2
     $                -13.D0/15.D0*G(291)**2)*ID(I,J)
          DUM1U2(I,J)=6.D0*THUYUD*ID(I,J)+4.D0*YUDHU(I,J)+2.D0
     $                *YDDHD(I,J)+(32.D0/3.D0*G(293)**2*G(323)+6.D0
     $                *G(292)**2*G(322)+26.D0/15.D0*G(291)**2*G(321))
     $                *ID(I,J)
          DUM1D1(I,J)=DTRACE(DUM1D11)*ID(I,J)+5.D0*YDDYD(I,J)+YUDYU(I,J)
     $                +(-16.D0/3.D0*G(293)**2-3.D0*G(292)**2-7.D0/15.D0
     $                *G(291)**2)*ID(I,J)
          DUM1D2(I,J)=DTRACE(DUM1D12)*ID(I,J)+4.D0*YDDHD(I,J)+2.D0
     $                *YUDHU(I,J)+(32.D0/3.D0*G(293)**2*G(323)+6.D0
     $                *G(292)**2*G(322)+14.D0/15.D0*G(291)**2*G(321))
     $                *ID(I,J)
          DUM1E1(I,J)=DTRACE(DUM1E11)*ID(I,J)+5.D0*YEDYE(I,J)+(-3.D0
     $                *G(292)**2-9.D0/5.D0*G(291)**2)*ID(I,J)
          DUM1E2(I,J)=DTRACE(DUM1E21)*ID(I,J)+4.D0*YEDHE(I,J)+(6.D0
     $                *G(292)**2*G(322)+18.D0/5.D0*G(291)**2*G(321))
     $                *ID(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2U1(I,J)=-3.D0*DTRACE(DUM2U11)*ID(I,J)-YDDYD(I,J)
     $                  *DTRACE(DUM2U12)-15.D0*YUDYU(I,J)*TYUDYU
     $                  -6.D0*YUDYUYUDYU(I,J)-2.D0*YDDYDYDDYD(I,J)
     $                  -4.D0*YDDYDYUDYU(I,J)+(16.D0*G(293)**2+4.D0/5.D0
     $                  *G(291)**2)*TYUDYU*ID(I,J)+12.D0*G(292)**2
     $                  *YUDYU(I,J)+2.D0/5.D0*G(291)**2*YDDYD(I,J)
     $                  +(-16.D0/9.D0*G(293)**4+8.D0*G(293)**2
     $                  *G(292)**2+136.D0/45.D0*G(293)**2*G(291)**2
     $                  +15.D0/2.D0*G(292)**4+G(292)**2*G(291)**2
     $                  +2743.D0/450.D0*G(291)**4)*ID(I,J)
            DUM2U2(I,J)=-6.D0*DTRACE(DUM2U21)*ID(I,J)-18.D0
     $                  *YUDYU(I,J)*THUYUD-YDDYD(I,J)*DTRACE(DUM2U22)
     $                  -12.D0*YUDHU(I,J)*TYUDYU-YDDHD(I,J)
     $                  *DTRACE(DUM2U23)-6.D0*YUDYUYUDHU(I,J)
     $                  -8.D0*YUDHUYUDYU(I,J)-4.D0*YDDYDYDDHD(I,J)-4.D0
     $                  *YDDHDYDDYD(I,J)-2.D0*YDDYDYUDHU(I,J)-4.D0
     $                  *YDDHDYUDYU(I,J)+(32.D0*G(293)**2+8.D0/5.D0
     $                  *G(291)**2)*THUYUD*ID(I,J)+(6.D0*G(292)**2
     $                  +6.D0/5.D0*G(291)**2)*YUDHU(I,J)+4.D0/5.D0
     $                  *G(291)**2*YDDHD(I,J)-(32.D0*G(293)**2*G(323)
     $                  +8.D0/5.D0*G(291)**2*G(321))*TYUDYU*ID(I,J)
     $                  -(12.D0*G(292)**2*G(322)+4.D0/5.D0*G(291)**2
     $                  *G(321))*YUDYU(I,J)-4.D0/5.D0*G(291)**2
     $                  *G(321)*YDDYD(I,J)+(64.D0/9.D0*G(293)**4
     $                  *G(323)-16.D0*G(293)**2*G(292)**2*(G(323)
     $                  +G(322))-272.D0/45.D0*G(293)**2*G(291)**2
     $                  *(G(323)+G(321))-30.D0*G(292)**4*G(322)-2.D0
     $                  *G(292)**2*G(291)**2*(G(322)+G(321))
     $                  -5486.D0/225.D0*G(291)**4*G(321))*ID(I,J)      
            DUM2D1(I,J)=-3.D0*DTRACE(DUM2D11)*ID(I,J)-3.D0*YUDYU(I,J)
     $                  *TYUDYU-5.D0*YDDYD(I,J)*DTRACE(DUM2D12)-6.D0
     $                  *YDDYDYDDYD(I,J)-2.D0*YUDYUYUDYU(I,J)-4.D0
     $                  *YUDYUYDDYD(I,J)+(16.D0*G(293)**2-2.D0/5.D0
     $                  *G(291)**2)*TYDDYD*ID(I,J)+6.D0/5.D0*G(291)**2
     $                  *TYEDYE*ID(I,J)+4.D0/5.D0*G(291)**2*YUDYU(I,J)
     $                  +(12.D0*G(292)**2+6.D0/5.D0*G(291)**2)
     $                  *YDDYD(I,J)+(-16.D0/9.D0*G(293)**4+8.D0
     $                  *G(293)**2*G(292)**2+8.D0/9.D0*G(293)**2
     $                  *G(291)**2+15.D0/2.D0*G(292)**4+G(292)**2
     $                  *G(291)**2+287.D0/90.D0*G(291)**4)*ID(I,J)
            DUM2D2(I,J)=-6.D0*DTRACE(DUM2D21)*ID(I,J)-6.D0*YUDYU(I,J)
     $                  *THUYUD-6.D0*YDDYD(I,J)*DTRACE(DUM2D22)-6.D0
     $                  *YUDHU(I,J)*TYUDYU-4.D0*YDDHD(I,J)
     $                  *DTRACE(DUM2D23)-6.D0*YDDYDYDDHD(I,J)-8.D0
     $                  *YDDHDYDDYD(I,J)-4.D0*YUDHUYUDYU(I,J)-4.D0
     $                  *YUDYUYUDHU(I,J)-4.D0*YUDHUYDDYD(I,J)-2.D0
     $                  *YUDYUYDDHD(I,J)+(32.D0*G(293)**2-4.D0/5.D0
     $                  *G(291)**2)*THDYDD*ID(I,J)+12.D0/5.D0*G(291)**2
     $                  *THEYED*ID(I,J)+8.D0/5.D0*G(291)**2*YUDHU(I,J)
     $                  +(6.D0*G(292)**2+6.D0/5.D0*G(291)**2)*YDDHD(I,J)
     $                  -(32.D0*G(293)**2*G(323)-4.D0/5.D0*G(291)**2
     $                  *G(321))*TYDDYD*ID(I,J)-12.D0/5.D0*G(291)**2
     $                  *G(321)*TYEDYE*ID(I,J)-(12.D0*G(292)**2*G(322)
     $                  +8.D0/5.D0*G(291)**2*G(321))*YDDYD(I,J)
     $                  -8.D0/5.D0*G(291)**2*G(321)*YUDYU(I,J)
     $                  +(64.D0/9.D0*G(293)**4*G(323)-16.D0*G(293)**2
     $                  *G(292)**2*(G(323)+G(322))-16.D0/9.D0
     $                  *G(293)**2*G(291)**2*(G(323)+G(321))-30.D0
     $                  *G(292)**4*G(322)-2.D0*G(292)**2*G(291)**2
     $                  *(G(322)+G(321))-574.D0/45.D0*G(291)**4*G(321))
     $                  *ID(I,J)
            DUM2E1(I,J)=-3.D0*DTRACE(DUM2E11)*ID(I,J)-5.D0*YEDYE(I,J)
     $                  *DTRACE(DUM2E12)-6.D0*YEDYEYEDYE(I,J)+(16.D0
     $                  *G(293)**2-2.D0/5.D0*G(291)**2)*TYDDYD*ID(I,J)
     $                  +6.D0/5.D0*G(291)**2*TYEDYE*ID(I,J)+(12.D0
     $                  *G(292)**2-6.D0/5.D0*G(291)**2)*YEDYE(I,J)
     $                  +(15.D0/2.D0*G(292)**4+9.D0/5.D0*G(292)**2
     $                  *G(291)**2+27.D0/2.D0*G(291)**4)*ID(I,J)
            DUM2E2(I,J)=-6.D0*DTRACE(DUM2E21)*ID(I,J)-4.D0*YEDHE(I,J)
     $                  *DTRACE(DUM2E22)-6.D0*YEDYE(I,J)*DTRACE(DUM2E23)
     $                  -6.D0*YEDYEYEDHE(I,J)-8.D0*YEDHEYEDYE(I,J)
     $                  +(32.D0*G(293)**2-4.D0/5.D0*G(291)**2)*THDYDD
     $                  *ID(I,J)+12.D0/5.D0*G(291)**2*THEYED*ID(I,J)
     $                  +(6.D0*G(292)**2+6.D0/5.D0*G(291)**2)*YEDHE(I,J)
     $                  -(32.D0*G(293)**2*G(323)-4.D0/5.D0*G(291)**2
     $                  *G(321))*TYDDYD*ID(I,J)-12.D0/5.D0*G(291)**2
     $                  *G(321)*TYEDYE*ID(I,J)-12.D0*G(292)**2*G(322)
     $                  *YEDYE(I,J)+(-30.D0*G(292)**4*G(322)-18.D0/5.D0
     $                  *G(292)**2*G(291)**2*(G(321)+G(322))-54.D0
     $                  *G(291)**4*G(321))*ID(I,J)
          ELSE
            DUM2U1(I,J)=0.D0
            DUM2U2(I,J)=0.D0
            DUM2D1(I,J)=0.D0
            DUM2D2(I,J)=0.D0
            DUM2E1(I,J)=0.D0
            DUM2E2(I,J)=0.D0
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the trilinears
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=DMATMUL(0,HU,DUM1U1,I,J)+DMATMUL(0,YU,DUM1U2,I,J)
          BETA1D(I,J)=DMATMUL(0,HD,DUM1D1,I,J)+DMATMUL(0,YD,DUM1D2,I,J)
          BETA1E(I,J)=DMATMUL(0,HE,DUM1E1,I,J)+DMATMUL(0,YE,DUM1E2,I,J)
          IF(SW2LP.EQ.1)THEN
            BETA2U(I,J)=DMATMUL(0,HU,DUM2U1,I,J)
     $                  +DMATMUL(0,YU,DUM2U2,I,J)
            BETA2D(I,J)=DMATMUL(0,HD,DUM2D1,I,J)
     $                  +DMATMUL(0,YD,DUM2D2,I,J)
            BETA2E(I,J)=DMATMUL(0,HE,DUM2E1,I,J)
     $                  +DMATMUL(0,YE,DUM2E2,I,J)
            B2HMU(I,J)=BETA2U(I,J)
            B2HMD(I,J)=BETA2D(I,J)
            B2HME(I,J)=BETA2E(I,J)
          END IF
!
!Calculate the differentials DH
!
          DH(1,I,J)=1.D0/16.D0/PI**2*BETA1U(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2U(I,J)
          DH(2,I,J)=1.D0/16.D0/PI**2*BETA1D(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2D(I,J)
          DH(3,I,J)=1.D0/16.D0/PI**2*BETA1E(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2E(I,J)
        END DO
      END DO
!
!Convert into form readable by RKSTP (converting to BT notation...)
!Need new do loop due to the (J,I) ordering
!
      DO I=1,3
        DO J=1,3
          F(323+(I-1)*3+J)=-DH(1,J,I)
          F(332+(I-1)*3+J)=-DH(2,J,I)
          F(341+(I-1)*3+J)=-DH(3,J,I)
        END DO
      END DO
!
!Now for the (mass)^2 parameters. DUMSP are dummy matrices for the
!calculation of S'
!
      IF(SW2LP.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            DUMSP1(I,J)=3.D0*G(351)*ID(I,J)+MQM(I,J)
            DUMSP2(I,J)=3.D0*G(352)*ID(I,J)-MQM(I,J)
            DUMSP3(I,J)=G(352)*ID(I,J)+MLM(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            DUMSP1(I,J)=0.D0
            DUMSP2(I,J)=0.D0
            DUMSP3(I,J)=0.D0
          END DO
        END DO
      END IF
!
!Dummy matrices for S, S' and SIGMAi
!
      DO I=1,3
        DO J=1,3
          DUMS(I,J)=MQM(I,J)-MLM(I,J)-2.D0*MUPM(I,J)+MDM(I,J)+MEM(I,J)
          IF(SW2LP.EQ.1)THEN
            DUMSP(I,J)=-DMATMUL(0,DUMSP1,YUDYU,I,J)+4.D0*YUDMUPMYU(I,J)
     $                 +DMATMUL(0,DUMSP2,YDDYD,I,J)-2.D0*YDDMDMYD(I,J)
     $                 +DMATMUL(0,DUMSP3,YEDYE,I,J)-2.D0*YEDMEMYE(I,J)
            DUMSIG1(I,J)=MQM(I,J)+3.D0*MLM(I,J)+8.D0*MUPM(I,J)
     $                   +2.D0*MDM(I,J)+6.D0*MEM(I,J)
            DUMSIG2(I,J)=3.D0*MQM(I,J)+MLM(I,J)
            DUMSIG3(I,J)=2.D0*MQM(I,J)+MUPM(I,J)+MDM(I,J)
          ELSE
            DUMSP(I,J)=0.D0
            DUMSIG1(I,J)=0.D0
            DUMSIG2(I,J)=0.D0
            DUMSIG3(I,J)=0.D0
          END IF
        END DO
      END DO
!
!These definitions make the equations shorter
!
      S=G(351)-G(352)+DTRACE(DUMS)
      IF(SW2LP.EQ.1)THEN
        SP=DTRACE(DUMSP)+(3.D0/2.D0*G(292)**2+3.D0/10.D0*G(291)**2)
     $     *(G(351)-G(352)-DTRACE(MLM))+(8.D0/3.D0*G(293)**2+3.D0/2.D0
     $     *G(292)**2+1.D0/30.D0*G(291)**2)*DTRACE(MQM)-(16.D0/3.D0
     $     *G(293)**2+16.D0/15.D0*G(291)**2)*DTRACE(MUPM)+(8.D0/3.D0
     $     *G(293)**2+2.D0/15.D0*G(291)**2)*DTRACE(MDM)+6.D0/5.D0
     $     *G(291)**2*DTRACE(MEM)
        SIG1=1.D0/5.D0*G(291)**2*(3.D0*(G(351)+G(352))
     $       +DTRACE(DUMSIG1))
        SIG2=G(292)**2*(G(351)+G(352)+DTRACE(DUMSIG2))
        SIG3=G(293)**2*DTRACE(DUMSIG3)
      ELSE
        SP=0.D0
        SIG1=0.D0
        SIG2=0.D0
        SIG3=0.D0
      END IF
!
!Now the dummy matrices for the higgs betas
!
      IF(SW2LP.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            DUM2HU11(I,J)=G(351)*ID(I,J)+MQM(I,J)
            DUM2HU12(I,J)=(G(351)+G(352))*ID(I,J)+MQM(I,J)
            DUM2HU21(I,J)=DUM2HU11(I,J)
            DUM2HD11(I,J)=G(352)*ID(I,J)+MQM(I,J)
            DUM2HD12(I,J)=DUM2HU12(I,J)
            DUM2HD13(I,J)=G(352)*ID(I,J)+MLM(I,J)
            DUM2HD21(I,J)=DUM2HD11(I,J)
            DUM2HD31(I,J)=DUM2HD13(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            DUM2HU11(I,J)=0.D0
            DUM2HU12(I,J)=0.D0
            DUM2HU21(I,J)=0.D0
            DUM2HD11(I,J)=0.D0
            DUM2HD12(I,J)=0.D0
            DUM2HD13(I,J)=0.D0
            DUM2HD21(I,J)=0.D0
            DUM2HD31(I,J)=0.D0
          END DO
        END DO
      END IF
      DO I=1,3
        DO J=1,3
          DUM1HU1(I,J)=G(351)*ID(I,J)+MQM(I,J)
          DUM1HD1(I,J)=G(352)*ID(I,J)+MQM(I,J)
          DUM1HD2(I,J)=G(352)*ID(I,J)+MLM(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2HU1(I,J)=6.D0*DMATMUL(0,DUM2HU11,YUDYUYUDYU,I,J)+6.D0
     $                   *YUDMUPMYUYUDYU(I,J)
     $                   +DMATMUL(0,DUM2HU12,YUDYUYDDYD,I,J)
     $                   +YUDMUPMYUYDDYD(I,J)+YUDYUMQMYDDYD(I,J)
     $                   +YUDYUYDDMDMYD(I,J)+6.D0*HUDHUYUDYU(I,J)+6.D0
     $                   *HUDYUYUDHU(I,J)+HDDHDYUDYU(I,J)
     $                   +YDDYDHUDHU(I,J)+HDDYDYUDHU(I,J)
     $                   +YDDHDHUDYU(I,J)
            DUM2HU2(I,J)=DMATMUL(0,DUM2HU21,YUDYU,I,J)+YUDMUPMYU(I,J)
     $                   +HUDHU(I,J)
            DUM2HD1(I,J)=6.D0*DMATMUL(0,DUM2HD11,YDDYDYDDYD,I,J)
     $                   +6.D0*YDDMDMYDYDDYD(I,J)
     $                   +DMATMUL(0,DUM2HD12,YUDYUYDDYD,I,J)
     $                   +YUDMUPMYUYDDYD(I,J)+YUDYUMQMYDDYD(I,J)
     $                   +YUDYUYDDMDMYD(I,J)+2.D0
     $                   *DMATMUL(0,DUM2HD13,YEDYEYEDYE,I,J)+2.D0
     $                   *YEDMEMYEYEDYE(I,J)+6.D0*HDDHDYDDYD(I,J)+6.D0
     $                   *HDDYDYDDHD(I,J)+HUDHUYDDYD(I,J)
     $                   +YUDYUHDDHD(I,J)+HUDYUYDDHD(I,J)
     $                   +YUDHUHDDYD(I,J)+2.D0*HEDHEYEDYE(I,J)+2.D0
     $                   *HEDYEYEDHE(I,J)
            DUM2HD2(I,J)=DMATMUL(0,DUM2HD21,YDDYD,I,J)+YDDMDMYD(I,J)
     $                   +HDDHD(I,J)
            DUM2HD3(I,J)=DMATMUL(0,DUM2HD31,YEDYE,I,J)+YEDMEMYE(I,J)
     $                   +HEDHE(I,J)
          ELSE
            DUM2HU1(I,J)=0.D0
            DUM2HU2(I,J)=0.D0
            DUM2HD1(I,J)=0.D0
            DUM2HD2(I,J)=0.D0
            DUM2HD3(I,J)=0.D0
          END IF
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          DUM1HU(I,J)=DMATMUL(0,DUM1HU1,YUDYU,I,J)+YUDMUPMYU(I,J)
     $                +HUDHU(I,J)
          DUM1HD(I,J)=6.D0*DMATMUL(0,DUM1HD1,YDDYD,I,J)+6.D0
     $                *YDDMDMYD(I,J)+2.D0*DMATMUL(0,DUM1HD2,YEDYE,I,J)
     $                +2.D0*YEDMEMYE(I,J)+6.D0*HDDHD(I,J)
     $                +2.D0*HEDHE(I,J)
        END DO
      END DO
!
!The higgs mass^2 beta parameters are not matrices:
!
      BETA1HU=6.D0*DTRACE(DUM1HU)-6.D0*G(292)**2*MODSQ(G(322))
     $        -6.D0/5.D0*G(291)**2*MODSQ(G(321))+3.D0/5.D0*G(291)**2*S
      BETA1HD=DTRACE(DUM1HD)-6.D0*G(292)**2*MODSQ(G(322))-6.D0/5.D0
     $        *G(291)**2*MODSQ(G(321))-3.D0/5.D0*G(291)**2*S
      IF(SW2LP.EQ.1)THEN
        BETA2HU=-6.D0*DTRACE(DUM2HU1)+(32.D0*G(293)**2+8.D0/5.D0
     $          *G(291)**2)*DTRACE(DUM2HU2)+32.D0*G(293)**2*(2.D0
     $          *MODSQ(G(323))*TYUDYU-CON(G(323))*TYUDHU-G(323)
     $          *THUDYU)+8.D0/5.D0*G(291)**2*(2.D0*MODSQ(G(321))
     $          *TYUDYU-CON(G(321))*TYUDHU-G(321)*THUDYU)+6.D0/5.D0
     $          *G(291)**2*SP+33.D0*G(292)**4*MODSQ(G(322))
     $          +18.D0/5.D0*G(292)**2*G(291)**2*(MODSQ(G(322))
     $          +MODSQ(G(321))+DRE(2,G(321),G(322)))+621.D0/25.D0
     $          *G(291)**4*MODSQ(G(321))+3.D0*G(292)**2*SIG2
     $          +3.D0/5.D0*G(291)**2*SIG1
        BETA2HD=-6.D0*DTRACE(DUM2HD1)+(32.D0*G(293)**2-4.D0/5.D0
     $          *G(291)**2)*DTRACE(DUM2HD2)+32.D0*G(293)**2*(2.D0
     $          *MODSQ(G(323))*TYDDYD-CON(G(323))*TYDDHD-G(323)
     $          *THDDYD)-4.D0/5.D0*G(291)**2*(2.D0*MODSQ(G(321))
     $          *TYDDYD-CON(G(321))*TYDDHD-G(321)*THDDYD)+12.D0/5.D0
     $          *G(291)**2*(DTRACE(DUM2HD3)+2.D0*MODSQ(G(321))*TYEDYE
     $          -G(321)*THEDYE-CON(G(321))*TYEDHE)-6.D0/5.D0
     $          *G(291)**2*SP+33.D0*G(292)**4*MODSQ(G(322))
     $          +18.D0/5.D0*G(292)**2*G(291)**2*(MODSQ(G(322))
     $          +MODSQ(G(321))+DRE(2,G(321),G(322)))+621.D0/25.D0
     $          *G(291)**4*MODSQ(G(321))+3.D0*G(292)**2*SIG2
     $          +3.D0/5.D0*G(291)**2*SIG1
         B2HUM=BETA2HU
         B2HDM=BETA2HD
      END IF
!
!Calculate the differentials F(351) and F(352)
!
      F(351)=1.D0/16.D0/PI**2*BETA1HU+1.D0/(16.D0*PI**2)**2*BETA2HU
      F(352)=1.D0/16.D0/PI**2*BETA1HD+1.D0/(16.D0*PI**2)**2*BETA2HD
!
!Now for the dummy matrices for the rest of the mass^2
!
      IF(SW2LP.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            DUM2Q31(I,J)=MQM(I,J)+4.D0*G(351)*ID(I,J)
            DUM2Q41(I,J)=MQM(I,J)+4.D0*G(352)*ID(I,J)
            DUM2L21(I,J)=MLM(I,J)+4.D0*G(352)*ID(I,J)
            DUM2U31(I,J)=MUPM(I,J)+4.D0*G(351)*ID(I,J)
            DUM2D31(I,J)=MDM(I,J)+4.D0*G(352)*ID(I,J)
            DUM2E21(I,J)=MEM(I,J)+4.D0*G(352)*ID(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            DUM2Q31(I,J)=0.D0
            DUM2Q41(I,J)=0.D0
            DUM2L21(I,J)=0.D0
            DUM2U31(I,J)=0.D0
            DUM2D31(I,J)=0.D0
            DUM2E21(I,J)=0.D0
          END DO
        END DO
      END IF
!
!This is the second round of dummy matrices
!
      DO I=1,3
        DO J=1,3
          DUM1Q1(I,J)=MQM(I,J)+2.D0*G(351)*ID(I,J)
          DUM1Q2(I,J)=MQM(I,J)+2.D0*G(352)*ID(I,J)
          DUM1Q3(I,J)=YUDYU(I,J)+YDDYD(I,J)
          DUM1L1(I,J)=MLM(I,J)+2.D0*G(352)*ID(I,J)
          DUM1U1(I,J)=2.D0*MUPM(I,J)+4.D0*G(351)*ID(I,J)
          DUM1D1(I,J)=2.D0*MDM(I,J)+4.D0*G(352)*ID(I,J)
          DUM1E1(I,J)=2.D0*MEM(I,J)+4.D0*G(352)*ID(I,J)
          IF(SW2LP.EQ.1)THEN
            DUM2Q1(I,J)=2.D0*DUM2Q31(I,J)
            DUM2Q2(I,J)=2.D0*DUM2Q41(I,J)
            DUM2Q3(I,J)=DMATMUL(0,DUM2Q31,YUDYU,I,J)+2.D0
     $                  *YUDMUPMYU(I,J)+YUDYUMQM(I,J)
            DUM2Q4(I,J)=DMATMUL(0,DUM2Q41,YDDYD,I,J)+2.D0
     $                  *YDDMDMYD(I,J)+YDDYDMQM(I,J)
            DUM2Q5(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2Q6(I,J)=MQMYUDYU(I,J)+YUDMUPMYU(I,J)
            DUM2Q7(I,J)=6.D0*MQMYDDYD(I,J)+6.D0*YDDMDMYD(I,J)+2.D0
     $                  *MLMYEDYE(I,J)+2.D0*YEDMEMYE(I,J)
            DUM2Q8(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2Q9(I,J)=6.D0*HDDHD(I,J)+2.D0*HEDHE(I,J)
            DUM2QA(I,J)=6.D0*YDDHD(I,J)+2.D0*YEDHE(I,J)
            DUM2QB(I,J)=6.D0*HDDYD(I,J)+2.D0*HEDYE(I,J)
            DUM2QC(I,J)=2.D0*DUM1Q1(I,J)
            DUM2QD(I,J)=MQM(I,J)+2.D0*G(352)*ID(I,J)
            DUM2L1(I,J)=2.D0*DUM2L21(I,J)
            DUM2L2(I,J)=DMATMUL(0,DUM2L21,YEDYE,I,J)+2.D0
     $                  *YEDMEMYE(I,J)+YEDYEMLM(I,J)
            DUM2L3(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2L4(I,J)=6.D0*MQMYDDYD(I,J)+6.D0*YDDMDMYD(I,J)+2.D0
     $                  *MLMYEDYE(I,J)+2.D0*YEDMEMYE(I,J)
            DUM2L5(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2L6(I,J)=6.D0*HDDHD(I,J)+2.D0*HEDHE(I,J)
            DUM2L7(I,J)=6.D0*YDDHD(I,J)+2.D0*YEDHE(I,J)
            DUM2L8(I,J)=6.D0*HDDYD(I,J)+2.D0*HEDYE(I,J)
            DUM2L9(I,J)=DUM1L1(I,J)
            DUM2U1(I,J)=2.D0*DUM2U31(I,J)
            DUM2U2(I,J)=2.D0*MUPM(I,J)+4.D0*(G(351)+G(352))*ID(I,J)
            DUM2U3(I,J)=DMATMUL(0,DUM2U31,YUYUD,I,J)+2.D0*YUMQMYUD(I,J)
     $                  +YUYUDMUPM(I,J)
            DUM2U4(I,J)=MQMYUDYU(I,J)+YUDMUPMYU(I,J)
            DUM2U5(I,J)=MUPM(I,J)+2.D0*G(351)*ID(I,J)
            DUM2D1(I,J)=2.D0*DUM2D31(I,J)
            DUM2D2(I,J)=2.D0*MDM(I,J)+4.D0*(G(351)+G(352))*ID(I,J)
            DUM2D3(I,J)=DMATMUL(0,DUM2D31,YDYDD,I,J)+2.D0*YDMQMYDD(I,J)
     $                  +YDYDDMDM(I,J)
            DUM2D4(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2D5(I,J)=3.D0*MQMYDDYD(I,J)+3.D0*YDDMDMYD(I,J)
     $                  +MLMYEDYE(I,J)+YEDMEMYE(I,J)
            DUM2D6(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2D7(I,J)=3.D0*HDDHD(I,J)+HEDHE(I,J)
            DUM2D8(I,J)=3.D0*HDDYD(I,J)+HEDYE(I,J)
            DUM2D9(I,J)=3.D0*YDDHD(I,J)+YEDHE(I,J)
            DUM2DA(I,J)=MDM(I,J)+2.D0*G(352)*ID(I,J)
            DUM2E1(I,J)=2.D0*DUM2E21(I,J)
            DUM2E2(I,J)=DMATMUL(0,DUM2E21,YEYED,I,J)+2.D0*YEMLMYED(I,J)
     $                  +YEYEDMEM(I,J)
            DUM2E3(I,J)=6.D0*YDDYD(I,J)+2.D0*YEDYE(I,J)
            DUM2E4(I,J)=3.D0*MQMYDDYD(I,J)+3.D0*YDDMDMYD(I,J)
     $                  +MLMYEDYE(I,J)+YEDMEMYE(I,J)
            DUM2E5(I,J)=3.D0*YDDYD(I,J)+YEDYE(I,J)
            DUM2E6(I,J)=3.D0*HDDHD(I,J)+HEDHE(I,J)
            DUM2E7(I,J)=3.D0*HDDYD(I,J)+HEDYE(I,J)
            DUM2E8(I,J)=3.D0*YDDHD(I,J)+YEDHE(I,J)
            DUM2E9(I,J)=MEM(I,J)+2.D0*G(352)*ID(I,J)
          ELSE
            DUM2Q1(I,J)=0.D0
            DUM2Q2(I,J)=0.D0
            DUM2Q3(I,J)=0.D0
            DUM2Q4(I,J)=0.D0
            DUM2Q5(I,J)=0.D0
            DUM2Q6(I,J)=0.D0
            DUM2Q7(I,J)=0.D0
            DUM2Q8(I,J)=0.D0
            DUM2Q9(I,J)=0.D0
            DUM2QA(I,J)=0.D0
            DUM2QB(I,J)=0.D0
            DUM2QC(I,J)=0.D0
            DUM2QD(I,J)=0.D0
            DUM2L1(I,J)=0.D0
            DUM2L2(I,J)=0.D0
            DUM2L3(I,J)=0.D0
            DUM2L4(I,J)=0.D0
            DUM2L5(I,J)=0.D0
            DUM2L6(I,J)=0.D0
            DUM2L7(I,J)=0.D0
            DUM2L8(I,J)=0.D0
            DUM2L9(I,J)=0.D0
            DUM2U1(I,J)=0.D0
            DUM2U2(I,J)=0.D0
            DUM2U3(I,J)=0.D0
            DUM2U4(I,J)=0.D0
            DUM2U5(I,J)=0.D0
            DUM2D1(I,J)=0.D0
            DUM2D2(I,J)=0.D0
            DUM2D3(I,J)=0.D0
            DUM2D4(I,J)=0.D0
            DUM2D5(I,J)=0.D0
            DUM2D6(I,J)=0.D0
            DUM2D7(I,J)=0.D0
            DUM2D8(I,J)=0.D0
            DUM2D9(I,J)=0.D0
            DUM2DA(I,J)=0.D0
            DUM2E1(I,J)=0.D0
            DUM2E2(I,J)=0.D0
            DUM2E3(I,J)=0.D0
            DUM2E4(I,J)=0.D0
            DUM2E5(I,J)=0.D0
            DUM2E6(I,J)=0.D0
            DUM2E7(I,J)=0.D0
            DUM2E8(I,J)=0.D0
            DUM2E9(I,J)=0.D0
          END IF
        END DO
      END DO
!
!Now calculate the beta functions for the rest of the mass^2
!
      DO I=1,3
        DO J=1,3
          BETA1MQ(I,J)=DMATMUL(0,DUM1Q1,YUDYU,I,J)
     $                 +DMATMUL(0,DUM1Q2,YDDYD,I,J)
     $                 +DMATMUL(0,DUM1Q3,MQM,I,J)+2.D0
     $                 *YUDMUPMYU(I,J)+2.D0*YDDMDMYD(I,J)+2.D0
     $                 *HUDHU(I,J)+2.D0*HDDHD(I,J)+(-32.D0/3.D0
     $                 *G(293)**2*MODSQ(G(323))-6.D0*G(292)**2
     $                 *MODSQ(G(322))-2.D0/15.D0*G(291)**2
     $                 *MODSQ(G(321))+1.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1ML(I,J)=DMATMUL(0,DUM1L1,YEDYE,I,J)+2.D0*YEDMEMYE(I,J)
     $                 +YEDYEMLM(I,J)+2.D0*HEDHE(I,J)+(-6.D0*G(292)**2
     $                 *MODSQ(G(322))-6.D0/5.D0*G(291)**2
     $                 *MODSQ(G(321))-3.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1MU(I,J)=DMATMUL(0,DUM1U1,YUYUD,I,J)+4.D0*YUMQMYUD(I,J)
     $                 +2.D0*YUYUDMUPM(I,J)+4.D0*HUHUD(I,J)+(-32.D0/3.D0
     $                 *G(293)**2*MODSQ(G(323))-32.D0/15.D0*G(291)**2
     $                 *MODSQ(G(321))-4.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1MD(I,J)=DMATMUL(0,DUM1D1,YDYDD,I,J)+4.D0*YDMQMYDD(I,J)
     $                 +2.D0*YDYDDMDM(I,J)+4.D0*HDHDD(I,J)+(-32.D0/3.D0
     $                 *G(293)**2*MODSQ(G(323))-8.D0/15.D0*G(291)**2
     $                 *MODSQ(G(321))+2.D0/5.D0*G(291)**2*S)*ID(I,J)
          BETA1ME(I,J)=DMATMUL(0,DUM1E1,YEYED,I,J)+4.D0*YEMLMYED(I,J)
     $                 +2.D0*YEYEDMEM(I,J)+4.D0*HEHED(I,J)+(-24.D0/5.D0
     $                 *G(291)**2*MODSQ(G(321))+6.D0/5.D0*G(291)**2*S)
     $                 *ID(I,J)
         IF(SW2LP.EQ.1)THEN
            BETA2MQ(I,J)=-DMATMUL(0,DUM2Q1,YUDYUYUDYU,I,J)-4.D0
     $                   *YUDMUPMYUYUDYU(I,J)-4.D0*YUDYUMQMYUDYU(I,J)
     $                   -4.D0*YUDYUYUDMUPMYU(I,J)
     $                   -2.D0*YUDYUYUDYUMQM(I,J)
     $                   -DMATMUL(0,DUM2Q2,YDDYDYDDYD,I,J)
     $                   -4.D0*YDDMDMYDYDDYD(I,J)
     $                   -4.D0*YDDYDMQMYDDYD(I,J)
     $                   -4.D0*YDDYDYDDMDMYD(I,J)
     $                   -2.D0*YDDYDYDDYDMQM(I,J)
     $                   -DUM2Q3(I,J)*3.D0*TYUDYU-DUM2Q4(I,J)
     $                   *DTRACE(DUM2Q5)-6.D0*YUDYU(I,J)*DTRACE(DUM2Q6)
     $                   -YDDYD(I,J)*DTRACE(DUM2Q7)-4.D0
     $                   *(YUDYUHUDHU(I,J)+HUDHUYUDYU(I,J)
     $                   +YUDHUHUDYU(I,J)+HUDYUYUDHU(I,J))-4.D0
     $                   *(YDDYDHDDHD(I,J)+HDDHDYDDYD(I,J)
     $                   +YDDHDHDDYD(I,J)+HDDYDYDDHD(I,J))-HUDHU(I,J)
     $                   *6.D0*TYUDYU-YUDYU(I,J)*6.D0*THUDHU-HUDYU(I,J)
     $                   *6.D0*TYUDHU-YUDHU(I,J)*6.D0*THUDYU-HDDHD(I,J)
     $                   *DTRACE(DUM2Q8)-YDDYD(I,J)*DTRACE(DUM2Q9)
     $                   -HDDYD(I,J)*DTRACE(DUM2QA)-YDDHD(I,J)
     $                   *DTRACE(DUM2QB)+2.D0/5.D0*G(291)**2
     $                   *(DMATMUL(0,DUM2QC,YUDYU,I,J)+4.D0
     $                   *YUDMUPMYU(I,J)+2.D0*YUDYUMQM(I,J)+4.D0
     $                   *HUDHU(I,J)-4.D0*G(321)*HUDYU(I,J)-4.D0
     $                   *CON(G(321))*YUDHU(I,J)+8.D0*MODSQ(G(321))
     $                   *YUDYU(I,J)+DMATMUL(0,DUM2QD,YDDYD,I,J)+2.D0
     $                   *YDDMDMYD(I,J)+YDDYDMQM(I,J)+2.D0*HDDHD(I,J)
     $                   -2.D0*G(321)*HDDYD(I,J)-2.D0*CON(G(321))
     $                   *YDDHD(I,J)+4.D0*MODSQ(G(321))*YDDYD(I,J))
     $                   +(2.D0/5.D0*G(291)**2*SP-128.D0/3.D0*G(293)**4
     $                   *MODSQ(G(323))+32.D0*G(293)**2*G(292)**2
     $                   *(MODSQ(G(323))+MODSQ(G(322))
     $                   +DRE(2,G(322),G(323)))+32.D0/45.D0*G(293)**2
     $                   *G(291)**2*(MODSQ(G(323))+MODSQ(G(321))
     $                   +DRE(2,G(321),G(323)))+33.D0*G(292)**4
     $                   *MODSQ(G(322))+2.D0/5.D0*G(292)**2
     $                   *G(291)**2*(MODSQ(G(322))+MODSQ(G(321))
     $                   +DRE(2,G(321),G(322)))+199.D0/75.D0*G(291)**4
     $                   *MODSQ(G(321))+16.D0/3.D0*G(293)**2*SIG3+3.D0
     $                   *G(292)**2*SIG2+1.D0/15.D0*G(291)**2*SIG1)
     $                   *ID(I,J)
            BETA2ML(I,J)=-DMATMUL(0,DUM2L1,YEDYEYEDYE,I,J)-4.D0
     $                   *YEDMEMYEYEDYE(I,J)-4.D0*YEDYEMLMYEDYE(I,J)
     $                   -4.D0*YEDYEYEDMEMYE(I,J)
     $                   -2.D0*YEDYEYEDYEMLM(I,J)
     $                   -DUM2L2(I,J)*DTRACE(DUM2L3)-YEDYE(I,J)
     $                   *DTRACE(DUM2L4)-4.D0*(YEDYEHEDHE(I,J)
     $                   +HEDHEYEDYE(I,J)+YEDHEHEDYE(I,J)
     $                   +HEDYEYEDHE(I,J))-HEDHE(I,J)*DTRACE(DUM2L5)
     $                   -YEDYE(I,J)*DTRACE(DUM2L6)-HEDYE(I,J)
     $                   *DTRACE(DUM2L7)-YEDHE(I,J)*DTRACE(DUM2L8)
     $                   +6.D0/5.D0*G(291)**2
     $                   *(DMATMUL(0,DUM2L9,YEDYE,I,J)
     $                   +2.D0*YEDMEMYE(I,J)+YEDYEMLM(I,J)
     $                   +2.D0*HEDHE(I,J)-2.D0*G(321)*HEDYE(I,J)
     $                   -2.D0*CON(G(321))*YEDHE(I,J)+4.D0
     $                   *MODSQ(G(321))*YEDYE(I,J))+(-6.D0/5.D0
     $                   *G(291)**2*SP+33.D0*G(292)**4*MODSQ(G(322))
     $                   +18.D0/5.D0*G(292)**2*G(291)**2
     $                   *(MODSQ(G(322))+MODSQ(G(321))
     $                   +DRE(2,G(321),G(322)))+621.D0/25.D0
     $                   *G(291)**4*MODSQ(G(321))+3.D0*G(292)**2*SIG2
     $                   +3.D0/5.D0*G(291)**2*SIG1)*ID(I,J)
            BETA2MU(I,J)=-DMATMUL(0,DUM2U1,YUYUDYUYUD,I,J)-4.D0
     $                   *YUMQMYUDYUYUD(I,J)-4.D0*YUYUDMUPMYUYUD(I,J)
     $                   -4.D0*YUYUDYUMQMYUD(I,J)
     $                   -2.D0*YUYUDYUYUDMUPM(I,J)
     $                   -DMATMUL(0,DUM2U2,YUYDDYDYUD,I,J)-4.D0
     $                   *YUMQMYDDYDYUD(I,J)-4.D0*YUYDDMDMYDYUD(I,J)
     $                   -4.D0*YUYDDYDMQMYUD(I,J)
     $                   -2.D0*YUYDDYDYUDMUPM(I,J)
     $                   -DUM2U3(I,J)*6.D0*TYUDYU-12.D0*YUYUD(I,J)
     $                   *DTRACE(DUM2U4)-4.D0*(HUHUDYUYUD(I,J)
     $                   +YUYUDHUHUD(I,J)+HUYUDYUHUD(I,J)
     $                   +YUHUDHUYUD(I,J))-4.D0*(HUHDDYDYUD(I,J)
     $                   +YUYDDHDHUD(I,J)+HUYDDYDHUD(I,J)
     $                   +YUHDDHDYUD(I,J))-12.D0*(HUHUD(I,J)*TYUDYU
     $                   +YUYUD(I,J)*THUDHU+HUYUD(I,J)*THUDYU+YUHUD(I,J)
     $                   *TYUDHU)+(6.D0*G(292)**2-2.D0/5.D0*G(291)**2)
     $                   *(DMATMUL(0,DUM2U5,YUYUD,I,J)
     $                   +2.D0*YUMQMYUD(I,J)+YUYUDMUPM(I,J)
     $                   +2.D0*HUHUD(I,J))+12.D0*G(292)**2
     $                   *(2.D0*MODSQ(G(322))*YUYUD(I,J)-CON(G(322))
     $                   *HUYUD(I,J)-G(322)*YUHUD(I,J))-4.D0/5.D0
     $                   *G(291)**2*(2.D0*MODSQ(G(321))*YUYUD(I,J)
     $                   -CON(G(321))*HUYUD(I,J)-G(321)*YUHUD(I,J))
     $                   +(-8.D0/5.D0*G(291)**2*SP-128.D0/3.D0*G(293)**4
     $                   *MODSQ(G(323))+512.D0/45.D0*G(293)**2
     $                   *G(291)**2*(MODSQ(G(323))+MODSQ(G(321))
     $                   +DRE(2,G(321),G(323)))+3424.D0/75.D0
     $                   *G(291)**4*MODSQ(G(321))+16.D0/3.D0*G(293)**2
     $                   *SIG3+16.D0/15.D0*G(291)**2*SIG1)*ID(I,J)
            BETA2MD(I,J)=-DMATMUL(0,DUM2D1,YDYDDYDYDD,I,J)-4.D0
     $                   *YDMQMYDDYDYDD(I,J)-4.D0*YDYDDMDMYDYDD(I,J)
     $                   -4.D0*YDYDDYDMQMYDD(I,J)
     $                   -2.D0*YDYDDYDYDDMDM(I,J)
     $                   -DMATMUL(0,DUM2D2,YDYUDYUYDD,I,J)-4.D0
     $                   *YDMQMYUDYUYDD(I,J)-4.D0*YDYUDMUPMYUYDD(I,J)
     $                   -4.D0*YDYUDYUMQMYDD(I,J)
     $                   -2.D0*YDYUDYUYDDMDM(I,J)
     $                   -DUM2D3(I,J)*DTRACE(DUM2D4)-4.D0*YDYDD(I,J)
     $                   *DTRACE(DUM2D5)-4.D0*(HDHDDYDYDD(I,J)
     $                   +YDYDDHDHDD(I,J)+HDYDDYDHDD(I,J)
     $                   +YDHDDHDYDD(I,J))-4.D0*(HDHUDYUYDD(I,J)
     $                   +YDYUDHUHDD(I,J)+HDYUDYUHDD(I,J)
     $                   +YDHUDHUYDD(I,J))-4.D0*HDHDD(I,J)
     $                   *DTRACE(DUM2D6)-4.D0*YDYDD(I,J)*DTRACE(DUM2D7)
     $                   -4.D0*HDYDD(I,J)*DTRACE(DUM2D8)-4.D0*YDHDD(I,J)
     $                   *DTRACE(DUM2D9)+(6.D0*G(292)**2+2.D0/5.D0
     $                   *G(291)**2)*(DMATMUL(0,DUM2DA,YDYDD,I,J)+2.D0
     $                   *YDMQMYDD(I,J)+YDYDDMDM(I,J)+2.D0*HDHDD(I,J))
     $                   +12.D0*G(292)**2*(2.D0*MODSQ(G(322))
     $                   *YDYDD(I,J)-CON(G(322))*HDYDD(I,J)-G(322)
     $                   *YDHDD(I,J))+4.D0/5.D0*G(291)**2*(2.D0
     $                   *MODSQ(G(321))*YDYDD(I,J)-CON(G(321))
     $                   *HDYDD(I,J)-G(321)*YDHDD(I,J))+(4.D0/5.D0
     $                   *G(291)**2*SP-128.D0/3.D0*G(293)**4
     $                   *MODSQ(G(323))+128.D0/45.D0*G(293)**2
     $                   *G(291)**2*(MODSQ(G(323))+MODSQ(G(321))
     $                   +DRE(2,G(321),G(323)))+808.D0/75.D0*G(291)**4
     $                   *MODSQ(G(321))+16.D0/3.D0*G(293)**2*SIG3
     $                   +4.D0/15.D0*G(291)**2*SIG1)*ID(I,J)
            BETA2ME(I,J)=-DMATMUL(0,DUM2E1,YEYEDYEYED,I,J)-4.D0
     $                   *YEMLMYEDYEYED(I,J)-4.D0*YEYEDMEMYEYED(I,J)
     $                   -4.D0*YEYEDYEMLMYED(I,J)
     $                   -2.D0*YEYEDYEYEDMEM(I,J)
     $                   -DUM2E2(I,J)*DTRACE(DUM2E3)-4.D0*YEYED(I,J)
     $                   *DTRACE(DUM2E4)-4.D0*(HEHEDYEYED(I,J)
     $                   +YEYEDHEHED(I,J)+HEYEDYEHED(I,J)
     $                   +YEHEDHEYED(I,J))-4.D0*HEHED(I,J)
     $                   *DTRACE(DUM2E5)-4.D0*YEYED(I,J)*DTRACE(DUM2E6)
     $                   -4.D0*HEYED(I,J)*DTRACE(DUM2E7)-4.D0*YEHED(I,J)
     $                   *DTRACE(DUM2E8)+(6.D0*G(292)**2-6.D0/5.D0
     $                   *G(291)**2)*(DMATMUL(0,DUM2E9,YEYED,I,J)+2.D0
     $                   *YEMLMYED(I,J)+YEYEDMEM(I,J)+2.D0*HEHED(I,J))
     $                   +12.D0*G(292)**2*(2.D0*MODSQ(G(322))
     $                   *YEYED(I,J)-CON(G(322))*HEYED(I,J)-G(322)
     $                   *YEHED(I,J))-12.D0/5.D0*G(291)**2*(2.D0
     $                   *MODSQ(G(321))*YEYED(I,J)-CON(G(321))
     $                   *HEYED(I,J)-G(321)*YEHED(I,J))+(12.D0/5.D0
     $                   *G(291)**2*SP+2808.D0/25.D0*G(291)**4
     $                   *MODSQ(G(321))+12.D0/5.D0*G(291)**2*SIG1)
     $                   *ID(I,J)
            B2MQM(I,J)=BETA2MQ(I,J)
            B2MLM(I,J)=BETA2ML(I,J)
            B2MUM(I,J)=BETA2MU(I,J)
            B2MDM(I,J)=BETA2MD(I,J)
            B2MEM(I,J)=BETA2ME(I,J)
          END  IF
!
!Calculate the differentials DM
!
          DM(1,I,J)=1.D0/16.D0/PI**2*BETA1MQ(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2MQ(I,J)
          DM(2,I,J)=1.D0/16.D0/PI**2*BETA1ML(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2ML(I,J)
          DM(3,I,J)=1.D0/16.D0/PI**2*BETA1MU(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2MU(I,J)
          DM(4,I,J)=1.D0/16.D0/PI**2*BETA1MD(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2MD(I,J)
          DM(5,I,J)=1.D0/16.D0/PI**2*BETA1ME(I,J)+1.D0/(16.D0*PI**2)**2
     $              *BETA2ME(I,J)
!
!Convert into form readable by RKSTP
!
          F(352+(I-1)*3+J)=DM(1,I,J)
          F(361+(I-1)*3+J)=DM(2,I,J)
          F(370+(I-1)*3+J)=DM(3,I,J)
          F(379+(I-1)*3+J)=DM(4,I,J)
          F(388+(I-1)*3+J)=DM(5,I,J)
        END DO
      END DO
!
!***********************************************************************
!END OF MSSM SECTION
!
!Now I can find the threshold gauge running with a change at m_H to SM
!
      DO I=1,3
        SUM=0.D0
        DO J=1,3
          IF(THHH.EQ.0)THEN
            SUM=SUM+B2LPSM(I,J)*G(J)**2
          ELSE
            SUM=SUM+B2LPM(I,J)*G(J)**2
          END IF
        END DO
        IF(THHH.EQ.0)THEN
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CSM(I,1)*TLYUDLYU+CSM(I,2)*TLYDDLYD
     $         +CSM(I,3)*TLYEDLYE)))
        ELSE
          F(I)=G(I)**3/16.D0/PI**2*(B1LP(I)+DBLE(SW2LP)/16.D0/PI**2
     $         *(SUM-(CM(I,1)*TYUDYU+CM(I,2)*TYDDYD+CM(I,3)*TYEDYE)))
        END IF
      END DO
!
!Next the full Yukawas
!
      DO I=1,3
        DO J=1,3
          DUMU1(I,J)=THSH*SQFTUQDFTUQ(I,J)
     $               +4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $               +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
          DUMU2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMD1(I,J)=THSH*SQFTDQDFTDQ(I,J)
     $               +1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $               +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
          DUMD2(I,J)=2.D0*THSH*SUFTUUFTUUD(I,J)
     $               +2.D0*THSH*SDFTDDFTDDD(I,J)
     $               +3.D0*THSW*SQGTQTGTQS(I,J)
     $               +1.D0/9.D0*THSB*SQGTPQTGTPQS(I,J)
     $               +16.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUME1(I,J)=THSH*SLFTELDFTEL(I,J)+THSB*SEGTPESGTPET(I,J)
          DUME2(I,J)=2.D0*THSH*SEFTEEFTEED(I,J)
     $               +3.D0*THSW*SLGTLTGTLS(I,J)+THSB*SLGTPLTGTPLS(I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
!
          B1U(I,J)=1.D0/2.D0*(3.D0*FUFUDFU(I,J)+FDFDDFU(I,J))
     $             +DMATMUL(0,FU,DUMU1,I,J)
     $             +1.D0/4.D0*DMATMUL(0,DUMU2,FU,I,J)
     $             +THSH*(-3.D0*THSW*G(204)
     $              *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $              *G(184)*SQGTPQTFTUQ(I,J))
     $             -4.D0/3.D0*THSB*THSH*G(184)
     $              *SUFTUUGTPUT(I,J)
     $             +FU(I,J)*3.D0*TFUDFU
     $             +1.D0/2.D0*THSH*FU(I,J)*(3.D0*THSW*MGTHUSQ
     $              +THSB*MGTPHUSQ)
     $             -FU(I,J)*(17.D0/20.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1D(I,J)=1.D0/2.D0*(3.D0*FDFDDFD(I,J)+FUFUDFD(I,J))
     $             +DMATMUL(0,FD,DUMD1,I,J)
     $             +1.D0/4.D0*DMATMUL(0,DUMD2,FD,I,J)
     $             +THSH*(-3.D0*THSW*G(205)
     $              *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $              *G(185)*SQGTPQTFTDQ(I,J))
     $             -2.D0/3.D0*THSB*THSH*G(185)
     $              *SDFTDDGTPDT(I,J)
     $             +FD(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FD(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FD(I,J)*(3.D0/12.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
!
          B1E(I,J)=3.D0/2.D0*FEFEDFE(I,J)
     $             +DMATMUL(0,FE,DUME1,I,J)
     $             +1.D0/4.D0*DMATMUL(0,DUME2,FE,I,J)
     $             +THSH*(-3.D0*THSW*G(205)
     $              *SLGTLTFTEL(I,J)+THSB
     $              *G(185)*SLGTPLTFTEL(I,J))
     $             -2.D0*THSB*THSH*G(185)
     $              *SEFTEEGTPET(I,J)
     $             +FE(I,J)*(3.D0*TFDDFD+TFEDFE)
     $             +1.D0/2.D0*THSH*FE(I,J)*(3.D0*THSW*MGTHDSQ
     $              +THSB*MGTPHDSQ)
     $             -FE(I,J)*(9.D0/4.D0*G(1)**2
     $              +9.D0/4.D0*G(2)**2)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
!
!Convert into form readable by RKSTP. The transpose in BETA2 takes
!account of the differences in notation.
!
          F(3+(I-1)*3+J)=(1.D0/16.D0/PI**2*B1U(I,J)
     $                   +1.D0/(16.D0*PI**2)**2*B2YMU(J,I))*THHH
          F(12+(I-1)*3+J)=(1.D0/16.D0/PI**2*B1D(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YMD(J,I))*THHH
          F(21+(I-1)*3+J)=(1.D0/16.D0/PI**2*B1E(I,J)
     $                    +1.D0/(16.D0*PI**2)**2*B2YME(J,I))*THHH
        END DO
      END DO
!
!The lambdas use the same dummy matrices. I will reuse the betas
!Only find the lambdas if we are below m_H
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            B1U(I,J)=0.D0
            B1D(I,J)=0.D0
            B1E(I,J)=0.D0
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            IF(SW2LP.EQ.1)THEN
              DUM2U(I,J)=3.D0/2.D0*LYUDLYU2(I,J)
     $                  -LYUDLYULYDDLYD(I,J)
     $                  -1.D0/4.D0*LYDDLYDLYUDLYU(I,J)
     $                  +11.D0/4.D0*LYDDLYD2(I,J)
     $                  +Y2*(5.D0/4.D0*LYDDLYD(I,J)-9.D0/4.D0
     $                  *LYUDLYU(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(429)**2*ID(I,J)-2.D0*G(429)*(3.D0
     $                  *LYUDLYU(I,J)+LYDDLYD(I,J))+(221.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYUDLYU(I,J)-(17.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYDDLYD(I,J)+Y4*ID(I,J)+((7.D0/150.D0
     $                  +2.D0/3.D0*NG)*G(1)**4-9.D0/20.D0*G(1)**2
     $                  *G(2)**2+19.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2D(I,J)=3.D0/2.D0*LYDDLYD2(I,J)
     $                  -LYDDLYDLYUDLYU(I,J)
     $                  -1.D0/4.D0*LYUDLYULYDDLYD(I,J)
     $                  +11.D0/4.D0*LYUDLYU2(I,J)
     $                  +Y2*(5.D0/4.D0*LYUDLYU(I,J)-9.D0/4.D0
     $                  *LYDDLYD(I,J))-CHI4*ID(I,J)+3.D0/2.D0
     $                  *G(429)**2*ID(I,J)-2.D0*G(429)*(3.D0
     $                  *LYDDLYD(I,J)+LYUDLYU(I,J))+(161.D0/80.D0
     $                  *G(1)**2+117.D0/16.D0*G(2)**2+20.D0
     $                  *G(3)**2)*LYDDLYD(I,J)-(77.D0/80.D0*G(1)**2
     $                  -27.D0/16.D0*G(2)**2+20.D0*G(3)**2)
     $                  *LYUDLYU(I,J)+Y4*ID(I,J)+(-(37.D0/300.D0
     $                  -4.D0/15.D0*NG)*G(1)**4-27.D0/20.D0*G(1)**2
     $                  *G(2)**2+31.D0/15.D0*G(1)**2*G(3)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4+9.D0*G(2)**2
     $                  *G(3)**2-(292.D0/3.D0-16.D0/3.D0*NG)
     $                  *G(3)**4)*ID(I,J)
              DUM2E(I,J)=3.D0/2.D0*LYEDLYE2(I,J)
     $                  -Y2*9.D0/4.D0*LYEDLYE(I,J)-CHI4*ID(I,J)
     $                  +3.D0/2.D0*G(429)**2*ID(I,J)-6.D0*G(429)
     $                  *LYEDLYE(I,J)+(441.D0/80.D0*G(1)**2
     $                  +117.D0/16.D0*G(2)**2)*LYEDLYE(I,J)
     $                  +Y4*ID(I,J)+((21.D0/100.D0+8.D0/5.D0*NG)
     $                  *G(1)**4+27.D0/20.D0*G(1)**2*G(2)**2
     $                  -(101.D0/8.D0-2.D0*NG)*G(2)**4)*ID(I,J)
            ELSE
              DUM2U(I,J)=0.D0
              DUM2D(I,J)=0.D0
              DUM2E(I,J)=0.D0
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            B1U(I,J)=1.D0/2.D0*(3.D0*LULUDLU(I,J)+LDLDDLU(I,J)
     $                -4.D0*LDLDDLU(I,J))
     $               +DMATMUL(0,LU,DUMU1,I,J)
     $               +1.D0/4.D0*DMATMUL(0,DUMU2,LU,I,J)
     $               +THSH*(-3.D0*THSW*G(289)
     $                *SQGTQTFTUQ(I,J)+1.D0/3.D0*THSB
     $                *G(287)*SQGTPQTFTUQ(I,J))
     $               -4.D0/3.D0*THSB*THSH*G(287)
     $                *SUFTUUGTPUT(I,J)
     $               +LU(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LU(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LU(I,J)*(17.D0/20.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1D(I,J)=1.D0/2.D0*(3.D0*LDLDDLD(I,J)+LULUDLD(I,J)
     $                -4.D0*LULUDLD(I,J))
     $               +DMATMUL(0,LD,DUMD1,I,J)
     $               +1.D0/4.D0*DMATMUL(0,DUMD2,LD,I,J)
     $               +THSH*(-3.D0*THSW*G(290)
     $                *SQGTQTFTDQ(I,J)-1.D0/3.D0*THSB
     $                *G(288)*SQGTPQTFTDQ(I,J))
     $               -2.D0/3.D0*THSB*THSH*G(288)
     $                *SDFTDDGTPDT(I,J)
     $               +LD(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LD(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LD(I,J)*(3.D0/12.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2+8.D0*G(3)**2)
            B1E(I,J)=3.D0/2.D0*LELEDLE(I,J)+DMATMUL(0,LE,DUME1,I,J)
     $               +1.D0/4.D0*DMATMUL(0,DUME2,LE,I,J)
     $               +THSH*(-3.D0*THSW*G(290)
     $                *SLGTLTFTEL(I,J)+THSB
     $                *G(288)*SLGTPLTFTEL(I,J))
     $               -2.D0*THSB*THSH*G(288)
     $                *SEFTEEGTPET(I,J)
     $               +LE(I,J)*(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)
     $               +1.D0/2.D0*THSH*LE(I,J)*(3.D0*THSW
     $                *(MSGTHUSQ+MCGTHDSQ)+THSB
     $                *(MSGTPHUSQ+MCGTPHDSQ))
     $               -LE(I,J)*(9.D0/4.D0*G(1)**2
     $                +9.D0/4.D0*G(2)**2)
!
            IF(SW2LP.EQ.1)THEN
              BETA2U(I,J)=DMATMUL(0,LYU,DUM2U,I,J)
              BETA2D(I,J)=DMATMUL(0,LYD,DUM2D,I,J)
              BETA2E(I,J)=DMATMUL(0,LYE,DUM2E,I,J)
            END IF
          END DO
        END DO
        DO I=1,3
          DO J=1,3
!
!Convert into form readable by RKSTP.
!
            F(111+(I-1)*3+J)=1.D0/16.D0/PI**2*B1U(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2U(J,I)
            F(120+(I-1)*3+J)=1.D0/16.D0/PI**2*B1D(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2D(J,I)
            F(129+(I-1)*3+J)=1.D0/16.D0/PI**2*B1E(I,J)
     $              +1.D0/(16.D0*PI**2)**2*BETA2E(J,I)
          END DO
        END DO
      END IF
!
!Next I am going to work out the gaugino terms, \mu and M_{1,2,3}
!and the running of B in MV notation
!
      DO I=1,3
        DO J=1,3
          DUMGRKMU1(I,J)=3.D0*SUFTUUFTUUD(I,J)+3.D0*SDFTDDFTDDD(I,J)
     $                   +SEFTEEFTEED(I,J)+3.D0*SQFTUQDFTUQ(I,J)
     $                   +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
!
          DUMM1(I,J)=1.D0/3.D0*SQGTPQDGTPQ(I,J)+SLGTPLDGTPL(I,J)
     $               +8.D0/3.D0*SUGTPUGTPUD(I,J)
     $               +2.D0/3.D0*SDGTPDGTPDD(I,J)+2.D0*SEGTPEGTPED(I,J)
!
          DUMM2(I,J)=3.D0*SQGTQDGTQ(I,J)+SLGTLDGTL(I,J)
!
          DUMM3(I,J)=2.D0*SQGTSQDGTSQ(I,J)+SUGTSUGTSUD(I,J)
     $               +SDGTSDGTSDD(I,J)
!     
        END DO
      END DO
      TDUMGRKMU=DTRACE(DUMGRKMU1)
      TDUMM1=DTRACE(DUMM1)
      TDUMM2=DTRACE(DUMM2)
      TDUMM3=DTRACE(DUMM3)
      SQTSUMTSFUSAUT=TSFMUL(THSQ,SUMTSFUSAUT)
      SQTSDMTSFDSADT=TSFMUL(THSQ,SDMTSFDSADT)
      SLTSEMTSFESAET=TSFMUL(THSL,SEMTSFESAET)
!
      IF(THHH.EQ.0)THEN
        BM(1)=G(31)*(TDUMM1+THSH*(MSGTPHUSQ+MCGTPHDSQ))
     $        -2.D0*THSH*(G(287)*G(108)*G(288)
     $                    +G(287)*G(108)*G(288))
!
        BM(2)=G(32)*(TDUMM2+THSH*(MSGTHUSQ+MCGTHDSQ))
     $        -2.D0*THSH*(G(289)*G(108)*G(290)
     $                    +G(289)*G(108)*G(290))
     $       -12.D0*G(32)*G(2)**2
!
        BM(3)=G(33)*TDUMM3-18.D0*G(33)*G(3)**2
!
        BGRKMU=1.D0/2.D0*G(108)*THSH*TDUMGRKMU
     $       +1.D0/4.D0*G(108)*THSH*((3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $        +(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ))
     $       -(3.D0*THSW*G(289)*(G(32)+(0.D0,1.D0)*G(600))*G(290)
     $        +THSB*G(287)*(G(31)+(0.D0,1.D0)*G(599))*G(288))
     $       -G(108)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
!
        BETA1B=0.D0
      ELSE
        BM(1)=G(31)*(TDUMM1+THSH*(MGTPHUSQ+MGTPHDSQ))
!
        BM(2)=G(32)*(TDUMM2+THSH*(MGTHUSQ+MGTHDSQ))
     $        -12.D0*G(32)*G(2)**2
!
        BM(3)=G(33)*TDUMM3-18.D0*G(33)*G(3)**2
!
        BGRKMU=1.D0/2.D0*G(108)*THSH*TDUMGRKMU
     $         +1.D0/4.D0*G(108)*THSH*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ
     $          +3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $         -G(108)*(9.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
        BETA1B=-1.D0/2.D0*(G(457)**2+3.D0*G(458)**2)*G(109)
     $         +6.D0*SQTSUMTSFUSAUT
     $         +6.D0*SQTSDMTSFDSADT
     $         +2.D0*SLTSEMTSFESAET
     $         -(2.D0*THSH*THSB*G(108)*G(184)*G(185)
     $          *(G(31)-(0.D0,1.D0)*G(599))+6.D0*THSH*THSW*G(108)
     $          *G(204)*G(205)*(G(32)-(0.D0,1.D0)*G(600)))
     $         -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(109)
     $         +1.D0/2.D0*(6.D0*TFUDFU+6.D0*TFDDFD+2.D0*TFEDFE
     $          +THSB*THSH*(MGTPHUSQ+MGTPHDSQ)+3.D0*THSW*THSH
     $          *(MGTHUSQ+MGTHDSQ))*G(109)
      END IF
!
!The RKSTP compatible derivatives are: (with a minus sign
!to convert the two loop mu running to BT notation)
!NB: There is some work here converting the MSSM MV notation
!    M_{1,2,3} 2-lp beta functions to BT notation.
!
      F(31)=THSB*(1.D0/16.D0/PI**2*BM(1)
     $      +1.D0/(16.D0*PI**2)**2*1.D0/2.D0*(B2M(1)+B2M(1)))
      F(32)=THSW*(1.D0/16.D0/PI**2*BM(2)
     $      +1.D0/(16.D0*PI**2)**2*1.D0/2.D0*(B2M(2)+B2M(2)))
      F(33)=THGL*(1.D0/16.D0/PI**2*BM(3)
     $      +1.D0/(16.D0*PI**2)**2*1.D0/2.D0*(B2M(3)+B2M(3)))
!
      F(108)=(1.D0/16.D0/PI**2*BGRKMU+1.D0/(16.D0*PI**2)**2*B2GRKMUM)
     $                                                           *THSH
      IF(THHH.EQ.1)THEN
        F(109)=1.D0/16.D0/PI**2*BETA1B-1.D0/(16.D0*PI**2)**2*BETA2BM
      END IF
!
!V_U and V_D - from PRD 49,4882 (1994)
!
      BETA1VU=3.D0/4.D0*(1.D0/5.D0*G(291)**2+G(292)**2)-3.D0*TYUDYU
      BETA1VD=3.D0/4.D0*(1.D0/5.D0*G(291)**2+G(292)**2)-3.D0*TYDDYD
     $        -TYEDYE
      IF(SW2LP.EQ.1)THEN
        BETA2VU=3.D0/4.D0*(3.D0*TYUDYUYUDYU+3.D0*TYUDYUYDDYD)-(19.D0
     $          /10.D0*G(291)**2+9.D0/2.D0*G(292)**2+20.D0*G(293)**2)
     $          *TYUDYU-(279.D0/800.D0+1803.D0/1600.D0*3.D0)*G(291)**4
     $          -(207.D0/32.D0+357.D0/64.D0*3.D0)*G(292)**4-(27.D0/80.D0
     $          +9.D0/80.D0*3.D0)*G(291)**2*G(292)**2
        BETA2VD=3.D0/4.D0*(3.D0*TYDDYDYDDYD+3.D0*TYDDYDYUDYU
     $          +TYEDYEYEDYE)-(2.D0/5.D0*G(291)**2+9.D0/2.D0*G(292)**2
     $          +20.D0*G(293)**2)*TYDDYD-(9.D0/5.D0*G(291)**2+3.D0/2.D0
     $          *G(292)**2)*TYEDYE-(279.D0/800.D0+1803.D0/1600.D0*3.D0)
     $          *G(291)**4-(207.D0/32.D0+357.D0/64.D0*3.D0)*G(292)**4
     $          -(27.D0/80.D0+9.D0/80.D0*3.D0)*G(291)**2*G(292)**2
      END IF
!
      F(110)=THHH*(G(110)*(1.D0/16.D0/PI**2*BETA1VU
     $             +1.D0/(16.D0*PI**2)**2*BETA2VU))
      F(111)=THHH*(G(111)*(1.D0/16.D0/PI**2*BETA1VD
     $             +1.D0/(16.D0*PI**2)**2*BETA2VD))
!
!Now the trilinear couplings - it should be ok to use the same
!dummy matrices, but the betas must be reset
!
      DO I=1,3
        DO J=1,3
          BETA1U(I,J)=0.D0
          BETA2U(I,J)=0.D0
          BETA1D(I,J)=0.D0
          BETA2D(I,J)=0.D0
          BETA1E(I,J)=0.D0
          BETA2E(I,J)=0.D0
        END DO
      END DO
!
!Some common terms
!
      DO I=1,3
        DO J=1,3
          DUM1U1(I,J)=8.D0/9.D0*THSB*GTPUTGTPUS(I,J)
     $               +8.D0/3.D0*THGL*GTSUTGTSUS(I,J)
     $               +2.D0*THSH*FTUUDFTUU(I,J)
          DUM1U2(I,J)=THSH*FTUQFTUQD(I,J)+THSH*FTDQFTDQD(I,J)
     $               +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $               +3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $               +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
!
          DUM1D1(I,J)=DUM1U2(I,J)
          DUM1D2(I,J)=2.D0/9.D0*THSB*GTPDTGTPDS(I,J)
     $               +8.D0/3.D0*THGL*GTSDTGTSDS(I,J)
     $               +2.D0*THSH*FTDDDFTDD(I,J)
!
          DUM1E1(I,J)=THSH*FTELFTELD(I,J)
     $               +1.D0/2.D0*THSB*GTPLSGTPLT(I,J)
     $               +3.D0/2.D0*THSW*GTLSGTLT(I,J)
          DUM1E2(I,J)=2.D0*THSB*GTPETGTPES(I,J)
     $               +2.D0*THSH*FTEEDFTEE(I,J)
        END DO
      END DO
!
!Running above m_H first
!
      IF(THHH.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            SUAUID(I,J)=SFMUL(THSU,AU,ID,I,J)
            SUAUFUQDFUQ(I,J)=SFMUL(THSU,AU,FUQDFUQ,I,J)
            SQIDAU(I,J)=SFMUL(THSQ,ID,AU,I,J)
            SQFUHUDAU(I,J)=SFMUL(THSQ,FUHUD,AU,I,J)
            SQFUURFUURDAU(I,J)=SFMUL(THSQ,FUURFUURD,AU,I,J)
            SDADFDQDFUQ(I,J)=SFMUL(THSD,AD,FDQDFUQ,I,J)
!
            SQIDAD(I,J)=SFMUL(THSQ,ID,AD,I,J)
            SQFDDRFDDRDAD(I,J)=SFMUL(THSQ,FDDRFDDRD,AD,I,J)
            SQFDHDDAD(I,J)=SFMUL(THSQ,FDHDD,AD,I,J)
            SLFEHDDAE(I,J)=SFMUL(THSL,FEHDD,AE,I,J)
            SDADID(I,J)=SFMUL(THSD,AD,ID,I,J)
            SDADFDQDFDQ(I,J)=SFMUL(THSD,AD,FDQDFDQ,I,J)
            SUAUFUQDFDQ(I,J)=SFMUL(THSU,AU,FUQDFDQ,I,J)
!
            SLIDAE(I,J)=SFMUL(THSL,ID,AE,I,J)
            SLFEERFEERDAE(I,J)=SFMUL(THSL,FEERFEERD,AE,I,J)
            SEAEID(I,J)=SFMUL(THSE,AE,ID,I,J)
            SEAEFELDFEL(I,J)=SFMUL(THSE,AE,FELDFEL,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SUSQIDAUID(I,J)=SFMUL(THSU,SQIDAU,ID,I,J)
            SDSQIDADID(I,J)=SFMUL(THSD,SQIDAD,ID,I,J)
            SESLIDAEID(I,J)=SFMUL(THSE,SLIDAE,ID,I,J)
          END DO
        END DO
        SUTSQFUHUDAU=TSFMUL(THSU,SQFUHUDAU)
        SDTSQFDHDDAD=TSFMUL(THSD,SQFDHDDAD)
        SETSLFEHDDAE=TSFMUL(THSE,SLFEHDDAE)
        DO I=1,3
          DO J=1,3
            BETA1U(I,J)=-2.D0/3.D0*G(457)**2*SUAUID(I,J)
     $                 +2.D0*SUAUFUQDFUQ(I,J)
     $                 -2.D0*(G(457)**2/9.D0+4.D0*G(459)**2/3.D0)
     $                  *SUSQIDAUID(I,J)
     $                 +6.D0*FUHU(I,J)*SUTSQFUHUDAU
     $                 +(G(457)**2/6.D0-3.D0*G(458)**2/2.D0)*SQIDAU(I,J)
     $                 +4.D0*SQFUURFUURDAU(I,J)
     $                 +2.D0*SDADFDQDFUQ(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*G(184)*GTPQSFTUU(I,J)
     $                  -4.D0/3.D0*GTPQSFUGTPUS(I,J)
     $                  -4.D0*THSH*FTUQGTPUS(I,J)*G(184))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSFUGTSUS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *G(204)*GTQSFTUU(I,J)
     $                 +SFMUL(THSU,AU,DUM1U1,I,J)
     $                 +(3.D0*TFUDFU+1.D0/2.D0*THSH*THSB*MGTPHUSQ
     $                  +3.D0/2.D0*THSH*THSW*MGTHUSQ)*AU(I,J)
     $                 +SFMUL(THSQ,DUM1U2,AU,I,J)
     $                 -(13.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*AU(I,J)
!
            BETA1D(I,J)=-(G(457)**2/6.D0+3.D0*G(458)**2/2.D0)
     $                  *SQIDAD(I,J)
     $                 +4.D0*SQFDDRFDDRDAD(I,J)
     $                 +2.D0*(G(457)**2/18.D0-4.D0*G(459)**2/3.D0)
     $                  *SDSQIDADID(I,J)
     $                 +6.D0*FDHD(I,J)*SDTSQFDHDDAD
     $                 +2.D0*FDHD(I,J)*SETSLFEHDDAE
     $                 -1.D0/3.D0*G(457)**2*SDADID(I,J)
     $                 +2.D0*SDADFDQDFDQ(I,J)
     $                 +2.D0*SUAUFUQDFDQ(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(-THSH*G(185)*GTPQSFTDD(I,J)
     $                  +2.D0/3.D0*GTPQSFDGTPDS(I,J)
     $                  -2.D0*THSH*FTDQGTPDS(I,J)*G(185))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSFDGTSDS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *G(205)*GTQSFTDD(I,J)
     $                 +SFMUL(THSQ,DUM1D1,AD,I,J)
     $                 +(3.D0*TFDDFD+TFEDFE
     $                  +1.D0/2.D0*THSH*THSB*MGTPHDSQ
     $                  +3.D0/2.D0*THSH*THSW*MGTHDSQ)*AD(I,J)
     $                 +SFMUL(THSD,AD,DUM1D2,I,J)
     $                 -(7.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*AD(I,J)
!
            BETA1E(I,J)=(G(457)**2/2.D0-3.D0*G(458)**2/2.D0)*SLIDAE(I,J)
     $                 +4.D0*SLFEERFEERDAE(I,J)
     $                 -G(457)**2*SESLIDAEID(I,J)
     $                 +2.D0*FEHD(I,J)*SETSLFEHDDAE
     $                 +6.D0*FEHD(I,J)*SDTSQFDHDDAD
     $                 -G(457)**2*SEAEID(I,J)
     $                 +2.D0*SEAEFELDFEL(I,J)
     $                 +2.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*G(185)*GTPLSFTEE(I,J)
     $                  -2.D0*GTPLSFEGTPES(I,J)
     $                  -2.D0*THSH*FTELGTPES(I,J)*G(185))
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *G(205)*GTLSFTEE(I,J)
     $                 +SFMUL(THSL,DUM1E1,AE,I,J)
     $                 +(3.D0*TFDDFD+TFEDFE
     $                  +1.D0/2.D0*THSH*THSB*MGTPHDSQ
     $                  +3.D0/2.D0*THSH*THSW*MGTHDSQ)*AE(I,J)
     $                 +SFMUL(THSE,AE,DUM1E2,I,J)
     $                 -(27.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $                  *AE(I,J)
          END DO
        END DO
!
!Convert into form readable by RKSTP (converting 2-lp to BT notation...)
!
        DO I=1,3
          DO J=1,3
            F(33+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BETA1U(I,J)
     $                      -1.D0/(16.D0*PI**2)**2*B2HMU(J,I))*THSU(J)
            F(42+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BETA1D(I,J)
     $                      -1.D0/(16.D0*PI**2)**2*B2HMD(J,I))*THSD(J)
            F(51+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BETA1E(I,J)
     $                      -1.D0/(16.D0*PI**2)**2*B2HME(J,I))*THSE(J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            SUTRIUID(I,J)=SFMUL(THSU,TRIU,ID,I,J)
            SUTRIUSFUQDSFUQ(I,J)=SFMUL(THSU,TRIU,SFUQDSFUQ,I,J)
            SQIDTRIU(I,J)=SFMUL(THSQ,ID,TRIU,I,J)
            SQFUHUDTRIU(I,J)=SFMUL(THSQ,FUHUD,TRIU,I,J)
            SQSFUURSFUURDTRIU(I,J)=SFMUL(THSQ,SFUURSFUURD,TRIU,I,J)
            SQCFDDRCFDDRDTRIU(I,J)=SFMUL(THSQ,CFDDRCFDDRD,TRIU,I,J)
!
            SQIDTRID(I,J)=SFMUL(THSQ,ID,TRID,I,J)
            SQSFUURSFUURDTRID(I,J)=SFMUL(THSQ,SFUURSFUURD,TRID,I,J)
            SQCFDDRCFDDRDTRID(I,J)=SFMUL(THSQ,CFDDRCFDDRD,TRID,I,J)
            SQFDHDDTRID(I,J)=SFMUL(THSQ,FDHDD,TRID,I,J)
            SLFEHDDTRIE(I,J)=SFMUL(THSL,FEHDD,TRIE,I,J)
            SDTRIDID(I,J)=SFMUL(THSD,TRID,ID,I,J)
            SDTRIDCFDQDCFDQ(I,J)=SFMUL(THSD,TRID,CFDQDCFDQ,I,J)
!
            SLIDTRIE(I,J)=SFMUL(THSL,ID,TRIE,I,J)
            SLCFEERCFEERDTRIE(I,J)=SFMUL(THSL,CFEERCFEERD,TRIE,I,J)
            SETRIEID(I,J)=SFMUL(THSE,TRIE,ID,I,J)
            SETRIECFELDCFEL(I,J)=SFMUL(THSE,TRIE,CFELDCFEL,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SUSQIDTRIUID(I,J)=SFMUL(THSU,SQIDTRIU,ID,I,J)
            SDSQIDTRIDID(I,J)=SFMUL(THSD,SQIDTRID,ID,I,J)
            SESLIDTRIEID(I,J)=SFMUL(THSE,SLIDTRIE,ID,I,J)
          END DO
        END DO
        SUTSQFUHUDTRIU=TSFMUL(THSU,SQFUHUDTRIU)
        SDTSQFDHDDTRID=TSFMUL(THSD,SQFDHDDTRID)
        SETSLFEHDDTRIE=TSFMUL(THSE,SLFEHDDTRIE)
        DO I=1,3
          DO J=1,3
            B1TRIU(I,J)=2.D0*G(541)/3.D0*SUTRIUID(I,J)
     $                 +2.D0*SUTRIUSFUQDSFUQ(I,J)
     $                 -2.D0*(G(457)**2/9.D0+4.D0*G(459)**2/3.D0)
     $                  *SUSQIDTRIUID(I,J)
     $                 +6.D0*FUHU(I,J)*SUTSQFUHUDTRIU
     $                 -2.D0*(G(541)/12.D0-3.D0*G(542)/4.D0)
     $                  *SQIDTRIU(I,J)
     $                 +4.D0*SQSFUURSFUURDTRIU(I,J)
     $                 -2.D0*SQCFDDRCFDDRDTRIU(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*G(287)*GTPQSFTUU(I,J)
     $                  -4.D0/3.D0*GTPQSLUGTPUS(I,J)
     $                  -4.D0*THSH*FTUQGTPUS(I,J)*G(287))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSLUGTSUS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *G(289)*GTQSFTUU(I,J)
     $                 +2.D0/3.D0*THSB*THSH*G(108)*G(288)
     $                  *(4.D0*FTUQGTPUS(I,J)-GTPQSFTUU(I,J))
     $                 +6.D0*THSH*THSW*G(108)*G(290)
     $                  *GTQSFTUU(I,J)
     $                 -4.D0*THSH*G(108)*FTDQLDDFTUU(I,J)
     $                 +SFMUL(THSU,TRIU,DUM1U1,I,J)
     $                 +(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)*TRIU(I,J)
     $                 +1.D0/2.D0*(THSB*MCGTPHDSQ+3.D0*THSW*MCGTHDSQ
     $                  +THSB*MSGTPHUSQ+3.D0*MSGTHUSQ)*THSH*TRIU(I,J)
     $                 +SFMUL(THSQ,DUM1U2,TRIU,I,J)
     $                 -(13.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*TRIU(I,J)
!
            B1TRID(I,J)=-2.D0*(G(541)/12.D0+3.D0*G(542)/4.D0)
     $                  *SQIDTRID(I,J)
     $                 -2.D0*SQSFUURSFUURDTRID(I,J)
     $                 +4.D0*SQCFDDRCFDDRDTRID(I,J)
     $                 +2.D0*(G(457)**2/18.D0-4.D0*G(459)**2/3.D0)
     $                  *SDSQIDTRIDID(I,J)
     $                 +6.D0*FDHD(I,J)*SDTSQFDHDDTRID
     $                 +2.D0*FDHD(I,J)*SETSLFEHDDTRIE
     $                 -G(541)/3.D0*SDTRIDID(I,J)
     $                 +2.D0*SDTRIDCFDQDCFDQ(I,J)
     $                 -4.D0*THSH*G(108)*FTUQLUDFTDD(I,J)
     $                 +2.D0/3.D0*THSB*THSH*G(108)*G(287)
     $                  *(2.D0*FTDQGTPDS(I,J)+GTPQSFTDD(I,J))
     $                 +6.D0*THSH*THSW*G(108)*G(289)
     $                  *GTQSFTDD(I,J)
     $                 +2.D0/3.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(-THSH*G(288)*GTPQSFTDD(I,J)
     $                  +2.D0/3.D0*GTPQSLDGTPDS(I,J)
     $                  -2.D0*THSH*FTDQGTPDS(I,J)*G(288))
     $                 -32.D0/3.D0*THGL*(G(33)-(0.D0,1.D0)*G(601))
     $                  *GTSQSLDGTSDS(I,J)
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *G(290)*GTQSFTDD(I,J)
     $                 +SFMUL(THSQ,DUM1D1,TRID,I,J)
     $                 +(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)*TRID(I,J)
     $                 +1.D0/2.D0*(THSB*MCGTPHDSQ+3.D0*THSW*MCGTHDSQ
     $                  +THSB*MSGTPHUSQ+3.D0*MSGTHUSQ)*THSH*TRID(I,J)
     $                 +SFMUL(THSD,TRID,DUM1D2,I,J)
     $                 -(7.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*TRID(I,J)
!
            B1TRIE(I,J)=2.D0*(G(541)/4.D0-3.D0*G(542)/4.D0)
     $                  *SLIDTRIE(I,J)
     $                 +4.D0*SLCFEERCFEERDTRIE(I,J)
     $                 -G(457)**2*SESLIDTRIEID(I,J)
     $                 +2.D0*FEHD(I,J)*SETSLFEHDDTRIE
     $                 +6.D0*FEHD(I,J)*SDTSQFDHDDTRID
     $                 -G(541)*SETRIEID(I,J)
     $                 +2.D0*SETRIECFELDCFEL(I,J)
     $                 +2.D0*THSB*THSH*G(108)*G(287)
     $                  *(2.D0*FTELGTPES(I,J)-GTPLSFTEE(I,J))
     $                 +6.D0*THSH*THSW*G(108)*G(289)
     $                  *GTLSFTEE(I,J)
     $                 +2.D0*THSB*(G(31)-(0.D0,1.D0)*G(599))
     $                  *(THSH*G(288)*GTPLSFTEE(I,J)
     $                  -2.D0*GTPLSLEGTPES(I,J)
     $                  -2.D0*THSH*FTELGTPES(I,J)*G(288))
     $                 -6.D0*THSW*THSH*(G(32)-(0.D0,1.D0)*G(600))
     $                  *G(290)*GTLSFTEE(I,J)
     $                 +SFMUL(THSL,DUM1E1,TRIE,I,J)
     $                 +(3.D0*TLUDLU+3.D0*TLDDLD+TLEDLE)*TRIE(I,J)
     $                 +1.D0/2.D0*(THSB*MCGTPHDSQ+3.D0*THSW*MCGTHDSQ
     $                  +THSB*MSGTPHUSQ+3.D0*MSGTHUSQ)*THSH*TRIE(I,J)
     $                 +SFMUL(THSE,TRIE,DUM1E2,I,J)
     $                 -(27.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $                  *TRIE(I,J)
          END DO
        END DO
!
!Convert into form readable by RKSTP (converting 2-lp to BT notation...)
!
        DO I=1,3
          DO J=1,3
            F(399+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*B1TRIU(I,J)
     $                       +1.D0/(16.D0*PI**2)**2*(-SINB
     $                       *B2HMU(J,I)-COSB*(G(398)
     $                       *B2YMU(J,I)+YU(J,I)*B2GRKMUM)))
     $                       *THSU(J)
            F(408+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*B1TRID(I,J)
     $                       +1.D0/(16.D0*PI**2)**2*(-COSB
     $                       *B2HMD(J,I)-SINB*(G(398)
     $                       *B2YMD(J,I)+YD(J,I)*B2GRKMUM)))
     $                       *THSD(J)
            F(417+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*B1TRIE(I,J)
     $                       +1.D0/(16.D0*PI**2)**2*(-COSB
     $                       *B2HME(J,I)-SINB*(G(398)
     $                       *B2YME(J,I)+YE(J,I)*B2GRKMUM)))
     $                       *THSE(J)
          END DO
        END DO
      END IF
!
!Now for the (mass)^2 parameters. DUMSP are dummy matrices for the
!calculation of S'
!
      DO I=1,3
        DO J=1,3
          SUIDMU(I,J)=SFMUL(THSU,ID,MUP,I,J)
          SQIDMQ(I,J)=SFMUL(THSQ,ID,MQ,I,J)
          SDIDMD(I,J)=SFMUL(THSD,ID,MD,I,J)
          SLIDML(I,J)=SFMUL(THSL,ID,ML,I,J)
          SEIDME(I,J)=SFMUL(THSE,ID,ME,I,J)
          SUFUHUSMU(I,J)=SFMUL(THSU,FUHUS,MUP,I,J)
          SDFDHDSMD(I,J)=SFMUL(THSD,FDHDS,MD,I,J)
          SQFUHUTMQ(I,J)=SFMUL(THSQ,FUHUT,MQ,I,J)
          SQFDHDTMQ(I,J)=SFMUL(THSQ,FDHDT,MQ,I,J)
          SEFEHDSME(I,J)=SFMUL(THSE,FEHDS,ME,I,J)
          SLFEHDTML(I,J)=SFMUL(THSL,FEHDT,ML,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          SQSQIDMQID(I,J)=SFMUL(THSQ,SQIDMQ,ID,I,J)
          SUSUFUHUSMUFUHUT(I,J)=SFMUL(THSU,SUFUHUSMU,FUHUT,I,J)
          SDSDFDHDSMDFDHDT(I,J)=SFMUL(THSD,SDFDHDSMD,FDHDT,I,J)
          SUAUSAUT(I,J)=SFMUL(THSU,AUS,AUT,I,J)
          DUM1Q1(I,J)=1.D0/18.D0*THSB*GTPQGTPQD(I,J)
     $               +3.D0/2.D0*THSW*GTQGTQD(I,J)
     $               +8.D0/3.D0*THGL*GTSQGTSQD(I,J)
     $               +THSH*FTUQSFTUQT(I,J)
     $               +THSH*FTDQSFTDQT(I,J)
!
          SUSUIDMUID(I,J)=SFMUL(THSU,SUIDMU,ID,I,J)
          SQSQFUHUTMQFUHUS(I,J)=SFMUL(THSQ,SQFUHUTMQ,FUHUS,I,J)
          SQAUTAUS(I,J)=SFMUL(THSQ,AUT,AUS,I,J)
          DUM1U1(I,J)=8.D0/9.D0*THSB*GTPUDGTPU(I,J)
     $               +8.D0/3.D0*THGL*GTSUDGTSU(I,J)
     $               +2.D0*THSH*FTUUTFTUUS(I,J)
!
          SDSDIDMDID(I,J)=SFMUL(THSD,SDIDMD,ID,I,J)
          SQSQFDHDTMQFDHDS(I,J)=SFMUL(THSQ,SQFDHDTMQ,FDHDS,I,J)
          SQADTADS(I,J)=SFMUL(THSQ,ADT,ADS,I,J)
          DUM1D1(I,J)=2.D0/9.D0*THSB*GTPDDGTPD(I,J)
     $               +8.D0/3.D0*THGL*GTSDDGTSD(I,J)
     $               +2.D0*THSH*FTDDTFTDDS(I,J)
!
          SLSLIDMLID(I,J)=SFMUL(THSL,SLIDML,ID,I,J)
          SESEFEHDSMEFEHDT(I,J)=SFMUL(THSE,SEFEHDSME,FEHDT,I,J)
          SDADSADT(I,J)=SFMUL(THSD,ADS,ADT,I,J)
          SEAESAET(I,J)=SFMUL(THSE,AES,AET,I,J)
          DUM1L1(I,J)=1.D0/2.D0*THSB*GTPLGTPLD(I,J)
     $               +3.D0/2.D0*THSW*GTLGTLD(I,J)
     $               +THSH*FTELSFTELT(I,J)
!
          SESEIDMEID(I,J)=SFMUL(THSE,SEIDME,ID,I,J)
          SLSLFEHDTMLFEHDS(I,J)=SFMUL(THSL,SLFEHDTML,FEHDS,I,J)
          SLAETAES(I,J)=SFMUL(THSL,AET,AES,I,J)
          DUM1E1(I,J)=2.D0*THSB*GTPEDGTPE(I,J)
     $               +2.D0*THSH*FTEETFTEES(I,J)
        END DO
      END DO
      SUTSUIDMU=TSFMUL(THSU,SUIDMU)
      SQTSQIDMQ=TSFMUL(THSQ,SQIDMQ)
      SDTSDIDMD=TSFMUL(THSD,SDIDMD)
      SLTSLIDML=TSFMUL(THSL,SLIDML)
      SETSEIDME=TSFMUL(THSE,SEIDME)
      SQTMQ=TSFMUL(THSQ,MQ)
      SUTMU=TSFMUL(THSU,MUP)
      SDTMD=TSFMUL(THSD,MD)
      SLTML=TSFMUL(THSL,ML)
      SETME=TSFMUL(THSE,ME)
!
      IF(THHH.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            SUFUQTFUQSMU(I,J)=SFMUL(THSU,FUQTFUQS,MUP,I,J)
            SQFUURSFUURTMQ(I,J)=SFMUL(THSQ,FUURSFUURT,MQ,I,J)
            SDMTSFDSMTSFDT(I,J)=SFMUL(THSD,MTSFDS,MTSFDT,I,J)
            SEMTSFESMTSFET(I,J)=SFMUL(THSE,MTSFES,MTSFET,I,J)
!
            SQFDDRSFDDRTMQ(I,J)=SFMUL(THSQ,FDDRSFDDRT,MQ,I,J)
            SDFDQTFDQSMD(I,J)=SFMUL(THSD,FDQTFDQS,MD,I,J)
            SLFEERSFEERTML(I,J)=SFMUL(THSL,FEERSFEERT,ML,I,J)
            SEFELTFELSME(I,J)=SFMUL(THSE,FELTFELS,ME,I,J)
            SUMTSFUSMTSFUT(I,J)=SFMUL(THSU,MTSFUS,MTSFUT,I,J)
          END DO
        END DO
        SUTSUFUQTFUQSMU=TSFMUL(THSU,SUFUQTFUQSMU)
        SQTSQFUURSFUURTMQ=TSFMUL(THSQ,SQFUURSFUURTMQ)
        SQTSUAUSAUT=TSFMUL(THSQ,SUAUSAUT)
        SQTSDMTSFDSMTSFDT=TSFMUL(THSQ,SDMTSFDSMTSFDT)
        SLTSEMTSFESMTSFET=TSFMUL(THSL,SEMTSFESMTSFET)
!
        SQTSQFDDRSFDDRTMQ=TSFMUL(THSQ,SQFDDRSFDDRTMQ)
        SDTSDFDQTFDQSMD=TSFMUL(THSD,SDFDQTFDQSMD)
        SLTSLFEERSFEERTML=TSFMUL(THSL,SLFEERSFEERTML)
        SETSEFELTFELSME=TSFMUL(THSE,SEFELTFELSME)
        SQTSUMTSFUSMTSFUT=TSFMUL(THSQ,SUMTSFUSMTSFUT)
        SQTSDADSADT=TSFMUL(THSQ,SDADSADT)
        SLTSEAESAET=TSFMUL(THSL,SEAESAET)
!
        BMHUPMT=3.D0/2.D0*(G(457)**2+G(458)**2)*G(61)
     $         -G(457)**2*G(62)
     $         -2.D0*G(457)**2*SUTSUIDMU
     $         +6.D0*SUTSUFUQTFUQSMU
     $         +G(457)**2*SQTSQIDMQ
     $         +6.D0*SQTSQFUURSFUURTMQ
     $         +G(457)**2*SDTSDIDMD
     $         -G(457)**2*SLTSLIDML
     $         +G(457)**2*SETSEIDME
     $         +6.D0*SQTSUAUSAUT
     $         +6.D0*SQTSDMTSFDSMTSFDT
     $         +2.D0*SLTSEMTSFESMTSFET
     $         -2.D0*THSH*MMUSQ*(THSB*MGTPHUSQ+3.D0*THSW*MGTHUSQ)
     $         -2.D0*THSH*(THSB*M1PM1PSQ*MGTPHUSQ+3.D0*THSW
     $          *M2PM2PSQ*MGTHUSQ)
     $         -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(61)
     $         +(6.D0*TFUDFU+THSB*THSH*MGTPHUSQ+3.D0*THSW*THSH
     $          *MGTHUSQ)*G(61)
!
        BMHDPMT=-G(457)**2*G(61)
     $         +3.D0/2.D0*(G(457)**2+G(458)**2)*G(62)
     $         +2.D0*G(457)**2*SUTSUIDMU
     $         -G(457)**2*SQTSQIDMQ
     $         +6.D0*SQTSQFDDRSFDDRTMQ
     $         -G(457)**2*SDTSDIDMD
     $         +6.D0*SDTSDFDQTFDQSMD
     $         +G(457)**2*SLTSLIDML
     $         +2.D0*SLTSLFEERSFEERTML
     $         -G(457)**2*SETSEIDME
     $         +2.D0*SETSEFELTFELSME
     $         +6.D0*SQTSUMTSFUSMTSFUT
     $         +6.D0*SQTSDADSADT
     $         +2.D0*SLTSEAESAET
     $         -2.D0*THSH*MMUSQ*(THSB*MGTPHDSQ+3.D0*THSW*MGTHDSQ)
     $         -2.D0*THSH*(THSB*M1PM1PSQ*MGTPHDSQ+3.D0*THSW
     $          *M2PM2PSQ*MGTHDSQ)
     $         -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(62)
     $         +(6.D0*TFDDFD+2.D0*TFEDFE+THSB*THSH*MGTPHDSQ
     $          +3.D0*THSW*THSH*MGTHDSQ)*G(62)
!
        F(61)=1.D0/16.D0/PI**2*BMHUPMT+1.D0/(16.D0*PI**2)**2
     $       *(B2HUM+G(398)*B2GRKMUM+G(398)*B2GRKMUM)
        F(62)=1.D0/16.D0/PI**2*BMHDPMT+1.D0/(16.D0*PI**2)**2
     $       *(B2HDM+G(398)*B2GRKMUM+G(398)*B2GRKMUM)
!
!Now the soft mass matrices
!
        DO I=1,3
          DO J=1,3
            SQMTSFUTMTSFUS(I,J)=SFMUL(THSQ,MTSFUT,MTSFUS,I,J)
            SQMTSFDTMTSFDS(I,J)=SFMUL(THSQ,MTSFDT,MTSFDS,I,J)
            SLMTSFETMTSFES(I,J)=SFMUL(THSL,MTSFET,MTSFES,I,J)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            BETA1MQ(I,J)=(1.D0/3.D0*G(457)**2*ID(I,J)+2.D0
     $                   *FUURSFUURT(I,J))*G(61)
     $                  +(-1.D0/3.D0*G(457)**2*ID(I,J)+2.D0
     $                   *FDDRSFDDRT(I,J))*G(62)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +(G(457)**2/18.D0+3.D0*G(458)**2/2.D0
     $                   +8.D0*G(459)**2/3.D0)*SQSQIDMQID(I,J)
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -1.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SUSUFUHUSMUFUHUT(I,J)
     $                  +2.D0*SDSDFDHDSMDFDHDT(I,J)
     $                  +2.D0*SUAUSAUT(I,J)
     $                  +2.D0*SUMTSFUSMTSFUT(I,J)
     $                  +2.D0*SDADSADT(I,J)
     $                  +2.D0*SDMTSFDSMTSFDT(I,J)
     $                  -2.D0/9.D0*THSB*M1PM1PSQ*GTPQGTPQD(I,J)
     $                  -6.D0*THSW*M2PM2PSQ*GTQGTQD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSQGTSQD(I,J)
     $                  -4.D0*THSH*MMUSQ*(FTUQSFTUQT(I,J)
     $                   +FTDQSFTDQT(I,J))
     $                  -(G(1)**2/10.D0+9.D0/2.D0*G(2)**2
     $                   +8.D0*G(3)**2)*MQ(I,J)
     $                  +SFMUL(THSQ,DUM1Q1,MQ,I,J)
     $                  +SFMUL(THSQ,MQ,DUM1Q1,I,J)
!
            BETA1MU(I,J)=(-4.D0/3.D0*G(457)**2*ID(I,J)+4.D0
     $                   *FUQTFUQS(I,J))*G(61)
     $                  +4.D0/3.D0*G(457)**2*ID(I,J)*G(62)
     $                  +8.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +8.D0/3.D0*(G(457)**2/3.D0+G(459)**2)
     $                   *SUSUIDMUID(I,J)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +4.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFUHUTMQFUHUS(I,J)
     $                  +4.D0*SQAUTAUS(I,J)
     $                  +4.D0*SQMTSFUTMTSFUS(I,J)
     $                  -32.D0/9.D0*THSB*M1PM1PSQ*GTPUDGTPU(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSUDGTSU(I,J)
     $                  -8.D0*THSH*MMUSQ*FTUUTFTUUS(I,J)
     $                  -(8.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MUP(I,J)
     $                  +SFMUL(THSU,DUM1U1,MUP,I,J)
     $                  +SFMUL(THSU,MUP,DUM1U1,I,J)
!
            BETA1MD(I,J)=2.D0/3.D0*G(457)**2*ID(I,J)*G(61)
     $                  +(-2.D0/3.D0*G(457)**2*ID(I,J)+4.D0
     $                   *FDQTFDQS(I,J))*G(62)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +2.D0/3.D0*(G(457)**2/3.D0+4.D0*G(459)**2)
     $                   *SDSDIDMDID(I,J)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFDHDTMQFDHDS(I,J)
     $                  +4.D0*SQADTADS(I,J)
     $                  +4.D0*SQMTSFDTMTSFDS(I,J)
     $                  -8.D0/9.D0*THSB*M1PM1PSQ*GTPDDGTPD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSDDGTSD(I,J)
     $                  -8.D0*THSH*MMUSQ*FTDDTFTDDS(I,J)
     $                  -(2.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MD(I,J)
     $                  +SFMUL(THSD,DUM1D1,MD,I,J)
     $                  +SFMUL(THSD,MD,DUM1D1,I,J)
!
            BETA1ML(I,J)=-G(457)**2*ID(I,J)*G(61)
     $                  +(G(457)**2*ID(I,J)+2.D0*FEERSFEERT(I,J))*G(62)
     $                  +2.D0*G(457)**2*ID(I,J)*SUTMU
     $                  -G(457)**2*ID(I,J)*SQTMQ
     $                  -G(457)**2*ID(I,J)*SDTMD
     $                  +G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/2.D0*(G(457)**2+3.D0*G(458)**2)
     $                   *SLSLIDMLID(I,J)
     $                  -G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SESEFEHDSMEFEHDT(I,J)
     $                  +2.D0*SEAESAET(I,J)
     $                  +2.D0*SEMTSFESMTSFET(I,J)
     $                  -2.D0*THSB*M1PM1PSQ*GTPLGTPLD(I,J)
     $                  -6.D0*THGL*M2PM2PSQ*GTLGTLD(I,J)
     $                  -4.D0*THSH*MMUSQ*FTELSFTELT(I,J)
     $                  -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*ML(I,J)
     $                  +SFMUL(THSL,DUM1L1,ML,I,J)
     $                  +SFMUL(THSL,ML,DUM1L1,I,J)
!
            BETA1ME(I,J)=2.D0*G(457)**2*ID(I,J)*G(61)
     $                  +(-2.D0*G(457)**2*ID(I,J)+4.D0
     $                   *FELTFELS(I,J))*G(62)
     $                  -4.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -2.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*G(457)**2*SESEIDMEID(I,J)
     $                  +4.D0*SLSLFEHDTMLFEHDS(I,J)
     $                  +4.D0*SLAETAES(I,J)
     $                  +4.D0*SLMTSFETMTSFES(I,J)
     $                  -8.D0*THSB*M1PM1PSQ*GTPEDGTPE(I,J)
     $                  -8.D0*THSH*MMUSQ*FTEETFTEES(I,J)
     $                  -18.D0*G(1)**2/5.D0*ME(I,J)
     $                  +SFMUL(THSE,DUM1E1,ME,I,J)
     $                  +SFMUL(THSE,ME,DUM1E1,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            F(62+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MQ(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MQM(I,J)
            F(71+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ML(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MLM(I,J)
            F(80+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MU(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MUM(I,J)
            F(89+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MD(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MDM(I,J)
            F(98+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ME(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MEM(I,J)
          END DO
        END DO
      ELSE
        DO I=1,3
          DO J=1,3
            SUSFUQTSFUQSMU(I,J)=SFMUL(THSU,SFUQTSFUQS,MUP,I,J)
            SQSFUURSSFUURTMQ(I,J)=SFMUL(THSQ,SFUURSSFUURT,MQ,I,J)
            SQCFDDRSCFDDRTMQ(I,J)=SFMUL(THSQ,CFDDRSCFDDRT,MQ,I,J)
            SDCFDQTCFDQSMD(I,J)=SFMUL(THSD,CFDQTCFDQS,MD,I,J)
            SLCFEERSCFEERTML(I,J)=SFMUL(THSL,CFEERSCFEERT,ML,I,J)
            SECFELTCFELSME(I,J)=SFMUL(THSE,CFELTCFELS,ME,I,J)
            SUTRIUTRIUD(I,J)=SFMUL(THSU,TRIU,TRIUD,I,J)
            SDTRIDTRIDD(I,J)=SFMUL(THSD,TRID,TRIDD,I,J)
            SETRIETRIED(I,J)=SFMUL(THSE,TRIE,TRIED,I,J)
          END DO
        END DO
        SUTSUSFUQTSFUQSMU=TSFMUL(THSU,SUSFUQTSFUQSMU)
        SQTSQSFUURSSFUURTMQ=TSFMUL(THSQ,SQSFUURSSFUURTMQ)
        SQTSQCFDDRSCFDDRTMQ=TSFMUL(THSQ,SQCFDDRSCFDDRTMQ)
        SDTSDCFDQTCFDQSMD=TSFMUL(THSD,SDCFDQTCFDQSMD)
        SLTSLCFEERSCFEERTML=TSFMUL(THSL,SLCFEERSCFEERTML)
        SETSECFELTCFELSME=TSFMUL(THSE,SECFELTCFELSME)
        SQTSUTRIUTRIUD=TSFMUL(THSQ,SUTRIUTRIUD)
        SQTSDTRIDTRIDD=TSFMUL(THSQ,SDTRIDTRIDD)
        SLTSETRIETRIED=TSFMUL(THSL,SETRIETRIED)
!      
        BMHUD=3.D0/2.D0*(G(543)**2+G(544)**2)*G(427)
     $       +2.D0*G(541)*SUTSUIDMU
     $       +6.D0*SUTSUSFUQTSFUQSMU
     $       -G(541)*SQTSQIDMQ
     $       +6.D0*SQTSQSFUURSSFUURTMQ
     $       +6.D0*SQTSQCFDDRSCFDDRTMQ
     $       -G(541)*SDTSDIDMD
     $       +6.D0*SDTSDCFDQTCFDQSMD
     $       +G(541)*SLTSLIDML
     $       +2.D0*SLTSLCFEERSCFEERTML
     $       -G(541)*SETSEIDME
     $       +2.D0*SETSECFELTCFELSME
     $       +6.D0*SQTSUTRIUTRIUD+6.D0*SQTSDTRIDTRIDD
     $        +2.D0*SLTSETRIETRIED
     $       -2.D0*THSH*MMUSQ*(THSB*(MSGTPHUSQ+MCGTPHDSQ)
     $        +3.D0*THSW*(MSGTHUSQ+MCGTHDSQ))
     $       -2.D0*THSH*(THSB*M1PM1PSQ*(MSGTPHUSQ+MCGTPHDSQ)
     $        +3.D0*THSW*M2PM2PSQ*(MSGTHUSQ+MCGTHDSQ))
     $       -1.D0/2.D0*(-4.D0*THSH*THSB*G(108)*G(287)*G(288)
     $        *(G(31)+(0.D0,1.D0)*G(599))-12.D0*THSH*THSW*G(108)
     $        *G(289)*G(290)*(G(32)+(0.D0,1.D0)*G(600)))
     $       -1.D0/2.D0*(-4.D0*THSH*THSB*G(108)*G(287)
     $        *G(288)*(G(31)-(0.D0,1.D0)*G(599))-12.D0*THSH*THSW
     $        *G(108)*G(289)*G(290)*(G(32)-(0.D0,1.D0)
     $        *G(600)))
     $       -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*G(427)
     $       +(6.D0*TLUDLU+THSB*THSH*MSGTPHUSQ+3.D0*THSW*THSH*MSGTHUSQ
     $        +6.D0*TLDDLD+2.D0*TLEDLE+THSB*THSH*MCGTPHDSQ
     $        +3.D0*THSW*THSH*MCGTHDSQ)*G(427)
!
        F(427)=1.D0/16.D0/PI**2*BMHUD+1.D0/(16.D0*PI**2)**2
     $        *(SINB**2*B2HUM+COSB**2*B2HDM+G(398)*B2GRKMUM
     $         +G(398)*B2GRKMUM+SINB*COSB*(BETA2B+BETA2B))
!
!Now the soft mass matrices
!
        DO I=1,3
          DO J=1,3
            SUTRIUSTRIUT(I,J)=SFMUL(THSU,TRIUS,TRIUT,I,J)
            SDTRIDSTRIDT(I,J)=SFMUL(THSD,TRIDS,TRIDT,I,J)
            SQTRIUTTRIUS(I,J)=SFMUL(THSQ,TRIUT,TRIUS,I,J)
            SQTRIDTTRIDS(I,J)=SFMUL(THSQ,TRIDT,TRIDS,I,J)
            SETRIESTRIET(I,J)=SFMUL(THSE,TRIES,TRIET,I,J)
            SLTRIETTRIES(I,J)=SFMUL(THSL,TRIET,TRIES,I,J)
          END DO
        END DO
!
        DO I=1,3
          DO J=1,3
            BETA1MQ(I,J)=(-1.D0/3.D0*G(541)*ID(I,J)+2.D0
     $                   *SFUURSSFUURT(I,J)+2.D0
     $                   *CFDDRSCFDDRT(I,J))*G(427)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +(G(457)**2/18.D0+3.D0*G(458)**2/2.D0
     $                   +8.D0*G(459)**2/3.D0)*SQSQIDMQID(I,J)
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -1.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SUSUFUHUSMUFUHUT(I,J)
     $                  +2.D0*SDSDFDHDSMDFDHDT(I,J)
     $                  +2.D0*SUTRIUSTRIUT(I,J)
     $                  +2.D0*SDTRIDSTRIDT(I,J)
     $                  -2.D0/9.D0*THSB*M1PM1PSQ*GTPQGTPQD(I,J)
     $                  -6.D0*THSW*M2PM2PSQ*GTQGTQD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSQGTSQD(I,J)
     $                  -4.D0*THSH*MMUSQ*(FTUQSFTUQT(I,J)
     $                   +FTDQSFTDQT(I,J))
     $                  -(G(1)**2/10.D0+9.D0/2.D0*G(2)**2
     $                   +8.D0*G(3)**2)*MQ(I,J)
     $                  +SFMUL(THSQ,DUM1Q1,MQ,I,J)
     $                  +SFMUL(THSQ,MQ,DUM1Q1,I,J)
!
            BETA1MU(I,J)=(4.D0/3.D0*G(541)*ID(I,J)+4.D0
     $                   *SFUQTSFUQS(I,J))*G(427)
     $                  +8.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +8.D0/3.D0*(G(457)**2/3.D0+G(459)**2)
     $                   *SUSUIDMUID(I,J)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +4.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFUHUTMQFUHUS(I,J)
     $                  +4.D0*SQTRIUTTRIUS(I,J)
     $                  -32.D0/9.D0*THSB*M1PM1PSQ*GTPUDGTPU(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSUDGTSU(I,J)
     $                  -8.D0*THSH*MMUSQ*FTUUTFTUUS(I,J)
     $                  -(8.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MUP(I,J)
     $                  +SFMUL(THSU,DUM1U1,MUP,I,J)
     $                  +SFMUL(THSU,MUP,DUM1U1,I,J)
!
            BETA1MD(I,J)=(-2.D0/3.D0*G(541)*ID(I,J)+4.D0
     $                   *CFDQTCFDQS(I,J))*G(427)
     $                  -4.D0/3.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SDTMD
     $                  +2.D0/3.D0*(G(457)**2/3.D0+4.D0*G(459)**2)
     $                   *SDSDIDMDID(I,J)
     $                  -2.D0/3.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0/3.D0*G(457)**2*ID(I,J)*SETME
     $                  +4.D0*SQSQFDHDTMQFDHDS(I,J)
     $                  +4.D0*SQTRIDTTRIDS(I,J)
     $                  -8.D0/9.D0*THSB*M1PM1PSQ*GTPDDGTPD(I,J)
     $                  -32.D0/3.D0*THGL*M3PM3PSQ*GTSDDGTSD(I,J)
     $                  -8.D0*THSH*MMUSQ*FTDDTFTDDS(I,J)
     $                  -(2.D0*G(1)**2/5.D0+8.D0*G(3)**2)*MD(I,J)
     $                  +SFMUL(THSD,DUM1D1,MD,I,J)
     $                  +SFMUL(THSD,MD,DUM1D1,I,J)
!
            BETA1ML(I,J)=(G(541)*ID(I,J)+2.D0*CFEERSCFEERT(I,J))
     $                   *G(427)
     $                  +2.D0*G(457)**2*ID(I,J)*SUTMU
     $                  -G(457)**2*ID(I,J)*SQTMQ
     $                  -G(457)**2*ID(I,J)*SDTMD
     $                  +G(457)**2*ID(I,J)*SLTML
     $                  +1.D0/2.D0*(G(457)**2+3.D0*G(458)**2)
     $                   *SLSLIDMLID(I,J)
     $                  -G(457)**2*ID(I,J)*SETME
     $                  +2.D0*SESEFEHDSMEFEHDT(I,J)
     $                  +2.D0*SETRIESTRIET(I,J)
     $                  -2.D0*THSB*M1PM1PSQ*GTPLGTPLD(I,J)
     $                  -6.D0*THGL*M2PM2PSQ*GTLGTLD(I,J)
     $                  -4.D0*THSH*MMUSQ*FTELSFTELT(I,J)
     $                  -(9.D0*G(1)**2/10.D0+9.D0*G(2)**2/2.D0)*ML(I,J)
     $                  +SFMUL(THSL,DUM1L1,ML,I,J)
     $                  +SFMUL(THSL,ML,DUM1L1,I,J)
!
            BETA1ME(I,J)=(-2.D0*G(541)*ID(I,J)
     $                   +4.D0*CFELTCFELS(I,J))
     $                   *G(427)
     $                  -4.D0*G(457)**2*ID(I,J)*SUTMU
     $                  +2.D0*G(457)**2*ID(I,J)*SQTMQ
     $                  +2.D0*G(457)**2*ID(I,J)*SDTMD
     $                  -2.D0*G(457)**2*ID(I,J)*SLTML
     $                  +2.D0*G(457)**2*ID(I,J)*SETME
     $                  +2.D0*G(457)**2*SESEIDMEID(I,J)
     $                  +4.D0*SLSLFEHDTMLFEHDS(I,J)
     $                  +4.D0*SLTRIETTRIES(I,J)
     $                  -8.D0*THSB*M1PM1PSQ*GTPEDGTPE(I,J)
     $                  -8.D0*THSH*MMUSQ*FTEETFTEES(I,J)
     $                  -18.D0*G(1)**2/5.D0*ME(I,J)
     $                  +SFMUL(THSE,DUM1E1,ME,I,J)
     $                  +SFMUL(THSE,ME,DUM1E1,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            F(62+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MQ(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MQM(I,J)
            F(71+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ML(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MLM(I,J)
            F(80+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MU(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MUM(I,J)
            F(89+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1MD(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MDM(I,J)
            F(98+(I-1)*3+J)=1.D0/16.D0/PI**2*BETA1ME(I,J)
     $                      +1.D0/(16.D0*PI**2)**2*B2MEM(I,J)
          END DO
        END DO
      END IF
!
      IF(THHH.EQ.0)THEN
!
!Finally we have the running of the Higgs Quartic Coupling and SM VEV.
!Programmed here is the MS-bar running. It therefore needs the MS-bar
!gauge and Yukawas.
!The gauge couplings and Yukawas needed to be converted to MS-bar using
!the Martin and Vaughn conversion in hep-ph/9308222.
!The following is after the conversion, so all Yukawas and Gauge
!couplings are still in the DR-bar scheme.
!
        BETALAM1=12*G(429)**2-(9.D0/5.D0*G(1)**2+9.D0*G(2)**2)
     $           *G(429)+9.D0/4.D0*(3.D0/25.D0*G(1)**4+2.D0/5.D0
     $           *G(1)**2*G(2)**2+G(2)**4)+4.D0*Y2*G(429)-4*H
        IF(SW2LP.EQ.1)THEN
          BETALAM2=-78.D0*G(429)**3+18.D0*(3.D0/5.D0*G(1)**2
     $            +3.D0*G(2)**2)*G(429)**2-((265.D0/8.D0-10*NG)
     $            *G(2)**4-117.D0/20.D0*G(1)**2*G(2)**2
     $            -9.D0/25.D0*(229.D0/24.D0+50.D0/9.D0*NG)*G(1)**4)
     $            *G(429)+(473.D0/8.D0-8.D0*NG)*G(2)**6-3.D0/5.D0
     $            *(121.D0/24.D0+8.D0/3.D0*NG)*G(1)**2*G(2)**4
     $            -9.D0/25.D0*(239.D0/24.D0+40.D0/9.D0*NG)
     $            *G(1)**4*G(2)**2-27.D0/125.D0*(59.D0/24.D0
     $            +40.D0/9.D0*NG)*G(1)**6+(-14.D0/5.D0*G(1)**2
     $            +18.D0*G(2)**2-128.D0*G(3)**2)*TLYUDLYU2
     $            +(34.D0/5.D0*G(1)**2+18.D0*G(2)**2-128.D0
     $            *G(3)**2)*TLYDDLYD2+(-42.D0/5.D0*G(1)**2
     $            +6.D0*G(2)**2)*TLYEDLYE2-3.D0/2.D0*G(2)**4
     $            *Y2+G(429)*((83.D0/10.D0*G(1)**2+27.D0/2.D0
     $            *G(2)**2+112.D0*G(3)**2)*TLYUDLYU+(-1.D0/10.D0
     $            *G(1)**2+27.D0/2.D0*G(2)**2+112.D0*G(3)**2)
     $            *TLYDDLYD+(93.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $            *TLYEDLYE)+3.D0/5.D0*G(1)**2*((-57.D0/10.D0
     $            *G(1)**2+21.D0*G(2)**2)*TLYUDLYU+(3.D0/2.D0
     $            *G(1)**2+9.D0*G(2)**2)*TLYDDLYD+(-15.D0/2.D0
     $            *G(1)**2+11.D0*G(2)**2)*TLYEDLYE)-24.D0
     $            *G(429)**2*Y2-G(429)*H+6.D0*G(429)
     $            *TLYUDLYULYDDLYD+20.D0*(3.D0*TLYUDLYU3+3.D0*TLYDDLYD3
     $            +TLYEDLYE3)-12.D0*TDUMLUD
        END IF
!
        F(429)=1.D0/(16.D0*PI**2)*BETALAM1+1.D0/(16.D0*PI**2)**2
     $         *BETALAM2
!
!Calculate the betas for the standard model vev.
!As with lambda this is the MS-bar running with DR-bar inputs except
!v and lambda
!
      BETAVEV1=9.D0/4.D0*(1.D0/5.D0*G(1)**2+G(2)**2)-Y2
      IF(SW2LP.EQ.1)THEN
        BETAVEV2=-3.D0/2.D0*G(429)**2-(83.D0/40.D0*G(1)**2+27.D0/8.D0
     $           *G(2)**2+28.D0*G(3)**2)*TYUDYU-(-1.D0/40.D0*G(1)**2
     $           +27.D0/8.D0*G(2)**2+28.D0*G(3)**2)*TYDDYD
     $           -(93.D0/40.D0*G(1)**2+9.D0/8.D0*G(2)**2)*TYEDYE+CHI4
     $           -27.D0/80.D0*G(1)**2*G(2)**2-(93.D0/800.D0+1.D0/2.D0
     $           *NG)*G(1)**4+(463.D0/32.D0-5.D0/2.D0*NG)*G(2)**4
      END IF
!
      F(428)=G(428)*(1.D0/(16.D0*PI**2)*BETAVEV1+1.D0/(16.D0*PI**2)**2
     $      *BETAVEV2)
!
      END IF
!
!Now the tilde terms are calculated.
!
!
      DO I=1,3
        DO J=1,3
          DUMTUQ1(I,J)=4.D0/9.D0*THSB*SUGTPUSGTPUT(I,J)
     $                 +4.D0/3.D0*THGL*SUGTSUSGTSUT(I,J)
     $                 +THSH*SQFTUQDFTUQ(I,J)
          DUMTUQ2(I,J)=SUFTUUFTUUD(I,J)+SQFTUQDFTUQ(I,J)
          IF(THHH.EQ.0)THEN
            DUMTUQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(289)
     $                   +1.D0/3.D0*THSB*GTPQS(I,J)*G(287)
          ELSE
            DUMTUQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(204)
     $                   +1.D0/3.D0*THSB*GTPQS(I,J)*G(184)
          END IF
          DUMTUQ4(I,J)=3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $                 +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $                 +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
          DUMTUQ5(I,J)=FTUQFTUQD(I,J)+FTDQFTDQD(I,J)
!
          DUMTDQ1(I,J)=1.D0/9.D0*THSB*SDGTPDSGTPDT(I,J)
     $                 +4.D0/3.D0*THGL*SDGTSDSGTSDT(I,J)
     $                 +THSH*SQFTDQDFTDQ(I,J)
          DUMTDQ2(I,J)=3.D0*SDFTDDFTDDD(I,J)+SEFTEEFTEED(I,J)
     $                 +3.D0*SQFTDQDFTDQ(I,J)+SLFTELDFTEL(I,J)
          IF(THHH.EQ.0)THEN
            DUMTDQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(290)
     $                   -1.D0/3.D0*THSB*GTPQS(I,J)*G(288)
          ELSE
            DUMTDQ3(I,J)=-3.D0*THSW*GTQS(I,J)*G(205)
     $                   -1.D0/3.D0*THSB*GTPQS(I,J)*G(185)
          END IF
          DUMTDQ4(I,J)=3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $                 +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $                 +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
          DUMTDQ5(I,J)=FTUQFTUQD(I,J)+FTDQFTDQD(I,J)
!
          DUMTEL1(I,J)=THSB*SEGTPESGTPET(I,J)
     $                 +THSH*SLFTELDFTEL(I,J)
          IF(THHH.EQ.0)THEN
            DUMTEL3(I,J)=-3.D0*THSW*GTLS(I,J)*G(290)
     $                   +THSB*GTPLS(I,J)*G(288)
          ELSE
            DUMTEL3(I,J)=-3.D0*THSW*GTLS(I,J)*G(205)
     $                   +THSB*GTPLS(I,J)*G(185)
          END IF
          DUMTEL4(I,J)=3.D0*THSW*GTLSGTLT(I,J)
     $                 +THSB*GTPLSGTPLT(I,J)
          DUMTEL5(I,J)=FTELFTELD(I,J)
!
          IF(THHH.EQ.0)THEN
            DUMTUU2(I,J)=LULUD(I,J)+LDLDD(I,J)
          ELSE
            DUMTUU2(I,J)=FUFUD(I,J)+FDFDD(I,J)
          END IF
          DUMTUU3(I,J)=THSH*(SUFTUUFTUUD(I,J)+SDFTDDFTDDD(I,J))
     $                 +3.D0/2.D0*THSW*SQGTQTGTQS(I,J)
     $                 +1.D0/18.D0*THSB*SQGTPQTGTPQS(I,J)
     $                 +8.D0/3.D0*THGL*SQGTSQTGTSQS(I,J)
          DUMTUU4(I,J)=8.D0/9.D0*THSB*GTPUTGTPUS(I,J)
     $                 +8.D0/3.D0*THGL*GTSUTGTSUS(I,J)
!
          DUMTDD2(I,J)=DUMTUU2(I,J)
          DUMTDD3(I,J)=DUMTUU3(I,J)
          DUMTDD4(I,J)=2.D0/9.D0*THSB*GTPDTGTPDS(I,J)
     $                 +8.D0/3.D0*THGL*GTSDTGTSDS(I,J)
!
          DUMTEE3(I,J)=THSH*SEFTEEFTEED(I,J)
     $                 +3.D0/2.D0*THSW*SLGTLTGTLS(I,J)
     $                 +1.D0/2.D0*THSB*SLGTPLTGTPLS(I,J)
        END DO
      END DO
      TDUMTUQ2=DTRACE(DUMTUQ2)
      TDUMTDQ2=DTRACE(DUMTDQ2)
      TDUMTEL2=TDUMTDQ2
      TDUMTUU1=TDUMTUQ2
      TDUMTDD1=TDUMTDQ2
      TDUMTEE1=TDUMTEL2
      DO I=1,3
        DO J=1,3
          IF(THHH.EQ.0)THEN
            BFTUQ(I,J)=FTUQLUDLU(I,J)
     $                 +DMATMUL(0,FTUQ,DUMTUQ1,I,J)
     $                 +3.D0/2.D0*THSH*FTUQ(I,J)*TDUMTUQ2
     $                 +1.D0/4.D0*THSH*FTUQ(I,J)
     $                  *(3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $                 +DMATMUL(0,DUMTUQ3,LU,I,J)
     $                 -4.D0/9.D0*THSB
     $                  *DMATMUL(0,GTPQS,SUFTUUGTPUT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *DMATMUL(0,GTSQS,SUFTUUGTSUT,I,J)
     $                 +SFMUL(THSQ,DUMTUQ4,FTUQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTUQ5,FTUQ,I,J)
     $                 -FTUQ(I,J)*(5.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDQ(I,J)=FTDQLDDLD(I,J)
     $                 +DMATMUL(0,FTDQ,DUMTDQ1,I,J)
     $                 +1.D0/2.D0*THSH*FTDQ(I,J)*TDUMTDQ2
     $                 +1.D0/4.D0*THSH*FTDQ(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +DMATMUL(0,DUMTDQ3,LD,I,J)
     $                 +2.D0/9.D0*THSB
     $                  *DMATMUL(0,GTPQS,SDFTDDGTPDT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *DMATMUL(0,GTSQS,SDFTDDGTSDT,I,J)
     $                 +SFMUL(THSQ,DUMTDQ4,FTDQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTDQ5,FTDQ,I,J)
     $                 -FTDQ(I,J)*(13.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEL(I,J)=FTELLEDLE(I,J)
     $                 +DMATMUL(0,FTEL,DUMTEL1,I,J)
     $                 +1.D0/2.D0*THSH*FTEL(I,J)*TDUMTEL2
     $                 +1.D0/4.D0*THSH*FTEL(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +DMATMUL(0,DUMTEL3,LE,I,J)
     $                 -2.D0*THSB*DMATMUL(0,GTPLS,SEFTEEGTPET,I,J)
     $                 +1.D0/2.D0*SFMUL(THSL,DUMTEL4,FTEL,I,J)
     $                 +THSH*SFMUL(THSL,DUMTEL5,FTEL,I,J)
     $                 -FTEL(I,J)*(9.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BFTUU(I,J)=3.D0/2.D0*THSH*FTUU(I,J)*TDUMTUU1
     $                 +1.D0/4.D0*THSH*FTUU(I,J)
     $                  *(3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTUU2,FTUU,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTUU3,FTUU,I,J)
     $                 -4.D0/9.D0*THSB*DMATMUL(0,SQGTPQTFTUQ,GTPUS,I,J)
     $                 -16.D0/3.D0*THGL*DMATMUL(0,SQGTSQTFTUQ,GTSUS,I,J)
     $                 -4.D0/3.D0*THSB*G(287)*LUGTPUS(I,J)
     $                 +2.D0*THSH*DMATMUL(0,SUFTUUFTUUD,FTUU,I,J)
     $                 +SFMUL(THSU,FTUU,DUMTUU4,I,J)
     $                 -FTUU(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDD(I,J)=1.D0/2.D0*THSH*FTDD(I,J)*TDUMTDD1
     $                 +1.D0/4.D0*THSH*FTDD(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTDD2,FTDD,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTDD3,FTDD,I,J)
     $                 +2.D0/9.D0*THSB*DMATMUL(0,SQGTPQTFTDQ,GTPDS,I,J)
     $                 -16.D0/3.D0*THGL*DMATMUL(0,SQGTSQTFTDQ,GTSDS,I,J)
     $                 -2.D0/3.D0*THSB*G(288)*LDGTPDS(I,J)
     $                 +2.D0*THSH*DMATMUL(0,SDFTDDFTDDD,FTDD,I,J)
     $                 +SFMUL(THSD,FTDD,DUMTDD4,I,J)
     $                 -FTDD(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEE(I,J)=1.D0/2.D0*THSH*FTEE(I,J)*TDUMTEE1
     $                 +1.D0/4.D0*THSH*FTEE(I,J)
     $                  *(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,LELED,FTEE,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTEE3,FTEE,I,J)
     $                 -2.D0*THSB*DMATMUL(0,SLGTPLTFTEL,GTPES,I,J)
     $                 -2.D0*THSB*G(288)*LEGTPES(I,J)
     $                 +2.D0*THSH*DMATMUL(0,SEFTEEFTEED,FTEE,I,J)
     $                 +2.D0*THSB*DMATMUL(0,SEFTEEGTPET,GTPES,I,J)
     $                 -FTEE(I,J)*(9.D0/10.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2)
          ELSE
            BFTUQ(I,J)=FTUQFUDFU(I,J)
     $                 +DMATMUL(0,FTUQ,DUMTUQ1,I,J)
     $                 +3.D0/2.D0*THSH*FTUQ(I,J)*TDUMTUQ2
     $                 +1.D0/4.D0*THSH*FTUQ(I,J)
     $                  *(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $                 +DMATMUL(0,DUMTUQ3,FU,I,J)
     $                 -4.D0/9.D0*THSB
     $                  *DMATMUL(0,GTPQS,SUFTUUGTPUT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *DMATMUL(0,GTSQS,SUFTUUGTSUT,I,J)
     $                 +SFMUL(THSQ,DUMTUQ4,FTUQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTUQ5,FTUQ,I,J)
     $                 -FTUQ(I,J)*(5.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDQ(I,J)=FTDQFDDFD(I,J)
     $                 +DMATMUL(0,FTDQ,DUMTDQ1,I,J)
     $                 +1.D0/2.D0*THSH*FTDQ(I,J)*TDUMTDQ2
     $                 +1.D0/4.D0*THSH*FTDQ(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +DMATMUL(0,DUMTDQ3,FD,I,J)
     $                 +2.D0/9.D0*THSB
     $                  *DMATMUL(0,GTPQS,SDFTDDGTPDT,I,J)
     $                 -16.D0/3.D0*THGL
     $                  *DMATMUL(0,GTSQS,SDFTDDGTSDT,I,J)
     $                 +SFMUL(THSQ,DUMTDQ4,FTDQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTDQ5,FTDQ,I,J)
     $                 -FTDQ(I,J)*(13.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEL(I,J)=FTELFEDFE(I,J)
     $                 +DMATMUL(0,FTEL,DUMTEL1,I,J)
     $                 +1.D0/2.D0*THSH*FTEL(I,J)*TDUMTEL2
     $                 +1.D0/4.D0*THSH*FTEL(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +DMATMUL(0,DUMTEL3,FE,I,J)
     $                 -2.D0*THSB*DMATMUL(0,GTPLS,SEFTEEGTPET,I,J)
     $                 +1.D0/2.D0*SFMUL(THSL,DUMTEL4,FTEL,I,J)
     $                 +THSH*SFMUL(THSL,DUMTEL5,FTEL,I,J)
     $                 -FTEL(I,J)*(9.D0/4.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BFTUU(I,J)=3.D0/2.D0*THSH*FTUU(I,J)*TDUMTUU1
     $                 +1.D0/4.D0*THSH*FTUU(I,J)
     $                  *(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTUU2,FTUU,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTUU3,FTUU,I,J)
     $                 -4.D0/9.D0*THSB*DMATMUL(0,SQGTPQTFTUQ,GTPUS,I,J)
     $                 -16.D0/3.D0*THGL*DMATMUL(0,SQGTSQTFTUQ,GTSUS,I,J)
     $                 -4.D0/3.D0*THSB*G(184)*FUGTPUS(I,J)
     $                 +2.D0*THSH*DMATMUL(0,SUFTUUFTUUD,FTUU,I,J)
     $                 +SFMUL(THSU,FTUU,DUMTUU4,I,J)
     $                 -FTUU(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTDD(I,J)=1.D0/2.D0*THSH*FTDD(I,J)*TDUMTDD1
     $                 +1.D0/4.D0*THSH*FTDD(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTDD2,FTDD,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTDD3,FTDD,I,J)
     $                 +2.D0/9.D0*THSB*DMATMUL(0,SQGTPQTFTDQ,GTPDS,I,J)
     $                 -16.D0/3.D0*THGL*DMATMUL(0,SQGTSQTFTDQ,GTSDS,I,J)
     $                 -2.D0/3.D0*THSB*G(185)*FDGTPDS(I,J)
     $                 +2.D0*THSH*DMATMUL(0,SDFTDDFTDDD,FTDD,I,J)
     $                 +SFMUL(THSD,FTDD,DUMTDD4,I,J)
     $                 -FTDD(I,J)*(1.D0/2.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2+4.D0*G(3)**2)
!
            BFTEE(I,J)=1.D0/2.D0*THSH*FTEE(I,J)*TDUMTEE1
     $                 +1.D0/4.D0*THSH*FTEE(I,J)
     $                  *(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,FEFED,FTEE,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,DUMTEE3,FTEE,I,J)
     $                 -2.D0*THSB*DMATMUL(0,SLGTPLTFTEL,GTPES,I,J)
     $                 -2.D0*THSB*G(185)*FEGTPES(I,J)
     $                 +2.D0*THSH*DMATMUL(0,SEFTEEFTEED,FTEE,I,J)
     $                 +2.D0*THSB*DMATMUL(0,SEFTEEGTPET,GTPES,I,J)
     $                 -FTEE(I,J)*(9.D0/10.D0*G(1)**2
     $                  +9.D0/2.D0*G(2)**2)
          END IF
!
!Convert into form readable by RKSTP. Two loop running is based
!on MV notation RGEs so I must take the transpose
!
          F(232+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BFTUQ(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMU(J,I))
          F(241+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BFTDQ(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMD(J,I))
          F(250+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BFTEL(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YME(J,I))
          F(259+(I-1)*3+J)=(1.D0/16.D0/PI**2*BFTUU(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMU(J,I))*THSU(J)
          F(268+(I-1)*3+J)=(1.D0/16.D0/PI**2*BFTDD(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YMD(J,I))*THSD(J)
          F(277+(I-1)*3+J)=(1.D0/16.D0/PI**2*BFTEE(I,J)
     $                     +1.D0/(16.D0*PI**2)**2*B2YME(J,I))*THSE(J)
        END DO
      END DO
!
!That's the yukawas done, now for the gauge couplings
!
      DO I=1,3
        DO J=1,3
          DUMTQ1(I,J)=3.D0*SQGTQDGTQ(I,J)+SLGTLDGTL(I,J)
          IF(THHH.EQ.0)THEN
            DUMTQ2(I,J)=LUSLUT(I,J)+LDSLDT(I,J)
          ELSE
            DUMTQ2(I,J)=FUSFUT(I,J)+FDSFDT(I,J)
          END IF
          DUMTQ3(I,J)=THSH*(SUFTUUSFTUUT(I,J)+SDFTDDSFTDDT(I,J))
     $                +3.D0/2.D0*THSW*SQGTQDGTQ(I,J)
     $                +1.D0/18.D0*THSB*SQGTPQDGTPQ(I,J)
     $                +8.D0/3.D0*THGL*SQGTSQDGTSQ(I,J)
          DUMTQ4(I,J)=3.D0*THSW*GTQGTQD(I,J)+1.D0/9.D0*THSB
     $                 *GTPQGTPQD(I,J)+16.D0/3.D0*THGL*GTSQGTSQD(I,J) 
          DUMTQ5(I,J)=FTUQSFTUQT(I,J)+FTDQSFTDQT(I,J)
!
          DUMTL2(I,J)=THSH*SEFTEESFTEET(I,J)
     $                +3.D0/2.D0*THSW*SLGTLDGTL(I,J)
     $                +1.D0/2.D0*THSB*SLGTPLDGTPL(I,J)
          DUMTL3(I,J)=3.D0*THSW*GTLGTLD(I,J)+THSB*GTPLGTPLD(I,J)
!
          DUMTPQ1(I,J)=1.D0/3.D0*SQGTPQDGTPQ(I,J)+SLGTPLDGTPL(I,J)
     $                 +8.D0/3.D0*SUGTPUGTPUD(I,J)
     $                 +2.D0/3.D0*SDGTPDGTPDD(I,J)
     $                 +2.D0*SEGTPEGTPED(I,J)
          DUMTPQ2(I,J)=DUMTQ2(I,J)
          DUMTPQ3(I,J)=DUMTQ3(I,J)
          DUMTPQ4(I,J)=DUMTQ4(I,J)
          DUMTPQ5(I,J)=DUMTQ5(I,J)
!
          DUMTPL2(I,J)=DUMTL2(I,J)
          DUMTPL3(I,J)=DUMTL3(I,J)
!
          DUMTPU1(I,J)=4.D0/9.D0*THSB*SUGTPUGTPUD(I,J)
     $                 +4.D0/3.D0*THGL*SUGTSUGTSUD(I,J)
     $                 +THSH*SQFTUQTFTUQS(I,J)
          DUMTPU3(I,J)=8.D0/9.D0*THSB*GTPUDGTPU(I,J)
     $                 +8.D0/3.D0*THGL*GTSUDGTSU(I,J)
!
          DUMTPD1(I,J)=1.D0/9.D0*THSB*SDGTPDGTPDD(I,J)
     $                 +4.D0/3.D0*THGL*SDGTSDGTSDD(I,J)
     $                 +THSH*SQFTDQTFTDQS(I,J)
          DUMTPD3(I,J)=2.D0/9.D0*THSB*GTPDDGTPD(I,J)
     $                 +8.D0/3.D0*THGL*GTSDDGTSD(I,J)
!
          DUMTPE1(I,J)=THSB*SEGTPEGTPED(I,J)+THSH*SLFTELTFTELS(I,J)
        END DO
      END DO
      TDUMTQ1=DTRACE(DUMTQ1)
      TDUMTL1=TDUMTQ1
      TDUMTPQ1=DTRACE(DUMTPQ1)
      TDUMTPL1=TDUMTPQ1
      TDUMTPU2=TDUMTPQ1
      TDUMTPD2=TDUMTPU2
      TDUMTPE2=TDUMTPD2
      DO I=1,3
        DO J=1,3
          IF(THHH.EQ.0)THEN
            BGTQ(I,J)=1.D0/2.D0*THSW*GTQ(I,J)*TDUMTQ1
     $                +1.D0/2.D0*THSH*THSW*GTQ(I,J)
     $                 *(MSGTHUSQ+MCGTHDSQ)
     $                +1.D0/2.D0*DMATMUL(0,GTQ,DUMTQ2,I,J)
     $                +1.D0/2.D0*DMATMUL(0,GTQ,DUMTQ3,I,J)
     $                -2.D0*THSH*(G(289)*FTUQSLUT(I,J)
     $                 +G(290)*FTDQSLDT(I,J))
     $                +1.D0/2.D0*SFMUL(THSQ,DUMTQ4,GTQ,I,J)
     $                +THSH*SFMUL(THSQ,DUMTQ5,GTQ,I,J)
     $                -GTQ(I,J)*(1.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTL(I,J)=1.D0/2.D0*THSW*GTL(I,J)*TDUMTL1
     $                +1.D0/2.D0*THSH*THSW*GTL(I,J)
     $                 *(MSGTHUSQ+MCGTHDSQ)
     $                +1.D0/2.D0*DMATMUL(0,GTL,LESLET,I,J)
     $                +1.D0/2.D0*DMATMUL(0,GTL,DUMTL2,I,J)
     $                -2.D0*THSH*G(290)*FTELSLET(I,J)
     $                +1.D0/2.D0*SFMUL(THSL,DUMTL3,GTL,I,J)
     $                +THSH*SFMUL(THSL,FTELSFTELT,GTL,I,J)
     $                -GTL(I,J)*(9.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2)
!
            BGTPQ(I,J)=1.D0/2.D0*THSB*GTPQ(I,J)*TDUMTPQ1
     $                 +1.D0/2.D0*THSH*THSB*GTPQ(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,GTPQ,DUMTPQ2,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTPQ,DUMTPQ3,I,J)
     $                 +4.D0*THSH*(-2.D0
     $                  *DMATMUL(0,FTUQS,SUGTPUFTUUT,I,J)
     $                  +DMATMUL(0,FTDQS,SDGTPDFTDDT,I,J))
     $                 +6.D0*THSH*(G(287)*FTUQSLUT(I,J)
     $                  -G(288)*FTDQSLDT(I,J))
     $                 +1.D0/2.D0*SFMUL(THSQ,DUMTPQ4,GTPQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTPQ5,GTPQ,I,J)
     $                 -GTPQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTPL(I,J)=1.D0/2.D0*THSB*GTPL(I,J)*TDUMTPL1
     $                 +1.D0/2.D0*THSH*THSB*GTPL(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,GTPL,LESLET,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTPL,DUMTPL2,I,J)
     $                 -4.D0*THSH*DMATMUL(0,FTELS,SEGTPEFTEET,I,J)
     $                 +2.D0*THSH*G(288)*FTELSLET(I,J)
     $                 +1.D0/2.D0*SFMUL(THSL,DUMTPL3,GTPL,I,J)
     $                 +THSH*SFMUL(THSL,FTELSFTELT,GTPL,I,J)
     $                 -GTPL(I,J)*(9.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BGTPU(I,J)=LUTLUSGTPU(I,J)
     $                 +DMATMUL(0,DUMTPU1,GTPU,I,J)
     $                 +1.D0/2.D0*THSB*GTPU(I,J)*TDUMTPU2
     $                 +1.D0/2.D0*THSB*THSH*GTPU(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 -3.D0*THSB*THSH*LUTFTUUS(I,J)*G(287)
     $                 -THSH*DMATMUL(0,SQFTUQTGTPQ,FTUUS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SUGTPUFTUUT,FTUUS,I,J)
     $                 +SFMUL(THSU,GTPU,DUMTPU3,I,J)
     $                 -GTPU(I,J)*(4.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPD(I,J)=LDTLDSGTPD(I,J)
     $                 +DMATMUL(0,DUMTPD1,GTPD,I,J)
     $                 +1.D0/2.D0*THSB*GTPD(I,J)*TDUMTPD2
     $                 +1.D0/2.D0*THSB*THSH*GTPD(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 -6.D0*THSB*THSH*LDTFTDDS(I,J)*G(288)
     $                 +2.D0*THSH*DMATMUL(0,SQFTDQTGTPQ,FTDDS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SDGTPDFTDDT,FTDDS,I,J)
     $                 +SFMUL(THSD,GTPD,DUMTPD3,I,J)
     $                 -GTPD(I,J)*(1.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPE(I,J)=LETLESGTPE(I,J)
     $                 +DMATMUL(0,DUMTPE1,GTPE,I,J)
     $                 +1.D0/2.D0*THSB*GTPE(I,J)*TDUMTPE2
     $                 +1.D0/2.D0*THSB*THSH*GTPE(I,J)
     $                  *(MSGTPHUSQ+MCGTPHDSQ)
     $                 -2.D0*THSB*THSH*LETFTEES(I,J)*G(288)
     $                 -2.D0*THSH*DMATMUL(0,SLFTELTGTPL,FTEES,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SEGTPEFTEET,FTEES,I,J)
     $                 +2.D0*THSB*DMATMUL(0,SEGTPEGTPED,GTPE,I,J)
     $                 -GTPE(I,J)*9.D0/5.D0*G(1)**2
          ELSE 
            BGTQ(I,J)=1.D0/2.D0*THSW*GTQ(I,J)*TDUMTQ1
     $                +1.D0/2.D0*THSH*THSW*GTQ(I,J)
     $                 *(MGTHUSQ+MGTHDSQ)
     $                +1.D0/2.D0*DMATMUL(0,GTQ,DUMTQ2,I,J)
     $                +1.D0/2.D0*DMATMUL(0,GTQ,DUMTQ3,I,J)
     $                -2.D0*THSH*(G(204)*FTUQSFUT(I,J)
     $                 +G(205)*FTDQSFDT(I,J))
     $                +1.D0/2.D0*SFMUL(THSQ,DUMTQ4,GTQ,I,J)
     $                +THSH*SFMUL(THSQ,DUMTQ5,GTQ,I,J)
     $                -GTQ(I,J)*(1.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTL(I,J)=1.D0/2.D0*THSW*GTL(I,J)*TDUMTL1
     $                +1.D0/2.D0*THSH*THSW*GTL(I,J)
     $                 *(MGTHUSQ+MGTHDSQ)
     $                +1.D0/2.D0*DMATMUL(0,GTL,FESFET,I,J)
     $                +1.D0/2.D0*DMATMUL(0,GTL,DUMTL2,I,J)
     $                -2.D0*THSH*G(205)*FTELSFET(I,J)
     $                +1.D0/2.D0*SFMUL(THSL,DUMTL3,GTL,I,J)
     $                +THSH*SFMUL(THSL,FTELSFTELT,GTL,I,J)
     $                -GTL(I,J)*(9.D0/20.D0*G(1)**2
     $                 +33.D0/4.D0*G(2)**2)
!
            BGTPQ(I,J)=1.D0/2.D0*THSB*GTPQ(I,J)*TDUMTPQ1
     $                 +1.D0/2.D0*THSH*THSB*GTPQ(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,GTPQ,DUMTPQ2,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTPQ,DUMTPQ3,I,J)
     $                 +4.D0*THSH*(-2.D0
     $                  *DMATMUL(0,FTUQS,SUGTPUFTUUT,I,J)
     $                  +DMATMUL(0,FTDQS,SDGTPDFTDDT,I,J))
     $                 +6.D0*THSH*(G(184)*FTUQSFUT(I,J)
     $                  -G(185)*FTDQSFDT(I,J))
     $                 +1.D0/2.D0*SFMUL(THSQ,DUMTPQ4,GTPQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTPQ5,GTPQ,I,J)
     $                 -GTPQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+4.D0*G(3)**2)
!
            BGTPL(I,J)=1.D0/2.D0*THSB*GTPL(I,J)*TDUMTPL1
     $                 +1.D0/2.D0*THSH*THSB*GTPL(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 +1.D0/2.D0*DMATMUL(0,GTPL,FESFET,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTPL,DUMTPL2,I,J)
     $                 -4.D0*THSH*DMATMUL(0,FTELS,SEGTPEFTEET,I,J)
     $                 +2.D0*THSH*G(185)*FTELSFET(I,J)
     $                 +1.D0/2.D0*SFMUL(THSL,DUMTPL3,GTPL,I,J)
     $                 +THSH*SFMUL(THSL,FTELSFTELT,GTPL,I,J)
     $                 -GTPL(I,J)*(9.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2)
!
            BGTPU(I,J)=FUTFUSGTPU(I,J)
     $                 +DMATMUL(0,DUMTPU1,GTPU,I,J)
     $                 +1.D0/2.D0*THSB*GTPU(I,J)*TDUMTPU2
     $                 +1.D0/2.D0*THSB*THSH*GTPU(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 -3.D0*THSB*THSH*FUTFTUUS(I,J)*G(184)
     $                 -THSH*DMATMUL(0,SQFTUQTGTPQ,FTUUS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SUGTPUFTUUT,FTUUS,I,J)
     $                 +SFMUL(THSU,GTPU,DUMTPU3,I,J)
     $                 -GTPU(I,J)*(4.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPD(I,J)=FDTFDSGTPD(I,J)
     $                 +DMATMUL(0,DUMTPD1,GTPD,I,J)
     $                 +1.D0/2.D0*THSB*GTPD(I,J)*TDUMTPD2
     $                 +1.D0/2.D0*THSB*THSH*GTPD(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 -6.D0*THSB*THSH*FDTFTDDS(I,J)*G(185)
     $                 +2.D0*THSH*DMATMUL(0,SQFTDQTGTPQ,FTDDS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SDGTPDFTDDT,FTDDS,I,J)
     $                 +SFMUL(THSD,GTPD,DUMTPD3,I,J)
     $                 -GTPD(I,J)*(1.D0/5.D0*G(1)**2+4.D0*G(3)**2)
!
            BGTPE(I,J)=FETFESGTPE(I,J)
     $                 +DMATMUL(0,DUMTPE1,GTPE,I,J)
     $                 +1.D0/2.D0*THSB*GTPE(I,J)*TDUMTPE2
     $                 +1.D0/2.D0*THSB*THSH*GTPE(I,J)
     $                  *(MGTPHUSQ+MGTPHDSQ)
     $                 -2.D0*THSB*THSH*FETFTEES(I,J)*G(185)
     $                 -2.D0*THSH*DMATMUL(0,SLFTELTGTPL,FTEES,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SEGTPEFTEET,FTEES,I,J)
     $                 +2.D0*THSB*DMATMUL(0,SEGTPEGTPED,GTPE,I,J)
     $                 -GTPE(I,J)*9.D0/5.D0*G(1)**2
          END IF
!
!Convert into form readable by RKSTP
!
          F(185+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BGTQ(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(2)*ID(I,J))
          F(194+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BGTL(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(2)*ID(I,J))
          F(138+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BGTPQ(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
          F(147+(I-1)*3+J)=THSL(I)*(1.D0/16.D0/PI**2*BGTPL(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
          F(156+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTPU(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
     $       *THSU(J)
          F(165+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTPD(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
     $       *THSD(J)
          F(174+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTPE(I,J)
     $       +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0)*ID(I,J))
     $       *THSE(J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          DUMTSQ1(I,J)=2.D0*SQGTSQDGTSQ(I,J)+SUGTSUGTSUD(I,J)
     $                 +SDGTSDGTSDD(I,J)
          IF(THHH.EQ.0)THEN
            DUMTSQ2(I,J)=LUSLUT(I,J)+LDSLDT(I,J)
          ELSE
            DUMTSQ2(I,J)=FUSFUT(I,J)+FDSFDT(I,J)
          END IF
          DUMTSQ3(I,J)=DUMTQ3(I,J)
          DUMTSQ4(I,J)=DUMTQ4(I,J)
          DUMTSQ5(I,J)=DUMTQ5(I,J)
!
          DUMTSU1(I,J)=DUMTPU1(I,J)
          DUMTSU3(I,J)=SUGTSUGTSUD(I,J)+SDGTSDGTSDD(I,J)
          DUMTSU4(I,J)=DUMTPU3(I,J)
!
          DUMTSD1(I,J)=DUMTPD1(I,J)
          DUMTSD4(I,J)=DUMTPD3(I,J)
        END DO
      END DO
      TDUMTSQ1=DTRACE(DUMTSQ1)
      TDUMTSU2=DTRACE(SQGTSQDGTSQ)
      TDUMTSU3=DTRACE(DUMTSU3)
      TDUMTSD2=TDUMTSU2
      TDUMTSD3=TDUMTSU3
      DO I=1,3
        DO J=1,3
          IF(THHH.EQ.0)THEN
            BGTSQ(I,J)=1.D0/2.D0*THGL*TDUMTSQ1*GTSQ(I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTSQ,DUMTSQ2,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTSQ,DUMTSQ3,I,J)
     $                 -2.D0*THSH*(DMATMUL(0,FTUQS,SUGTSUFTUUT,I,J)
     $                  +DMATMUL(0,FTDQS,SDGTSDFTDDT,I,J))
     $                 +1.D0/2.D0*SFMUL(THSQ,DUMTSQ4,GTSQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTSQ5,GTSQ,I,J)
     $                 -GTSQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+13.D0*G(3)**2)
!
            BGTSU(I,J)=LUTLUSGTSU(I,J)
     $                 +DMATMUL(0,DUMTSU1,GTSU,I,J)
     $                 +THGL*GTSU(I,J)*TDUMTSU2
     $                 +1.D0/2.D0*THGL*GTSU(I,J)*TDUMTSU3
     $                 -4.D0*THSH*DMATMUL(0,SQFTUQTGTSQ,FTUUS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SUGTSUFTUUT,FTUUS,I,J)
     $                 +SFMUL(THSU,GTSU,DUMTSU4,I,J)
     $                 -GTSU(I,J)*(4.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)
!
            BGTSD(I,J)=LDTLDSGTSD(I,J)
     $                 +DMATMUL(0,DUMTSD1,GTSD,I,J)
     $                 +THGL*GTSD(I,J)*TDUMTSD2
     $                 +1.D0/2.D0*THGL*GTSD(I,J)*TDUMTSD3
     $                 -4.D0*THSH*DMATMUL(0,SQFTDQTGTSQ,FTDDS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SDGTSDFTDDT,FTDDS,I,J)
     $                 +SFMUL(THSD,GTSD,DUMTSD4,I,J)
     $                 -GTSD(I,J)*(1.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)
          ELSE
            BGTSQ(I,J)=1.D0/2.D0*THGL*TDUMTSQ1*GTSQ(I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTSQ,DUMTSQ2,I,J)
     $                 +1.D0/2.D0*DMATMUL(0,GTSQ,DUMTSQ3,I,J)
     $                 -2.D0*THSH*(DMATMUL(0,FTUQS,SUGTSUFTUUT,I,J)
     $                  +DMATMUL(0,FTDQS,SDGTSDFTDDT,I,J))
     $                 +1.D0/2.D0*SFMUL(THSQ,DUMTSQ4,GTSQ,I,J)
     $                 +THSH*SFMUL(THSQ,DUMTSQ5,GTSQ,I,J)
     $                 -GTSQ(I,J)*(1.D0/20.D0*G(1)**2
     $                  +9.D0/4.D0*G(2)**2+13.D0*G(3)**2)
!
            BGTSU(I,J)=FUTFUSGTSU(I,J)
     $                 +DMATMUL(0,DUMTSU1,GTSU,I,J)
     $                 +THGL*GTSU(I,J)*TDUMTSU2
     $                 +1.D0/2.D0*THGL*GTSU(I,J)*TDUMTSU3
     $                 -4.D0*THSH*DMATMUL(0,SQFTUQTGTSQ,FTUUS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SUGTSUFTUUT,FTUUS,I,J)
     $                 +SFMUL(THSU,GTSU,DUMTSU4,I,J)
     $                 -GTSU(I,J)*(4.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)
!
            BGTSD(I,J)=FDTFDSGTSD(I,J)
     $                 +DMATMUL(0,DUMTSD1,GTSD,I,J)
     $                 +THGL*GTSD(I,J)*TDUMTSD2
     $                 +1.D0/2.D0*THGL*GTSD(I,J)*TDUMTSD3
     $                 -4.D0*THSH*DMATMUL(0,SQFTDQTGTSQ,FTDDS,I,J)
     $                 +2.D0*THSH*DMATMUL(0,SDGTSDFTDDT,FTDDS,I,J)
     $                 +SFMUL(THSD,GTSD,DUMTSD4,I,J)
     $                 -GTSD(I,J)*(1.D0/5.D0*G(1)**2
     $                  +13.D0*G(3)**2)

          END IF
!
!Convert into form readable by RKSTP
!
          F(205+(I-1)*3+J)=THSQ(I)*(1.D0/16.D0/PI**2*BGTSQ(I,J)
     $           +1.D0/(16.D0*PI**2)**2*B2GM(3)*ID(I,J))
          F(214+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTSU(I,J)
     $           +1.D0/(16.D0*PI**2)**2*B2GM(3)*ID(I,J))*THSU(J)
          F(223+(I-1)*3+J)=(1.D0/16.D0/PI**2*BGTSD(I,J)
     $           +1.D0/(16.D0*PI**2)**2*B2GM(3)*ID(I,J))*THSD(J)
        END DO
      END DO
!
!The Higgs tilde gauge couplings are calculated twice,
!once for below m_H when s and c are absorbed into the definition.
!Some terms are valid in both regimes and do not need to be
!calculated twice.
!
      DO I=1,3
        DO J=1,3
          DUMGTHU1(I,J)=DMATMUL(0,SQFTUQTGTQ,FUS,I,J)
          DUMGTHD1(I,J)=3.D0*DMATMUL(0,SQFTDQTGTQ,FDS,I,J)
     $                  +DMATMUL(0,SLFTELTGTL,FES,I,J)
          DUMGTPHU1(I,J)=DMATMUL(0,SQFTUQTGTPQ,FUS,I,J)
          DUMGTPHU2(I,J)=DMATMUL(0,SUGTPUFTUUT,FUS,I,J)
          DUMGTPHD1(I,J)=-DMATMUL(0,SQFTDQTGTPQ,FDS,I,J)
     $                   +DMATMUL(0,SLFTELTGTPL,FES,I,J)
          DUMGTPHD2(I,J)=DMATMUL(0,SDGTPDFTDDT,FDS,I,J)
     $                   +DMATMUL(0,SEGTPEFTEET,FES,I,J)
        END DO
      END DO
      TSQGTQDGTQ=DTRACE(SQGTQDGTQ)
      TSLGTLDGTL=DTRACE(SLGTLDGTL)
      TSUFTUUFTUUD=DTRACE(SUFTUUFTUUD)
      TSQFTUQDFTUQ=DTRACE(SQFTUQDFTUQ)
      TDUMGTHU1=DTRACE(DUMGTHU1)
      TFUSFUT=DTRACE(FUSFUT)
      TFDSFDT=DTRACE(FDSFDT)
      TFESFET=DTRACE(FESFET)
      TSDFTDDFTDDD=DTRACE(SDFTDDFTDDD)
      TSEFTEEFTEED=DTRACE(SEFTEEFTEED)
      TSQFTDQDFTDQ=DTRACE(SQFTDQDFTDQ)
      TSLFTELDFTEL=DTRACE(SLFTELDFTEL)
      TDUMGTHD1=DTRACE(DUMGTHD1)
      TSQGTPQDGTPQ=DTRACE(SQGTPQDGTPQ)
      TSLGTPLDGTPL=DTRACE(SLGTPLDGTPL)
      TSUGTPUGTPUD=DTRACE(SUGTPUGTPUD)
      TSDGTPDGTPDD=DTRACE(SDGTPDGTPDD)
      TSEGTPEGTPED=DTRACE(SEGTPEGTPED)
      TDUMGTPHU1=DTRACE(DUMGTPHU1)
      TDUMGTPHU2=DTRACE(DUMGTPHU2)
      TDUMGTPHD1=DTRACE(DUMGTPHD1)
      TDUMGTPHD2=DTRACE(DUMGTPHD2)
!
      BGTHU=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(204)
     $      +1.D0/2.D0*THSH*THSW*(MGTHUSQ+MGTHDSQ)*G(204)
     $      +1.D0/2.D0*THSH*(3.D0*G(204)*TSUFTUUFTUUD
     $       +3.D0*G(204)*TSQFTUQDFTUQ)
     $      +1.D0/4.D0*THSH*G(204)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $      -6.D0*TDUMGTHU1
     $      +G(204)*3.D0*TFUSFUT
     $      +1.D0/2.D0*THSH*G(204)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $      -G(204)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
      BGTHD=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(205)
     $      +1.D0/2.D0*THSH*THSW*(MGTHUSQ+MGTHDSQ)*G(205)
     $      +1.D0/2.D0*THSH*(3.D0*G(205)*TSDFTDDFTDDD
     $       +G(205)*TSEFTEEFTEED+3.D0*G(205)*TSQFTDQDFTDQ
     $       +G(205)*TSLFTELDFTEL)
     $      +1.D0/4.D0*THSH*G(205)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $      -2.D0*TDUMGTHD1
     $      +G(205)*(3.D0*TFDSFDT+TFESFET)
     $      +1.D0/2.D0*THSH*G(205)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $      -G(205)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
      BGTPHU=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $        +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $        +2.D0*TSEGTPEGTPED)*G(184)
     $       +1.D0/2.D0*THSH*THSB*(MGTPHUSQ+MGTPHDSQ)*G(184)
     $       +1.D0/2.D0*THSH*(3.D0*G(184)*TSUFTUUFTUUD
     $        +3.D0*G(184)*TSQFTUQDFTUQ)
     $       +1.D0/4.D0*THSH*G(184)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $       +2.D0*TDUMGTPHU1-8.D0*TDUMGTPHU2
     $       +G(184)*3.D0*TFUSFUT
     $       +1.D0/2.D0*THSH*G(184)*(3.D0*THSW*MGTHUSQ+THSB*MGTPHUSQ)
     $       -G(184)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
      BGTPHD=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $        +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $        +2.D0*TSEGTPEGTPED)*G(185)
     $       +1.D0/2.D0*THSH*THSB*(MGTPHUSQ+MGTPHDSQ)*G(185)
     $       +1.D0/2.D0*THSH*(3.D0*G(185)*TSDFTDDFTDDD
     $        +G(185)*TSEFTEEFTEED+3.D0*G(185)*TSQFTDQDFTDQ
     $        +G(185)*TSLFTELDFTEL)
     $       +1.D0/4.D0*THSH*G(185)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       +2.D0*TDUMGTPHD1-4.D0*TDUMGTPHD2
     $       +G(185)*(3.D0*TFDSFDT+TFESFET)
     $       +1.D0/2.D0*THSH*G(185)*(3.D0*THSW*MGTHDSQ+THSB*MGTPHDSQ)
     $       -G(185)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
!Convert into form readable by RKSTP
!
      F(184)=THHH*(1.D0/16.D0/PI**2*BGTPHU
     $           +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0))
      F(185)=THHH*(1.D0/16.D0/PI**2*BGTPHD
     $           +1.D0/(16.D0*PI**2)**2*B2GM(1)*DSQRT(3.D0/5.D0))
      F(204)=THHH*(1.D0/16.D0/PI**2*BGTHU
     $           +1.D0/(16.D0*PI**2)**2*B2GM(2))
      F(205)=THHH*(1.D0/16.D0/PI**2*BGTHD
     $           +1.D0/(16.D0*PI**2)**2*B2GM(2))
!
      IF(THHH.EQ.0)THEN
        DO I=1,3
          DO J=1,3
            DUMGTHU1(I,J)=DMATMUL(0,SQFTUQTGTQ,LUS,I,J)
            DUMGTHD1(I,J)=3.D0*DMATMUL(0,SQFTDQTGTQ,LDS,I,J)
     $                    +DMATMUL(0,SLFTELTGTL,LES,I,J)
            DUMGTPHU1(I,J)=DMATMUL(0,SQFTUQTGTPQ,LUS,I,J)
            DUMGTPHU2(I,J)=DMATMUL(0,SUGTPUFTUUT,LUS,I,J)
            DUMGTPHD1(I,J)=-DMATMUL(0,SQFTDQTGTPQ,LDS,I,J)
     $                     +DMATMUL(0,SLFTELTGTPL,LES,I,J)
            DUMGTPHD2(I,J)=DMATMUL(0,SDGTPDFTDDT,LDS,I,J)
     $                     +DMATMUL(0,SEGTPEFTEET,LES,I,J)
          END DO
        END DO
        TDUMGTHU1=DTRACE(DUMGTHU1)
        TLUSLUT=DTRACE(LUSLUT)
        TLDSLDT=DTRACE(LDSLDT)
        TLESLET=DTRACE(LESLET)
        TDUMGTHD1=DTRACE(DUMGTHD1)
        TDUMGTPHU1=DTRACE(DUMGTPHU1)
        TDUMGTPHU2=DTRACE(DUMGTPHU2)
        TDUMGTPHD1=DTRACE(DUMGTPHD1)
        TDUMGTPHD2=DTRACE(DUMGTPHD2)
!
        BSGTHU=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(289)
     $        +1.D0/2.D0*THSH*THSW*(MSGTHUSQ+MCGTHDSQ)*G(289)
     $        +1.D0/2.D0*THSH*(3.D0*G(289)*TSUFTUUFTUUD
     $         +3.D0*G(289)*TSQFTUQDFTUQ)
     $        +1.D0/4.D0*THSH*G(289)
     $         *(3.D0*THSW*MSGTHUSQ+THSB*MSGTPHUSQ)
     $        -6.D0*TDUMGTHU1
****
     $        +THSB*THSH*G(290)*G(288)*G(287)
     $        -THSW*THSH*MCGTHDSQ*G(289)
****
     $        +G(289)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $        +1.D0/2.D0*THSH*G(289)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $         +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $        -G(289)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
        BCGTHD=1.D0/2.D0*THSW*(3.D0*TSQGTQDGTQ+TSLGTLDGTL)*G(290)
     $        +1.D0/2.D0*THSH*THSW*(MSGTHUSQ+MCGTHDSQ)*G(290)
     $        +1.D0/2.D0*THSH*(3.D0*G(290)*TSDFTDDFTDDD
     $         +G(290)*TSEFTEEFTEED+3.D0*G(290)*TSQFTDQDFTDQ
     $         +G(290)*TSLFTELDFTEL)
     $        +1.D0/4.D0*THSH*G(290)*(3.D0*THSW*MCGTHDSQ+THSB*MCGTPHDSQ)
     $        -2.D0*TDUMGTHD1
****
     $        +THSB*THSH*G(289)*G(287)*G(288)
     $        -THSW*THSH*MSGTHUSQ*G(290)
****
     $        +G(290)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $        +1.D0/2.D0*THSH*G(290)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $         +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $        -G(290)*(9.D0/20.D0*G(1)**2+33.D0/4.D0*G(2)**2)
!
        BSGTPHU=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $          +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $          +2.D0*TSEGTPEGTPED)*G(287)
     $         +1.D0/2.D0*THSH*THSB*(MSGTPHUSQ+MCGTPHDSQ)*G(287)
     $         +1.D0/2.D0*THSH*(3.D0*G(287)*TSUFTUUFTUUD
     $          +3.D0*G(287)*TSQFTUQDFTUQ)
     $         +1.D0/4.D0*THSH*G(287)*(3.D0*THSW*MSGTHUSQ
     $          +THSB*MSGTPHUSQ)
     $         +2.D0*TDUMGTPHU1-8.D0*TDUMGTPHU2
****
     $         +THSB*THSH*MCGTPHDSQ*G(287)
     $         +3.D0*THSW*THSH*G(288)*G(290)*G(289)
****
     $         +G(287)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $         +1.D0/2.D0*THSH*G(287)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $          +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $         -G(287)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
        BCGTPHD=2.D0/4.D0*THSB*(1.D0/3.D0*TSQGTPQDGTPQ+TSLGTPLDGTPL
     $          +8.D0/3.D0*TSUGTPUGTPUD+2.D0/3.D0*TSDGTPDGTPDD
     $          +2.D0*TSEGTPEGTPED)*G(288)
     $         +1.D0/2.D0*THSH*THSB*(MSGTPHUSQ+MCGTPHDSQ)*G(288)
     $         +1.D0/2.D0*THSH*(3.D0*G(288)*TSDFTDDFTDDD
     $          +G(288)*TSEFTEEFTEED+3.D0*G(288)*TSQFTDQDFTDQ
     $          +G(288)*TSLFTELDFTEL)
     $         +1.D0/4.D0*THSH*G(288)*(3.D0*THSW*MCGTHDSQ
     $          +THSB*MCGTPHDSQ)
     $         +2.D0*TDUMGTPHD1-4.D0*TDUMGTPHD2
****
     $         +THSB*THSH*MSGTPHUSQ*G(288)
     $         +3.D0*THSW*THSH*G(287)*G(289)*G(290)
****
     $         +G(288)*(3.D0*TLUSLUT+3.D0*TLDSLDT+TLESLET)
     $         +1.D0/2.D0*THSH*G(288)*(3.D0*THSW*(MSGTHUSQ+MCGTHDSQ)
     $          +THSB*(MSGTPHUSQ+MCGTPHDSQ))
     $         -G(288)*(9.D0/20.D0*G(1)**2+9.D0/4.D0*G(2)**2)
!
!The two loop part of the running will be the same as the
!ordinary g running multiplied by an s (or corresponding c)
!
!
!Convert into form readable by RKSTP
!
        F(287)=1.D0/16.D0/PI**2*BSGTPHU
     $             +1.D0/(16.D0*PI**2)**2*SINB*B2GM(1)*DSQRT(3.D0/5.D0)
        F(288)=1.D0/16.D0/PI**2*BCGTPHD
     $             +1.D0/(16.D0*PI**2)**2*COSB*B2GM(1)*DSQRT(3.D0/5.D0)
        F(289)=1.D0/16.D0/PI**2*BSGTHU
     $             +1.D0/(16.D0*PI**2)**2*SINB*B2GM(2)
        F(290)=1.D0/16.D0/PI**2*BCGTHD
     $             +1.D0/(16.D0*PI**2)**2*COSB*B2GM(2)
      END IF
!
!Now I work out the \tilde{\mu}^* f^{h_u}_u terms. These only exist
!above m_H
!
      IF(THHH.EQ.1)THEN
        DO I=1,3
          DO J=1,3
            SQIDMTSFU(I,J)=SFMUL(THSQ,ID,MTSFU,I,J)
            SQIDMTSFD(I,J)=SFMUL(THSQ,ID,MTSFD,I,J)
            SQFUHUDMTSFU(I,J)=SFMUL(THSQ,FUHUD,MTSFU,I,J)
            SQFDHDDMTSFD(I,J)=SFMUL(THSQ,FDHDD,MTSFD,I,J)
            SLIDMTSFE(I,J)=SFMUL(THSL,ID,MTSFE,I,J)
            SLFEHDDMTSFE(I,J)=SFMUL(THSL,FEHDD,MTSFE,I,J)
            SUMTSFUID(I,J)=SFMUL(THSU,MTSFU,ID,I,J)
            SUMTSFUFUQDFDQ(I,J)=SFMUL(THSU,MTSFU,FUQDFDQ,I,J)
            SDMTSFDID(I,J)=SFMUL(THSD,MTSFD,ID,I,J)
            SDMTSFDFDQDFUQ(I,J)=SFMUL(THSD,MTSFD,FDQDFUQ,I,J)
            SEMTSFEID(I,J)=SFMUL(THSE,MTSFE,ID,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            SUSQIDMTSFUID(I,J)=SFMUL(THSU,SQIDMTSFU,ID,I,J)
            SDSQIDMTSFDID(I,J)=SFMUL(THSD,SQIDMTSFD,ID,I,J)
            SESLIDMTSFEID(I,J)=SFMUL(THSE,SLIDMTSFE,ID,I,J)
          END DO
        END DO
        SUTSQFUHUDMTSFU=TSFMUL(THSU,SQFUHUDMTSFU)
        SDTSQFDHDDMTSFD=TSFMUL(THSD,SQFDHDDMTSFD)
        SETSLFEHDDMTSFE=TSFMUL(THSE,SLFEHDDMTSFE)
        DO I=1,3
          DO J=1,3
            DUMMTSFU1(I,J)=(G(457)**2/6.D0-3.D0*G(458)**2/2.D0)*ID(I,J)
     $                    +2.D0*FDDRFDDRD(I,J)
            DUMMTSFU2(I,J)=8.D0/9.D0*THSB*GTPUTGTPUS(I,J)
     $                    +8.D0/3.D0*THGL*GTSUTGTSUS(I,J)
     $                    +2.D0*THSH*FTUUDFTUU(I,J)
            DUMMTSFU3(I,J)=THSH*FTUQFTUQD(I,J)+THSH*FTDQFTDQD(I,J)
     $                    +1.D0/18.D0*THSB*GTPQSGTPQT(I,J)
     $                    +3.D0/2.D0*THSW*GTQSGTQT(I,J)
     $                    +8.D0/3.D0*THGL*GTSQSGTSQT(I,J)
!
            DUMMTSFD1(I,J)=(G(457)**2/6.D0+3.D0*G(458)**2/2.D0)*ID(I,J)
     $                    -2.D0*FUURFUURD(I,J)
            DUMMTSFD2(I,J)=2.D0/9.D0*THSB*GTPDTGTPDS(I,J)
     $                    +8.D0/3.D0*THGL*GTSDTGTSDS(I,J)
     $                    +2.D0*THSH*FTDDDFTDD(I,J)
            DUMMTSFD3(I,J)=DUMMTSFU3(I,J)
!
            DUMMTSFE1(I,J)=(G(457)**2/2.D0-3.D0*G(458)**2/2.D0)*ID(I,J)
            DUMMTSFE2(I,J)=2.D0*THSB*GTPETGTPES(I,J)
     $                    +2.D0*THSH*FTEEDFTEE(I,J)
            DUMMTSFE3(I,J)=THSH*FTELFTELD(I,J)
     $                    +1.D0/2.D0*THSB*GTPLSGTPLT(I,J)
     $                    +3.D0/2.D0*THSW*GTLSGTLT(I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            BMTSFU(I,J)=2.D0/3.D0*G(457)**2*SUMTSFUID(I,J)
     $                 -2.D0*(G(457)**2/9.D0+4.D0*G(459)**2/3.D0)
     $                  *SUSQIDMTSFUID(I,J)
     $                 +6.D0*FUHU(I,J)*SUTSQFUHUDMTSFU
     $                 -SFMUL(THSQ,DUMMTSFU1,MTSFU,I,J)
     $                 -2.D0*SDMTSFDFDQDFUQ(I,J)
     $                 -2.D0/3.D0*THSB*THSH*G(108)*G(185)
     $                  *(4.D0*FTUQGTPUS(I,J)-GTPQSFTUU(I,J))
     $                 -6.D0*THSH*THSW*G(108)*G(205)
     $                  *GTQSFTUU(I,J)
     $                 +4.D0*THSH*G(108)*FTDQFDDFTUU(I,J)
     $                 +SFMUL(THSU,MTSFU,DUMMTSFU2,I,J)
     $                 +(3.D0*TFDDFD+TFEDFE+1.D0/2.D0*THSB*THSH*MGTPHDSQ
     $                  +3.D0/2.D0*THSW*THSH*MGTHDSQ)*MTSFU(I,J)
     $                 +SFMUL(THSQ,DUMMTSFU3,MTSFU,I,J)
     $                 -(13.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*MTSFU(I,J)
!
            BMTSFD(I,J)=1.D0/3.D0*G(457)**2*SDMTSFDID(I,J)
     $                 +2.D0*(G(457)**2/18.D0-4.D0*G(459)**2/3.D0)
     $                  *SDSQIDMTSFDID(I,J)
     $                 +6.D0*FDHD(I,J)*SDTSQFDHDDMTSFD
     $                 +2.D0*FDHD(I,J)*SETSLFEHDDMTSFE
     $                 +SFMUL(THSQ,DUMMTSFD1,MTSFD,I,J)
     $                 -2.D0*SUMTSFUFUQDFDQ(I,J)
     $                 -2.D0/3.D0*THSB*THSH*G(108)*G(184)
     $                  *(2.D0*FTDQGTPDS(I,J)+GTPQSFTDD(I,J))
     $                 -6.D0*THSH*THSW*G(108)*G(204)
     $                  *GTQSFTDD(I,J)
     $                 +4.D0*THSH*G(108)*FTUQFUDFTDD(I,J)
     $                 +SFMUL(THSD,MTSFD,DUMMTSFD2,I,J)
     $                 +(3.D0*TFUDFU+1.D0/2.D0*THSB*THSH*MGTPHUSQ
     $                  +3.D0/2.D0*THSW*THSH*MGTHUSQ)*MTSFD(I,J)
     $                 +SFMUL(THSQ,DUMMTSFD3,MTSFD,I,J)
     $                 -(7.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2
     $                  +8.D0*G(3)**2)*MTSFD(I,J)
!
            BMTSFE(I,J)=G(457)**2*SEMTSFEID(I,J)
     $                 -G(457)**2*SESLIDMTSFEID(I,J)
     $                 +2.D0*FEHD(I,J)*SETSLFEHDDMTSFE
     $                 +6.D0*FEHD(I,J)*SDTSQFDHDDMTSFD
     $                 -SFMUL(THSL,DUMMTSFE1,MTSFE,I,J)
     $                 -2.D0*THSB*THSH*G(108)*G(184)
     $                  *(2.D0*FTELGTPES(I,J)-GTPLSFTEE(I,J))
     $                 -6.D0*THSH*THSW*G(108)*G(204)
     $                  *GTLSFTEE(I,J)
     $                 +SFMUL(THSE,MTSFE,DUMMTSFE2,I,J)
     $                 +(3.D0*TFUDFU+1.D0/2.D0*THSB*THSH*MGTPHUSQ
     $                  +3.D0/2.D0*THSW*THSH*MGTHUSQ)*MTSFE(I,J)
     $                 +SFMUL(THSL,DUMMTSFE3,MTSFE,I,J)
     $                 -(27.D0/10.D0*G(1)**2+9.D0/2.D0*G(2)**2)
     $                  *MTSFE(I,J)
!
!Convert into form readable by RKSTP
!
            F(429+(I-1)*3+J)=1.D0/16.D0/PI**2*BMTSFU(I,J)
     $           +1.D0/(16.D0*PI**2)**2*(G(398)*B2YMU(J,I)
     $                                    +YU(J,I)*B2GRKMUM)
            F(438+(I-1)*3+J)=1.D0/16.D0/PI**2*BMTSFD(I,J)
     $           +1.D0/(16.D0*PI**2)**2*(G(398)*B2YMD(J,I)
     $                                    +YD(J,I)*B2GRKMUM)
            F(447+(I-1)*3+J)=1.D0/16.D0/PI**2*BMTSFE(I,J)
     $           +1.D0/(16.D0*PI**2)**2*(G(398)*B2YME(J,I)
     $                                    +YE(J,I)*B2GRKMUM)

          END DO
        END DO
      END IF
!
!
!Next I am going to work out the primed gaugino terms, M'_{1,2,3}
!
!
      IF(THHH.EQ.0)THEN
        BMP(1)=G(599)*(TDUMM1+THSH*(MSGTPHUSQ+MCGTPHDSQ))
     $        -2.D0*(0.D0,1.D0)*THSH*(G(287)*G(108)*G(288)
     $                            -G(287)*G(108)*G(288))
!
        BMP(2)=G(600)*(TDUMM2+THSH*(MSGTHUSQ+MCGTHDSQ))
     $        -2.D0*(0.D0,1.D0)*THSH*(G(289)*G(108)*G(290)
     $                            -G(289)*G(108)*G(290))
     $       -12.D0*G(600)*G(2)**2
!
        BMP(3)=G(601)*TDUMM3-18.D0*G(601)*G(3)**2
      ELSE
        BMP(1)=G(599)*(TDUMM1+THSH*(MGTPHUSQ+MGTPHDSQ))
!
        BMP(2)=G(600)*(TDUMM2+THSH*(MGTHUSQ+MGTHDSQ))
     $        -12.D0*G(600)*G(2)**2
!
        BMP(3)=G(601)*TDUMM3-18.D0*G(601)*G(3)**2
      END IF
!
!NB: There is some work here converting the MSSM MV notation
!    M_{1,2,3} 2-lp beta functions to BT notation.
!The RKSTP compatible derivatives are:
!
      F(599)=THSB*(1.D0/16.D0/PI**2*BMP(1)
     $         +1.D0/(16.D0*PI**2)**2*(0.D0,1.D0)/2.D0
     $                                *(B2M(1)-B2M(1)))
      F(600)=THSW*(1.D0/16.D0/PI**2*BMP(2)
     $         +1.D0/(16.D0*PI**2)**2*(0.D0,1.D0)/2.D0
     $                                *(B2M(2)-B2M(2)))
      F(601)=THGL*(1.D0/16.D0/PI**2*BMP(3)
     $         +1.D0/(16.D0*PI**2)**2*(0.D0,1.D0)/2.D0
     $                                *(B2M(3)-B2M(3)))
!
!Rotate F to original basis
!
      DO I=1,601
        FCURR(I)=F(I)
      END DO
!
      IF(NSQ+NSU+NSD+NSL+NSE.LT.15)THEN
!
!If necessary, rotate the Fs back into the old current basis.
!
        CALL DROTSQBACK(F,FCURR)
      END IF
!
!Set the OLDN for the next step
!
      OLDNSQ=NSQ
      OLDNSU=NSU
      OLDNSD=NSD
      OLDNSL=NSL
      OLDNSE=NSE
!
      RETURN
      END
!
!
      SUBROUTINE DROTSQBACK(GROT,GUNROT)
!
!Purpose: To rotate the couplings back from the squark mass
!         basis to the old current basis
!
      IMPLICIT NONE
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE PRECISION GUNROT(601),GROT(601),DMATMUL
!
      DOUBLE PRECISION FU(3,3),FD(3,3),FE(3,3),AU(3,3),AD(3,3),AE(3,3)
      DOUBLE PRECISION MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE PRECISION FUM(3,3),FDM(3,3),FEM(3,3)
      DOUBLE PRECISION AUM(3,3),ADM(3,3),AEM(3,3)
      DOUBLE PRECISION MQM(3,3),MLM(3,3),MUPM(3,3),MDM(3,3),MEM(3,3)
      DOUBLE PRECISION LU(3,3),LD(3,3),LE(3,3)
      DOUBLE PRECISION GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3)
      DOUBLE PRECISION GTPE(3,3),GTQ(3,3),GTL(3,3),GTSQ(3,3)
      DOUBLE PRECISION GTSU(3,3),GTSD(3,3),FTUQ(3,3),FTDQ(3,3)
      DOUBLE PRECISION FTEL(3,3),FTUU(3,3),FTDD(3,3),FTEE(3,3)
      DOUBLE PRECISION TRIU(3,3),TRID(3,3),TRIE(3,3)
      DOUBLE PRECISION MTSFU(3,3),MTSFD(3,3),MTSFE(3,3)
!
      DOUBLE PRECISION FUT(3,3),FDT(3,3),FET(3,3)
      DOUBLE PRECISION AUT(3,3),ADT(3,3),AET(3,3)
      DOUBLE PRECISION MQT(3,3),MLT(3,3),MUPT(3,3),MDT(3,3),MET(3,3)
      DOUBLE PRECISION FUMT(3,3),FDMT(3,3),FEMT(3,3)
      DOUBLE PRECISION AUMT(3,3),ADMT(3,3),AEMT(3,3)
      DOUBLE PRECISION MQMT(3,3),MLMT(3,3),MUPMT(3,3)
      DOUBLE PRECISION MDMT(3,3),MEMT(3,3)
      DOUBLE PRECISION LUT(3,3),LDT(3,3),LET(3,3)
      DOUBLE PRECISION GTPQT(3,3),GTPLT(3,3),GTPUT(3,3),GTPDT(3,3)
      DOUBLE PRECISION GTPET(3,3),GTQT(3,3),GTLT(3,3),GTSQT(3,3)
      DOUBLE PRECISION GTSUT(3,3),GTSDT(3,3),FTUQT(3,3),FTDQT(3,3)
      DOUBLE PRECISION FTELT(3,3),FTUUT(3,3),FTDDT(3,3),FTEET(3,3)
      DOUBLE PRECISION TRIUT(3,3),TRIDT(3,3),TRIET(3,3)
      DOUBLE PRECISION MTSFUT(3,3),MTSFDT(3,3),MTSFET(3,3)
!
      DOUBLE PRECISION DRQ(3,3),DRUP(3,3),DRD(3,3),DRL(3,3),DRE(3,3)
      DOUBLE PRECISION RQT(3,3),RUPT(3,3),RDT(3,3),RLT(3,3),RET(3,3)
!
      INTEGER I,J
!
      DO I=1,601
        GUNROT(I)=GROT(I)
      END DO
      DO I=1,3
        DO J=1,3
          DRQ(I,J)=DBLE(RQTOT(I,J))
          DRUP(I,J)=DBLE(RUPTOT(I,J))
          DRD(I,J)=DBLE(RDTOT(I,J))
          DRL(I,J)=DBLE(RLTOT(I,J))
          DRE(I,J)=DBLE(RETOT(I,J))
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=GROT(3+(I-1)*3+J)
          FD(I,J)=GROT(12+(I-1)*3+J)
          FE(I,J)=GROT(21+(I-1)*3+J)
          AU(I,J)=GROT(33+(I-1)*3+J)
          AD(I,J)=GROT(42+(I-1)*3+J)
          AE(I,J)=GROT(51+(I-1)*3+J)
          MQ(I,J)=GROT(62+(I-1)*3+J)
          ML(I,J)=GROT(71+(I-1)*3+J)
          MUP(I,J)=GROT(80+(I-1)*3+J)
          MD(I,J)=GROT(89+(I-1)*3+J)
          ME(I,J)=GROT(98+(I-1)*3+J)
!
          FUM(I,J)=GROT(293+(I-1)*3+J)
          FDM(I,J)=GROT(302+(I-1)*3+J)
          FEM(I,J)=GROT(311+(I-1)*3+J)
          AUM(I,J)=GROT(323+(I-1)*3+J)
          ADM(I,J)=GROT(332+(I-1)*3+J)
          AEM(I,J)=GROT(341+(I-1)*3+J)
          MQM(I,J)=GROT(352+(I-1)*3+J)
          MLM(I,J)=GROT(361+(I-1)*3+J)
          MUPM(I,J)=GROT(370+(I-1)*3+J)
          MDM(I,J)=GROT(379+(I-1)*3+J)
          MEM(I,J)=GROT(388+(I-1)*3+J)
!
          LU(I,J)=GROT(111+(I-1)*3+J)
          LD(I,J)=GROT(120+(I-1)*3+J)
          LE(I,J)=GROT(129+(I-1)*3+J)
!
          GTPQ(I,J)=GROT(138+(I-1)*3+J)
          GTPL(I,J)=GROT(147+(I-1)*3+J)
          GTPU(I,J)=GROT(156+(I-1)*3+J)
          GTPD(I,J)=GROT(165+(I-1)*3+J)
          GTPE(I,J)=GROT(174+(I-1)*3+J)
          GTQ(I,J)=GROT(185+(I-1)*3+J)
          GTL(I,J)=GROT(194+(I-1)*3+J)
          GTSQ(I,J)=GROT(205+(I-1)*3+J)
          GTSU(I,J)=GROT(214+(I-1)*3+J)
          GTSD(I,J)=GROT(223+(I-1)*3+J)
          FTUQ(I,J)=GROT(232+(I-1)*3+J)
          FTDQ(I,J)=GROT(241+(I-1)*3+J)
          FTEL(I,J)=GROT(250+(I-1)*3+J)
          FTUU(I,J)=GROT(259+(I-1)*3+J)
          FTDD(I,J)=GROT(268+(I-1)*3+J)
          FTEE(I,J)=GROT(277+(I-1)*3+J)
!
          TRIU(I,J)=GROT(399+(I-1)*3+J)
          TRID(I,J)=GROT(408+(I-1)*3+J)
          TRIE(I,J)=GROT(417+(I-1)*3+J)
          MTSFU(I,J)=GROT(429+(I-1)*3+J)
          MTSFD(I,J)=GROT(438+(I-1)*3+J)
          MTSFE(I,J)=GROT(447+(I-1)*3+J)
!
          RQT(I,J)=DRQ(J,I)
          RUPT(I,J)=DRUP(J,I)
          RDT(I,J)=DRD(J,I)
          RLT(I,J)=DRL(J,I)
          RET(I,J)=DRE(J,I)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUT(I,J)=DMATMUL(0,FU,RUPT,I,J)
          FDT(I,J)=DMATMUL(0,FD,RDT,I,J)
          FET(I,J)=DMATMUL(0,FE,RET,I,J)
          AUT(I,J)=DMATMUL(0,AU,RUPT,I,J)
          ADT(I,J)=DMATMUL(0,AD,RDT,I,J)
          AET(I,J)=DMATMUL(0,AE,RET,I,J)
          MQT(I,J)=DMATMUL(2,MQ,DRQ,I,J)
          MLT(I,J)=DMATMUL(2,ML,DRL,I,J)
          MUPT(I,J)=DMATMUL(2,MUP,DRUP,I,J)
          MDT(I,J)=DMATMUL(2,MD,DRD,I,J)
          MET(I,J)=DMATMUL(2,ME,DRE,I,J)
!
          FUMT(I,J)=DMATMUL(0,FUM,RUPT,I,J)
          FDMT(I,J)=DMATMUL(0,FDM,RDT,I,J)
          FEMT(I,J)=DMATMUL(0,FEM,RET,I,J)
          AUMT(I,J)=DMATMUL(0,AUM,RUPT,I,J)
          ADMT(I,J)=DMATMUL(0,ADM,RDT,I,J)
          AEMT(I,J)=DMATMUL(0,AEM,RET,I,J)
          MQMT(I,J)=DMATMUL(2,MQM,DRQ,I,J)
          MLMT(I,J)=DMATMUL(2,MLM,DRL,I,J)
          MUPMT(I,J)=DMATMUL(2,MUPM,DRUP,I,J)
          MDMT(I,J)=DMATMUL(2,MDM,DRD,I,J)
          MEMT(I,J)=DMATMUL(2,MEM,DRE,I,J)
!
          LUT(I,J)=DMATMUL(0,LU,RUPT,I,J)
          LDT(I,J)=DMATMUL(0,LD,RDT,I,J)
          LET(I,J)=DMATMUL(0,LE,RET,I,J)
!
          GTPQT(I,J)=DMATMUL(2,GTPQ,DRQ,I,J)
          GTPLT(I,J)=DMATMUL(2,GTPL,DRL,I,J)
          GTPUT(I,J)=DMATMUL(2,GTPU,DRUP,I,J)
          GTPDT(I,J)=DMATMUL(2,GTPD,DRD,I,J)
          GTPET(I,J)=DMATMUL(2,GTPE,DRE,I,J)
          GTQT(I,J)=DMATMUL(2,GTQ,DRQ,I,J)
          GTLT(I,J)=DMATMUL(2,GTL,DRL,I,J)
          GTSQT(I,J)=DMATMUL(2,GTSQ,DRQ,I,J)
          GTSUT(I,J)=DMATMUL(2,GTSU,DRUP,I,J)
          GTSDT(I,J)=DMATMUL(2,GTSD,DRD,I,J)
          FTUQT(I,J)=DMATMUL(0,FTUQ,RUPT,I,J)
          FTDQT(I,J)=DMATMUL(0,FTDQ,RDT,I,J)
          FTELT(I,J)=DMATMUL(0,FTEL,RET,I,J)
          FTUUT(I,J)=DMATMUL(0,FTUU,RUPT,I,J)
          FTDDT(I,J)=DMATMUL(0,FTDD,RDT,I,J)
          FTEET(I,J)=DMATMUL(0,FTEE,RET,I,J)
!
          TRIUT(I,J)=DMATMUL(0,TRIU,RUPT,I,J)
          TRIDT(I,J)=DMATMUL(0,TRID,RDT,I,J)
          TRIET(I,J)=DMATMUL(0,TRIE,RET,I,J)
          MTSFUT(I,J)=DMATMUL(0,MTSFU,RUPT,I,J)
          MTSFDT(I,J)=DMATMUL(0,MTSFD,RDT,I,J)
          MTSFET(I,J)=DMATMUL(0,MTSFE,RET,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GUNROT(3+(I-1)*3+J)=DMATMUL(1,RQT,FUT,I,J)
          GUNROT(12+(I-1)*3+J)=DMATMUL(1,RQT,FDT,I,J)
          GUNROT(21+(I-1)*3+J)=DMATMUL(1,RLT,FET,I,J)
          GUNROT(33+(I-1)*3+J)=DMATMUL(1,RQT,AUT,I,J)
          GUNROT(42+(I-1)*3+J)=DMATMUL(1,RQT,ADT,I,J)
          GUNROT(51+(I-1)*3+J)=DMATMUL(1,RLT,AET,I,J)
          GUNROT(62+(I-1)*3+J)=DMATMUL(0,DRQ,MQT,I,J)
          GUNROT(71+(I-1)*3+J)=DMATMUL(0,DRL,MLT,I,J)
          GUNROT(80+(I-1)*3+J)=DMATMUL(0,DRUP,MUPT,I,J)
          GUNROT(89+(I-1)*3+J)=DMATMUL(0,DRD,MDT,I,J)
          GUNROT(98+(I-1)*3+J)=DMATMUL(0,DRE,MET,I,J)
!
          GUNROT(293+(I-1)*3+J)=DMATMUL(1,RQT,FUMT,I,J)
          GUNROT(302+(I-1)*3+J)=DMATMUL(1,RQT,FDMT,I,J)
          GUNROT(311+(I-1)*3+J)=DMATMUL(1,RLT,FEMT,I,J)
          GUNROT(323+(I-1)*3+J)=DMATMUL(1,RQT,AUMT,I,J)
          GUNROT(332+(I-1)*3+J)=DMATMUL(1,RQT,ADMT,I,J)
          GUNROT(341+(I-1)*3+J)=DMATMUL(1,RLT,AEMT,I,J)
          GUNROT(352+(I-1)*3+J)=DMATMUL(0,DRQ,MQMT,I,J)
          GUNROT(361+(I-1)*3+J)=DMATMUL(0,DRL,MLMT,I,J)
          GUNROT(370+(I-1)*3+J)=DMATMUL(0,DRUP,MUPMT,I,J)
          GUNROT(379+(I-1)*3+J)=DMATMUL(0,DRD,MDMT,I,J)
          GUNROT(388+(I-1)*3+J)=DMATMUL(0,DRE,MEMT,I,J)
!
          GUNROT(111+(I-1)*3+J)=DMATMUL(1,RQT,LUT,I,J)
          GUNROT(120+(I-1)*3+J)=DMATMUL(1,RQT,LDT,I,J)
          GUNROT(129+(I-1)*3+J)=DMATMUL(1,RLT,LET,I,J)
!
          GUNROT(138+(I-1)*3+J)=DMATMUL(0,DRQ,GTPQT,I,J)
          GUNROT(147+(I-1)*3+J)=DMATMUL(0,DRL,GTPLT,I,J)
          GUNROT(156+(I-1)*3+J)=DMATMUL(0,DRUP,GTPUT,I,J)
          GUNROT(165+(I-1)*3+J)=DMATMUL(0,DRD,GTPDT,I,J)
          GUNROT(174+(I-1)*3+J)=DMATMUL(0,DRE,GTPET,I,J)
          GUNROT(185+(I-1)*3+J)=DMATMUL(0,DRQ,GTQT,I,J)
          GUNROT(194+(I-1)*3+J)=DMATMUL(0,DRL,GTLT,I,J)
          GUNROT(205+(I-1)*3+J)=DMATMUL(0,DRQ,GTSQT,I,J)
          GUNROT(214+(I-1)*3+J)=DMATMUL(0,DRUP,GTSUT,I,J)
          GUNROT(223+(I-1)*3+J)=DMATMUL(0,DRD,GTSDT,I,J)
          GUNROT(232+(I-1)*3+J)=DMATMUL(1,RQT,FTUQT,I,J)
          GUNROT(241+(I-1)*3+J)=DMATMUL(1,RQT,FTDQT,I,J)
          GUNROT(250+(I-1)*3+J)=DMATMUL(1,RLT,FTELT,I,J)
          GUNROT(259+(I-1)*3+J)=DMATMUL(1,RQT,FTUUT,I,J)
          GUNROT(268+(I-1)*3+J)=DMATMUL(1,RQT,FTDDT,I,J)
          GUNROT(277+(I-1)*3+J)=DMATMUL(1,RLT,FTEET,I,J)
!
          GUNROT(399+(I-1)*3+J)=DMATMUL(1,RQT,TRIUT,I,J)
          GUNROT(408+(I-1)*3+J)=DMATMUL(1,RQT,TRIDT,I,J)
          GUNROT(417+(I-1)*3+J)=DMATMUL(1,RLT,TRIET,I,J)
          GUNROT(429+(I-1)*3+J)=DMATMUL(1,RQT,MTSFUT,I,J)
          GUNROT(438+(I-1)*3+J)=DMATMUL(1,RQT,MTSFDT,I,J)
          GUNROT(447+(I-1)*3+J)=DMATMUL(1,RLT,MTSFET,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE DROTSQ(GUNROT,GROT)
!
!Purpose: To rotate from the current basis to the squark mass basis
!
      IMPLICIT NONE
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE PRECISION GUNROT(601),GROT(601),DMATMUL
!
      DOUBLE PRECISION FU(3,3),FD(3,3),FE(3,3),AU(3,3),AD(3,3),AE(3,3)
      DOUBLE PRECISION MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE PRECISION FUM(3,3),FDM(3,3),FEM(3,3)
      DOUBLE PRECISION AUM(3,3),ADM(3,3),AEM(3,3)
      DOUBLE PRECISION MQM(3,3),MLM(3,3),MUPM(3,3),MDM(3,3),MEM(3,3)
      DOUBLE PRECISION LU(3,3),LD(3,3),LE(3,3)
      DOUBLE PRECISION GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3)
      DOUBLE PRECISION GTPE(3,3),GTQ(3,3),GTL(3,3),GTSQ(3,3)
      DOUBLE PRECISION GTSU(3,3),GTSD(3,3),FTUQ(3,3),FTDQ(3,3)
      DOUBLE PRECISION FTEL(3,3),FTUU(3,3),FTDD(3,3),FTEE(3,3)
      DOUBLE PRECISION TRIU(3,3),TRID(3,3),TRIE(3,3)
      DOUBLE PRECISION MTSFU(3,3),MTSFD(3,3),MTSFE(3,3)
!
      DOUBLE PRECISION FUT(3,3),FDT(3,3),FET(3,3)
      DOUBLE PRECISION AUT(3,3),ADT(3,3),AET(3,3)
      DOUBLE PRECISION MQT(3,3),MLT(3,3),MUPT(3,3),MDT(3,3),MET(3,3)
      DOUBLE PRECISION FUMT(3,3),FDMT(3,3),FEMT(3,3)
      DOUBLE PRECISION AUMT(3,3),ADMT(3,3),AEMT(3,3)
      DOUBLE PRECISION MQMT(3,3),MLMT(3,3),MUPMT(3,3)
      DOUBLE PRECISION MDMT(3,3),MEMT(3,3)
      DOUBLE PRECISION LUT(3,3),LDT(3,3),LET(3,3)
      DOUBLE PRECISION GTPQT(3,3),GTPLT(3,3),GTPUT(3,3),GTPDT(3,3)
      DOUBLE PRECISION GTPET(3,3),GTQT(3,3),GTLT(3,3),GTSQT(3,3)
      DOUBLE PRECISION GTSUT(3,3),GTSDT(3,3),FTUQT(3,3),FTDQT(3,3)
      DOUBLE PRECISION FTELT(3,3),FTUUT(3,3),FTDDT(3,3),FTEET(3,3)
      DOUBLE PRECISION TRIUT(3,3),TRIDT(3,3),TRIET(3,3)
      DOUBLE PRECISION MTSFUT(3,3),MTSFDT(3,3),MTSFET(3,3)
!
      DOUBLE PRECISION DRQ(3,3),DRUP(3,3),DRD(3,3),DRL(3,3),DRE(3,3)
      DOUBLE PRECISION RQT(3,3),RUPT(3,3),RDT(3,3),RLT(3,3),RET(3,3)
!
      INTEGER I,J
!
      DO I=1,601
        GROT(I)=GUNROT(I)
      END DO
      DO I=1,3
        DO J=1,3
          DRQ(I,J)=DBLE(RQTOT(I,J))
          DRUP(I,J)=DBLE(RUPTOT(I,J))
          DRD(I,J)=DBLE(RDTOT(I,J))
          DRL(I,J)=DBLE(RLTOT(I,J))
          DRE(I,J)=DBLE(RETOT(I,J))
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=GUNROT(3+(I-1)*3+J)
          FD(I,J)=GUNROT(12+(I-1)*3+J)
          FE(I,J)=GUNROT(21+(I-1)*3+J)
          AU(I,J)=GUNROT(33+(I-1)*3+J)
          AD(I,J)=GUNROT(42+(I-1)*3+J)
          AE(I,J)=GUNROT(51+(I-1)*3+J)
          MQ(I,J)=GUNROT(62+(I-1)*3+J)
          ML(I,J)=GUNROT(71+(I-1)*3+J)
          MUP(I,J)=GUNROT(80+(I-1)*3+J)
          MD(I,J)=GUNROT(89+(I-1)*3+J)
          ME(I,J)=GUNROT(98+(I-1)*3+J)
!
          FUM(I,J)=GUNROT(293+(I-1)*3+J)
          FDM(I,J)=GUNROT(302+(I-1)*3+J)
          FEM(I,J)=GUNROT(311+(I-1)*3+J)
          AUM(I,J)=GUNROT(323+(I-1)*3+J)
          ADM(I,J)=GUNROT(332+(I-1)*3+J)
          AEM(I,J)=GUNROT(341+(I-1)*3+J)
          MQM(I,J)=GUNROT(352+(I-1)*3+J)
          MLM(I,J)=GUNROT(361+(I-1)*3+J)
          MUPM(I,J)=GUNROT(370+(I-1)*3+J)
          MDM(I,J)=GUNROT(379+(I-1)*3+J)
          MEM(I,J)=GUNROT(388+(I-1)*3+J)
!
          LU(I,J)=GUNROT(111+(I-1)*3+J)
          LD(I,J)=GUNROT(120+(I-1)*3+J)
          LE(I,J)=GUNROT(129+(I-1)*3+J)
!
          GTPQ(I,J)=GUNROT(138+(I-1)*3+J)
          GTPL(I,J)=GUNROT(147+(I-1)*3+J)
          GTPU(I,J)=GUNROT(156+(I-1)*3+J)
          GTPD(I,J)=GUNROT(165+(I-1)*3+J)
          GTPE(I,J)=GUNROT(174+(I-1)*3+J)
          GTQ(I,J)=GUNROT(185+(I-1)*3+J)
          GTL(I,J)=GUNROT(194+(I-1)*3+J)
          GTSQ(I,J)=GUNROT(205+(I-1)*3+J)
          GTSU(I,J)=GUNROT(214+(I-1)*3+J)
          GTSD(I,J)=GUNROT(223+(I-1)*3+J)
          FTUQ(I,J)=GUNROT(232+(I-1)*3+J)
          FTDQ(I,J)=GUNROT(241+(I-1)*3+J)
          FTEL(I,J)=GUNROT(250+(I-1)*3+J)
          FTUU(I,J)=GUNROT(259+(I-1)*3+J)
          FTDD(I,J)=GUNROT(268+(I-1)*3+J)
          FTEE(I,J)=GUNROT(277+(I-1)*3+J)
!
          TRIU(I,J)=GUNROT(399+(I-1)*3+J)
          TRID(I,J)=GUNROT(408+(I-1)*3+J)
          TRIE(I,J)=GUNROT(417+(I-1)*3+J)
          MTSFU(I,J)=GUNROT(429+(I-1)*3+J)
          MTSFD(I,J)=GUNROT(438+(I-1)*3+J)
          MTSFE(I,J)=GUNROT(447+(I-1)*3+J)
!
          RQT(I,J)=DRQ(J,I)
          RUPT(I,J)=DRUP(J,I)
          RDT(I,J)=DRD(J,I)
          RLT(I,J)=DRL(J,I)
          RET(I,J)=DRE(J,I)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUT(I,J)=DMATMUL(2,FU,RUPT,I,J)
          FDT(I,J)=DMATMUL(2,FD,RDT,I,J)
          FET(I,J)=DMATMUL(2,FE,RET,I,J)
          AUT(I,J)=DMATMUL(2,AU,RUPT,I,J)
          ADT(I,J)=DMATMUL(2,AD,RDT,I,J)
          AET(I,J)=DMATMUL(2,AE,RET,I,J)
          MQT(I,J)=DMATMUL(0,MQ,DRQ,I,J)
          MLT(I,J)=DMATMUL(0,ML,DRL,I,J)
          MUPT(I,J)=DMATMUL(0,MUP,DRUP,I,J)
          MDT(I,J)=DMATMUL(0,MD,DRD,I,J)
          MET(I,J)=DMATMUL(0,ME,DRE,I,J)
!
          FUMT(I,J)=DMATMUL(2,FUM,RUPT,I,J)
          FDMT(I,J)=DMATMUL(2,FDM,RDT,I,J)
          FEMT(I,J)=DMATMUL(2,FEM,RET,I,J)
          AUMT(I,J)=DMATMUL(2,AUM,RUPT,I,J)
          ADMT(I,J)=DMATMUL(2,ADM,RDT,I,J)
          AEMT(I,J)=DMATMUL(2,AEM,RET,I,J)
          MQMT(I,J)=DMATMUL(0,MQM,DRQ,I,J)
          MLMT(I,J)=DMATMUL(0,MLM,DRL,I,J)
          MUPMT(I,J)=DMATMUL(0,MUPM,DRUP,I,J)
          MDMT(I,J)=DMATMUL(0,MDM,DRD,I,J)
          MEMT(I,J)=DMATMUL(0,MEM,DRE,I,J)
!
          LUT(I,J)=DMATMUL(2,LU,RUPT,I,J)
          LDT(I,J)=DMATMUL(2,LD,RDT,I,J)
          LET(I,J)=DMATMUL(2,LE,RET,I,J)
!
          GTPQT(I,J)=DMATMUL(0,GTPQ,DRQ,I,J)
          GTPLT(I,J)=DMATMUL(0,GTPL,DRL,I,J)
          GTPUT(I,J)=DMATMUL(0,GTPU,DRUP,I,J)
          GTPDT(I,J)=DMATMUL(0,GTPD,DRD,I,J)
          GTPET(I,J)=DMATMUL(0,GTPE,DRE,I,J)
          GTQT(I,J)=DMATMUL(0,GTQ,DRQ,I,J)
          GTLT(I,J)=DMATMUL(0,GTL,DRL,I,J)
          GTSQT(I,J)=DMATMUL(0,GTSQ,DRQ,I,J)
          GTSUT(I,J)=DMATMUL(0,GTSU,DRUP,I,J)
          GTSDT(I,J)=DMATMUL(0,GTSD,DRD,I,J)
          FTUQT(I,J)=DMATMUL(2,FTUQ,RUPT,I,J)
          FTDQT(I,J)=DMATMUL(2,FTDQ,RDT,I,J)
          FTELT(I,J)=DMATMUL(2,FTEL,RET,I,J)
          FTUUT(I,J)=DMATMUL(2,FTUU,RUPT,I,J)
          FTDDT(I,J)=DMATMUL(2,FTDD,RDT,I,J)
          FTEET(I,J)=DMATMUL(2,FTEE,RET,I,J)
!
          TRIUT(I,J)=DMATMUL(2,TRIU,RUPT,I,J)
          TRIDT(I,J)=DMATMUL(2,TRID,RDT,I,J)
          TRIET(I,J)=DMATMUL(2,TRIE,RET,I,J)
          MTSFUT(I,J)=DMATMUL(2,MTSFU,RUPT,I,J)
          MTSFDT(I,J)=DMATMUL(2,MTSFD,RDT,I,J)
          MTSFET(I,J)=DMATMUL(2,MTSFE,RET,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GROT(3+(I-1)*3+J)=DMATMUL(0,RQT,FUT,I,J)
          GROT(12+(I-1)*3+J)=DMATMUL(0,RQT,FDT,I,J)
          GROT(21+(I-1)*3+J)=DMATMUL(0,RLT,FET,I,J)
          GROT(33+(I-1)*3+J)=DMATMUL(0,RQT,AUT,I,J)
          GROT(42+(I-1)*3+J)=DMATMUL(0,RQT,ADT,I,J)
          GROT(51+(I-1)*3+J)=DMATMUL(0,RLT,AET,I,J)
          GROT(62+(I-1)*3+J)=DMATMUL(1,DRQ,MQT,I,J)
          GROT(71+(I-1)*3+J)=DMATMUL(1,DRL,MLT,I,J)
          GROT(80+(I-1)*3+J)=DMATMUL(1,DRUP,MUPT,I,J)
          GROT(89+(I-1)*3+J)=DMATMUL(1,DRD,MDT,I,J)
          GROT(98+(I-1)*3+J)=DMATMUL(1,DRE,MET,I,J)
!
          GROT(293+(I-1)*3+J)=DMATMUL(0,RQT,FUMT,I,J)
          GROT(302+(I-1)*3+J)=DMATMUL(0,RQT,FDMT,I,J)
          GROT(311+(I-1)*3+J)=DMATMUL(0,RLT,FEMT,I,J)
          GROT(323+(I-1)*3+J)=DMATMUL(0,RQT,AUMT,I,J)
          GROT(332+(I-1)*3+J)=DMATMUL(0,RQT,ADMT,I,J)
          GROT(341+(I-1)*3+J)=DMATMUL(0,RLT,AEMT,I,J)
          GROT(352+(I-1)*3+J)=DMATMUL(1,DRQ,MQMT,I,J)
          GROT(361+(I-1)*3+J)=DMATMUL(1,DRL,MLMT,I,J)
          GROT(370+(I-1)*3+J)=DMATMUL(1,DRUP,MUPMT,I,J)
          GROT(379+(I-1)*3+J)=DMATMUL(1,DRD,MDMT,I,J)
          GROT(388+(I-1)*3+J)=DMATMUL(1,DRE,MEMT,I,J)
!
          GROT(111+(I-1)*3+J)=DMATMUL(0,RQT,LUT,I,J)
          GROT(120+(I-1)*3+J)=DMATMUL(0,RQT,LDT,I,J)
          GROT(129+(I-1)*3+J)=DMATMUL(0,RLT,LET,I,J)
!
          GROT(138+(I-1)*3+J)=DMATMUL(1,DRQ,GTPQT,I,J)
          GROT(147+(I-1)*3+J)=DMATMUL(1,DRL,GTPLT,I,J)
          GROT(156+(I-1)*3+J)=DMATMUL(1,DRUP,GTPUT,I,J)
          GROT(165+(I-1)*3+J)=DMATMUL(1,DRD,GTPDT,I,J)
          GROT(174+(I-1)*3+J)=DMATMUL(1,DRE,GTPET,I,J)
          GROT(185+(I-1)*3+J)=DMATMUL(1,DRQ,GTQT,I,J)
          GROT(194+(I-1)*3+J)=DMATMUL(1,DRL,GTLT,I,J)
          GROT(205+(I-1)*3+J)=DMATMUL(1,DRQ,GTSQT,I,J)
          GROT(214+(I-1)*3+J)=DMATMUL(1,DRUP,GTSUT,I,J)
          GROT(223+(I-1)*3+J)=DMATMUL(1,DRD,GTSDT,I,J)
          GROT(232+(I-1)*3+J)=DMATMUL(0,RQT,FTUQT,I,J)
          GROT(241+(I-1)*3+J)=DMATMUL(0,RQT,FTDQT,I,J)
          GROT(250+(I-1)*3+J)=DMATMUL(0,RLT,FTELT,I,J)
          GROT(259+(I-1)*3+J)=DMATMUL(0,RQT,FTUUT,I,J)
          GROT(268+(I-1)*3+J)=DMATMUL(0,RQT,FTDDT,I,J)
          GROT(277+(I-1)*3+J)=DMATMUL(0,RLT,FTEET,I,J)
!
          GROT(399+(I-1)*3+J)=DMATMUL(0,RQT,TRIUT,I,J)
          GROT(408+(I-1)*3+J)=DMATMUL(0,RQT,TRIDT,I,J)
          GROT(417+(I-1)*3+J)=DMATMUL(0,RLT,TRIET,I,J)
          GROT(429+(I-1)*3+J)=DMATMUL(0,RQT,MTSFUT,I,J)
          GROT(438+(I-1)*3+J)=DMATMUL(0,RQT,MTSFDT,I,J)
          GROT(447+(I-1)*3+J)=DMATMUL(0,RLT,MTSFET,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE DSMMA(G,Q,SVLQ)
!
!Purpose: Find the down squark mass matrix, sort the eigenvectors and
!         eigenvalues and then (if required) write them to a file
!         called sqm2d.dat
!
      IMPLICIT NONE
!
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,STRADD,SQM2D
      DOUBLE COMPLEX B1EVE(6,6),B1EVA(6),DSQM(6,6),COSTHB,SINTHB
      DOUBLE PRECISION MSB1,MSB2
      DOUBLE PRECISION SUM,Q
      DOUBLE COMPLEX G(601),CDSQM(6,6),TEMP(7)
      DOUBLE COMPLEX EVERTMP(6,6),CWORK1(99),CWORK2(12)
      INTEGER I,J,K,CIERR,SVLQ
!
      DO I=1,3
        DO J=1,3
          B1EVE(I,J)=(0.D0,0.D0)
        END DO
        B1EVA(I)=(0.D0,0.D0)
      END DO
      COSTHB=(0.D0,0.D0)
      SINTHB=(0.D0,0.D0)
      MSB1=0.D0
      MSB2=0.D0
!
      CALL DOWNSQM(G,Q,DSQM)
      DO I=1,6
        DO J=1,6
          CDSQM(I,J)=DSQM(I,J)
        END DO
      END DO
      SQM2D=STRADD(FNRGE,'.sqm2d')
      OPEN(25,FILE=SQM2D,STATUS='UNKNOWN',FORM='FORMATTED')
!
      WRITE(25,*)
      WRITE(25,*)'THE 6X6 DOWN-TYPE MATRIX AT THE SCALE: ',Q
      WRITE(25,*)'IN THE CHOSEN BASIS, IS:'
!
      WRITE(25,*)
      WRITE(25,15)'d_l','s_l','b_l','d_r','s_r','b_r'
      WRITE(25,*)
!
      WRITE(25,11)'d_l',DSQM(1,1),DSQM(1,2),DSQM(1,3),
     $            DSQM(1,4),DSQM(1,5),DSQM(1,6)
      WRITE(25,11)'s_l',DSQM(2,1),DSQM(2,2),DSQM(2,3),
     $            DSQM(2,4),DSQM(2,5),DSQM(2,6)
      WRITE(25,11)'b_l',DSQM(3,1),DSQM(3,2),DSQM(3,3),
     $            DSQM(3,4),DSQM(3,5),DSQM(3,6)
      WRITE(25,11)'d_r',DSQM(4,1),DSQM(4,2),DSQM(4,3),
     $            DSQM(4,4),DSQM(4,5),DSQM(4,6)
      WRITE(25,11)'s_r',DSQM(5,1),DSQM(5,2),DSQM(5,3),
     $            DSQM(5,4),DSQM(5,5),DSQM(5,6)
      WRITE(25,11)'b_r',DSQM(6,1),DSQM(6,2),DSQM(6,3),
     $            DSQM(6,4),DSQM(6,5),DSQM(6,6)
!
      WRITE(25,*)
!
      IF(SVLQ.EQ.0)THEN
!
!Diagonalise the matrix using the LAPACK routine ZGEEV.
!
        CALL ZGEEV('V','V',6,CDSQM,6,B1EVA,B1EVE,6,EVERTMP,6,CWORK1,99
     $             ,CWORK2,CIERR)
!
        DO J=1,6
          SUM=0.D0
          DO I=1,6
            SUM=SUM+DBLE(B1EVE(I,J)*CONJG(B1EVE(I,J)))
          END DO
          DO I=1,6
            B1EVE(I,J)=B1EVE(I,J)/DSQRT(SUM)
          END DO
        END DO
!
!Find the eigenvectors with the most amount of sbottom
!
        DO I=1,5
          DO J=I+1,6
            IF(DSQRT(DBLE(B1EVE(3,I)*CONJG(B1EVE(3,I)))+DBLE(B1EVE(6,I)
     $              *CONJG(B1EVE(6,I)))).LT.
     $         DSQRT(DBLE(B1EVE(3,J)*CONJG(B1EVE(3,J)))+DBLE(B1EVE(6,J)
     $              *CONJG(B1EVE(6,J)))))THEN
              DO K=1,6
                TEMP(K)=B1EVE(K,I)
                B1EVE(K,I)=B1EVE(K,J)
                B1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=B1EVA(I)
              B1EVA(I)=B1EVA(J)
              B1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Now find next two with most amount of scalar strange
!
        DO I=3,5
          DO J=I+1,6
            IF(DSQRT(DBLE(B1EVE(2,I)*CONJG(B1EVE(2,I)))+DBLE(B1EVE(5,I)
     $              *CONJG(B1EVE(5,I)))).LT.
     $         DSQRT(DBLE(B1EVE(2,J)*CONJG(B1EVE(2,J)))+DBLE(B1EVE(5,J)
     $              *CONJG(B1EVE(5,J)))))THEN
              DO K=1,6
                TEMP(K)=B1EVE(K,I)
                B1EVE(K,I)=B1EVE(K,J)
                B1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=B1EVA(I)
              B1EVA(I)=B1EVA(J)
              B1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Sort the first two eigenvectors by mass
!
        IF(ABS(B1EVA(1)).GT.ABS(B1EVA(2)))THEN
          DO I=1,6
            TEMP(I)=B1EVE(I,1)
            B1EVE(I,1)=B1EVE(I,2)
            B1EVE(I,2)=TEMP(I)
          END DO
          TEMP(7)=B1EVA(1)
          B1EVA(1)=B1EVA(2)
          B1EVA(2)=TEMP(7)
        END IF
!
!Sort the next two by mass
!
        IF(ABS(B1EVA(3)).GT.ABS(B1EVA(4)))THEN
          DO I=1,6
            TEMP(I)=B1EVE(I,3)
            B1EVE(I,3)=B1EVE(I,4)
            B1EVE(I,4)=TEMP(I)
          END DO
          TEMP(7)=B1EVA(3)
          B1EVA(3)=B1EVA(4)
          B1EVA(4)=TEMP(7)
        END IF
!
!And the final two
!
        IF(ABS(B1EVA(5)).GT.ABS(B1EVA(6)))THEN
          DO I=1,6
            TEMP(I)=B1EVE(I,5)
            B1EVE(I,5)=B1EVE(I,6)
            B1EVE(I,6)=TEMP(I)
          END DO
          TEMP(7)=B1EVA(5)
          B1EVA(5)=B1EVA(6)
          B1EVA(6)=TEMP(7)
        END IF
!
!Note about sin{\theta_b}: If I define sin{\theta_b} the same way as
!BT, then -sin{\theta_b}=conjg(b1eve(6,1))
!
        COSTHB=CONJG(B1EVE(3,1))
        SINTHB=-CONJG(B1EVE(6,1))
        MSB1=DSQRT(ABS(B1EVA(1)))
        MSB2=DSQRT(ABS(B1EVA(2)))
!
        WRITE(25,*)'EIGENVECTORS ARE:'
        WRITE(25,*)
        WRITE(25,15)'b_1','b_2','s_1','s_2','d_1','d_2'
        WRITE(25,*)
!     
        WRITE(25,11)'d_l',B1EVE(1,1),B1EVE(1,2),B1EVE(1,3),
     $              B1EVE(1,4),B1EVE(1,5),B1EVE(1,6)
        WRITE(25,11)'s_l',B1EVE(2,1),B1EVE(2,2),B1EVE(2,3),
     $              B1EVE(2,4),B1EVE(2,5),B1EVE(2,6)
        WRITE(25,11)'b_l',B1EVE(3,1),B1EVE(3,2),B1EVE(3,3),
     $              B1EVE(3,4),B1EVE(3,5),B1EVE(3,6)
        WRITE(25,11)'d_r',B1EVE(4,1),B1EVE(4,2),B1EVE(4,3),
     $              B1EVE(4,4),B1EVE(4,5),B1EVE(4,6)
        WRITE(25,11)'s_r',B1EVE(5,1),B1EVE(5,2),B1EVE(5,3),
     $              B1EVE(5,4),B1EVE(5,5),B1EVE(5,6)
        WRITE(25,11)'b_r',B1EVE(6,1),B1EVE(6,2),B1EVE(6,3),
     $              B1EVE(6,4),B1EVE(6,5),B1EVE(6,6)
!
        WRITE(25,*)
        WRITE(25,*)'EIGENVALUES ARE:'
        WRITE(25,16)B1EVA(1),B1EVA(2),B1EVA(3),B1EVA(4),B1EVA(5)
     $             ,B1EVA(6)
        WRITE(25,*)
        WRITE(25,*)'MASSES ARE:'
        WRITE(25,17)DSQRT(DSQRT(DBLE(B1EVA(1)*CONJG(B1EVA(1)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(2)*CONJG(B1EVA(2)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(3)*CONJG(B1EVA(3)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(4)*CONJG(B1EVA(4)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(5)*CONJG(B1EVA(5)))))
     $             ,DSQRT(DSQRT(DBLE(B1EVA(6)*CONJG(B1EVA(6)))))
      END IF
!
!
 11   FORMAT(SP,1P,A,3X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 15   FORMAT(15X,A,19X,A,19X,A,19X,A,19X,A,19X,A)
 16   FORMAT(SP,1P,6X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 17   FORMAT(10X,F11.4,11X,F11.4,11X,F11.4,11X,F11.4,11X
     $      ,F11.4,11X,F11.4)
      WRITE(25,*)
      CLOSE(25)
!
      RETURN
      END
!
      FUNCTION DTRACE(A)
!
!Purpose: To compute the trace of matrix A
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION DTRACE
      INTEGER I
!
      DTRACE=0.D0
      DO I=1,3
        DTRACE=DTRACE+A(I,I)
      END DO
!
      RETURN
      END
!
      SUBROUTINE GOUT601(G,QEND,PRCH,OUTFILE)
!
!Purpose: Writes all the results in a readable format
!         PRCH is a switch to alter the output.
!             =1 -> output with no tilde-couplings or lambdas,
!                   most appropriate for high scale output
!             =2 -> output with tilde-couplings and lambdas,
!                   no regular Yukawas,
!                   most appropriate for low scale output
!             =3 -> everything,
!                   only appropriate at m_H
!
      IMPLICIT NONE
!
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,STRADD,GUTOUT,WKOUT
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      DOUBLE PRECISION QEND
      DOUBLE COMPLEX G(601)
      CHARACTER*20 OUTFILE
      INTEGER I,SW,PRCH
!
!Now write out the results
!
      OPEN(15,FILE=OUTFILE,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(15,*)
      WRITE(15,*)'Q IS',QEND,' COMP IS',COMP
      WRITE(15,*)
!
      WRITE(15,*)'COUPLINGS:'
      WRITE(15,10)G(1),G(2),G(3)
      WRITE(15,*)
!
      IF(PRCH.NE.2)THEN
        WRITE(15,*)'f_U:'
        DO I=1,3
          WRITE(15,10)G(4+3*(I-1)),G(5+3*(I-1)),G(6+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'f_D:'
        DO I=1,3
          WRITE(15,10)G(13+3*(I-1)),G(14+3*(I-1)),G(15+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'f_E:'
        DO I=1,3
          WRITE(15,10)G(22+3*(I-1)),G(23+3*(I-1)),G(24+3*(I-1))
        END DO
        WRITE(15,*)
      END IF
!
      WRITE(15,*)'GAUGINO MASSES M:'
      WRITE(15,10)G(31),G(32),G(33)
      WRITE(15,*)
!
      IF(COMP.EQ.1)THEN
        WRITE(15,*)'GAUGINO MASSES M'':'
        WRITE(15,10)G(599),G(600),G(601)
        WRITE(15,*)
      END IF
!
      WRITE(15,*)'a_U:'
      DO I=1,3
        WRITE(15,10)G(34+3*(I-1)),G(35+3*(I-1)),G(36+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'a_D:'
      DO I=1,3
        WRITE(15,10)G(43+3*(I-1)),G(44+3*(I-1)),G(45+3*(I-1)) 
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'a_E:'
      DO I=1,3
        WRITE(15,10)G(52+3*(I-1)),G(53+3*(I-1)),G(54+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'SQUARED HIGGS MASSES:'
      WRITE(15,5)G(61)-G(108)**2,G(62)-G(108)**2
      WRITE(15,*)
!
      WRITE(15,*)'M_Q^2:'
      DO I=1,3
        WRITE(15,10)G(63+3*(I-1)),G(64+3*(I-1)),G(65+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_L^2:'
      DO I=1,3
        WRITE(15,10)G(72+3*(I-1)),G(73+3*(I-1)),G(74+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_U^2:'
      DO I=1,3
        WRITE(15,10)G(81+3*(I-1)),G(82+3*(I-1)),G(83+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_D^2:'
      DO I=1,3
        WRITE(15,10)G(90+3*(I-1)),G(91+3*(I-1)),G(92+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'M_E^2:'
      DO I=1,3
        WRITE(15,10)G(99+3*(I-1)),G(100+3*(I-1)),G(101+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MU AND B:'
      WRITE(15,5)G(108),G(109)
      WRITE(15,*)
      WRITE(15,*)'V_U AND V_D:'
      WRITE(15,5)G(110),G(111)
      WRITE(15,*)
!
      IF(PRCH.NE.1)THEN
        WRITE(15,*)'lambda_U:'
        DO I=1,3
          WRITE(15,10)G(112+3*(I-1)),G(113+3*(I-1)),G(114+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'lambda_D:'
        DO I=1,3
          WRITE(15,10)G(121+3*(I-1)),G(122+3*(I-1)),G(123+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'lambda_E:'
        DO I=1,3
          WRITE(15,10)G(130+3*(I-1)),G(131+3*(I-1)),G(132+3*(I-1))
        END DO
        WRITE(15,*)
  
        WRITE(15,*)'GTP_Q:'
        DO I=1,3
          WRITE(15,10)G(139+3*(I-1)),G(140+3*(I-1)),G(141+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_L:'
        DO I=1,3
          WRITE(15,10)G(148+3*(I-1)),G(149+3*(I-1)),G(150+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_U:'
        DO I=1,3
          WRITE(15,10)G(157+3*(I-1)),G(158+3*(I-1)),G(159+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_D:'
        DO I=1,3
          WRITE(15,10)G(166+3*(I-1)),G(167+3*(I-1)),G(168+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTP_E:'
        DO I=1,3
          WRITE(15,10)G(175+3*(I-1)),G(176+3*(I-1)),G(177+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTPH_U AND GTPH_D:'
        WRITE(15,5)G(184),G(185)
        WRITE(15,*)
!
        WRITE(15,*)'GT_Q:'
        DO I=1,3
          WRITE(15,10)G(186+3*(I-1)),G(187+3*(I-1)),G(188+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GT_L:'
        DO I=1,3
          WRITE(15,10)G(195+3*(I-1)),G(196+3*(I-1)),G(197+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTH_U AND GTH_D:'
        WRITE(15,5)G(204),G(205)
        WRITE(15,*)
!
        WRITE(15,*)'GTS_Q:'
        DO I=1,3
          WRITE(15,10)G(206+3*(I-1)),G(207+3*(I-1)),G(208+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTS_U:'
        DO I=1,3
          WRITE(15,10)G(215+3*(I-1)),G(216+3*(I-1)),G(217+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'GTS_D:'
        DO I=1,3
          WRITE(15,10)G(224+3*(I-1)),G(225+3*(I-1)),G(226+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTQ_U:'
        DO I=1,3
          WRITE(15,10)G(233+3*(I-1)),G(234+3*(I-1)),G(235+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTQ_D:'
        DO I=1,3
          WRITE(15,10)G(242+3*(I-1)),G(243+3*(I-1)),G(244+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTL_E:'
        DO I=1,3
          WRITE(15,10)G(251+3*(I-1)),G(252+3*(I-1)),G(253+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTU_U:'
        DO I=1,3
          WRITE(15,10)G(260+3*(I-1)),G(261+3*(I-1)),G(262+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTD_D:'
        DO I=1,3
          WRITE(15,10)G(269+3*(I-1)),G(270+3*(I-1)),G(271+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'FTE_E:'
        DO I=1,3
          WRITE(15,10)G(278+3*(I-1)),G(279+3*(I-1)),G(280+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'sGTPH_U AND cGTPH_D:'
        WRITE(15,5)G(287),G(288)
        WRITE(15,*)
!
        WRITE(15,*)'sGTH_U AND cGTH_D:'
        WRITE(15,5)G(289),G(290)
        WRITE(15,*)
        WRITE(15,*)
      END IF
!
      WRITE(15,*)'MSSM COUPLINGS:'
      WRITE(15,10)G(291),G(292),G(293)
      WRITE(15,*)
!
      WRITE(15,*)'MSSM f_U:'
      DO I=1,3
        WRITE(15,10)G(294+3*(I-1)),G(295+3*(I-1)),G(296+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM f_D:'
      DO I=1,3
        WRITE(15,10)G(303+3*(I-1)),G(304+3*(I-1)),G(305+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM f_E:'
      DO I=1,3
        WRITE(15,10)G(312+3*(I-1)),G(313+3*(I-1)),G(314+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM GAUGINO MASSES:'
      WRITE(15,10)G(321),G(322),G(323)
      WRITE(15,*)
!
      WRITE(15,*)'MSSM a_U:'
      DO I=1,3
        WRITE(15,10)G(324+3*(I-1)),G(325+3*(I-1)),G(326+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM a_D:'
      DO I=1,3
        WRITE(15,10)G(333+3*(I-1)),G(334+3*(I-1)),G(335+3*(I-1)) 
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM a_E:'
      DO I=1,3
        WRITE(15,10)G(342+3*(I-1)),G(343+3*(I-1)),G(344+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM SQUARED HIGGS MASSES:'
      WRITE(15,5)G(351),G(352)
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_Q^2:'
      DO I=1,3
        WRITE(15,10)G(353+3*(I-1)),G(354+3*(I-1)),G(355+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_L^2:'
      DO I=1,3
        WRITE(15,10)G(362+3*(I-1)),G(363+3*(I-1)),G(364+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_U^2:'
      DO I=1,3
        WRITE(15,10)G(371+3*(I-1)),G(372+3*(I-1)),G(373+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_D^2:'
      DO I=1,3
        WRITE(15,10)G(380+3*(I-1)),G(381+3*(I-1)),G(382+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM M_E^2:'
      DO I=1,3
        WRITE(15,10)G(389+3*(I-1)),G(390+3*(I-1)),G(391+3*(I-1))
      END DO
      WRITE(15,*)
!
      WRITE(15,*)'MSSM MU AND B:'
      WRITE(15,5)G(398),G(399)
      WRITE(15,*)
!
      IF(PRCH.NE.1)THEN
        WRITE(15,*)'TRI_U:'
        DO I=1,3
          WRITE(15,10)G(400+3*(I-1)),G(401+3*(I-1)),G(402+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'TRI_D:'
        DO I=1,3
          WRITE(15,10)G(409+3*(I-1)),G(410+3*(I-1)),G(411+3*(I-1)) 
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'TRI_E:'
        DO I=1,3
          WRITE(15,10)G(418+3*(I-1)),G(419+3*(I-1)),G(420+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'M_HUD:'
        WRITE(15,7)G(427)
        WRITE(15,*)
!
        WRITE(15,*)'SM VEV AND LAMBDA_SM:'
        WRITE(15,5)G(428),G(429)
        WRITE(15,*)
!
        WRITE(15,*)'MTSF_U:'
        DO I=1,3
          WRITE(15,10)G(430+3*(I-1)),G(431+3*(I-1)),G(432+3*(I-1))
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'MTSF_D:'
        DO I=1,3
          WRITE(15,10)G(439+3*(I-1)),G(440+3*(I-1)),G(441+3*(I-1)) 
        END DO
        WRITE(15,*)
!
        WRITE(15,*)'MTSF_E:'
        DO I=1,3
          WRITE(15,10)G(448+3*(I-1)),G(449+3*(I-1)),G(450+3*(I-1))
        END DO
        WRITE(15,*)
      END IF
!
      CLOSE(15)
!
 5    FORMAT(SP,1P,4X,'(',D11.4,',',D11.4,')',2X,'(',D11.4,',',D11.4,')'
     $       ) !Used to be 11.4
 7    FORMAT(SP,1P,4X,'(',D11.4,',',D11.4,')')
 10   FORMAT(SP,1P,4X,'(',D11.4,',',D11.4,')',2X,'(',D11.4,',',D11.4,')'
     $      ,2X,'(',D11.4,',',D11.4,')')
!
      RETURN
      END
!
      SUBROUTINE HERMTEST(IN,OUT)
!
!Purpose: Test if complex matrix IN(3,3) is hermitian.
!         Return 1 if not, 0 otherwise.
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX IN(3,3)
      INTEGER OUT
!
      OUT=0
      IF(DIMAG(IN(1,1)).NE.0.D0)OUT=1
      IF(DIMAG(IN(2,2)).NE.0.D0)OUT=1
      IF(DIMAG(IN(3,3)).NE.0.D0)OUT=1
      IF(IN(1,2).NE.CONJG(IN(2,1)))OUT=1
      IF(IN(1,3).NE.CONJG(IN(3,1)))OUT=1
      IF(IN(2,3).NE.CONJG(IN(3,2)))OUT=1
!
      RETURN
      END
!
      SUBROUTINE HIGHIN(TIME)
!
!Purpose: To put the high scale conditions into G(601) ready
!         for running back down
!
!The switch 'TIME' is used to tell if this is the first call to
!HIGHIN. If so, the MSSM mu parameter is also set here.
!In future calls it will be set along with b at m_H.
!
!NOTE ABOUT INPUT OF SOFT MATRICES:
!
!All soft matrices are input in the current basis.
!To ensure that it is possible to have good control over additional
!flavour violating inputs, it is possible to set the coefficients
!for three independent matrices which have been carefully chosen 
!so that they do not introduce new flavour structure. In addition
!to the coefficients of these three matrices, there is also the
!possibility for introducing an entirely general matrix, which should
!be chosen with care.
!
!a_{ude}(GUT)=f_{ude}[A_{{ude}}\1+W_{ude}f^\dagger_{ude}f_{ude}
!         +X_{ude}f^\dagger_{ude}f_{ude}f^\dagger_{ude}f_{ude}]+Z_{ude}
!
!m^2_{UDE}(GUT)=m^2_{{UDE}0}[\1+R_{UDE}f^T_{ude}f^*_{ude}
!                 +S_{UDE}f^T_{ude}f^*_{ude}f^T_{ude}f^*_{ude}]+T_{UDE}
!
!m^2_{QL}(GUT)=m^2_{{QL}0}\1+T_{QL}
!
!where \1 is the unit matrix, denoted in the code by CID(I,J).
!
      IMPLICIT NONE
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/HISUSY/CM1,CM2,CM3,CM1P,CM2P,CM3P
     $             ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $             ,CTQ,CTU,CTD,CTL,CTE
     $             ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $             ,CZU,CZD,CZE,CMHU,CMHD,CU,CD,CE
      DOUBLE COMPLEX CM1,CM2,CM3,CM1P,CM2P,CM3P
     $              ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $              ,CTQ(3,3),CTU(3,3),CTD(3,3),CTL(3,3),CTE(3,3)
     $              ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $              ,CZU(3,3),CZD(3,3),CZE(3,3),CMHU,CMHD
      INTEGER CU,CD,CE
      SAVE/HISUSY/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
      COMMON /RGEFNM/ FNRGE
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3),FUS(3,3),FDS(3,3),FES(3,3)
     $              ,AU(3,3),AD(3,3),AE(3,3)
     $              ,MQSQ(3,3),MLSQ(3,3),MUPSQ(3,3),MDSQ(3,3),MESQ(3,3)
     $              ,FUDFU(3,3),FDDFD(3,3),FEDFE(3,3)
     $              ,FUDFUFUDFU(3,3),FDDFDFDDFD(3,3),FEDFEFEDFE(3,3)
     $              ,FUTFUS(3,3),FDTFDS(3,3),FETFES(3,3)
     $              ,FUTFUSFUTFUS(3,3),FDTFDSFDTFDS(3,3)
     $              ,FETFESFETFES(3,3)
      DOUBLE COMPLEX DUMAU(3,3),DUMAD(3,3),DUMAE(3,3),FUSFUT(3,3)
     $              ,FDSFDT(3,3)
      DOUBLE COMPLEX CID(3,3),CMATMUL,GTMP(601)
!
      DOUBLE PRECISION MUSQ,TU,TD
      INTEGER I,J,K,TIME
      CHARACTER*128 FILENAME,STRADD,FNRGE
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!Set the initial values for the old number of matter sfermions
!and their rotations.
!
      OLDNSQ=3
      OLDNSU=3
      OLDNSD=3
      OLDNSL=3
      OLDNSE=3
      DO I=1,3
        THSQ(I)=1
        THSU(I)=1
        THSD(I)=1
        THSL(I)=1
        THSE(I)=1
      END DO
      DO I=1,3
        DO J=1,3
          RQTOT(I,J)=CID(I,J)
          RUPTOT(I,J)=CID(I,J)
          RDTOT(I,J)=CID(I,J)
          RLTOT(I,J)=CID(I,J)
          RETOT(I,J)=CID(I,J)
          DO K=1,2
            RQSAV(K,I,J)=CID(I,J)
            RUPSAV(K,I,J)=CID(I,J)
            RDSAV(K,I,J)=CID(I,J)
            RLSAV(K,I,J)=CID(I,J)
            RESAV(K,I,J)=CID(I,J)            
          END DO
        END DO
      END DO
!
!Convert the G's to Matrices
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          FE(I,J)=G(21+(I-1)*3+J)
        END DO
      END DO
!
!If we are using mSUGRA boundary conditions
!
      IF(SUG.EQ.1)THEN
!
        DO I=1,3
          DO J=1,3
            AU(I,J)=DCMPLX(A0)*FU(I,J)
            AD(I,J)=DCMPLX(A0)*FD(I,J)
            AE(I,J)=DCMPLX(A0)*FE(I,J)
            MQSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MLSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MUPSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MDSQ(I,J)=DCMPLX(M0**2)*CID(I,J)
            MESQ(I,J)=DCMPLX(M0**2)*CID(I,J)
          END DO
        END DO
!
        G(31)=DCMPLX(M12)
        G(32)=DCMPLX(M12)
        G(33)=DCMPLX(M12)
        G(599)=(0.D0,0.D0)
        G(600)=(0.D0,0.D0)
        G(601)=(0.D0,0.D0)
        IF(TIME.EQ.1)THEN
          G(61)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
          G(62)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
        ELSE
!
!(G(61)+G(62)) WAS FIXED BY EWSB AT M_SUSY
!ABS(G(108)) WAS FIXED AT HIGHEST SUSY THRESHOLD
!
          G(108)=ABS(G(108))*EXP((0.D0,1.D0)*PHASEMU)
          G(61)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
          G(62)=DCMPLX(M0**2)+G(108)*CONJG(G(108))
        END IF
        G(351)=DCMPLX(M0**2)
        G(352)=DCMPLX(M0**2)
!
!Or use the non-mSUGRA matrices entered by the user
!
      ELSE
!
        DO I=1,3
          DO J=1,3
            FUS(I,J)=CONJG(FU(I,J))
            FDS(I,J)=CONJG(FD(I,J))
            FES(I,J)=CONJG(FE(I,J))
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            FUDFU(I,J)=CMATMUL(1,FU,FU,I,J)
            FDDFD(I,J)=CMATMUL(1,FD,FD,I,J)
            FEDFE(I,J)=CMATMUL(1,FE,FE,I,J)
            FUTFUS(I,J)=CMATMUL(1,FUS,FUS,I,J)
            FDTFDS(I,J)=CMATMUL(1,FDS,FDS,I,J)
            FETFES(I,J)=CMATMUL(1,FES,FES,I,J)
            FUSFUT(I,J)=CMATMUL(2,FUS,FUS,I,J)
            FDSFDT(I,J)=CMATMUL(2,FDS,FDS,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            FUTFUSFUTFUS(I,J)=CMATMUL(0,FUTFUS,FUTFUS,I,J)
            FDTFDSFDTFDS(I,J)=CMATMUL(0,FDTFDS,FDTFDS,I,J)
            FETFESFETFES(I,J)=CMATMUL(0,FETFES,FETFES,I,J)
            FUDFUFUDFU(I,J)=CMATMUL(0,FUDFU,FUDFU,I,J)
            FDDFDFDDFD(I,J)=CMATMUL(0,FDDFD,FDDFD,I,J)
            FEDFEFEDFE(I,J)=CMATMUL(0,FEDFE,FEDFE,I,J)
          END DO
        END DO
!
!Here we use the relations in the comment at the start of this
!subroutine
!
        DO I=1,3
          DO J=1,3
            DUMAU(I,J)=CAU0*CID(I,J)+CWU*FUDFU(I,J)+CXU*FUDFUFUDFU(I,J)
            DUMAD(I,J)=CAD0*CID(I,J)+CWD*FDDFD(I,J)+CXD*FDDFDFDDFD(I,J)
            DUMAE(I,J)=CAE0*CID(I,J)+CWE*FEDFE(I,J)+CXE*FEDFEFEDFE(I,J)
          END DO
        END DO
!
!The next two lines can be hardwired to allow for entry of a T_Q
!which is proportional to the Yukawas, a la Minimal Flavour Violation
!See arXiv:0810.5765.
!
        TU=0.D0
        TD=0.D0
        IF(TU.NE.0.D0.OR.TD.NE.0.D0)THEN
          WRITE(*,*)'NON-ZERO MFV CONDITIONS HARDWIRED AT GUT'
          WRITE(*,*)'TU=',TU
          WRITE(*,*)'TD=',TD
        END IF
        DO I=1,3
          DO J=1,3
            AU(I,J)=CMATMUL(0,FU,DUMAU,I,J)+CZU(I,J)
            AD(I,J)=CMATMUL(0,FD,DUMAD,I,J)+CZD(I,J)
            AE(I,J)=CMATMUL(0,FE,DUMAE,I,J)+CZE(I,J)
!
            MQSQ(I,J)=CMQ0**2*CID(I,J)+CTQ(I,J)
            MQSQ(I,J)=CMQ0**2*(CID(I,J)+TU*FUSFUT(I,J)
     $                        +TD*FDSFDT(I,J))+CTQ(I,J)
            MUPSQ(I,J)=CMU0**2*(DBLE(CU)*CID(I,J)+CRU*FUTFUS(I,J)
     $                         +CSU*FUTFUSFUTFUS(I,J))+CTU(I,J)
            MDSQ(I,J)=CMD0**2*(DBLE(CD)*CID(I,J)+CRD*FDTFDS(I,J)
     $                         +CSD*FDTFDSFDTFDS(I,J))+CTD(I,J)
            MLSQ(I,J)=CML0**2*CID(I,J)+CTL(I,J)
            MESQ(I,J)=CME0**2*(DBLE(CE)*CID(I,J)+CRE*FETFES(I,J)
     $                         +CSE*FETFESFETFES(I,J))+CTE(I,J)
          END DO
        END DO
!
        G(31)=CM1
        G(32)=CM2
        G(33)=CM3
        G(599)=CM1P
        G(600)=CM2P
        G(601)=CM3P
        IF(TIME.EQ.1)THEN
          G(61)=CMHU**2+G(108)*CONJG(G(108))
          G(62)=CMHD**2+G(108)*CONJG(G(108))
        ELSE
!
!(G(61)+G(62)) WAS FIXED BY EWSB AT M_SUSY
!ABS(G(108)) WAS FIXED AT HIGHEST SUSY THRESHOLD
!
          G(108)=ABS(G(108))*EXP((0.D0,1.D0)*PHASEMU)
          G(61)=CMHU**2+G(108)*CONJG(G(108))
          G(62)=CMHD**2+G(108)*CONJG(G(108))
        END IF
        G(351)=CMHU**2
        G(352)=CMHD**2
      END IF
!
!Now put back into the g's
!
      DO I=1,3
        DO J=1,3
          G(3+(I-1)*3+J)=FU(I,J)
          G(12+(I-1)*3+J)=FD(I,J)
          G(21+(I-1)*3+J)=FE(I,J)
          G(33+(I-1)*3+J)=AU(I,J)
          G(42+(I-1)*3+J)=AD(I,J)
          G(51+(I-1)*3+J)=AE(I,J)
          G(62+(I-1)*3+J)=MQSQ(I,J)
          G(71+(I-1)*3+J)=MLSQ(I,J)
          G(80+(I-1)*3+J)=MUPSQ(I,J)
          G(89+(I-1)*3+J)=MDSQ(I,J)
          G(98+(I-1)*3+J)=MESQ(I,J)
        END DO
      END DO
!
!Next, I need to set the starting values for the tilde couplings and
!MSSM terms
!
      DO I=1,3
        DO J=1,3
          G(138+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(147+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(156+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(165+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(174+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(185+(I-1)*3+J)=G(2)*CID(I,J)
          G(194+(I-1)*3+J)=G(2)*CID(I,J)
          G(205+(I-1)*3+J)=G(3)*CID(I,J)
          G(214+(I-1)*3+J)=G(3)*CID(I,J)
          G(223+(I-1)*3+J)=G(3)*CID(I,J)
          G(232+(I-1)*3+J)=G(3+(I-1)*3+J)
          G(241+(I-1)*3+J)=G(12+(I-1)*3+J)
          G(250+(I-1)*3+J)=G(21+(I-1)*3+J)
          G(259+(I-1)*3+J)=G(3+(I-1)*3+J)
          G(268+(I-1)*3+J)=G(12+(I-1)*3+J)
          G(277+(I-1)*3+J)=G(21+(I-1)*3+J)
          G(293+(I-1)*3+J)=G(3+(I-1)*3+J)
          G(302+(I-1)*3+J)=G(12+(I-1)*3+J)
          G(311+(I-1)*3+J)=G(21+(I-1)*3+J)
          G(323+(I-1)*3+J)=G(33+(I-1)*3+J)
          G(332+(I-1)*3+J)=G(42+(I-1)*3+J)
          G(341+(I-1)*3+J)=G(51+(I-1)*3+J)
          G(352+(I-1)*3+J)=G(62+(I-1)*3+J)
          G(361+(I-1)*3+J)=G(71+(I-1)*3+J)
          G(370+(I-1)*3+J)=G(80+(I-1)*3+J)
          G(379+(I-1)*3+J)=G(89+(I-1)*3+J)
          G(388+(I-1)*3+J)=G(98+(I-1)*3+J)
          G(429+(I-1)*3+J)=CONJG(G(108))*G(3+(I-1)*3+J)
          G(438+(I-1)*3+J)=CONJG(G(108))*G(12+(I-1)*3+J)
          G(447+(I-1)*3+J)=CONJG(G(108))*G(21+(I-1)*3+J)
        END DO
      END DO
      G(184)=DSQRT(3.D0/5.D0)*G(1)
      G(185)=DSQRT(3.D0/5.D0)*G(1)
      G(204)=G(2)
      G(205)=G(2)
      G(291)=G(1)
      G(292)=G(2)
      G(293)=G(3)
      G(321)=G(31)-(0.D0,1.D0)*G(599)
      G(322)=G(32)-(0.D0,1.D0)*G(600)
      G(323)=G(33)-(0.D0,1.D0)*G(601)
      IF(TIME.EQ.1)THEN
        G(398)=G(108)
      END IF
!
!Rotate the the basis where the Yukawas are diagonal at m_t and
!write the results to a file.
!
      DO I=1,601 !Save the original couplings
        GTMP(I)=G(I)
      END DO
      CALL ROTBACK(0)
      FILENAME=STRADD(FNRGE,'.gtout')
      CALL GOUT601(G,MHIGH,1,FILENAME)
!
!Replace the rotated couplings by the old saved ones
!This removes any inaccuracies introduced by the rotation.
!
      DO I=1,601
        G(I)=GTMP(I)
      END DO
!
      RETURN
      END
!
      SUBROUTINE KMIN
!
!Purpose: Contains the entries of the KM matrix.
!         From PDG, using central values
!
      IMPLICIT NONE
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      DOUBLE PRECISION S12,C12,S13,C13,S23,C23,D13
      INTEGER I
!
      S12=0.2243D0
      S23=0.0413D0
      S13=0.0037D0
      IF(COMP.EQ.1)THEN
        D13=1.047198D0 !This is the phase in the KM
      ELSE
        D13=0.D0
      END IF
!
      C12=DCOS(DASIN(S12))
      C23=DCOS(DASIN(S23))
      C13=DCOS(DASIN(S13))
!
      KM(1,1)=DCMPLX(C12*C13)
      KM(1,2)=DCMPLX(S12*C13)
      KM(1,3)=S13*DCMPLX(DCOS(D13),-DSIN(D13))
      KM(2,1)=-DCMPLX(S12*C23)-C12*S23*S13*DCMPLX(DCOS(D13),DSIN(D13))
      KM(2,2)=DCMPLX(C12*C23)-S12*S23*S13*DCMPLX(DCOS(D13),DSIN(D13))
      KM(2,3)=DCMPLX(S23*C13)
      KM(3,1)=DCMPLX(S12*S23)-C12*C23*S13*DCMPLX(DCOS(D13),DSIN(D13))
      KM(3,2)=-DCMPLX(C12*S23)-S12*C23*S13*DCMPLX(DCOS(D13),DSIN(D13))
      KM(3,3)=DCMPLX(C23*C13)
!
      RETURN
      END
!
      SUBROUTINE MASS
!
!Purpose: Input the low scale masses of the quarks and leptons
!         in the order u,c,d,s,e,mu and thus derive the yukawas.
!         From paper by Fusaoka and Koide.
!
      IMPLICIT NONE
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      DOUBLE PRECISION MATQ(2,3,3),MATL(3,3)
      INTEGER I,J
!
!Use same (Pierce) prescription as SUGRA.F for vev at MZ
!Eq. 19 NPB491 p.3 (1997)
!Remove log since it is presumably the bulk of the SUSY 
!contribution if all sparticles are decoupled at M_Z
!
!      VSMMZ=DCMPLX((248.6D0+0.9D0*LOG(RGEMS/MZ))/DSQRT(2.D0))
      VSMMZ=DCMPLX(248.6D0/DSQRT(2.D0))
!
      MWEAK(1)=.00233D0
      MWEAK(2)=.677D0
      MWEAK(3)=.00469D0
      MWEAK(4)=.0934D0
      MWEAK(5)=0.4873D-3
      MWEAK(6)=102.87D-3
!
!The values of f_b, f_t and f_tau are derived from RGEs
!
      RETURN
      END
!
      SUBROUTINE MASSSQM(GRUN)
!
!Purpose: To construct the quark mass basis soft mass matrices.
!
      IMPLICIT NONE
!
!Needed to save the values of mq and mu for use in decay calculations
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE COMPLEX EQ(3,3),EU(3,3),ED(3,3),EL(3,3),EE(3,3)
      DOUBLE COMPLEX MQMASS(3,3),MUMASS(3,3),MDMASS(3,3)
      DOUBLE COMPLEX MLMASS(3,3),MEMASS(3,3)
      DOUBLE COMPLEX MQTMP(3,3),MUTMP(3,3),MDTMP(3,3)
      DOUBLE COMPLEX MLTMP(3,3),METMP(3,3)
      DOUBLE COMPLEX MQCURR(3,3),MUCURR(3,3),MDCURR(3,3)
      DOUBLE COMPLEX MLCURR(3,3),MECURR(3,3)
      DOUBLE COMPLEX SUM,RQTMP(3,3),RUTMP(3,3),RDTMP(3,3)
      DOUBLE COMPLEX RLTMP(3,3),RETMP(3,3)
      DOUBLE COMPLEX GRUN(601),CMATMUL,CID(3,3),VLQ(3,3)
      DOUBLE COMPLEX EVAQ(3),EVAU(3),EVAD(3),EVAL(3),EVAE(3)
      INTEGER I,J,K
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
      DO I=1,3
        DO J=1,3
          EQ(I,J)=CID(I,J)
          EU(I,J)=CID(I,J)
          ED(I,J)=CID(I,J)
          EL(I,J)=CID(I,J)
          EE(I,J)=CID(I,J)
          MQMASS(I,J)=CID(I,J)
          MUMASS(I,J)=CID(I,J)
          MDMASS(I,J)=CID(I,J)
          MLMASS(I,J)=CID(I,J)
          MEMASS(I,J)=CID(I,J)
          MQTMP(I,J)=CID(I,J)
          MUTMP(I,J)=CID(I,J)
          MDTMP(I,J)=CID(I,J)
          MLTMP(I,J)=CID(I,J)
          METMP(I,J)=CID(I,J)
          RQTMP(I,J)=CID(I,J)
          RUTMP(I,J)=CID(I,J)
          RDTMP(I,J)=CID(I,J)
          RLTMP(I,J)=CID(I,J)
          RETMP(I,J)=CID(I,J)
        END DO
        EVAQ(I)=(0.D0,0.D0)
        EVAU(I)=(0.D0,0.D0)
        EVAD(I)=(0.D0,0.D0)
        EVAL(I)=(0.D0,0.D0)
        EVAE(I)=(0.D0,0.D0)
      END DO
!
!Use the familiar rotation for the doublet.
!
      DO I=1,3
        DO J=1,3
          IF(SVLQ.EQ.1)THEN
            VLQ(I,J)=VLU(I,J)
          ELSE
            VLQ(I,J)=VLD(I,J)
          END IF
        END DO
      END DO
!
!Use OLD* for the number of active squarks since we still haven't
!taken a step at this new scale
!
!
!If there are two active squarks, I need to find the combined rotation
!first from the squark mass basis to the basis where the squarks are
!diagonal at the heaviest decoupling and then from this basis to
!the original current basis. If the second squark has also decoupled
!this combined rotation is saved elsewhere.
!
      IF(OLDNSQ.LT.3.OR.OLDNSU.LT.3.OR.OLDNSD.LT.3
     $                     .OR.OLDNSL.LT.3.OR.OLDNSE.LT.3)THEN
        DO I=1,3
          DO J=1,3
            IF(OLDNSQ.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RQSAV(2,I,K)*MQVE(K,J)
              END DO
              RQTMP(I,J)=SUM
            END IF
            IF(OLDNSU.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RUPSAV(2,I,K)*MUPVE(K,J)
              END DO
              RUTMP(I,J)=SUM
            END IF          
            IF(OLDNSD.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RDSAV(2,I,K)*MDVE(K,J)
              END DO
              RDTMP(I,J)=SUM
            END IF          
            IF(OLDNSL.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RLSAV(2,I,K)*MLVE(K,J)
              END DO
              RLTMP(I,J)=SUM
            END IF          
            IF(OLDNSE.EQ.2)THEN
              SUM=(0.D0,0.D0)
              DO K=1,3
                SUM=SUM+RESAV(2,I,K)*MEVE(K,J)
              END DO
              RETMP(I,J)=SUM
            END IF          
          END DO
        END DO
!
!Work out the eigenvectors in terms of the current basis
!states. If there has been no decoupling, then we can just use
!the current basis soft mass matrices as they are.
!
!Also in this loop, work out the mass matrix in the squark mass
!basis.
!
        DO I=1,3
          IF(OLDNSQ.LT.3)THEN
            IF(OLDNSQ.EQ.2)THEN
              DO J=1,3
                EQ(I,J)=RQTMP(I,J)
              END DO
              IF(I.LT.3)EVAQ(I)=MQVA(I)
            ELSE IF(OLDNSQ.LT.2)THEN
              DO J=1,3
                EQ(I,J)=RQSAV(1,I,J)
              END DO
              IF(OLDNSQ.EQ.1.AND.I.EQ.1)EVAQ(I)=MQVA(I)
              IF(OLDNSQ.EQ.0.AND.I.EQ.1)EVAQ(I)=MQSAV(I,4,I)
              IF(I.EQ.2)EVAQ(I)=MQSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAQ(I)=MQSAV(I,4,I)
          END IF
!
          IF(OLDNSU.LT.3)THEN
            IF(OLDNSU.EQ.2)THEN
              DO J=1,3
                EU(I,J)=RUTMP(I,J)
              END DO
              IF(I.LT.3)EVAU(I)=MUPVA(I)
            ELSE IF(OLDNSU.LT.2)THEN
              DO J=1,3
                EU(I,J)=RUPSAV(1,I,J)
              END DO
              IF(OLDNSU.EQ.1.AND.I.EQ.1)EVAU(I)=MUPVA(I)
              IF(OLDNSU.EQ.0.AND.I.EQ.1)EVAU(I)=MUPSAV(I,4,I)
              IF(I.EQ.2)EVAU(I)=MUPSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAU(I)=MUPSAV(I,4,I)
          END IF
!
          IF(OLDNSD.LT.3)THEN
            IF(OLDNSD.EQ.2)THEN
              DO J=1,3
                ED(I,J)=RDTMP(I,J)
              END DO
              IF(I.LT.3)EVAD(I)=MDVA(I)
            ELSE IF(OLDNSD.LT.2)THEN
              DO J=1,3
                ED(I,J)=RDSAV(1,I,J)
              END DO
              IF(OLDNSD.EQ.1.AND.I.EQ.1)EVAD(I)=MDVA(I)
              IF(OLDNSD.EQ.0.AND.I.EQ.1)EVAD(I)=MDSAV(I,4,I)
              IF(I.EQ.2)EVAD(I)=MDSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAD(I)=MDSAV(I,4,I)
          END IF
!
          IF(OLDNSL.LT.3)THEN
            IF(OLDNSL.EQ.2)THEN
              DO J=1,3
                EL(I,J)=RLTMP(I,J)
              END DO
              IF(I.LT.3)EVAL(I)=MLVA(I)
            ELSE IF(OLDNSL.LT.2)THEN
              DO J=1,3
                EL(I,J)=RLSAV(1,I,J)
              END DO
              IF(OLDNSL.EQ.1.AND.I.EQ.1)EVAL(I)=MLVA(I)
              IF(OLDNSL.EQ.0.AND.I.EQ.1)EVAL(I)=MLSAV(I,4,I)
              IF(I.EQ.2)EVAL(I)=MLSAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAL(I)=MLSAV(I,4,I)
          END IF
!
          IF(OLDNSE.LT.3)THEN
            IF(OLDNSE.EQ.2)THEN
              DO J=1,3
                EE(I,J)=RETMP(I,J)
              END DO
              IF(I.LT.3)EVAE(I)=MEVA(I)
            ELSE IF(OLDNSE.LT.2)THEN
              DO J=1,3
                EE(I,J)=RESAV(1,I,J)
              END DO
              IF(OLDNSE.EQ.1.AND.I.EQ.1)EVAE(I)=MEVA(I)
              IF(OLDNSE.EQ.0.AND.I.EQ.1)EVAE(I)=MESAV(I,4,I)
              IF(I.EQ.2)EVAE(I)=MESAV(I,4,I)
            END IF
            IF(I.EQ.3)EVAE(I)=MESAV(I,4,I)
          END IF
        END DO
!
!Fix the rotation matrices using ORTH as before so that 
!the rotation error is in the least damaging entries.
!
        IF(OLDNSQ.LT.3)CALL ORTH(EVAQ,EQ,OFFMAXQ,OFFMAXQVAL,VLQ,0)
        IF(OLDNSU.LT.3)CALL ORTH(EVAU,EU,OFFMAXU,OFFMAXUVAL,VRU,0)
        IF(OLDNSD.LT.3)CALL ORTH(EVAD,ED,OFFMAXD,OFFMAXDVAL,VRD,0)
        IF(OLDNSL.LT.3)CALL ORTH(EVAL,EL,OFFMAXL,OFFMAXLVAL,CID,0)
        IF(OLDNSE.LT.3)CALL ORTH(EVAE,EE,OFFMAXE,OFFMAXEVAL,CID,0)
!
        DO I=1,3
          MQMASS(I,I)=EVAQ(I)
          MUMASS(I,I)=EVAU(I)
          MDMASS(I,I)=EVAD(I)
          MLMASS(I,I)=EVAL(I)
          MEMASS(I,I)=EVAE(I)
        END DO
!
!I must now rotate from the squark mass basis to
!the original current basis and then to the quark mass basis.
!
        DO I=1,3
          DO J=1,3
            MQTMP(I,J)=CMATMUL(2,MQMASS,EQ,I,J)
            MUTMP(I,J)=CMATMUL(2,MUMASS,EU,I,J)
            MDTMP(I,J)=CMATMUL(2,MDMASS,ED,I,J)
            MLTMP(I,J)=CMATMUL(2,MLMASS,EL,I,J)
            METMP(I,J)=CMATMUL(2,MEMASS,EE,I,J)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            MQCURR(I,J)=CMATMUL(0,EQ,MQTMP,I,J)
            MUCURR(I,J)=CMATMUL(0,EU,MUTMP,I,J)
            MDCURR(I,J)=CMATMUL(0,ED,MDTMP,I,J)
            MLCURR(I,J)=CMATMUL(0,EL,MLTMP,I,J)
            MECURR(I,J)=CMATMUL(0,EE,METMP,I,J)
          END DO
        END DO
      END IF
!
      IF(OLDNSQ.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MQCURR(I,J)=GRUN(62+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSU.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MUCURR(I,J)=GRUN(80+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSD.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MDCURR(I,J)=GRUN(89+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSL.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MLCURR(I,J)=GRUN(71+(I-1)*3+J)
          END DO
        END DO
      END IF
      IF(OLDNSE.EQ.3)THEN
        DO I=1,3
          DO J=1,3
            MECURR(I,J)=GRUN(98+(I-1)*3+J)
          END DO
        END DO
      END IF
!
      DO I=1,3
        DO J=1,3
          MQTMP(I,J)=CMATMUL(0,MQCURR,VLQ,I,J)
          MUTMP(I,J)=CMATMUL(0,MUCURR,VRU,I,J)
          MDTMP(I,J)=CMATMUL(0,MDCURR,VRD,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MQQMASS(I,J)=CMATMUL(1,VLQ,MQTMP,I,J)
          MUQMASS(I,J)=CMATMUL(1,VRU,MUTMP,I,J)
          MDQMASS(I,J)=CMATMUL(1,VRD,MDTMP,I,J)
          MLQMASS(I,J)=MLCURR(I,J)
          MEQMASS(I,J)=MECURR(I,J)
        END DO
      END DO
!
   51 FORMAT(1X,A3,3X,D27.20,3x,D27.20)
   52 FORMAT(1X,I1,I1,3X,D27.20,3x,D27.20)
   53 FORMAT(1X,D27.20,3x,D27.20)
!
      RETURN
      END
!
      FUNCTION DMATMUL(DAG,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B
!         If DAG=1 compute dagger of A
!         If DAG=2 compute dagger of B
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3),B(3,3)
      DOUBLE PRECISION DMATMUL,DAGA(3,3),DAGB(3,3)
      INTEGER LOOP,DAG,I,J
!
      DMATMUL=0.D0
      IF(DAG.EQ.1)THEN
        CALL DAGGER(A,DAGA)
        DO LOOP=1,3
          DMATMUL=DMATMUL+DAGA(I,LOOP)*B(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.2)THEN
        CALL DAGGER(B,DAGB)
        DO LOOP=1,3
          DMATMUL=DMATMUL+A(I,LOOP)*DAGB(LOOP,J)
        END DO
      ELSE IF(DAG.EQ.0)THEN
        DO LOOP=1,3
          DMATMUL=DMATMUL+A(I,LOOP)*B(LOOP,J)
        END DO
      ELSE
        WRITE(*,*)'WRONG DAG IN DMATMUL',DAG
      END IF
!
      RETURN
      END
!
      FUNCTION MODSQ(A)
!
!Purpose: Calculates the mod squared of a value. This will only
!         be a necessary function if A is complex
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A,MODSQ
!
      MODSQ=A**2
!
      RETURN
      END
!
      SUBROUTINE ORTH(VA,VE,QBLARG,QBLARGVAL,QBROT,P)
!
!Purpose: To make sure that the eigenvectors are orthogonal in the
!         case that the eigenvectors are close to one another.
!         Resets one entry in one eigenvector. Chooses entry
!         to reset by finding the largest multiple in the
!         orthogonality sum.
!         QBLARG is the location of the largest off-diagonal entry
!                of M in quark basis.
!         QBLARGVAL is the value of the largest off-diagonal entry
!         QBROT is rotation from current to quark basis.
!
!The diagonalisation routine is not perfectly accurate. Just as the
!Vs come with intrinsic error in unitarity of 10^-16, the diagonalisation
!has an intrinsic deviation from unitarity, which is particulary
!acute if the eigenvectors are close, when the error can be
!as large as 10^-11. This is a problem when we are trying to calculate
!small numbers, so we must shift the error into the most convenient
!location. ORTH will use information about the eigenvectors, eigenvalues,
!and current->quark mass basis rotation to shift this error so that it
!affects the largest of the off-diagonal entries of the squark mass matrix.
!By shifting the error to affect the largest off-diagonal entry, we ensure
!that it has the smallest possible effect. For example, if we use mSUGRA
!and m_0=300, we find that the effect on the off-diagonal entries is
!of the order of 10^-6 (~10^-11*300^2). If the off-diagonal entry in the
!squared squark mass matrix is larger than 10^-6, the numerical fluctuations
!are swamped by the actual value. See Appendix A of PRD 79, (2009) 035004
!
      IMPLICIT NONE
!
      COMMON/ORTHWARNING/ORTHFLAG
      INTEGER ORTHFLAG
      SAVE/ORTHWARNING/
!
      DOUBLE COMPLEX VE(3,3),VA(3),QBROT(3,3),COMBROT(3,3)
      DOUBLE COMPLEX SUM,CMATMUL,QBLARGVAL,MAXCOMBOFF
      DOUBLE PRECISION MULT(3)
      INTEGER I,J,K,A,B,QBLARG(2),FIXE(2),FIXQB(2),P
!
      IF(P.EQ.1)WRITE(*,*)'BEFORE'
      DO I=1,3
        DO J=1,3
          COMBROT(I,J)=CMATMUL(1,VE,QBROT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          IF(P.EQ.1)WRITE(*,*)I,J,ABS(CMATMUL(1,COMBROT,COMBROT,I,J))
        END DO
      END DO
      IF(P.EQ.1)WRITE(*,*)
!
!Find two closest eigenvalues
!
      DO I=1,2
        IF(ABS(VA(1)-VA(2)).LT.ABS(VA(1)-VA(3)).AND.
     $     ABS(VA(1)-VA(2)).LT.ABS(VA(2)-VA(3)))THEN
          FIXE(1)=1
          FIXE(2)=2
        ELSE IF(ABS(VA(2)-VA(3)).LT.ABS(VA(1)-VA(3)))THEN
          FIXE(1)=2
          FIXE(2)=3
        ELSE
          FIXE(1)=1
          FIXE(2)=3
        END IF
      END DO
!
!Find the diagonal quark basis entry which corresponds to these 
!eigenvalues.
!The first index on combrot is the squark mass basis and the second
!is the quark mass basis.
!
      DO J=1,2
        IF(ABS(COMBROT(FIXE(J),1)).GT.ABS(COMBROT(FIXE(J),2)).AND.
     $     ABS(COMBROT(FIXE(J),1)).GT.ABS(COMBROT(FIXE(J),3)))THEN
          FIXQB(J)=1
        ELSE IF(ABS(COMBROT(FIXE(J),2)).GT.ABS(COMBROT(FIXE(J),3)))THEN
          FIXQB(J)=2
        ELSE
          FIXQB(J)=3
        END IF
      END DO
      IF(FIXQB(1).EQ.FIXQB(2))THEN
        WRITE(*,*)'ERROR IN FIXQB IN ORTH - EIGENVALUES ARE'
        WRITE(*,*)'PROBABLY VERY CLOSE. CHECK THE SQUARK'
        WRITE(*,*)'DIAGONALIZATION'
!
!Note: An error here is probably due combrot being far away from
!      the unit matrix. In this case, there is no clear separation
!      of entries in one basis to be identified with the eigenvalues.
!      This will probably not be a problem, since in this case it
!      is unlikely that the eigenvalues will be particularly close
!      and the fix in this subroutine will not be necessary.
!
        RETURN
!
      END IF
!
!Now fix the eigenvectors specified by FIXE to be orthogonal by
!changing one entry in the eigenvector. The eigenvector to alter
!corresponds to the contributor in the largest off-diagonal entry
!in the quark mass basis.
!
      IF((QBLARG(1)+QBLARG(2)).EQ.(FIXQB(1)+FIXQB(2)))THEN
!Don't fix anything if the error is already in the largest
!off-diagonal element
        A=0
        B=0
      ELSE IF(QBLARG(1).EQ.FIXQB(1).OR.QBLARG(2).EQ.FIXQB(1))THEN
!Fix the eigenvector corresponding to FIXQB(1)
        A=1 !A is entry of FIXE that is being fixed
        B=2
      ELSE
!Fix the eigenvector corresponding to FIXQB(2)
        A=2
        B=1
      END IF
!
      IF(A.NE.B)THEN
        DO K=1,3
          MULT(K)=ABS(CONJG(VE(K,FIXE(B)))*VE(K,FIXE(A)))
        END DO
        IF(MULT(1).GT.MULT(2).AND.MULT(1).GT.MULT(3))THEN
          VE(1,FIXE(A))=-(VE(2,FIXE(A))*CONJG(VE(2,FIXE(B)))
     $                   +VE(3,FIXE(A))*CONJG(VE(3,FIXE(B))))
     $                                          /CONJG(VE(1,FIXE(B)))
        ELSE IF(MULT(2).GT.MULT(1).AND.MULT(2).GT.MULT(3))THEN
          VE(2,FIXE(A))=-(VE(1,FIXE(A))*CONJG(VE(1,FIXE(B)))
     $                   +VE(3,FIXE(A))*CONJG(VE(3,FIXE(B))))
     $                                          /CONJG(VE(2,FIXE(B)))
        ELSE IF(MULT(3).GT.MULT(1).AND.MULT(3).GT.MULT(2))THEN
          VE(3,FIXE(A))=-(VE(1,FIXE(A))*CONJG(VE(1,FIXE(B)))
     $                   +VE(2,FIXE(A))*CONJG(VE(2,FIXE(B))))
     $                                          /CONJG(VE(3,FIXE(B)))
        ELSE IF(MULT(1).NE.MULT(2).AND.MULT(1).NE.MULT(3)
     $                            .AND.MULT(2).NE.MULT(3))THEN
          WRITE(*,*)'ERROR IN MULT OF ORTH'
        END IF 
      END IF
      IF(P.EQ.1)WRITE(*,*)'AFTER'
      DO I=1,3
        DO J=1,3
          COMBROT(I,J)=CMATMUL(1,VE,QBROT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          IF(P.EQ.1)WRITE(*,*)I,J,ABS(CMATMUL(1,COMBROT,COMBROT,I,J))
        END DO
      END DO
!
!Find largest off-diagonal value of COMBROT^dagger x COMBROT for
!for error checking purposes.
!
      MAXCOMBOFF=(0.D0,0.D0)
      DO I=1,3
        DO J=1,3
          IF(I.NE.J)THEN
            IF(ABS(MAXCOMBOFF).LT.ABS(CMATMUL(1,COMBROT,COMBROT,I,J)))
     $                         MAXCOMBOFF=CMATMUL(1,COMBROT,COMBROT,I,J)
          END IF
        END DO
      END DO
!
!If the smallest eigenvalue multiplied by the largest off-diagonal
!entry in COMBROT is greater than the largest off-diagonal entry of
!the squark mass matrix in chosen quark 'mass' basis then add 1 to
!ORTHFLAG. See the main routine for RGEFLAV.
!
      IF(ABS(VA(1)*MAXCOMBOFF).GT.ABS(QBLARGVAL))ORTHFLAG=ORTHFLAG+1
!
      RETURN
      END
!
      SUBROUTINE OUTCOUP(Q,UD)
!
!Purpose: To output each matrix to a different file, in either
!         the 'mass' or current basis
!
!         UD=1 if we are running up, =0 if we are running down.
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      DOUBLE COMPLEX GTMP(601)
      DOUBLE PRECISION Q,COSB,SINB
      INTEGER I,J,SWROT,UD
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!Switch which allows the user to turn off the rotation.
!
      SWROT=1
!
!Remember the old gs
!
      DO I=1,601
        GTMP(I)=G(I)
      END DO
!
!Perform the rotation to the mass basis if nec.
!
!      IF(SWROT.EQ.1)CALL ROTBACK(0)
      IF(SWROT.EQ.1)CALL ROTBACK(1) !This line is used to go back to the 
                                    !true mass basis.
!
!Insert the soft masses
!
      DO I=1,3
        DO J=1,3
          G(62+(I-1)*3+J)=MQQMASS(I,J)
          G(71+(I-1)*3+J)=MLQMASS(I,J)
          G(80+(I-1)*3+J)=MUQMASS(I,J)
          G(89+(I-1)*3+J)=MDQMASS(I,J)
          G(98+(I-1)*3+J)=MEQMASS(I,J)
        END DO
      END DO
!
!Make some adjustments to facilitate comparisons
!
      G(1)=G(1)
      DO I=4,12
        G(I+108)=G(I+108)/SINB
      END DO
      DO I=13,30
        G(I+108)=G(I+108)/COSB
      END DO
      G(287)=G(287)/SINB
      G(289)=G(289)/SINB
      G(288)=G(288)/COSB
      G(290)=G(290)/COSB
      IF(Q.LE.QNH-EPS)THEN
        DO I=1,27
          IF(I.LT.10)THEN
            G(I+323)=SINB*G(I+323)-COSB*CONJG(G(398))*G(I+293)
          ELSE
            G(I+323)=COSB*G(I+323)-SINB*CONJG(G(398))*G(I+293)
          END IF
        END DO
      END IF
      G(351)=G(351)+ABS(G(398))**2
      G(352)=G(352)+ABS(G(398))**2
      DO I=1,601
        IF(ABS(G(I)).LT.1.D-30)THEN
          G(I)=(1.D-30,1.D-30)
        END IF
      END DO
!
!Print each set of terms out - absolute squares first
!
      WRITE(68,56)Q,ABS(G(1)),DBLE(G(1))
      WRITE(67,56)Q,ABS(G(2)),DBLE(G(2))
      WRITE(66,56)Q,ABS(G(3)),DBLE(G(3))
      WRITE(62,56)Q,ABS(G(31)),DBLE(G(31))
      WRITE(52,55)Q,ABS(G(31)),ABS(G(32)),ABS(G(32))/ABS(G(31))
     $             ,ABS(G(321)),ABS(G(322)),ABS(G(322))/ABS(G(321))
     $             ,abs(g(2))**2/(abs(g(1))/dsqrt(3.d0/5.d0))**2
     $             ,4.D0*ABS(G(2))**2*ABS(G(108))*sinb*cosb/
     $                              (12.D0*ABS(G(32))*ABS(G(2))**2)
      WRITE(61,56)Q,ABS(G(32)),DBLE(G(32))
      WRITE(60,56)Q,ABS(G(33)),DBLE(G(33))
      WRITE(55,56)Q,ABS(G(110)),G(110)
      WRITE(54,56)Q,ABS(G(111)),G(111)
      IF((Q.GE.QNH-ABS(EPS).AND.UD.EQ.0).OR.
     $                  (Q.GE.QNH+ABS(EPS).AND.UD.EQ.1))THEN
        I=4
        WRITE(63,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=13
        WRITE(64,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=22
        WRITE(65,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=34
        WRITE(71,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=43
        WRITE(72,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=52
        WRITE(73,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=61
        WRITE(74,56)Q,ABS(G(I)),G(I)
        I=62
        WRITE(75,56)Q,ABS(G(I)),G(I)
        I=430
        WRITE(83,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=439
        WRITE(84,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=448
        WRITE(85,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=109
        WRITE(82,56)Q,ABS(G(I)),G(I)
      ELSE
        I=112
        WRITE(63,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=121
        WRITE(64,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=130
        WRITE(65,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=400
        WRITE(86,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=409
        WRITE(87,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=418
        WRITE(88,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $  ,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $  ,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
        I=427
        WRITE(89,56)Q,ABS(G(I)),G(I)
      END IF
      I=63
      WRITE(76,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=72
      WRITE(79,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=81
      WRITE(77,57)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
     $,abs(2.d0*(g(81)-g(85))/(g(81)+g(85)))
      IF((G(89)-G(85)).NE.0)WRITE(53,56)Q,ABS(G(86)/(G(89)-G(85)))
     $,G(86)/(G(89)-G(85))
      I=90
      WRITE(78,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=99
      WRITE(80,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=108
      WRITE(81,56)Q,ABS(G(I)),G(I)
      I=139
      WRITE(59,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=157
      WRITE(58,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=186
      WRITE(57,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=233
      WRITE(56,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=324
      WRITE(90,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=333
      WRITE(91,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=342
      WRITE(92,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=351
      WRITE(93,55)Q,ABS(G(I)),G(I)
      I=352
      WRITE(94,55)Q,ABS(G(I)),G(I)
      I=353
      WRITE(95,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=362
      WRITE(98,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=371
      WRITE(96,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=380
      WRITE(97,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=389
      WRITE(99,55)Q,ABS(G(I)),ABS(G(I+1)),ABS(G(I+2)),ABS(G(I+3))
     $,ABS(G(I+4)),ABS(G(I+5)),ABS(G(I+6)),ABS(G(I+7)),ABS(G(I+8))
     $,G(I),G(I+1),G(I+2),G(I+3),G(I+4),G(I+5),G(I+6),G(I+7),G(I+8)
      I=398
      WRITE(70,56)Q,ABS(G(I)),G(I)
      I=399
      WRITE(69,56)Q,ABS(G(I)),G(I)
!
!Return the gs to their original values
!
      DO I=1,601
        G(I)=GTMP(I)
      END DO
!
   51 FORMAT(1X,A3,3X,D27.20)
   55 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
   57 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X
     $,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
   56 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
!
      RETURN
      END
!
      SUBROUTINE REMSF(G)
!
!Purpose: To remove the entry in the matter sfermion mass matrix
!         which corresponds to the decoupled particle, and
!         return new entries for current basis and mass basis g's
!
      IMPLICIT NONE
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE COMPLEX G(601),CMATMUL
      DOUBLE COMPLEX MQ(3,3),MUP(3,3),MD(3,3),ML(3,3),ME(3,3)
      DOUBLE COMPLEX MQT(3,3),MUPT(3,3),MDT(3,3),MLT(3,3),MET(3,3)
      DOUBLE COMPLEX MQMASS(3,3),MUPMASS(3,3),MDMASS(3,3)
      DOUBLE COMPLEX MLMASS(3,3),MEMASS(3,3)
      INTEGER I,J,NSQ,NSU,NSD,NSL,NSE
!
      NSQ=0
      NSU=0
      NSD=0
      NSL=0
      NSE=0
      DO I=1,3
        NSQ=NSQ+THSQ(I)
        NSU=NSU+THSU(I)
        NSD=NSD+THSD(I)
        NSL=NSL+THSL(I)
        NSE=NSE+THSE(I)
      END DO
!
!THE NEXT LINE SHOULD BE UNNECESSARY
!It just checks to see that we are running down, not up.
!
      IF(OLDNSQ.LT.NSQ.OR.OLDNSU.LT.NSU.OR.OLDNSD.LT.NSD.OR.
     $   OLDNSL.LT.NSL.OR.OLDNSE.LT.NSE)RETURN
!
!Rotate the relevant matrices
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=G(62+(I-1)*3+J)
          ML(I,J)=G(71+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
          ME(I,J)=G(98+(I-1)*3+J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MQT(I,J)=CMATMUL(0,MQ,RQTOT,I,J)
          MLT(I,J)=CMATMUL(0,ML,RLTOT,I,J)
          MUPT(I,J)=CMATMUL(0,MUP,RUPTOT,I,J)
          MDT(I,J)=CMATMUL(0,MD,RDTOT,I,J)
          MET(I,J)=CMATMUL(0,ME,RETOT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MQMASS(I,J)=CMATMUL(1,RQTOT,MQT,I,J)
          MLMASS(I,J)=CMATMUL(1,RLTOT,MLT,I,J)
          MUPMASS(I,J)=CMATMUL(1,RUPTOT,MUPT,I,J)
          MDMASS(I,J)=CMATMUL(1,RDTOT,MDT,I,J)
          MEMASS(I,J)=CMATMUL(1,RETOT,MET,I,J)
        END DO
      END DO
!
!Remove the entries from the mass basis matrices
!
      DO I=1,3
        DO J=1,3
          IF(NSQ.LT.3.AND.OLDNSQ-NSQ.NE.0)THEN
            IF(I.EQ.J)THEN
              MQMASS(I,J)=MQMASS(I,J)*THSQ(I)
            ELSE
              MQMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSU.LT.3.AND.OLDNSU-NSU.NE.0)THEN
            IF(I.EQ.J)THEN
              MUPMASS(I,I)=MUPMASS(I,J)*THSU(I)
            ELSE
              MUPMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSD.LT.3.AND.OLDNSD-NSD.NE.0)THEN
            IF(I.EQ.J)THEN
              MDMASS(I,J)=MDMASS(I,J)*THSD(I)
            ELSE
              MDMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSL.LT.3.AND.OLDNSL-NSL.NE.0)THEN
            IF(I.EQ.J)THEN
              MLMASS(I,J)=MLMASS(I,J)*THSL(I)
            ELSE
              MLMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
          IF(NSE.LT.3.AND.OLDNSE-NSE.NE.0)THEN
            IF(I.EQ.J)THEN
              MEMASS(I,J)=MEMASS(I,J)*THSE(I)
            ELSE
              MEMASS(I,J)=(0.D0,0.D0)
            END IF
          END IF
        END DO
      END DO
!
!Finally, rotate these matrices back
!
      DO I=1,3
        DO J=1,3
          MQT(I,J)=CMATMUL(2,MQMASS,RQTOT,I,J)
          MLT(I,J)=CMATMUL(2,MLMASS,RLTOT,I,J)
          MUPT(I,J)=CMATMUL(2,MUPMASS,RUPTOT,I,J)
          MDT(I,J)=CMATMUL(2,MDMASS,RDTOT,I,J)
          MET(I,J)=CMATMUL(2,MEMASS,RETOT,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          G(62+(I-1)*3+J)=CMATMUL(0,RQTOT,MQT,I,J)
          G(71+(I-1)*3+J)=CMATMUL(0,RLTOT,MLT,I,J)
          G(80+(I-1)*3+J)=CMATMUL(0,RUPTOT,MUPT,I,J)
          G(89+(I-1)*3+J)=CMATMUL(0,RDTOT,MDT,I,J)
          G(98+(I-1)*3+J)=CMATMUL(0,RETOT,MET,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE RGEFLAV(RM0,RM12,RA0,RTANB,RSIGNMU,RMT)
!
!Purpose: To run the RGEs and derive values for all the parameters at
!         M_t/m_H
!
!Hard wired constants:
!        KM(3,3) - the mixing matrix
!        MWEAK(6) - the low scale masses
!
!Inputs - see inrge.dat:
!        High scale parameters
!        MHIGH - Value of high scale if not GUT scale
!        SUG - Whether to use mSUGRA
!        UNI - Whether to run up to unified gauge couplings 
!        ACC - Accuracy of the RGEs (ie. one or two loops)
!        COMP - Whether to compute the complex RGEs
!        PHASEMU - phase for \mu at the EWSB scale
!        SVLQ - Basis in which to output results
!        VLU,VRU,VLD,VRD - Quark rotation matrices
!        g(1-3),f_b,f_tau,lambda,vev (at M_z) and f_t(M_t)
!
!Some of the inputs are passed since RGEFLAV is used as a subroutine
!
!Requires:
!        EISPACK - cg.f and dependencies
!        RGEs - crkstp.f, crgeutils.f, csmrgedr.f, crge215.o, crge601.o
!               drkstp.f, drgeutils.f, dsmrgedr.f, drge215.o, drge601.o
!        INPUTS - rgeread.f,inrge.dat
!        OTHER - sugefffl.f,sqsix.f
!
!Output - Writes to a file (weakout.dat) the full complement of
!         couplings and masses at m_H. These are written
!         in the basis where one of either the up- or down-type
!         Yukawa matrices are diagonal, depending on user input.
!
!         If SQSIX is called, there is also output to the files sqm2u.dat
!         and sqm2d.dat which contian the (6x6) up-type and down-type
!         squark mass matrices and (when in the appropriate basis)
!         corresponding eigenvectors and eigenvalues.
!
!         In addition, there are commented lines in:
!            UPMZMHIGH
!            DOWNMHIGHMZ
!            UPMMHIGH2
!            DECRUN (located in sqsix.f)
!
!         that can be reinstated to print to separate files the running
!         of various parameters. The lines include OPEN and CLOSE
!         commands as well as calls to OUTCOUP, UOUTCOUP and DECOUTCOUP
!         NB: These commands require the presence of a subdirectory
!             called 'out', which must be created before running the
!             executable.
!
!HOW TO SEPARATE OUT THE GUT SCALE FROM THE FLAVOUR SCALE:
!
!   If at some later time it is required that the GUT scale can
!   be separated from the flavour scale, it may be easiest
!   to implement by slightly altering DOWNMSCOND. Simply
!   run all the way to the GUT scale, but at each run down,
!   check if the scale has passed the required flavour scale.
!   When the running passes this scale, call the new subroutine
!   that would contain the flavour BCs.
!
!
!NOTE: For the time being, RGEFLAV must be called with
!      equivalent mSUGRA inputs. This is not a necessity,
!      but at present mSUGRA is the default GUT scale input.
!
      IMPLICIT NONE
!
!     SUGRA parameters passed by ISAJET
!
      COMMON /RGEFNM/ FNRGE
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
!     All 601 RGE terms. Definitions are in *rge601.f
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
!     Used to tell RGEs that we want to run to
!     two loops
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
!     NU is the number of active up-like quarks and SMDR2LP tells
!     csmrgedr.f to run to two loops
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
!     Constants read in by RGEREAD from a file
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
!Warning flag from the accuracy subroutine: ORTH
!
      COMMON/ORTHWARNING/ORTHFLAG
      INTEGER ORTHFLAG
      SAVE/ORTHWARNING/
!
      INTEGER I,SWSQSIX,NSTEP,FINALIT,THCNG
      DOUBLE PRECISION QEND,PI
      REAL RM0,RM12,RA0,RTANB,RSIGNMU,RMT
      CHARACTER*128 FILENAME,STRADD,FNRGE
!
      PI=4.D0*DATAN(1.D0)
!
!Convert real inputs to double precision.
!I change some signs here since RGEFLAV uses book notation
!and ISAJET inputs are in MV notation.
!
      M0=DBLE(RM0)
      M12=DBLE(RM12)
      A0=-DBLE(RA0)
      TANB=DBLE(RTANB)
      SIGNMU=-DBLE(RSIGNMU)
      PHASEMU=PI*(1.D0-SIGNMU)/2.D0
      MT=RMT
!
      SWSQSIX=1
!
!May need to be increased to allow for more accurate convergence
!of iterations.
!
      FINALIT=16
      THCNG=0
!
!Next the rest of the inputs from INRGE.DAT
!AS PART OF RGEREAD, READ IN THE SWITCHES ACC AND COMP
!
      CALL RGEREAD
!
      SW2LP=ACC   !These two lines set the accuracy of the MSSM
      SMDR2LP=ACC !and SM running
!
      DO I=1,601
        G(I)=(0.D0,0.D0)
      END DO
!
!The light weak-scale quark masses are derived at Mz
!
      CALL MASS
!
!We are ready to run up.
!
      CALL UPMZMHIGH
!
!Now insert the high scale conditions read in earlier
!
      CALL HIGHIN(1)
!
!Run back down to the weak scale and iterate
!
      DO I=1,FINALIT
!
        ORTHFLAG=0
!
!May need to be increased to allow for more accurate running.
!
        NSTEP=INT(100.D0*1.6D0**I)
        WRITE(*,*)'ITERATION NUMBER ',I
!
!Sets when the iteration of the threshold locations ends. Could be
!set based upon a test of convergence. We continue iterating
!after fixing the thresholds since if the thresholds are still
!changing the Yukawas will not converge on a solution where
!they are diagonal at M_t.
!
        IF(I.GT.5)THEN
          NSTEP=INT(100.D0*1.6D0**5)
          IF(I.GT.10)THCNG=1
        END IF
        CALL DOWNMHIGHMZ(QEND,NSTEP,THCNG)
        CALL UPMZMHIGH2(NSTEP)
        CALL HIGHIN(0)
!
!If ORTHFLAG is too large, there is a high probability that there
!have been accuracy issues in the rotation between the squark
!mass basis and the output basis. ORTHFLAG obtains contributions
!to its value whenever ORTH is unlikely to have been able to shift
!the error in the diagonalisation into an entry that is much larger
!than the size of the error.
!Print a warning
!
        IF(ORTHFLAG.GT.INT(NSTEP/10))WRITE(*,*)'WARNING FROM ORTHFLAG'
!
      END DO
!
!One final run back down
!
      IF(SWSQSIX.EQ.0)THEN
        CALL DOWNMHIGHMZ(QEND,NSTEP,1)
        CALL ROTBACK(0)
!
!Write the results to a file
!
        FILENAME=STRADD(FNRGE,'.wkout')
        CALL GOUT601(G,QEND,3,FILENAME)
!
      ELSE
!
!If required, the next lines are the call to the decay calculation prog.
!
        IF(SVLQ.EQ.0)WRITE(*,*)'DOWN BASIS CHOSEN - NO STOP CALCULATION'
        CALL DOWNMHIGHMH(QEND,NSTEP)
        CALL ROTBACK(0)
!
!Write the results to a file
!
        FILENAME=STRADD(FNRGE,'.wkout')
        CALL GOUT601(G,QEND,3,FILENAME)
!
        CALL SQSIX(G,QEND)
      END IF
!
      RETURN
      END
!
      SUBROUTINE RGEREAD
!
!Purpose: To read in the inputs needed by RGEFLAV and populate the
!         common blocks.
!
!If the user decides to enter complex variables, they must be entered as
!"(x,y)" where z=x+iy.
!
      IMPLICIT NONE
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
!Contains the unitary matrices that are used in the rotation
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,FNINPUT,STRADD
!
!High scale SUSY inputs
!
      COMMON/HISUSY/CM1,CM2,CM3,CM1P,CM2P,CM3P
     $             ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $             ,CTQ,CTU,CTD,CTL,CTE
     $             ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $             ,CZU,CZD,CZE,CMHU,CMHD,CU,CD,CE
      DOUBLE COMPLEX CM1,CM2,CM3,CM1P,CM2P,CM3P
     $              ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $              ,CTQ(3,3),CTU(3,3),CTD(3,3),CTL(3,3),CTE(3,3)
     $              ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $              ,CZU(3,3),CZD(3,3),CZE(3,3),CMHU,CMHD
      INTEGER CU,CD,CE
      SAVE/HISUSY/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      DOUBLE COMPLEX COMBM1,COMBM2,COMBM3,DAGKM(3,3),CMATMUL
      DOUBLE PRECISION ALP,BET,GAM,DEL,ID(3,3)
      DOUBLE PRECISION DM1,DM2,DM3
     $                ,DMQ0,DMU0,DMD0,DML0,DME0,DRU,DRD,DRE,DSU,DSD,DSE
     $                ,DTQ(3,3),DTU(3,3),DTD(3,3),DTL(3,3),DTE(3,3)
     $                ,DAU0,DAD0,DAE0,DWU,DWD,DWE,DXU,DXD,DXE
     $                ,DZU(3,3),DZD(3,3),DZE(3,3),DMHU,DMHD
      INTEGER I,J,VRID,VLRIN,VLKM,TMP,ERRNO,ERRNOMR,HERMBAD
      CHARACTER COMPTEST
!
      DATA ID(1,1)/1.D0/,ID(1,2)/0.D0/,ID(1,3)/0.D0/
      DATA ID(2,1)/0.D0/,ID(2,2)/1.D0/,ID(2,3)/0.D0/
      DATA ID(3,1)/0.D0/,ID(3,2)/0.D0/,ID(3,3)/1.D0/
!
      ERRNO=0
      ERRNOMR=0
!
      DO I=1,3
        DO J=1,3
          DTQ(I,J)=0.D0
          DTU(I,J)=0.D0
          DTD(I,J)=0.D0
          DTL(I,J)=0.D0
          DTE(I,J)=0.D0
          DZU(I,J)=0.D0
          DZD(I,J)=0.D0
          DZE(I,J)=0.D0
        END DO
      END DO
!
!M_Z inputs. Start with ALPHAEM and SIN2THW from the PDG
!Pole masses from PDG rather than running masses.
!Updated from PDG 2008. For ALPHAEM, see Sec.10: EW Model and
!Constraints on New Physics
!
      ALPHAEM=1.D0/127.925D0
      ALPHASMSB=.1176D0
      XWMSB=.23119D0
      MZ=91.1876D0!Used in later subroutines, could be passed
                  !from an earlier programme
      MW=80.403D0
!
!Now read from the input file
!
      FNINPUT=STRADD(FNRGE,'.rgein')
      OPEN(12,FILE=FNINPUT,STATUS='UNKNOWN')
      READ(12,*)
      READ(12,*)ACC
      IF(ACC.EQ.0)WRITE(*,*)'USING ONE LOOP RGES'
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)COMP
      IF(COMP.EQ.0)WRITE(*,*)'USING REAL KM AND RUNNING'
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)TMP
      IF(TMP.EQ.0.OR.COMP.EQ.0)THEN
        DO I=1,3
          READ(12,*)
        END DO
      ELSE
        DO I=1,2
          READ(12,*)
        END DO
        READ(12,*)PHASEMU
      END IF
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)ISADEC
!
!Now we know COMP, we can call the KM routine
!
      CALL KMIN
      CALL CDAGGER(KM,DAGKM)
!
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)SUG
!
      IF(SUG.EQ.1)THEN
        UNI=1
        DO I=1,169
          READ(12,*)
        END DO
      ELSE
!
!Set all the inputs to follow to their mSUGRA defaults.
!It is only if the user chooses non-mSURA conditions that
!these defaults won't be used.
!
        UNI=1
        CM1=DCMPLX(M12,0.D0)
        CM2=DCMPLX(M12,0.D0)
        CM3=DCMPLX(M12,0.D0)
        CM1P=(0.D0,0.D0)
        CM2P=(0.D0,0.D0)
        CM3P=(0.D0,0.D0)
        CMHU=DCMPLX(M0)
        CMHD=DCMPLX(M0)
        CMQ0=DCMPLX(M0)
        CMU0=DCMPLX(M0)
        CMD0=DCMPLX(M0)
        CML0=DCMPLX(M0)
        CME0=DCMPLX(M0)
        CU=1
        CD=1
        CE=1
        CRU=(0.D0,0.D0)
        CRD=(0.D0,0.D0)
        CRE=(0.D0,0.D0)
        CSU=(0.D0,0.D0)
        CSD=(0.D0,0.D0)
        CSE=(0.D0,0.D0)
        CAU0=DCMPLX(A0)
        CAD0=DCMPLX(A0)
        CAE0=DCMPLX(A0)
        CWU=(0.D0,0.D0)
        CWD=(0.D0,0.D0)
        CWE=(0.D0,0.D0)
        CXU=(0.D0,0.D0)
        CXD=(0.D0,0.D0)
        CXE=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            CTQ(I,J)=(0.D0,0.D0)
            CTU(I,J)=(0.D0,0.D0)
            CTD(I,J)=(0.D0,0.D0)
            CTL(I,J)=(0.D0,0.D0)
            CTE(I,J)=(0.D0,0.D0)
            CZU(I,J)=(0.D0,0.D0)
            CZD(I,J)=(0.D0,0.D0)
            CZE(I,J)=(0.D0,0.D0)
          END DO
        END DO
!
        DO I=1,4
          READ(12,*)
        END DO
        READ(12,*)TMP
        IF(TMP.EQ.0)THEN
          UNI=0
          DO I=1,2
            READ(12,*)
          END DO
          READ(12,*)MHIGH
        ELSE
          DO I=1,3
            READ(12,*)
          END DO
        END IF
!
!First the complex read statements. Each variable can either be the
!same as the mSUGRA values or a user defined value
!
        IF(COMP.EQ.1)THEN
          DO I=1,2
            READ(12,*)
          END DO
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
             READ(12,*)
            END DO
!
!The following three lines check if the user has input the numbers
!in the correct form. If not, there is an error message and the
!programme stops
!
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=1
              GOTO 50
            END IF
            BACKSPACE 12 !NOW GO BACK TO READ THE LINE AGAIN
            READ(12,*)COMBM1,COMBM2,COMBM3
            CM1=DCMPLX(DBLE(COMBM1),0.D0)
            CM2=DCMPLX(DBLE(COMBM2),0.D0)
            CM3=DCMPLX(DBLE(COMBM3),0.D0)
            CM1P=DCMPLX(DIMAG(COMBM1),0.D0)
            CM2P=DCMPLX(DIMAG(COMBM2),0.D0)
            CM3P=DCMPLX(DIMAG(COMBM3),0.D0)
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=2
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CMHU,CMHD
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=3
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMQ0
            CMQ0=DCMPLX(DMQ0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=3
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTQ(I,1),CTQ(I,2),CTQ(I,3)
            END DO
            CALL HERMTEST(CTQ,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_Q IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMU0
            CMU0=DCMPLX(DMU0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CU
            IF(CU.NE.1.AND.CU.NE.0)THEN
              WRITE(*,*)'ERROR READING CU'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRU
            CRU=DCMPLX(DRU)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSU
            CSU=DCMPLX(DSU)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=4
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTU(I,1),CTU(I,2),CTU(I,3)
            END DO
            CALL HERMTEST(CTU,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_U IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMD0
            CMD0=DCMPLX(DMD0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CD
            IF(CD.NE.1.AND.CD.NE.0)THEN
              WRITE(*,*)'ERROR READING CD'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRD
            CRD=DCMPLX(DRD)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSD
            CSD=DCMPLX(DSD)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=5
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTD(I,1),CTD(I,2),CTD(I,3)
            END DO
            CALL HERMTEST(CTD,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_D IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=6
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DML0
            CML0=DCMPLX(DML0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=6
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTL(I,1),CTL(I,2),CTL(I,3)
            END DO
            CALL HERMTEST(CTL,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_L IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DME0
            CME0=DCMPLX(DME0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CE
            IF(CE.NE.1.AND.CE.NE.0)THEN
              WRITE(*,*)'ERROR READING CE'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRE
            CRE=DCMPLX(DRE)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSE
            CSE=DCMPLX(DSE)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=7
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTE(I,1),CTE(I,2),CTE(I,3)
            END DO
            CALL HERMTEST(CTE,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_E IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=8
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CAU0
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=8
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CWU
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=8
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CXU
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=8
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CZU(I,1),CZU(I,2),CZU(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=9
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CAD0
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=9
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CWD
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=9
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CXD
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=9
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CZD(I,1),CZD(I,2),CZD(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=10
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CAE0
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=10
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CWE
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=10
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CXE
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=10
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CZE(I,1),CZE(I,2),CZE(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
!Now the real read statements. All variables are converted to complex
!at the end
!
        ELSE
          DO I=1,2
            READ(12,*)
          END DO
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
             READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DM1,DM2,DM3
            CM1=DCMPLX(DM1,0.D0)
            CM2=DCMPLX(DM2,0.D0)
            CM3=DCMPLX(DM3,0.D0)
            CM1P=(0.D0,0.D0)
            CM2P=(0.D0,0.D0)
            CM3P=(0.D0,0.D0)
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMHU,DMHD
            CMHU=DCMPLX(DMHU,0.D0)
            CMHD=DCMPLX(DMHD,0.D0)
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=3
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMQ0
            CMQ0=DCMPLX(DMQ0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=3
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTQ(I,1),DTQ(I,2),DTQ(I,3)
            END DO
            CALL SYMMTEST(DTQ,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_Q IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMU0
            CMU0=DCMPLX(DMU0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CU
            IF(CU.NE.1.AND.CU.NE.0)THEN
              WRITE(*,*)'ERROR READING CU'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRU
            CRU=DCMPLX(DRU)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSU
            CSU=DCMPLX(DSU)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=4
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTU(I,1),DTU(I,2),DTU(I,3)
            END DO
            CALL SYMMTEST(DTU,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_U IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMD0
            CMD0=DCMPLX(DMD0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CD
            IF(CD.NE.1.AND.CD.NE.0)THEN
              WRITE(*,*)'ERROR READING CD'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRD
            CRD=DCMPLX(DRD)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSD
            CSD=DCMPLX(DSD)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=5
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTD(I,1),DTD(I,2),DTD(I,3)
            END DO
            CALL SYMMTEST(DTD,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_D IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=6
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DML0
            CML0=DCMPLX(DML0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=6
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTL(I,1),DTL(I,2),DTL(I,3)
            END DO
            CALL SYMMTEST(DTL,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_L IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DME0
            CME0=DCMPLX(DME0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CE
            IF(CE.NE.1.AND.CE.NE.0)THEN
              WRITE(*,*)'ERROR READING CE'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRE
            CRE=DCMPLX(DRE)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSE
            CSE=DCMPLX(DSE)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=7
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTE(I,1),DTE(I,2),DTE(I,3)
            END DO
            CALL SYMMTEST(DTE,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_E IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=8
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DAU0
            CAU0=DCMPLX(DAU0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=8
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DWU
            CWU=DCMPLX(DWU)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=8
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DXU
            CXU=DCMPLX(DXU)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=8
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DZU(I,1),DZU(I,2),DZU(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=9
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DAD0
            CAD0=DCMPLX(DAD0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=9
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DWD
            CWD=DCMPLX(DWD)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=9
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DXD
            CXD=DCMPLX(DXD)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=9
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DZD(I,1),DZD(I,2),DZD(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=10
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DAE0
            CAE0=DCMPLX(DAE0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=10
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DWE
            CWE=DCMPLX(DWE)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=10
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DXE
            CXE=DCMPLX(DXE)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=10
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DZE(I,1),DZE(I,2),DZE(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
          DO I=1,3
            DO J=1,3
              CTQ(I,J)=DCMPLX(DTQ(I,J))
              CTU(I,J)=DCMPLX(DTU(I,J))
              CTD(I,J)=DCMPLX(DTD(I,J))
              CTL(I,J)=DCMPLX(DTL(I,J))
              CTE(I,J)=DCMPLX(DTE(I,J))
              CZU(I,J)=DCMPLX(DZU(I,J))
              CZD(I,J)=DCMPLX(DZD(I,J))
              CZE(I,J)=DCMPLX(DZE(I,J))
            END DO
          END DO
        END IF
      END IF
!
!Now deal with the rotation matrices
!
      READ(12,*)SVLQ !SVLQ=1 for up quark diagonal basis
      IF(SVLQ.NE.0.AND.SVLQ.NE.1)THEN
        WRITE(*,*)'ERROR: UNKNOWN CHOICE FOR QUARK DOUBLET ROTATION'
        WRITE(*,*)'OUTPUT WILL BE IN UP QUARK DIAGONAL BASIS'
        SVLQ=1
      END IF
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)VLRIN !Does user input Vs?
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)VLKM !Is V^U_L=KM or KM^dagger?
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)VRID !Are V_Rs identity?
!
!Deal with user input Vs
!
      IF(VLRIN.EQ.1)THEN
        DO I=1,5
          READ(12,*)
        END DO
        READ(12,*)ALP,BET,GAM,DEL
        IF(COMP.EQ.0)DEL=0.D0 !This line makes sure that the V''s have
                              !no phase if the user chose real KM
        CALL VGEN(ALP,BET,GAM,DEL,VLU)
        DO I=1,2
          READ(12,*)
        END DO
        READ(12,*)ALP,BET,GAM,DEL
        IF(COMP.EQ.0)DEL=0.D0
        CALL VGEN(ALP,BET,GAM,DEL,VRU)
        DO I=1,2
          READ(12,*)
        END DO
        READ(12,*)ALP,BET,GAM,DEL
        IF(COMP.EQ.0)DEL=0.D0
        CALL VGEN(ALP,BET,GAM,DEL,VRD)
        DO I=1,3
          DO J=1,3
            VLD(I,J)=CMATMUL(0,VLU,KM,I,J)
          END DO
        END DO
!
!Set the Vs if they are not input by the user
!
      ELSE
!
!For the Left-handed matrices, just V^U_L is an input.
!
        IF(VLKM.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              VLU(I,J)=DCMPLX(ID(I,J))
            END DO
          END DO
        ELSE
          DO I=1,3
            DO J=1,3
              VLU(I,J)=KM(I,J)
            END DO
          END DO
        END IF
!
!V^D_L is set so that (V^U_L)^dagger * V^D_L = KM
!
        DO I=1,3
          DO J=1,3
            VLD(I,J)=CMATMUL(0,VLU,KM,I,J)
          END DO
        END DO
!
!The V_Rs are arbitrary. If VRID=0 then the KM is used to set them
!
        IF(VRID.EQ.1)THEN
          DO I=1,3
             DO J=1,3
              VRU(I,J)=DCMPLX(ID(I,J))
              VRD(I,J)=DCMPLX(ID(I,J))
            END DO
          END DO
        ELSE
          DO I=1,3
            DO J=1,3
              VRU(I,J)=DAGKM(I,J)
              VRD(I,J)=KM(I,J)
            END DO
          END DO
        END IF
      END IF
!
      RETURN
!
  50  CONTINUE
      WRITE(*,*)'REAL NUMBERS READ WHERE COMPLEX NUMBERS EXPECTED'
      GOTO 52
  51  CONTINUE
      WRITE(*,*)'COMPLEX NUMBERS READ WHERE REAL NUMBERS EXPECTED'
  52  CONTINUE
      IF(ERRNO.EQ.1)THEN
        WRITE(*,*)'ERROR WAS IN M1,M2,M3'
      ELSE IF(ERRNO.EQ.2)THEN
        WRITE(*,*)'ERROR WAS IN M_{H_U}, M_{H_D}'
      ELSE IF(ERRNO.EQ.3)THEN
        WRITE(*,*)'ERROR WAS IN M_Q'
      ELSE IF(ERRNO.EQ.4)THEN
        WRITE(*,*)'ERROR WAS IN M_U'
      ELSE IF(ERRNO.EQ.5)THEN
        WRITE(*,*)'ERROR WAS IN M_D'
      ELSE IF(ERRNO.EQ.6)THEN
        WRITE(*,*)'ERROR WAS IN M_L'
      ELSE IF(ERRNO.EQ.7)THEN
        WRITE(*,*)'ERROR WAS IN M_E'
      ELSE IF(ERRNO.EQ.8)THEN
        WRITE(*,*)'ERROR WAS IN A_U'
      ELSE IF(ERRNO.EQ.9)THEN
        WRITE(*,*)'ERROR WAS IN A_D'
      ELSE IF(ERRNO.EQ.10)THEN
        WRITE(*,*)'ERROR WAS IN A_E'
      END IF
      IF(ERRNOMR.EQ.1)THEN
        WRITE(*,*)'M_0 WAS COMPLEX'
      ELSE IF(ERRNOMR.EQ.2)THEN
        WRITE(*,*)'R WAS COMPLEX'
      ELSE IF(ERRNOMR.EQ.3)THEN
        WRITE(*,*)'S WAS COMPLEX'
      ELSE IF(ERRNOMR.EQ.4)THEN
        WRITE(*,*)'T WAS OF THE WRONG KIND'
      END IF
      STOP99
!
      END
!
      SUBROUTINE ROTATE215(G215)
!
!Purpose: Rotate all the matrices from the quark mass
!         basis to the current basis using the 
!         previously derived V's
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),LUP(3,3),LDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3),LUDUM(3,3),LDDUM(3,3)
!
      DOUBLE COMPLEX VLUT(3,3),VRUT(3,3),VLDT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL,G215(215)
      INTEGER I,J
!
      DO I=1,3
        DO J=1,3
          VLUT(I,J)=VLU(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDT(I,J)=VLD(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO   
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=G215(3+(I-1)*3+J)
          FDP(I,J)=G215(12+(I-1)*3+J)
          LUP(I,J)=G215(33+(I-1)*3+J)
          LDP(I,J)=G215(42+(I-1)*3+J)
        END DO
      END DO
!
!Now rotate the matrices back
!
      DO I=1,3 
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(0,LUP,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(0,LDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=CMATMUL(0,VLUT,FUDUM,I,J)
          LUP(I,J)=CMATMUL(0,VLUT,LUDUM,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G215(3+(I-1)*3+J)=CMATMUL(1,VLUT,FUDUM,I,J)
          G215(12+(I-1)*3+J)=CMATMUL(1,VLDT,FDDUM,I,J)
          G215(33+(I-1)*3+J)=CMATMUL(1,VLUT,LUDUM,I,J)
          G215(42+(I-1)*3+J)=CMATMUL(1,VLDT,LDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTATE(SMASS)
!
!Purpose: Rotate the parameters to the original current
!         basis
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),MQP(3,3),MUP(3,3),MDP(3,3)
      DOUBLE COMPLEX AUP(3,3),ADP(3,3),LUP(3,3),LDP(3,3)
      DOUBLE COMPLEX GTPQP(3,3),GTPUP(3,3),GTPDP(3,3),GTQP(3,3)
      DOUBLE COMPLEX GTSQP(3,3),GTSUP(3,3),GTSDP(3,3)
      DOUBLE COMPLEX FTUQP(3,3),FTDQP(3,3),FTUUP(3,3),FTDDP(3,3)
      DOUBLE COMPLEX FMUP(3,3),FMDP(3,3),MQMP(3,3),MUMP(3,3),MDMP(3,3)
      DOUBLE COMPLEX AUMP(3,3),ADMP(3,3),TRIUP(3,3),TRIDP(3,3)
      DOUBLE COMPLEX MTSFUP(3,3),MTSFDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX MQDUM(3,3),MUDUM(3,3),MDDUM(3,3)
      DOUBLE COMPLEX AUDUM(3,3),ADDUM(3,3),LUDUM(3,3),LDDUM(3,3)
      DOUBLE COMPLEX GTPQDUM(3,3),GTPUDUM(3,3),GTPDDUM(3,3),GTQDUM(3,3)
      DOUBLE COMPLEX GTSQDUM(3,3),GTSUDUM(3,3),GTSDDUM(3,3)
      DOUBLE COMPLEX FTUQDUM(3,3),FTDQDUM(3,3),FTUUDUM(3,3),FTDDDUM(3,3)
      DOUBLE COMPLEX FMUDUM(3,3),FMDDUM(3,3)
      DOUBLE COMPLEX MQMDUM(3,3),MUMDUM(3,3),MDMDUM(3,3)
      DOUBLE COMPLEX AUMDUM(3,3),ADMDUM(3,3),TRIUDUM(3,3),TRIDDUM(3,3)
      DOUBLE COMPLEX MTSFUDUM(3,3),MTSFDDUM(3,3)
!
      DOUBLE COMPLEX VLUQ(3,3),VLDQ(3,3)
      DOUBLE COMPLEX VLUQT(3,3),VRUT(3,3),VLDQT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL      
      INTEGER I,J,SMASS
!
!Set the SU(2) doublet rotation
!If SMASS=1 we rotate from the quark mass basis
!Otherwise we rotate between different current bases.
!
      DO I=1,3
        DO J=1,3
          IF(SMASS.EQ.1)THEN
            VLUQ(I,J)=VLU(I,J)
            VLDQ(I,J)=VLD(I,J)
          ELSE
            IF(SVLQ.EQ.1)THEN
              VLUQ(I,J)=VLU(I,J)
              VLDQ(I,J)=VLU(I,J)
            ELSE
              VLUQ(I,J)=VLD(I,J)
              VLDQ(I,J)=VLD(I,J)
            END IF
          END IF
        END DO
      END DO
!
!Convert the G's to matrices
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=G(3+(I-1)*3+J)
          FDP(I,J)=G(12+(I-1)*3+J)
          MQP(I,J)=G(62+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MDP(I,J)=G(89+(I-1)*3+J)
          AUP(I,J)=G(33+(I-1)*3+J)
          ADP(I,J)=G(42+(I-1)*3+J)
!
          LUP(I,J)=G(111+(I-1)*3+J)
          LDP(I,J)=G(120+(I-1)*3+J)
!
          GTPQP(I,J)=G(138+(I-1)*3+J)
          GTPUP(I,J)=G(156+(I-1)*3+J)
          GTPDP(I,J)=G(165+(I-1)*3+J)
          GTQP(I,J)=G(185+(I-1)*3+J)
          GTSQP(I,J)=G(205+(I-1)*3+J)
          GTSUP(I,J)=G(214+(I-1)*3+J)
          GTSDP(I,J)=G(223+(I-1)*3+J)
          FTUQP(I,J)=G(232+(I-1)*3+J)
          FTDQP(I,J)=G(241+(I-1)*3+J)
          FTUUP(I,J)=G(259+(I-1)*3+J)
          FTDDP(I,J)=G(268+(I-1)*3+J)
!
          FMUP(I,J)=G(293+(I-1)*3+J)
          FMDP(I,J)=G(302+(I-1)*3+J)
          MQMP(I,J)=G(352+(I-1)*3+J)
          MUMP(I,J)=G(370+(I-1)*3+J)
          MDMP(I,J)=G(379+(I-1)*3+J)
          AUMP(I,J)=G(323+(I-1)*3+J)
          ADMP(I,J)=G(332+(I-1)*3+J)
!
          TRIUP(I,J)=G(399+(I-1)*3+J)
          TRIDP(I,J)=G(408+(I-1)*3+J)
          MTSFUP(I,J)=G(429+(I-1)*3+J)
          MTSFDP(I,J)=G(438+(I-1)*3+J)
        END DO
      END DO
!
!Calculate the transposes
!
      DO I=1,3
        DO J=1,3
          VLUQT(I,J)=VLUQ(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDQT(I,J)=VLDQ(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO
!
!Now rotate the matrices
!
      DO I=1,3
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
          MQDUM(I,J)=CMATMUL(2,MQP,VLUQ,I,J)
          MUDUM(I,J)=CMATMUL(2,MUP,VRU,I,J)
          MDDUM(I,J)=CMATMUL(2,MDP,VRD,I,J)
          AUDUM(I,J)=CMATMUL(0,AUP,VRUT,I,J)
          ADDUM(I,J)=CMATMUL(0,ADP,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(0,LUP,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(0,LDP,VRDT,I,J)
          GTPQDUM(I,J)=CMATMUL(2,GTPQP,VLUQ,I,J)
          GTPUDUM(I,J)=CMATMUL(2,GTPUP,VRU,I,J)
          GTPDDUM(I,J)=CMATMUL(2,GTPDP,VRD,I,J)
          GTQDUM(I,J)=CMATMUL(2,GTQP,VLUQ,I,J)
          GTSQDUM(I,J)=CMATMUL(2,GTSQP,VLUQ,I,J)
          GTSUDUM(I,J)=CMATMUL(2,GTSUP,VRU,I,J)
          GTSDDUM(I,J)=CMATMUL(2,GTSDP,VRD,I,J)
          FTUQDUM(I,J)=CMATMUL(0,FTUQP,VRUT,I,J)
          FTDQDUM(I,J)=CMATMUL(0,FTDQP,VRDT,I,J)
          FTUUDUM(I,J)=CMATMUL(0,FTUUP,VRUT,I,J)
          FTDDDUM(I,J)=CMATMUL(0,FTDDP,VRDT,I,J)
          FMUDUM(I,J)=CMATMUL(0,FMUP,VRUT,I,J)
          FMDDUM(I,J)=CMATMUL(0,FMDP,VRDT,I,J)
          MQMDUM(I,J)=CMATMUL(2,MQMP,VLUQ,I,J)
          MUMDUM(I,J)=CMATMUL(2,MUMP,VRU,I,J)
          MDMDUM(I,J)=CMATMUL(2,MDMP,VRD,I,J)
          AUMDUM(I,J)=CMATMUL(0,AUMP,VRUT,I,J)
          ADMDUM(I,J)=CMATMUL(0,ADMP,VRDT,I,J)
          TRIUDUM(I,J)=CMATMUL(0,TRIUP,VRUT,I,J)
          TRIDDUM(I,J)=CMATMUL(0,TRIDP,VRDT,I,J)
          MTSFUDUM(I,J)=CMATMUL(0,MTSFUP,VRUT,I,J)
          MTSFDDUM(I,J)=CMATMUL(0,MTSFDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G(3+(I-1)*3+J)=CMATMUL(1,VLUQT,FUDUM,I,J)
          G(12+(I-1)*3+J)=CMATMUL(1,VLDQT,FDDUM,I,J)
          G(62+(I-1)*3+J)=CMATMUL(0,VLUQ,MQDUM,I,J)
          G(80+(I-1)*3+J)=CMATMUL(0,VRU,MUDUM,I,J)
          G(89+(I-1)*3+J)=CMATMUL(0,VRD,MDDUM,I,J)
          G(33+(I-1)*3+J)=CMATMUL(1,VLUQT,AUDUM,I,J)
          G(42+(I-1)*3+J)=CMATMUL(1,VLDQT,ADDUM,I,J)
!
          G(111+(I-1)*3+J)=CMATMUL(1,VLUQT,LUDUM,I,J)
          G(120+(I-1)*3+J)=CMATMUL(1,VLDQT,LDDUM,I,J)
          G(138+(I-1)*3+J)=CMATMUL(0,VLUQ,GTPQDUM,I,J)
          G(156+(I-1)*3+J)=CMATMUL(0,VRU,GTPUDUM,I,J)
          G(165+(I-1)*3+J)=CMATMUL(0,VRD,GTPDDUM,I,J)
          G(185+(I-1)*3+J)=CMATMUL(0,VLUQ,GTQDUM,I,J)
          G(205+(I-1)*3+J)=CMATMUL(0,VLUQ,GTSQDUM,I,J)
          G(214+(I-1)*3+J)=CMATMUL(0,VRU,GTSUDUM,I,J)
          G(223+(I-1)*3+J)=CMATMUL(0,VRD,GTSDDUM,I,J)
          G(232+(I-1)*3+J)=CMATMUL(1,VLUQT,FTUQDUM,I,J)
          G(241+(I-1)*3+J)=CMATMUL(1,VLDQT,FTDQDUM,I,J)
          G(259+(I-1)*3+J)=CMATMUL(1,VLUQT,FTUUDUM,I,J)
          G(268+(I-1)*3+J)=CMATMUL(1,VLDQT,FTDDDUM,I,J)
!
          G(293+(I-1)*3+J)=CMATMUL(1,VLUQT,FMUDUM,I,J)
          G(302+(I-1)*3+J)=CMATMUL(1,VLDQT,FMDDUM,I,J)
          G(352+(I-1)*3+J)=CMATMUL(0,VLUQ,MQMDUM,I,J)
          G(370+(I-1)*3+J)=CMATMUL(0,VRU,MUMDUM,I,J)
          G(379+(I-1)*3+J)=CMATMUL(0,VRD,MDMDUM,I,J)
          G(323+(I-1)*3+J)=CMATMUL(1,VLUQT,AUMDUM,I,J)
          G(332+(I-1)*3+J)=CMATMUL(1,VLDQT,ADMDUM,I,J)
!
          G(399+(I-1)*3+J)=CMATMUL(1,VLUQT,TRIUDUM,I,J)
          G(408+(I-1)*3+J)=CMATMUL(1,VLDQT,TRIDDUM,I,J)
          G(429+(I-1)*3+J)=CMATMUL(1,VLUQT,MTSFUDUM,I,J)
          G(438+(I-1)*3+J)=CMATMUL(1,VLDQT,MTSFDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTATESM(GROT)
!
!Purpose: Rotate the Yukawa matrices from the basis in which the
!         chosen Yukawa matrix is diagonal to the general current
!         basis.
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FUP(3,3),FDP(3,3)
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX GROT(32)
!
      DOUBLE COMPLEX VLUT(3,3),VRUT(3,3),VLDT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J
!
!Convert the G's to matrices
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=GROT(3+(I-1)*3+J)
          FDP(I,J)=GROT(12+(I-1)*3+J)
        END DO
      END DO
!
!Calculate the transposes
!
      DO I=1,3
        DO J=1,3
          VLUT(I,J)=VLU(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDT(I,J)=VLD(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO
!
!Now rotate the matrices
!
      DO I=1,3
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GROT(3+(I-1)*3+J)=CMATMUL(1,VLUT,FUDUM,I,J)
          GROT(12+(I-1)*3+J)=CMATMUL(1,VLDT,FDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTBACK215(G215)
!
!Purpose: Rotate all the matrices back to the quark mass
!         basis using the previously derived V's
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),LU(3,3),LD(3,3)
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),LUP(3,3),LDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3),LUDUM(3,3),LDDUM(3,3)
!
      DOUBLE COMPLEX VLUT(3,3),VRUT(3,3),VLDT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL,G215(215)
      INTEGER I,J
!
      DO I=1,3
        DO J=1,3
          VLUT(I,J)=VLU(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDT(I,J)=VLD(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO   
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G215(3+(I-1)*3+J)
          FD(I,J)=G215(12+(I-1)*3+J)
          LU(I,J)=G215(33+(I-1)*3+J)
          LD(I,J)=G215(42+(I-1)*3+J)
        END DO
      END DO
!
!Now rotate the matrices back
!
      DO I=1,3 
        DO J=1,3
          FUDUM(I,J)=CMATMUL(2,FU,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(2,FD,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(2,LU,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(2,LD,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=CMATMUL(0,VLUT,FUDUM,I,J)
          FDP(I,J)=CMATMUL(0,VLDT,FDDUM,I,J)
          LUP(I,J)=CMATMUL(0,VLUT,LUDUM,I,J)
          LDP(I,J)=CMATMUL(0,VLDT,LDDUM,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G215(3+(I-1)*3+J)=FUP(I,J)
          G215(12+(I-1)*3+J)=FDP(I,J)
          G215(33+(I-1)*3+J)=LUP(I,J)
          G215(42+(I-1)*3+J)=LDP(I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTBACK(SMASS)
!
!Purpose: Rotate all the matrices from the original current basis.
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),MQ(3,3),MU(3,3),MD(3,3)
      DOUBLE COMPLEX AU(3,3),AD(3,3),LU(3,3),LD(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPU(3,3),GTPD(3,3),GTQ(3,3)
      DOUBLE COMPLEX GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE COMPLEX FTUQ(3,3),FTDQ(3,3),FTUU(3,3),FTDD(3,3)
      DOUBLE COMPLEX FMU(3,3),FMD(3,3),MQM(3,3),MUM(3,3),MDM(3,3)
      DOUBLE COMPLEX AUM(3,3),ADM(3,3)
      DOUBLE COMPLEX TRIU(3,3),TRID(3,3),MTSFU(3,3),MTSFD(3,3)
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),MQP(3,3),MUP(3,3),MDP(3,3)
      DOUBLE COMPLEX AUP(3,3),ADP(3,3),LUP(3,3),LDP(3,3)
      DOUBLE COMPLEX GTPQP(3,3),GTPUP(3,3),GTPDP(3,3),GTQP(3,3)
      DOUBLE COMPLEX GTSQP(3,3),GTSUP(3,3),GTSDP(3,3)
      DOUBLE COMPLEX FTUQP(3,3),FTDQP(3,3),FTUUP(3,3),FTDDP(3,3)
      DOUBLE COMPLEX FMUP(3,3),FMDP(3,3),MQMP(3,3),MUMP(3,3),MDMP(3,3)
      DOUBLE COMPLEX AUMP(3,3),ADMP(3,3)
      DOUBLE COMPLEX TRIUP(3,3),TRIDP(3,3),MTSFUP(3,3),MTSFDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX MQDUM(3,3),MUDUM(3,3),MDDUM(3,3)
      DOUBLE COMPLEX AUDUM(3,3),ADDUM(3,3),LUDUM(3,3),LDDUM(3,3)
      DOUBLE COMPLEX GTPQDUM(3,3),GTPUDUM(3,3),GTPDDUM(3,3),GTQDUM(3,3)
      DOUBLE COMPLEX GTSQDUM(3,3),GTSUDUM(3,3),GTSDDUM(3,3)
      DOUBLE COMPLEX FTUQDUM(3,3),FTDQDUM(3,3),FTUUDUM(3,3),FTDDDUM(3,3)
      DOUBLE COMPLEX FMUDUM(3,3),FMDDUM(3,3)
      DOUBLE COMPLEX MQMDUM(3,3),MUMDUM(3,3),MDMDUM(3,3)
      DOUBLE COMPLEX AUMDUM(3,3),ADMDUM(3,3)
      DOUBLE COMPLEX TRIUDUM(3,3),TRIDDUM(3,3),MTSFUDUM(3,3)
      DOUBLE COMPLEX MTSFDDUM(3,3)
!
      DOUBLE COMPLEX VLUQ(3,3),VLDQ(3,3)
      DOUBLE COMPLEX VLUQT(3,3),VRUT(3,3),VLDQT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J,SMASS
!
!Set the SU(2) doublet rotation
!If SMASS=1 we rotate to the quark mass basis
!Otherwise we rotate between different current bases.
!
      DO I=1,3
        DO J=1,3
          IF(SMASS.EQ.1)THEN
            VLUQ(I,J)=VLU(I,J)
            VLDQ(I,J)=VLD(I,J)
          ELSE
            IF(SVLQ.EQ.1)THEN
              VLUQ(I,J)=VLU(I,J)
              VLDQ(I,J)=VLU(I,J)
            ELSE
              VLUQ(I,J)=VLD(I,J)
              VLDQ(I,J)=VLD(I,J)
            END IF
          END IF
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          VLUQT(I,J)=VLUQ(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDQT(I,J)=VLDQ(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO   
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=G(3+(I-1)*3+J)
          FD(I,J)=G(12+(I-1)*3+J)
          AU(I,J)=G(33+(I-1)*3+J)
          AD(I,J)=G(42+(I-1)*3+J)
          MQ(I,J)=G(62+(I-1)*3+J)
          MU(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
!
          LU(I,J)=G(111+(I-1)*3+J)
          LD(I,J)=G(120+(I-1)*3+J)
!
          GTPQ(I,J)=G(138+(I-1)*3+J)
          GTPU(I,J)=G(156+(I-1)*3+J)
          GTPD(I,J)=G(165+(I-1)*3+J)
          GTQ(I,J)=G(185+(I-1)*3+J)
          GTSQ(I,J)=G(205+(I-1)*3+J)
          GTSU(I,J)=G(214+(I-1)*3+J)
          GTSD(I,J)=G(223+(I-1)*3+J)
          FTUQ(I,J)=G(232+(I-1)*3+J)
          FTDQ(I,J)=G(241+(I-1)*3+J)
          FTUU(I,J)=G(259+(I-1)*3+J)
          FTDD(I,J)=G(268+(I-1)*3+J)
!
          FMU(I,J)=G(293+(I-1)*3+J)
          FMD(I,J)=G(302+(I-1)*3+J)
          AUM(I,J)=G(323+(I-1)*3+J)
          ADM(I,J)=G(332+(I-1)*3+J)
          MQM(I,J)=G(352+(I-1)*3+J)
          MUM(I,J)=G(370+(I-1)*3+J)
          MDM(I,J)=G(379+(I-1)*3+J)
!
          TRIU(I,J)=G(399+(I-1)*3+J)
          TRID(I,J)=G(408+(I-1)*3+J)
          MTSFU(I,J)=G(429+(I-1)*3+J)
          MTSFD(I,J)=G(438+(I-1)*3+J)
        END DO
      END DO
!
!Now rotate the matrices back
!
      DO I=1,3 
        DO J=1,3
          FUDUM(I,J)=CMATMUL(2,FU,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(2,FD,VRDT,I,J)
          MQDUM(I,J)=CMATMUL(0,MQ,VLUQ,I,J)
          MUDUM(I,J)=CMATMUL(0,MU,VRU,I,J)
          MDDUM(I,J)=CMATMUL(0,MD,VRD,I,J)
          AUDUM(I,J)=CMATMUL(2,AU,VRUT,I,J)
          ADDUM(I,J)=CMATMUL(2,AD,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(2,LU,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(2,LD,VRDT,I,J)
          GTPQDUM(I,J)=CMATMUL(0,GTPQ,VLUQ,I,J)
          GTPUDUM(I,J)=CMATMUL(0,GTPU,VRU,I,J)
          GTPDDUM(I,J)=CMATMUL(0,GTPD,VRD,I,J)
          GTQDUM(I,J)=CMATMUL(0,GTQ,VLUQ,I,J)
          GTSQDUM(I,J)=CMATMUL(0,GTSQ,VLUQ,I,J)
          GTSUDUM(I,J)=CMATMUL(0,GTSU,VRU,I,J)
          GTSDDUM(I,J)=CMATMUL(0,GTSD,VRD,I,J)
          FTUQDUM(I,J)=CMATMUL(2,FTUQ,VRUT,I,J)
          FTDQDUM(I,J)=CMATMUL(2,FTDQ,VRDT,I,J)
          FTUUDUM(I,J)=CMATMUL(2,FTUU,VRUT,I,J)
          FTDDDUM(I,J)=CMATMUL(2,FTDD,VRDT,I,J)
          FMUDUM(I,J)=CMATMUL(2,FMU,VRUT,I,J)
          FMDDUM(I,J)=CMATMUL(2,FMD,VRDT,I,J)
          MQMDUM(I,J)=CMATMUL(0,MQM,VLUQ,I,J)
          MUMDUM(I,J)=CMATMUL(0,MUM,VRU,I,J)
          MDMDUM(I,J)=CMATMUL(0,MDM,VRD,I,J)
          AUMDUM(I,J)=CMATMUL(2,AUM,VRUT,I,J)
          ADMDUM(I,J)=CMATMUL(2,ADM,VRDT,I,J)
          TRIUDUM(I,J)=CMATMUL(2,TRIU,VRUT,I,J)
          TRIDDUM(I,J)=CMATMUL(2,TRID,VRDT,I,J)
          MTSFUDUM(I,J)=CMATMUL(2,MTSFU,VRUT,I,J)
          MTSFDDUM(I,J)=CMATMUL(2,MTSFD,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=CMATMUL(0,VLUQT,FUDUM,I,J)
          FDP(I,J)=CMATMUL(0,VLDQT,FDDUM,I,J)
          MQP(I,J)=CMATMUL(1,VLUQ,MQDUM,I,J)
          MUP(I,J)=CMATMUL(1,VRU,MUDUM,I,J)
          MDP(I,J)=CMATMUL(1,VRD,MDDUM,I,J)
          AUP(I,J)=CMATMUL(0,VLUQT,AUDUM,I,J)
          ADP(I,J)=CMATMUL(0,VLDQT,ADDUM,I,J)
          LUP(I,J)=CMATMUL(0,VLUQT,LUDUM,I,J)
          LDP(I,J)=CMATMUL(0,VLDQT,LDDUM,I,J)
          GTPQP(I,J)=CMATMUL(1,VLUQ,GTPQDUM,I,J)
          GTPUP(I,J)=CMATMUL(1,VRU,GTPUDUM,I,J)
          GTPDP(I,J)=CMATMUL(1,VRD,GTPDDUM,I,J)
          GTQP(I,J)=CMATMUL(1,VLUQ,GTQDUM,I,J)
          GTSQP(I,J)=CMATMUL(1,VLUQ,GTSQDUM,I,J)
          GTSUP(I,J)=CMATMUL(1,VRU,GTSUDUM,I,J)
          GTSDP(I,J)=CMATMUL(1,VRD,GTSDDUM,I,J)
          FTUQP(I,J)=CMATMUL(0,VLUQT,FTUQDUM,I,J)
          FTDQP(I,J)=CMATMUL(0,VLDQT,FTDQDUM,I,J)
          FTUUP(I,J)=CMATMUL(0,VLUQT,FTUUDUM,I,J)
          FTDDP(I,J)=CMATMUL(0,VLDQT,FTDDDUM,I,J)
          FMUP(I,J)=CMATMUL(0,VLUQT,FMUDUM,I,J)
          FMDP(I,J)=CMATMUL(0,VLDQT,FMDDUM,I,J)
          MQMP(I,J)=CMATMUL(1,VLUQ,MQMDUM,I,J)
          MUMP(I,J)=CMATMUL(1,VRU,MUMDUM,I,J)
          MDMP(I,J)=CMATMUL(1,VRD,MDMDUM,I,J)
          AUMP(I,J)=CMATMUL(0,VLUQT,AUMDUM,I,J)
          ADMP(I,J)=CMATMUL(0,VLDQT,ADMDUM,I,J)
          TRIUP(I,J)=CMATMUL(0,VLUQT,TRIUDUM,I,J)
          TRIDP(I,J)=CMATMUL(0,VLDQT,TRIDDUM,I,J)
          MTSFUP(I,J)=CMATMUL(0,VLUQT,MTSFUDUM,I,J)
          MTSFDP(I,J)=CMATMUL(0,VLDQT,MTSFDDUM,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          G(3+(I-1)*3+J)=FUP(I,J)
          G(12+(I-1)*3+J)=FDP(I,J)
          G(62+(I-1)*3+J)=MQP(I,J)
          G(80+(I-1)*3+J)=MUP(I,J)
          G(89+(I-1)*3+J)=MDP(I,J)
          G(33+(I-1)*3+J)=AUP(I,J)
          G(42+(I-1)*3+J)=ADP(I,J)
!
          G(111+(I-1)*3+J)=LUP(I,J)
          G(120+(I-1)*3+J)=LDP(I,J)
!
          G(138+(I-1)*3+J)=GTPQP(I,J)
          G(156+(I-1)*3+J)=GTPUP(I,J)
          G(165+(I-1)*3+J)=GTPDP(I,J)
          G(185+(I-1)*3+J)=GTQP(I,J)
          G(205+(I-1)*3+J)=GTSQP(I,J)
          G(214+(I-1)*3+J)=GTSUP(I,J)
          G(223+(I-1)*3+J)=GTSDP(I,J)
          G(232+(I-1)*3+J)=FTUQP(I,J)
          G(241+(I-1)*3+J)=FTDQP(I,J)
          G(259+(I-1)*3+J)=FTUUP(I,J)
          G(268+(I-1)*3+J)=FTDDP(I,J)
!
          G(293+(I-1)*3+J)=FMUP(I,J)
          G(302+(I-1)*3+J)=FMDP(I,J)
          G(352+(I-1)*3+J)=MQMP(I,J)
          G(370+(I-1)*3+J)=MUMP(I,J)
          G(379+(I-1)*3+J)=MDMP(I,J)
          G(323+(I-1)*3+J)=AUMP(I,J)
          G(332+(I-1)*3+J)=ADMP(I,J)
!
          G(399+(I-1)*3+J)=TRIUP(I,J)
          G(408+(I-1)*3+J)=TRIDP(I,J)
          G(429+(I-1)*3+J)=MTSFUP(I,J)
          G(438+(I-1)*3+J)=MTSFDP(I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTBACKSM(GROT)
!
!Purpose: Rotate the Yukawa matrices from the general current
!         basis to the chosen 'mass' basis
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FUP(3,3),FDP(3,3)
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX GROT(32)
!
      DOUBLE COMPLEX VLUT(3,3),VRUT(3,3),VLDT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J
!
!Convert the G's to matrices
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=GROT(3+(I-1)*3+J)
          FDP(I,J)=GROT(12+(I-1)*3+J)
        END DO
      END DO
!
!Calculate the transposes
!
      DO I=1,3
        DO J=1,3
          VLUT(I,J)=VLU(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDT(I,J)=VLD(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO
!
!Now rotate the matrices
!
      DO I=1,3
        DO J=1,3
          FUDUM(I,J)=CMATMUL(2,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(2,FDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GROT(3+(I-1)*3+J)=CMATMUL(0,VLUT,FUDUM,I,J)
          GROT(12+(I-1)*3+J)=CMATMUL(0,VLDT,FDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTSQBACK(GROT,GUNROT)
!
!Purpose: To rotate the couplings back from the squark mass
!         basis to the old current basis
!
      IMPLICIT NONE
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE COMPLEX GUNROT(601),GROT(601),CMATMUL
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3),AU(3,3),AD(3,3),AE(3,3)
      DOUBLE COMPLEX MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE COMPLEX FUM(3,3),FDM(3,3),FEM(3,3)
      DOUBLE COMPLEX AUM(3,3),ADM(3,3),AEM(3,3)
      DOUBLE COMPLEX MQM(3,3),MLM(3,3),MUPM(3,3),MDM(3,3),MEM(3,3)
      DOUBLE COMPLEX LU(3,3),LD(3,3),LE(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3)
      DOUBLE COMPLEX GTPE(3,3),GTQ(3,3),GTL(3,3),GTSQ(3,3)
      DOUBLE COMPLEX GTSU(3,3),GTSD(3,3),FTUQ(3,3),FTDQ(3,3)
      DOUBLE COMPLEX FTEL(3,3),FTUU(3,3),FTDD(3,3),FTEE(3,3)
      DOUBLE COMPLEX TRIU(3,3),TRID(3,3),TRIE(3,3)
      DOUBLE COMPLEX MTSFU(3,3),MTSFD(3,3),MTSFE(3,3)
!
      DOUBLE COMPLEX FUT(3,3),FDT(3,3),FET(3,3)
      DOUBLE COMPLEX AUT(3,3),ADT(3,3),AET(3,3)
      DOUBLE COMPLEX MQT(3,3),MLT(3,3),MUPT(3,3),MDT(3,3),MET(3,3)
      DOUBLE COMPLEX FUMT(3,3),FDMT(3,3),FEMT(3,3)
      DOUBLE COMPLEX AUMT(3,3),ADMT(3,3),AEMT(3,3)
      DOUBLE COMPLEX MQMT(3,3),MLMT(3,3),MUPMT(3,3),MDMT(3,3),MEMT(3,3)
      DOUBLE COMPLEX LUT(3,3),LDT(3,3),LET(3,3)
      DOUBLE COMPLEX GTPQT(3,3),GTPLT(3,3),GTPUT(3,3),GTPDT(3,3)
      DOUBLE COMPLEX GTPET(3,3),GTQT(3,3),GTLT(3,3),GTSQT(3,3)
      DOUBLE COMPLEX GTSUT(3,3),GTSDT(3,3),FTUQT(3,3),FTDQT(3,3)
      DOUBLE COMPLEX FTELT(3,3),FTUUT(3,3),FTDDT(3,3),FTEET(3,3)
      DOUBLE COMPLEX TRIUT(3,3),TRIDT(3,3),TRIET(3,3)
      DOUBLE COMPLEX MTSFUT(3,3),MTSFDT(3,3),MTSFET(3,3)
!
      DOUBLE COMPLEX RQT(3,3),RUPT(3,3),RDT(3,3),RLT(3,3),RET(3,3)
!
      INTEGER I,J
!
      DO I=1,601
        GUNROT(I)=GROT(I)
      END DO
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=GROT(3+(I-1)*3+J)
          FD(I,J)=GROT(12+(I-1)*3+J)
          FE(I,J)=GROT(21+(I-1)*3+J)
          AU(I,J)=GROT(33+(I-1)*3+J)
          AD(I,J)=GROT(42+(I-1)*3+J)
          AE(I,J)=GROT(51+(I-1)*3+J)
          MQ(I,J)=GROT(62+(I-1)*3+J)
          ML(I,J)=GROT(71+(I-1)*3+J)
          MUP(I,J)=GROT(80+(I-1)*3+J)
          MD(I,J)=GROT(89+(I-1)*3+J)
          ME(I,J)=GROT(98+(I-1)*3+J)
!
          FUM(I,J)=GROT(293+(I-1)*3+J)
          FDM(I,J)=GROT(302+(I-1)*3+J)
          FEM(I,J)=GROT(311+(I-1)*3+J)
          AUM(I,J)=GROT(323+(I-1)*3+J)
          ADM(I,J)=GROT(332+(I-1)*3+J)
          AEM(I,J)=GROT(341+(I-1)*3+J)
          MQM(I,J)=GROT(352+(I-1)*3+J)
          MLM(I,J)=GROT(361+(I-1)*3+J)
          MUPM(I,J)=GROT(370+(I-1)*3+J)
          MDM(I,J)=GROT(379+(I-1)*3+J)
          MEM(I,J)=GROT(388+(I-1)*3+J)
!
          LU(I,J)=GROT(111+(I-1)*3+J)
          LD(I,J)=GROT(120+(I-1)*3+J)
          LE(I,J)=GROT(129+(I-1)*3+J)
!
          GTPQ(I,J)=GROT(138+(I-1)*3+J)
          GTPL(I,J)=GROT(147+(I-1)*3+J)
          GTPU(I,J)=GROT(156+(I-1)*3+J)
          GTPD(I,J)=GROT(165+(I-1)*3+J)
          GTPE(I,J)=GROT(174+(I-1)*3+J)
          GTQ(I,J)=GROT(185+(I-1)*3+J)
          GTL(I,J)=GROT(194+(I-1)*3+J)
          GTSQ(I,J)=GROT(205+(I-1)*3+J)
          GTSU(I,J)=GROT(214+(I-1)*3+J)
          GTSD(I,J)=GROT(223+(I-1)*3+J)
          FTUQ(I,J)=GROT(232+(I-1)*3+J)
          FTDQ(I,J)=GROT(241+(I-1)*3+J)
          FTEL(I,J)=GROT(250+(I-1)*3+J)
          FTUU(I,J)=GROT(259+(I-1)*3+J)
          FTDD(I,J)=GROT(268+(I-1)*3+J)
          FTEE(I,J)=GROT(277+(I-1)*3+J)
!
          TRIU(I,J)=GROT(399+(I-1)*3+J)
          TRID(I,J)=GROT(408+(I-1)*3+J)
          TRIE(I,J)=GROT(417+(I-1)*3+J)
          MTSFU(I,J)=GROT(429+(I-1)*3+J)
          MTSFD(I,J)=GROT(438+(I-1)*3+J)
          MTSFE(I,J)=GROT(447+(I-1)*3+J)
!
          RQT(I,J)=RQTOT(J,I)
          RUPT(I,J)=RUPTOT(J,I)
          RDT(I,J)=RDTOT(J,I)
          RLT(I,J)=RLTOT(J,I)
          RET(I,J)=RETOT(J,I)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUT(I,J)=CMATMUL(0,FU,RUPT,I,J)
          FDT(I,J)=CMATMUL(0,FD,RDT,I,J)
          FET(I,J)=CMATMUL(0,FE,RET,I,J)
          AUT(I,J)=CMATMUL(0,AU,RUPT,I,J)
          ADT(I,J)=CMATMUL(0,AD,RDT,I,J)
          AET(I,J)=CMATMUL(0,AE,RET,I,J)
          MQT(I,J)=CMATMUL(2,MQ,RQTOT,I,J)
          MLT(I,J)=CMATMUL(2,ML,RLTOT,I,J)
          MUPT(I,J)=CMATMUL(2,MUP,RUPTOT,I,J)
          MDT(I,J)=CMATMUL(2,MD,RDTOT,I,J)
          MET(I,J)=CMATMUL(2,ME,RETOT,I,J)
!
          FUMT(I,J)=CMATMUL(0,FUM,RUPT,I,J)
          FDMT(I,J)=CMATMUL(0,FDM,RDT,I,J)
          FEMT(I,J)=CMATMUL(0,FEM,RET,I,J)
          AUMT(I,J)=CMATMUL(0,AUM,RUPT,I,J)
          ADMT(I,J)=CMATMUL(0,ADM,RDT,I,J)
          AEMT(I,J)=CMATMUL(0,AEM,RET,I,J)
          MQMT(I,J)=CMATMUL(2,MQM,RQTOT,I,J)
          MLMT(I,J)=CMATMUL(2,MLM,RLTOT,I,J)
          MUPMT(I,J)=CMATMUL(2,MUPM,RUPTOT,I,J)
          MDMT(I,J)=CMATMUL(2,MDM,RDTOT,I,J)
          MEMT(I,J)=CMATMUL(2,MEM,RETOT,I,J)
!
          LUT(I,J)=CMATMUL(0,LU,RUPT,I,J)
          LDT(I,J)=CMATMUL(0,LD,RDT,I,J)
          LET(I,J)=CMATMUL(0,LE,RET,I,J)
!
          GTPQT(I,J)=CMATMUL(2,GTPQ,RQTOT,I,J)
          GTPLT(I,J)=CMATMUL(2,GTPL,RLTOT,I,J)
          GTPUT(I,J)=CMATMUL(2,GTPU,RUPTOT,I,J)
          GTPDT(I,J)=CMATMUL(2,GTPD,RDTOT,I,J)
          GTPET(I,J)=CMATMUL(2,GTPE,RETOT,I,J)
          GTQT(I,J)=CMATMUL(2,GTQ,RQTOT,I,J)
          GTLT(I,J)=CMATMUL(2,GTL,RLTOT,I,J)
          GTSQT(I,J)=CMATMUL(2,GTSQ,RQTOT,I,J)
          GTSUT(I,J)=CMATMUL(2,GTSU,RUPTOT,I,J)
          GTSDT(I,J)=CMATMUL(2,GTSD,RDTOT,I,J)
          FTUQT(I,J)=CMATMUL(0,FTUQ,RUPT,I,J)
          FTDQT(I,J)=CMATMUL(0,FTDQ,RDT,I,J)
          FTELT(I,J)=CMATMUL(0,FTEL,RET,I,J)
          FTUUT(I,J)=CMATMUL(0,FTUU,RUPT,I,J)
          FTDDT(I,J)=CMATMUL(0,FTDD,RDT,I,J)
          FTEET(I,J)=CMATMUL(0,FTEE,RET,I,J)
!
          TRIUT(I,J)=CMATMUL(0,TRIU,RUPT,I,J)
          TRIDT(I,J)=CMATMUL(0,TRID,RDT,I,J)
          TRIET(I,J)=CMATMUL(0,TRIE,RET,I,J)
          MTSFUT(I,J)=CMATMUL(0,MTSFU,RUPT,I,J)
          MTSFDT(I,J)=CMATMUL(0,MTSFD,RDT,I,J)
          MTSFET(I,J)=CMATMUL(0,MTSFE,RET,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GUNROT(3+(I-1)*3+J)=CMATMUL(1,RQT,FUT,I,J)
          GUNROT(12+(I-1)*3+J)=CMATMUL(1,RQT,FDT,I,J)
          GUNROT(21+(I-1)*3+J)=CMATMUL(1,RLT,FET,I,J)
          GUNROT(33+(I-1)*3+J)=CMATMUL(1,RQT,AUT,I,J)
          GUNROT(42+(I-1)*3+J)=CMATMUL(1,RQT,ADT,I,J)
          GUNROT(51+(I-1)*3+J)=CMATMUL(1,RLT,AET,I,J)
          GUNROT(62+(I-1)*3+J)=CMATMUL(0,RQTOT,MQT,I,J)
          GUNROT(71+(I-1)*3+J)=CMATMUL(0,RLTOT,MLT,I,J)
          GUNROT(80+(I-1)*3+J)=CMATMUL(0,RUPTOT,MUPT,I,J)
          GUNROT(89+(I-1)*3+J)=CMATMUL(0,RDTOT,MDT,I,J)
          GUNROT(98+(I-1)*3+J)=CMATMUL(0,RETOT,MET,I,J)
!
          GUNROT(293+(I-1)*3+J)=CMATMUL(1,RQT,FUMT,I,J)
          GUNROT(302+(I-1)*3+J)=CMATMUL(1,RQT,FDMT,I,J)
          GUNROT(311+(I-1)*3+J)=CMATMUL(1,RLT,FEMT,I,J)
          GUNROT(323+(I-1)*3+J)=CMATMUL(1,RQT,AUMT,I,J)
          GUNROT(332+(I-1)*3+J)=CMATMUL(1,RQT,ADMT,I,J)
          GUNROT(341+(I-1)*3+J)=CMATMUL(1,RLT,AEMT,I,J)
          GUNROT(352+(I-1)*3+J)=CMATMUL(0,RQTOT,MQMT,I,J)
          GUNROT(361+(I-1)*3+J)=CMATMUL(0,RLTOT,MLMT,I,J)
          GUNROT(370+(I-1)*3+J)=CMATMUL(0,RUPTOT,MUPMT,I,J)
          GUNROT(379+(I-1)*3+J)=CMATMUL(0,RDTOT,MDMT,I,J)
          GUNROT(388+(I-1)*3+J)=CMATMUL(0,RETOT,MEMT,I,J)
!
          GUNROT(111+(I-1)*3+J)=CMATMUL(1,RQT,LUT,I,J)
          GUNROT(120+(I-1)*3+J)=CMATMUL(1,RQT,LDT,I,J)
          GUNROT(129+(I-1)*3+J)=CMATMUL(1,RLT,LET,I,J)
!
          GUNROT(138+(I-1)*3+J)=CMATMUL(0,RQTOT,GTPQT,I,J)
          GUNROT(147+(I-1)*3+J)=CMATMUL(0,RLTOT,GTPLT,I,J)
          GUNROT(156+(I-1)*3+J)=CMATMUL(0,RUPTOT,GTPUT,I,J)
          GUNROT(165+(I-1)*3+J)=CMATMUL(0,RDTOT,GTPDT,I,J)
          GUNROT(174+(I-1)*3+J)=CMATMUL(0,RETOT,GTPET,I,J)
          GUNROT(185+(I-1)*3+J)=CMATMUL(0,RQTOT,GTQT,I,J)
          GUNROT(194+(I-1)*3+J)=CMATMUL(0,RLTOT,GTLT,I,J)
          GUNROT(205+(I-1)*3+J)=CMATMUL(0,RQTOT,GTSQT,I,J)
          GUNROT(214+(I-1)*3+J)=CMATMUL(0,RUPTOT,GTSUT,I,J)
          GUNROT(223+(I-1)*3+J)=CMATMUL(0,RDTOT,GTSDT,I,J)
          GUNROT(232+(I-1)*3+J)=CMATMUL(1,RQT,FTUQT,I,J)
          GUNROT(241+(I-1)*3+J)=CMATMUL(1,RQT,FTDQT,I,J)
          GUNROT(250+(I-1)*3+J)=CMATMUL(1,RLT,FTELT,I,J)
          GUNROT(259+(I-1)*3+J)=CMATMUL(1,RQT,FTUUT,I,J)
          GUNROT(268+(I-1)*3+J)=CMATMUL(1,RQT,FTDDT,I,J)
          GUNROT(277+(I-1)*3+J)=CMATMUL(1,RLT,FTEET,I,J)
!
          GUNROT(399+(I-1)*3+J)=CMATMUL(1,RQT,TRIUT,I,J)
          GUNROT(408+(I-1)*3+J)=CMATMUL(1,RQT,TRIDT,I,J)
          GUNROT(417+(I-1)*3+J)=CMATMUL(1,RLT,TRIET,I,J)
          GUNROT(429+(I-1)*3+J)=CMATMUL(1,RQT,MTSFUT,I,J)
          GUNROT(438+(I-1)*3+J)=CMATMUL(1,RQT,MTSFDT,I,J)
          GUNROT(447+(I-1)*3+J)=CMATMUL(1,RLT,MTSFET,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE ROTSQ(GUNROT,GROT)
!
!Purpose: To rotate from the current basis to the squark mass basis
!
      IMPLICIT NONE
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE COMPLEX GUNROT(601),GROT(601),CMATMUL
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),FE(3,3),AU(3,3),AD(3,3),AE(3,3)
      DOUBLE COMPLEX MQ(3,3),ML(3,3),MUP(3,3),MD(3,3),ME(3,3)
      DOUBLE COMPLEX FUM(3,3),FDM(3,3),FEM(3,3)
      DOUBLE COMPLEX AUM(3,3),ADM(3,3),AEM(3,3)
      DOUBLE COMPLEX MQM(3,3),MLM(3,3),MUPM(3,3),MDM(3,3),MEM(3,3)
      DOUBLE COMPLEX LU(3,3),LD(3,3),LE(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPL(3,3),GTPU(3,3),GTPD(3,3)
      DOUBLE COMPLEX GTPE(3,3),GTQ(3,3),GTL(3,3),GTSQ(3,3)
      DOUBLE COMPLEX GTSU(3,3),GTSD(3,3),FTUQ(3,3),FTDQ(3,3)
      DOUBLE COMPLEX FTEL(3,3),FTUU(3,3),FTDD(3,3),FTEE(3,3)
      DOUBLE COMPLEX TRIU(3,3),TRID(3,3),TRIE(3,3)
      DOUBLE COMPLEX MTSFU(3,3),MTSFD(3,3),MTSFE(3,3)
!
      DOUBLE COMPLEX FUT(3,3),FDT(3,3),FET(3,3)
      DOUBLE COMPLEX AUT(3,3),ADT(3,3),AET(3,3)
      DOUBLE COMPLEX MQT(3,3),MLT(3,3),MUPT(3,3),MDT(3,3),MET(3,3)
      DOUBLE COMPLEX FUMT(3,3),FDMT(3,3),FEMT(3,3)
      DOUBLE COMPLEX AUMT(3,3),ADMT(3,3),AEMT(3,3)
      DOUBLE COMPLEX MQMT(3,3),MLMT(3,3),MUPMT(3,3),MDMT(3,3),MEMT(3,3)
      DOUBLE COMPLEX LUT(3,3),LDT(3,3),LET(3,3)
      DOUBLE COMPLEX GTPQT(3,3),GTPLT(3,3),GTPUT(3,3),GTPDT(3,3)
      DOUBLE COMPLEX GTPET(3,3),GTQT(3,3),GTLT(3,3),GTSQT(3,3)
      DOUBLE COMPLEX GTSUT(3,3),GTSDT(3,3),FTUQT(3,3),FTDQT(3,3)
      DOUBLE COMPLEX FTELT(3,3),FTUUT(3,3),FTDDT(3,3),FTEET(3,3)
      DOUBLE COMPLEX TRIUT(3,3),TRIDT(3,3),TRIET(3,3)
      DOUBLE COMPLEX MTSFUT(3,3),MTSFDT(3,3),MTSFET(3,3)
!
      DOUBLE COMPLEX RQT(3,3),RUPT(3,3),RDT(3,3),RLT(3,3),RET(3,3)
!
      INTEGER I,J
!
      DO I=1,601
        GROT(I)=GUNROT(I)
      END DO
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=GUNROT(3+(I-1)*3+J)
          FD(I,J)=GUNROT(12+(I-1)*3+J)
          FE(I,J)=GUNROT(21+(I-1)*3+J)
          AU(I,J)=GUNROT(33+(I-1)*3+J)
          AD(I,J)=GUNROT(42+(I-1)*3+J)
          AE(I,J)=GUNROT(51+(I-1)*3+J)
          MQ(I,J)=GUNROT(62+(I-1)*3+J)
          ML(I,J)=GUNROT(71+(I-1)*3+J)
          MUP(I,J)=GUNROT(80+(I-1)*3+J)
          MD(I,J)=GUNROT(89+(I-1)*3+J)
          ME(I,J)=GUNROT(98+(I-1)*3+J)
!
          FUM(I,J)=GUNROT(293+(I-1)*3+J)
          FDM(I,J)=GUNROT(302+(I-1)*3+J)
          FEM(I,J)=GUNROT(311+(I-1)*3+J)
          AUM(I,J)=GUNROT(323+(I-1)*3+J)
          ADM(I,J)=GUNROT(332+(I-1)*3+J)
          AEM(I,J)=GUNROT(341+(I-1)*3+J)
          MQM(I,J)=GUNROT(352+(I-1)*3+J)
          MLM(I,J)=GUNROT(361+(I-1)*3+J)
          MUPM(I,J)=GUNROT(370+(I-1)*3+J)
          MDM(I,J)=GUNROT(379+(I-1)*3+J)
          MEM(I,J)=GUNROT(388+(I-1)*3+J)
!
          LU(I,J)=GUNROT(111+(I-1)*3+J)
          LD(I,J)=GUNROT(120+(I-1)*3+J)
          LE(I,J)=GUNROT(129+(I-1)*3+J)
!
          GTPQ(I,J)=GUNROT(138+(I-1)*3+J)
          GTPL(I,J)=GUNROT(147+(I-1)*3+J)
          GTPU(I,J)=GUNROT(156+(I-1)*3+J)
          GTPD(I,J)=GUNROT(165+(I-1)*3+J)
          GTPE(I,J)=GUNROT(174+(I-1)*3+J)
          GTQ(I,J)=GUNROT(185+(I-1)*3+J)
          GTL(I,J)=GUNROT(194+(I-1)*3+J)
          GTSQ(I,J)=GUNROT(205+(I-1)*3+J)
          GTSU(I,J)=GUNROT(214+(I-1)*3+J)
          GTSD(I,J)=GUNROT(223+(I-1)*3+J)
          FTUQ(I,J)=GUNROT(232+(I-1)*3+J)
          FTDQ(I,J)=GUNROT(241+(I-1)*3+J)
          FTEL(I,J)=GUNROT(250+(I-1)*3+J)
          FTUU(I,J)=GUNROT(259+(I-1)*3+J)
          FTDD(I,J)=GUNROT(268+(I-1)*3+J)
          FTEE(I,J)=GUNROT(277+(I-1)*3+J)
!
          TRIU(I,J)=GUNROT(399+(I-1)*3+J)
          TRID(I,J)=GUNROT(408+(I-1)*3+J)
          TRIE(I,J)=GUNROT(417+(I-1)*3+J)
          MTSFU(I,J)=GUNROT(429+(I-1)*3+J)
          MTSFD(I,J)=GUNROT(438+(I-1)*3+J)
          MTSFE(I,J)=GUNROT(447+(I-1)*3+J)
!
          RQT(I,J)=RQTOT(J,I)
          RUPT(I,J)=RUPTOT(J,I)
          RDT(I,J)=RDTOT(J,I)
          RLT(I,J)=RLTOT(J,I)
          RET(I,J)=RETOT(J,I)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUT(I,J)=CMATMUL(2,FU,RUPT,I,J)
          FDT(I,J)=CMATMUL(2,FD,RDT,I,J)
          FET(I,J)=CMATMUL(2,FE,RET,I,J)
          AUT(I,J)=CMATMUL(2,AU,RUPT,I,J)
          ADT(I,J)=CMATMUL(2,AD,RDT,I,J)
          AET(I,J)=CMATMUL(2,AE,RET,I,J)
          MQT(I,J)=CMATMUL(0,MQ,RQTOT,I,J)
          MLT(I,J)=CMATMUL(0,ML,RLTOT,I,J)
          MUPT(I,J)=CMATMUL(0,MUP,RUPTOT,I,J)
          MDT(I,J)=CMATMUL(0,MD,RDTOT,I,J)
          MET(I,J)=CMATMUL(0,ME,RETOT,I,J)
!
          FUMT(I,J)=CMATMUL(2,FUM,RUPT,I,J)
          FDMT(I,J)=CMATMUL(2,FDM,RDT,I,J)
          FEMT(I,J)=CMATMUL(2,FEM,RET,I,J)
          AUMT(I,J)=CMATMUL(2,AUM,RUPT,I,J)
          ADMT(I,J)=CMATMUL(2,ADM,RDT,I,J)
          AEMT(I,J)=CMATMUL(2,AEM,RET,I,J)
          MQMT(I,J)=CMATMUL(0,MQM,RQTOT,I,J)
          MLMT(I,J)=CMATMUL(0,MLM,RLTOT,I,J)
          MUPMT(I,J)=CMATMUL(0,MUPM,RUPTOT,I,J)
          MDMT(I,J)=CMATMUL(0,MDM,RDTOT,I,J)
          MEMT(I,J)=CMATMUL(0,MEM,RETOT,I,J)
!
          LUT(I,J)=CMATMUL(2,LU,RUPT,I,J)
          LDT(I,J)=CMATMUL(2,LD,RDT,I,J)
          LET(I,J)=CMATMUL(2,LE,RET,I,J)
!
          GTPQT(I,J)=CMATMUL(0,GTPQ,RQTOT,I,J)
          GTPLT(I,J)=CMATMUL(0,GTPL,RLTOT,I,J)
          GTPUT(I,J)=CMATMUL(0,GTPU,RUPTOT,I,J)
          GTPDT(I,J)=CMATMUL(0,GTPD,RDTOT,I,J)
          GTPET(I,J)=CMATMUL(0,GTPE,RETOT,I,J)
          GTQT(I,J)=CMATMUL(0,GTQ,RQTOT,I,J)
          GTLT(I,J)=CMATMUL(0,GTL,RLTOT,I,J)
          GTSQT(I,J)=CMATMUL(0,GTSQ,RQTOT,I,J)
          GTSUT(I,J)=CMATMUL(0,GTSU,RUPTOT,I,J)
          GTSDT(I,J)=CMATMUL(0,GTSD,RDTOT,I,J)
          FTUQT(I,J)=CMATMUL(2,FTUQ,RUPT,I,J)
          FTDQT(I,J)=CMATMUL(2,FTDQ,RDT,I,J)
          FTELT(I,J)=CMATMUL(2,FTEL,RET,I,J)
          FTUUT(I,J)=CMATMUL(2,FTUU,RUPT,I,J)
          FTDDT(I,J)=CMATMUL(2,FTDD,RDT,I,J)
          FTEET(I,J)=CMATMUL(2,FTEE,RET,I,J)
!
          TRIUT(I,J)=CMATMUL(2,TRIU,RUPT,I,J)
          TRIDT(I,J)=CMATMUL(2,TRID,RDT,I,J)
          TRIET(I,J)=CMATMUL(2,TRIE,RET,I,J)
          MTSFUT(I,J)=CMATMUL(2,MTSFU,RUPT,I,J)
          MTSFDT(I,J)=CMATMUL(2,MTSFD,RDT,I,J)
          MTSFET(I,J)=CMATMUL(2,MTSFE,RET,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GROT(3+(I-1)*3+J)=CMATMUL(0,RQT,FUT,I,J)
          GROT(12+(I-1)*3+J)=CMATMUL(0,RQT,FDT,I,J)
          GROT(21+(I-1)*3+J)=CMATMUL(0,RLT,FET,I,J)
          GROT(33+(I-1)*3+J)=CMATMUL(0,RQT,AUT,I,J)
          GROT(42+(I-1)*3+J)=CMATMUL(0,RQT,ADT,I,J)
          GROT(51+(I-1)*3+J)=CMATMUL(0,RLT,AET,I,J)
          GROT(62+(I-1)*3+J)=CMATMUL(1,RQTOT,MQT,I,J)
          GROT(71+(I-1)*3+J)=CMATMUL(1,RLTOT,MLT,I,J)
          GROT(80+(I-1)*3+J)=CMATMUL(1,RUPTOT,MUPT,I,J)
          GROT(89+(I-1)*3+J)=CMATMUL(1,RDTOT,MDT,I,J)
          GROT(98+(I-1)*3+J)=CMATMUL(1,RETOT,MET,I,J)
!
          GROT(293+(I-1)*3+J)=CMATMUL(0,RQT,FUMT,I,J)
          GROT(302+(I-1)*3+J)=CMATMUL(0,RQT,FDMT,I,J)
          GROT(311+(I-1)*3+J)=CMATMUL(0,RLT,FEMT,I,J)
          GROT(323+(I-1)*3+J)=CMATMUL(0,RQT,AUMT,I,J)
          GROT(332+(I-1)*3+J)=CMATMUL(0,RQT,ADMT,I,J)
          GROT(341+(I-1)*3+J)=CMATMUL(0,RLT,AEMT,I,J)
          GROT(352+(I-1)*3+J)=CMATMUL(1,RQTOT,MQMT,I,J)
          GROT(361+(I-1)*3+J)=CMATMUL(1,RLTOT,MLMT,I,J)
          GROT(370+(I-1)*3+J)=CMATMUL(1,RUPTOT,MUPMT,I,J)
          GROT(379+(I-1)*3+J)=CMATMUL(1,RDTOT,MDMT,I,J)
          GROT(388+(I-1)*3+J)=CMATMUL(1,RETOT,MEMT,I,J)
!
          GROT(111+(I-1)*3+J)=CMATMUL(0,RQT,LUT,I,J)
          GROT(120+(I-1)*3+J)=CMATMUL(0,RQT,LDT,I,J)
          GROT(129+(I-1)*3+J)=CMATMUL(0,RLT,LET,I,J)
!
          GROT(138+(I-1)*3+J)=CMATMUL(1,RQTOT,GTPQT,I,J)
          GROT(147+(I-1)*3+J)=CMATMUL(1,RLTOT,GTPLT,I,J)
          GROT(156+(I-1)*3+J)=CMATMUL(1,RUPTOT,GTPUT,I,J)
          GROT(165+(I-1)*3+J)=CMATMUL(1,RDTOT,GTPDT,I,J)
          GROT(174+(I-1)*3+J)=CMATMUL(1,RETOT,GTPET,I,J)
          GROT(185+(I-1)*3+J)=CMATMUL(1,RQTOT,GTQT,I,J)
          GROT(194+(I-1)*3+J)=CMATMUL(1,RLTOT,GTLT,I,J)
          GROT(205+(I-1)*3+J)=CMATMUL(1,RQTOT,GTSQT,I,J)
          GROT(214+(I-1)*3+J)=CMATMUL(1,RUPTOT,GTSUT,I,J)
          GROT(223+(I-1)*3+J)=CMATMUL(1,RDTOT,GTSDT,I,J)
          GROT(232+(I-1)*3+J)=CMATMUL(0,RQT,FTUQT,I,J)
          GROT(241+(I-1)*3+J)=CMATMUL(0,RQT,FTDQT,I,J)
          GROT(250+(I-1)*3+J)=CMATMUL(0,RLT,FTELT,I,J)
          GROT(259+(I-1)*3+J)=CMATMUL(0,RQT,FTUUT,I,J)
          GROT(268+(I-1)*3+J)=CMATMUL(0,RQT,FTDDT,I,J)
          GROT(277+(I-1)*3+J)=CMATMUL(0,RLT,FTEET,I,J)
!
          GROT(399+(I-1)*3+J)=CMATMUL(0,RQT,TRIUT,I,J)
          GROT(408+(I-1)*3+J)=CMATMUL(0,RQT,TRIDT,I,J)
          GROT(417+(I-1)*3+J)=CMATMUL(0,RLT,TRIET,I,J)
          GROT(429+(I-1)*3+J)=CMATMUL(0,RQT,MTSFUT,I,J)
          GROT(438+(I-1)*3+J)=CMATMUL(0,RQT,MTSFDT,I,J)
          GROT(447+(I-1)*3+J)=CMATMUL(0,RLT,MTSFET,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      FUNCTION SFMUL(TH,A,B,I,J)
!
!Purpose: To multiply the two matrices A and B but taking
!         into account restrictions on the multiple from
!         thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A, B input as (3x3) complex
!         I, J integer
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3),B(3,3),SFMUL
      INTEGER LOOP,TH(3),I,J
!
      SFMUL=0.D0
      DO LOOP=1,3
        SFMUL=SFMUL+A(I,LOOP)*B(LOOP,J)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
!
      SUBROUTINE SORTTH
!
!Purpose: Sorts the thresholds into ascending order
!
      IMPLICIT NONE
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE PRECISION TEMP
      INTEGER I,J,K
!
      DO I=1,3
        QTHSORT(I)=QTHQL(I)
        QTHSORT(I+3)=QTHUR(I)
        QTHSORT(I+6)=QTHDR(I)
        QTHSORT(I+9)=QTHLL(I)
        QTHSORT(I+12)=QTHER(I)
      END DO
      QTHSORT(16)=QNSH
      QTHSORT(17)=QNSG
      QTHSORT(18)=QNH
      QTHSORT(19)=QTHSB
      QTHSORT(20)=QTHSW
      LOCMH=18
      DO J=2,20
        TEMP=QTHSORT(J)
        DO I=J-1,1,-1
          IF(QTHSORT(I).LT.TEMP)GOTO 10
          QTHSORT(I+1)=QTHSORT(I)
!
          IF(I.EQ.LOCMH)LOCMH=I+1
!
        END DO
        I=0
  10    QTHSORT(I+1)=TEMP
!
        IF(J.EQ.18)LOCMH=I+1
!
      END DO
!
      RETURN
      END
!
      SUBROUTINE SORTZG(VA,VE,OLDNACTIVE)
!
!Purpose: To sort the eigenvectors into order of increasing
!         eigenvalue after normalising them.
!         If the heaviest sfermion has decoupled, place
!         the smallest eigenvalue in the place where the largest
!         would be.
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX VA(3),VE(3,3),SWAPVA,SWAPVE(3)
      DOUBLE PRECISION SUM,REZ,IMZ,THETA
      INTEGER I,J,K,OLDNACTIVE
!
      DO J=1,3
        SUM=0.D0
        DO I=1,3
          SUM=SUM+DBLE(VE(I,J)*CONJG(VE(I,J)))
        END DO
        DO I=1,3
          VE(I,J)=VE(I,J)/DSQRT(SUM)
        END DO
      END DO
!
!Fix the largest entry to be real and positive.
!
      DO I=1,3
        IF(ABS(VE(1,I)).GT.ABS(VE(2,I))
     $              .AND.ABS(VE(1,I)).GT.ABS(VE(3,I)))THEN
          REZ=DBLE(VE(1,I))
          IMZ=DIMAG(VE(1,I))
        ELSE IF(ABS(VE(2,I)).GT.ABS(VE(3,I)))THEN
          REZ=DBLE(VE(2,I))
          IMZ=DIMAG(VE(2,I))
        ELSE
          REZ=DBLE(VE(3,I))
          IMZ=DIMAG(VE(3,I))
        END IF
        THETA=DATAN2(IMZ,REZ)
        DO J=1,3
          VE(J,I)=VE(J,I)*DCMPLX(DCOS(THETA),-DSIN(THETA))
        END DO
      END DO
!
!Now sort the eigenvectors
!
      DO I=1,2
        DO J=I+1,3
          IF(DBLE(VA(I)).GT.DBLE(VA(J)))THEN
            DO K=1,3
              SWAPVE(K)=VE(K,J)
              VE(K,J)=VE(K,I)
              VE(K,I)=SWAPVE(K)
            END DO
            SWAPVA=VA(J)
            VA(J)=VA(I)
            VA(I)=SWAPVA
          END IF
        END DO
      END DO
!
!Redo the sorting if we do not have 3 active sfermions
!Need to be careful about negative eigenvalues, which are
!OK as long as the value of the eigenvalue is positive
!at the weak scale.
!
      IF(OLDNACTIVE.EQ.2)THEN
        IF(DBLE(VA(1)).LT.0.AND.DBLE(VA(2)).GE.0)THEN
          DO K=1,3
            SWAPVE(K)=VE(K,2)
            VE(K,2)=VE(K,3)
            VE(K,3)=SWAPVE(K)
          END DO
          SWAPVA=VA(2)
          VA(2)=VA(3)
          VA(3)=SWAPVA
        ELSE IF(DBLE(VA(1)).GE.0)THEN
          DO I=2,1,-1
            DO K=1,3
              SWAPVE(K)=VE(K,I)
              VE(K,I)=VE(K,3)
              VE(K,3)=SWAPVE(K)
            END DO
            SWAPVA=VA(I)
            VA(I)=VA(3)
            VA(3)=SWAPVA
          END DO
        END IF
      END IF
      IF(OLDNACTIVE.EQ.1.AND.DBLE(VA(1)).GE.0)THEN
        DO I=2,1,-1
          DO K=1,3
            SWAPVE(K)=VE(K,I)
            VE(K,I)=VE(K,I+1)
            VE(K,I+1)=SWAPVE(K)
          END DO
          SWAPVA=VA(I)
          VA(I)=VA(I+1)
          VA(I+1)=SWAPVA
        END DO
      END IF
!
      RETURN
      END
!
      SUBROUTINE SQDIAG(GIN)
!
!Purpose: To find the matrices which will rotate to the sfermion
!         mass basis by diagonalising.
!
!         Contains a number of conditional statements to account
!         for the decoupling of squarks.
!
      IMPLICIT NONE
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON /SQROT/ RQTOT,RUPTOT,RDTOT,RLTOT,RETOT
     $               ,RQSAV,RUPSAV,RDSAV,RLSAV,RESAV
     $               ,OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      DOUBLE COMPLEX RQTOT(3,3),RUPTOT(3,3),RDTOT(3,3)
      DOUBLE COMPLEX RLTOT(3,3),RETOT(3,3)
      DOUBLE COMPLEX RQSAV(2,3,3),RUPSAV(2,3,3),RDSAV(2,3,3)
      DOUBLE COMPLEX RLSAV(2,3,3),RESAV(2,3,3)
      INTEGER OLDNSQ,OLDNSU,OLDNSD,OLDNSL,OLDNSE
      SAVE /SQROT/
!
      DOUBLE COMPLEX GIN(601),G(601)
      DOUBLE COMPLEX MQ(3,3),MUP(3,3),MD(3,3),ML(3,3),ME(3,3)
      DOUBLE COMPLEX MQ2(2,2),MQVA2(2),MQVE2(2,2)
      DOUBLE COMPLEX MUP2(2,2),MUPVA2(2),MUPVE2(2,2)
      DOUBLE COMPLEX MD2(2,2),MDVA2(2),MDVE2(2,2)
      DOUBLE COMPLEX ML2(2,2),MLVA2(2),MLVE2(2,2)
      DOUBLE COMPLEX ME2(2,2),MEVA2(2),MEVE2(2,2)
      DOUBLE COMPLEX EVERTMP(3,3),CWORK1(99),CWORK2(6)
      DOUBLE COMPLEX EVERTMP2(2,2),CWORK12(99),CWORK22(4)
      INTEGER I,J,DIAGERR,CIERR1,CIERR2,CIERR3,CIERR4,CIERR5
!
      CIERR1=0
      CIERR2=0
      CIERR3=0
      CIERR4=0
      CIERR5=0
      DIAGERR=0
      DO I=1,3
        MQVA(I)=(0.D0,0.D0)
        MUPVA(I)=(0.D0,0.D0)
        MDVA(I)=(0.D0,0.D0)
        MLVA(I)=(0.D0,0.D0)
        MEVA(I)=(0.D0,0.D0)
      END DO
!
!First rotate to the correct basis for the scale we are
!considering.
!
      CALL ROTSQ(GIN,G)
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=G(62+(I-1)*3+J)
          MUP(I,J)=G(80+(I-1)*3+J)
          MD(I,J)=G(89+(I-1)*3+J)
          ML(I,J)=G(71+(I-1)*3+J)
          ME(I,J)=G(98+(I-1)*3+J)
        END DO
      END DO
!
!Separate out the (2x2) entries for use if the heaviest sfermion
!has already decoupled.
!
      DO I=1,2
        DO J=1,2
          IF(OLDNSQ.EQ.2)THEN
            MQ2(I,J)=MQ(I,J)
          END IF
          IF(OLDNSU.EQ.2)THEN
            MUP2(I,J)=MUP(I,J)
          END IF
          IF(OLDNSD.EQ.2)THEN
            MD2(I,J)=MD(I,J)
          END IF
          IF(OLDNSL.EQ.2)THEN
            ML2(I,J)=ML(I,J)
          END IF
          IF(OLDNSE.EQ.2)THEN
            ME2(I,J)=ME(I,J)
          END IF
        END DO
      END DO
!
!Find eigenvectors and eigenvalues. If there are only two active
!squarks in each group, then only diagonalise the (2x2).
!
      IF(OLDNSQ.EQ.3)THEN
        CALL ZGEEV('V','V',3,MQ,3,MQVA,MQVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR1)
      ELSE IF(OLDNSQ.EQ.2)THEN
        CALL ZGEEV('V','V',2,MQ2,2,MQVA2,MQVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR1)
      END IF
      IF(OLDNSU.EQ.3)THEN
        CALL ZGEEV('V','V',3,MUP,3,MUPVA,MUPVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR2)
      ELSE IF(OLDNSU.EQ.2)THEN
        CALL ZGEEV('V','V',2,MUP2,2,MUPVA2,MUPVE2,2,EVERTMP2,2,CWORK12
     $             ,99,CWORK22,CIERR2)
      END IF
      IF(OLDNSD.EQ.3)THEN
        CALL ZGEEV('V','V',3,MD,3,MDVA,MDVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR3)
      ELSE IF(OLDNSD.EQ.2)THEN
        CALL ZGEEV('V','V',2,MD2,2,MDVA2,MDVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR3)
      END IF
      IF(OLDNSL.EQ.3)THEN
        CALL ZGEEV('V','V',3,ML,3,MLVA,MLVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR4)
      ELSE IF(OLDNSL.EQ.2)THEN
        CALL ZGEEV('V','V',2,ML2,2,MLVA2,MLVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR4)
      END IF
      IF(OLDNSE.EQ.3)THEN
        CALL ZGEEV('V','V',3,ME,3,MEVA,MEVE,3,EVERTMP,3,CWORK1,99
     $             ,CWORK2,CIERR5)
      ELSE IF(OLDNSE.EQ.2)THEN
        CALL ZGEEV('V','V',2,ME2,2,MEVA2,MEVE2,2,EVERTMP2,2,CWORK12,99
     $             ,CWORK22,CIERR5)
      END IF
!
!If only two active sfermions, set the third eigenvalue to zero
!and the (3,3) entry in the rotation to (1.d0,0.d0)
!
      DO I=1,3
        DO J=1,3
          IF(OLDNSQ.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MQVE(I,J)=MQVE2(I,J)
            ELSE
              MQVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MQVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSQ.LT.2)THEN
            MQVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MQVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSU.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MUPVE(I,J)=MUPVE2(I,J)
            ELSE
              MUPVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MUPVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSU.LT.2)THEN
            MUPVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MUPVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSD.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MDVE(I,J)=MDVE2(I,J)
            ELSE
              MDVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MDVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSD.LT.2)THEN
            MDVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MDVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSL.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MLVE(I,J)=MLVE2(I,J)
            ELSE
              MLVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MLVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSL.LT.2)THEN
            MLVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MLVE(I,J)=(1.D0,0.D0)
          END IF
          IF(OLDNSE.EQ.2)THEN
            IF(I.LT.3.AND.J.LT.3)THEN
              MEVE(I,J)=MEVE2(I,J)
            ELSE
              MEVE(I,J)=(0.D0,0.D0)
              IF(I.EQ.J)MEVE(3,3)=(1.D0,0.D0)
            END IF
          ELSE IF(OLDNSE.LT.2)THEN
            MEVE(I,J)=(0.D0,0.D0)
            IF(I.EQ.J)MEVE(I,J)=(1.D0,0.D0)
          END IF
        END DO
        IF(OLDNSQ.EQ.2)THEN
          IF(I.LT.3)THEN
            MQVA(I)=MQVA2(I)
            MQVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSU.EQ.2)THEN
          IF(I.LT.3)THEN
           MUPVA(I)=MUPVA2(I)
           MUPVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSD.EQ.2)THEN
          IF(I.LT.3)THEN
            MDVA(I)=MDVA2(I)
            MDVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSL.EQ.2)THEN
          IF(I.LT.3)THEN
            MLVA(I)=MLVA2(I)
            MLVA(3)=(0.D0,0.D0)
          END IF
        END IF
        IF(OLDNSE.EQ.2)THEN
          IF(I.LT.3)THEN
            MEVA(I)=MEVA2(I)
            MEVA(3)=(0.D0,0.D0)
          END IF
        END IF
      END DO
      IF(OLDNSQ.EQ.1)THEN
        MQVA(1)=MQ(1,1)
        DO I=2,3
          MQVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSU.EQ.1)THEN
        MUPVA(1)=MUP(1,1)
        DO I=2,3
          MUPVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSD.EQ.1)THEN
        MDVA(1)=MD(1,1)
        DO I=2,3
          MDVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSL.EQ.1)THEN
        MLVA(1)=ML(1,1)
        DO I=2,3
          MLVA(I)=(0.D0,0.D0)
        END DO
      END IF
      IF(OLDNSE.EQ.1)THEN
        MEVA(1)=ME(1,1)
        DO I=2,3
          MEVA(I)=(0.D0,0.D0)
        END DO
      END IF
!
!Sort the eigenvectors into order and normalise them
!
      CALL SORTZG(MQVA,MQVE,OLDNSQ)
      CALL SORTZG(MUPVA,MUPVE,OLDNSU)
      CALL SORTZG(MDVA,MDVE,OLDNSD)
      CALL SORTZG(MLVA,MLVE,OLDNSL)
      CALL SORTZG(MEVA,MEVE,OLDNSE)
!
!Error checking
!
      DO I=1,3
        IF(MQVE(1,I).EQ.(0.D0,0.D0).AND.MQVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MQVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MUPVE(1,I).EQ.(0.D0,0.D0).AND.MUPVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MUPVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MDVE(1,I).EQ.(0.D0,0.D0).AND.MDVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MDVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MLVE(1,I).EQ.(0.D0,0.D0).AND.MLVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MLVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
        IF(MEVE(1,I).EQ.(0.D0,0.D0).AND.MEVE(2,I).EQ.(0.D0,0.D0)
     $    .AND.MEVE(3,I).EQ.(0.D0,0.D0))DIAGERR=1
      END DO
      IF(DIAGERR.EQ.1)WRITE(*,*)
     $           'ERROR IN SQDIAG'
!
      IF(CIERR1.NE.0.OR.CIERR2.NE.0.OR.CIERR3.NE.0.OR.CIERR4.NE.0
     $  .OR.CIERR5.NE.0)WRITE(*,*)'ERROR WITH DIAGONALISATION IN SQDIAG'
!
      RETURN
      END
!
      SUBROUTINE SQSIX(GIN,QEND)
!
!Purpose: Receives the weak scale parameters.
!         Computes both the up-type and down-type mass matrices,
!         and (if SVLQ=1) subsequently the flavour changing stop
!         decay rate.
!         Also carries out a calculation based upon the one-step
!         formula by Hikasa-Kobayashi for comparison with the full
!         RGE method.
!
!         Must be used in conjunction with isajet since it uses
!         isajet results for neutralino eigenvectors.
!
!         The couplings are run to the scale of the lightest decoupling
!         chosen from both left- and right-handed squarks suitable for
!         the basis chosen with SVLQ, before the mass matrices are
!         calculated.
!
      IMPLICIT NONE
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX GIN(601),GDEC(601)
      DOUBLE PRECISION Q,QEND
      INTEGER I
!
      Q=QEND !QEND is the scale at which RGEFLAV called SQSIX
!
!Isajet was called with the opposite signs of mu and A_0
!since it also flips M_0. To make sure everything is
!consistent, I calculate the rest using the isajet
!numerical inputs, ie I swap the signs of M_{1,2,3}, A_0
!and mu.
!
      DO I=31,60
        GIN(I)=-GIN(I)
        GIN(I+290)=-GIN(I+290)
        IF(I.LT.58)GIN(I+369)=-GIN(I+369)
        IF(I.LT.58)GIN(I+399)=-GIN(I+399)
        IF(I.LT.34)GIN(I+568)=-GIN(I+568)
      END DO
      GIN(108)=-GIN(108)
      GIN(398)=-GIN(398)
!
      DO I=1,601
        GDEC(I)=GIN(I)
      END DO
!
!Save the couplings at m_H for use by the one-step calculation.
!GIN is passed to SQSIX in the current basis specified by SVLQ.
!
      IF(ABS(Q-QNH).LT.1.D-5)THEN
        DO I=1,601
          GHIK(I)=GIN(I)
        END DO
      END IF
!
!Input is g at m_H. Run to place where we will calculate the
!squark mass matrices, Q_DEC.
!
      CALL DECRUN(GIN,GDEC,Q,SVLQ)
!
!Calculate the (6x6) matrix for both the up-type and down-type squarks
!
      CALL USMMA(GDEC,Q,SVLQ)
      CALL DSMMA(GDEC,Q,SVLQ)
!
!Finally calculate the flavour changing decay rate, stop1-->c Neutralino
!
      IF(SVLQ.EQ.1)THEN
        CALL ST1CNEU(GDEC,Q)
      END IF
!
      RETURN
      END
!
      SUBROUTINE STROTATE(GMASS,GCURR,SMASS)
!
!Purpose: Rotate the parameters to the original current
!         basis
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX GMASS(601),GCURR(601)
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),MQP(3,3),MUP(3,3),MDP(3,3)
      DOUBLE COMPLEX AUP(3,3),ADP(3,3),LUP(3,3),LDP(3,3)
      DOUBLE COMPLEX GTPQP(3,3),GTPUP(3,3),GTPDP(3,3),GTQP(3,3)
      DOUBLE COMPLEX GTSQP(3,3),GTSUP(3,3),GTSDP(3,3)
      DOUBLE COMPLEX FTUQP(3,3),FTDQP(3,3),FTUUP(3,3),FTDDP(3,3)
      DOUBLE COMPLEX FMUP(3,3),FMDP(3,3),MQMP(3,3),MUMP(3,3),MDMP(3,3)
      DOUBLE COMPLEX AUMP(3,3),ADMP(3,3),TRIUP(3,3),TRIDP(3,3)
      DOUBLE COMPLEX MTSFUP(3,3),MTSFDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX MQDUM(3,3),MUDUM(3,3),MDDUM(3,3)
      DOUBLE COMPLEX AUDUM(3,3),ADDUM(3,3),LUDUM(3,3),LDDUM(3,3)
      DOUBLE COMPLEX GTPQDUM(3,3),GTPUDUM(3,3),GTPDDUM(3,3),GTQDUM(3,3)
      DOUBLE COMPLEX GTSQDUM(3,3),GTSUDUM(3,3),GTSDDUM(3,3)
      DOUBLE COMPLEX FTUQDUM(3,3),FTDQDUM(3,3),FTUUDUM(3,3),FTDDDUM(3,3)
      DOUBLE COMPLEX FMUDUM(3,3),FMDDUM(3,3)
      DOUBLE COMPLEX MQMDUM(3,3),MUMDUM(3,3),MDMDUM(3,3)
      DOUBLE COMPLEX AUMDUM(3,3),ADMDUM(3,3),TRIUDUM(3,3),TRIDDUM(3,3)
      DOUBLE COMPLEX MTSFUDUM(3,3),MTSFDDUM(3,3)
!
      DOUBLE COMPLEX VLUQ(3,3),VLDQ(3,3)
      DOUBLE COMPLEX VLUQT(3,3),VRUT(3,3),VLDQT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J,SMASS
!
!Set the SU(2) doublet rotation
!If SMASS=1 we rotate from the quark mass basis
!Otherwise we rotate between different current bases.
!
      DO I=1,3
        DO J=1,3
          IF(SMASS.EQ.1)THEN
            VLUQ(I,J)=VLU(I,J)
            VLDQ(I,J)=VLD(I,J)
          ELSE
            IF(SVLQ.EQ.1)THEN
              VLUQ(I,J)=VLU(I,J)
              VLDQ(I,J)=VLU(I,J)
            ELSE
              VLUQ(I,J)=VLD(I,J)
              VLDQ(I,J)=VLD(I,J)
            END IF
          END IF
        END DO
      END DO
!
!Convert the G's to matrices
!
      DO I=1,601
        GCURR(I)=GMASS(I)
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=GMASS(3+(I-1)*3+J)
          FDP(I,J)=GMASS(12+(I-1)*3+J)
          MQP(I,J)=GMASS(62+(I-1)*3+J)
          MUP(I,J)=GMASS(80+(I-1)*3+J)
          MDP(I,J)=GMASS(89+(I-1)*3+J)
          AUP(I,J)=GMASS(33+(I-1)*3+J)
          ADP(I,J)=GMASS(42+(I-1)*3+J)
!
          LUP(I,J)=GMASS(111+(I-1)*3+J)
          LDP(I,J)=GMASS(120+(I-1)*3+J)
!
          GTPQP(I,J)=GMASS(138+(I-1)*3+J)
          GTPUP(I,J)=GMASS(156+(I-1)*3+J)
          GTPDP(I,J)=GMASS(165+(I-1)*3+J)
          GTQP(I,J)=GMASS(185+(I-1)*3+J)
          GTSQP(I,J)=GMASS(205+(I-1)*3+J)
          GTSUP(I,J)=GMASS(214+(I-1)*3+J)
          GTSDP(I,J)=GMASS(223+(I-1)*3+J)
          FTUQP(I,J)=GMASS(232+(I-1)*3+J)
          FTDQP(I,J)=GMASS(241+(I-1)*3+J)
          FTUUP(I,J)=GMASS(259+(I-1)*3+J)
          FTDDP(I,J)=GMASS(268+(I-1)*3+J)
!
          FMUP(I,J)=GMASS(293+(I-1)*3+J)
          FMDP(I,J)=GMASS(302+(I-1)*3+J)
          MQMP(I,J)=GMASS(352+(I-1)*3+J)
          MUMP(I,J)=GMASS(370+(I-1)*3+J)
          MDMP(I,J)=GMASS(379+(I-1)*3+J)
          AUMP(I,J)=GMASS(323+(I-1)*3+J)
          ADMP(I,J)=GMASS(332+(I-1)*3+J)
!
          TRIUP(I,J)=GMASS(399+(I-1)*3+J)
          TRIDP(I,J)=GMASS(408+(I-1)*3+J)
          MTSFUP(I,J)=GMASS(429+(I-1)*3+J)
          MTSFDP(I,J)=GMASS(438+(I-1)*3+J)
        END DO
      END DO
!
!Calculate the transposes
!
      DO I=1,3
        DO J=1,3
          VLUQT(I,J)=VLUQ(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDQT(I,J)=VLDQ(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO
!
!Now rotate the matrices
!
      DO I=1,3
        DO J=1,3
          FUDUM(I,J)=CMATMUL(0,FUP,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(0,FDP,VRDT,I,J)
          MQDUM(I,J)=CMATMUL(2,MQP,VLUQ,I,J)
          MUDUM(I,J)=CMATMUL(2,MUP,VRU,I,J)
          MDDUM(I,J)=CMATMUL(2,MDP,VRD,I,J)
          AUDUM(I,J)=CMATMUL(0,AUP,VRUT,I,J)
          ADDUM(I,J)=CMATMUL(0,ADP,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(0,LUP,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(0,LDP,VRDT,I,J)
          GTPQDUM(I,J)=CMATMUL(2,GTPQP,VLUQ,I,J)
          GTPUDUM(I,J)=CMATMUL(2,GTPUP,VRU,I,J)
          GTPDDUM(I,J)=CMATMUL(2,GTPDP,VRD,I,J)
          GTQDUM(I,J)=CMATMUL(2,GTQP,VLUQ,I,J)
          GTSQDUM(I,J)=CMATMUL(2,GTSQP,VLUQ,I,J)
          GTSUDUM(I,J)=CMATMUL(2,GTSUP,VRU,I,J)
          GTSDDUM(I,J)=CMATMUL(2,GTSDP,VRD,I,J)
          FTUQDUM(I,J)=CMATMUL(0,FTUQP,VRUT,I,J)
          FTDQDUM(I,J)=CMATMUL(0,FTDQP,VRDT,I,J)
          FTUUDUM(I,J)=CMATMUL(0,FTUUP,VRUT,I,J)
          FTDDDUM(I,J)=CMATMUL(0,FTDDP,VRDT,I,J)
          FMUDUM(I,J)=CMATMUL(0,FMUP,VRUT,I,J)
          FMDDUM(I,J)=CMATMUL(0,FMDP,VRDT,I,J)
          MQMDUM(I,J)=CMATMUL(2,MQMP,VLUQ,I,J)
          MUMDUM(I,J)=CMATMUL(2,MUMP,VRU,I,J)
          MDMDUM(I,J)=CMATMUL(2,MDMP,VRD,I,J)
          AUMDUM(I,J)=CMATMUL(0,AUMP,VRUT,I,J)
          ADMDUM(I,J)=CMATMUL(0,ADMP,VRDT,I,J)
          TRIUDUM(I,J)=CMATMUL(0,TRIUP,VRUT,I,J)
          TRIDDUM(I,J)=CMATMUL(0,TRIDP,VRDT,I,J)
          MTSFUDUM(I,J)=CMATMUL(0,MTSFUP,VRUT,I,J)
          MTSFDDUM(I,J)=CMATMUL(0,MTSFDP,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GCURR(3+(I-1)*3+J)=CMATMUL(1,VLUQT,FUDUM,I,J)
          GCURR(12+(I-1)*3+J)=CMATMUL(1,VLDQT,FDDUM,I,J)
          GCURR(62+(I-1)*3+J)=CMATMUL(0,VLUQ,MQDUM,I,J)
          GCURR(80+(I-1)*3+J)=CMATMUL(0,VRU,MUDUM,I,J)
          GCURR(89+(I-1)*3+J)=CMATMUL(0,VRD,MDDUM,I,J)
          GCURR(33+(I-1)*3+J)=CMATMUL(1,VLUQT,AUDUM,I,J)
          GCURR(42+(I-1)*3+J)=CMATMUL(1,VLDQT,ADDUM,I,J)
!
          GCURR(111+(I-1)*3+J)=CMATMUL(1,VLUQT,LUDUM,I,J)
          GCURR(120+(I-1)*3+J)=CMATMUL(1,VLDQT,LDDUM,I,J)
          GCURR(138+(I-1)*3+J)=CMATMUL(0,VLUQ,GTPQDUM,I,J)
          GCURR(156+(I-1)*3+J)=CMATMUL(0,VRU,GTPUDUM,I,J)
          GCURR(165+(I-1)*3+J)=CMATMUL(0,VRD,GTPDDUM,I,J)
          GCURR(185+(I-1)*3+J)=CMATMUL(0,VLUQ,GTQDUM,I,J)
          GCURR(205+(I-1)*3+J)=CMATMUL(0,VLUQ,GTSQDUM,I,J)
          GCURR(214+(I-1)*3+J)=CMATMUL(0,VRU,GTSUDUM,I,J)
          GCURR(223+(I-1)*3+J)=CMATMUL(0,VRD,GTSDDUM,I,J)
          GCURR(232+(I-1)*3+J)=CMATMUL(1,VLUQT,FTUQDUM,I,J)
          GCURR(241+(I-1)*3+J)=CMATMUL(1,VLDQT,FTDQDUM,I,J)
          GCURR(259+(I-1)*3+J)=CMATMUL(1,VLUQT,FTUUDUM,I,J)
          GCURR(268+(I-1)*3+J)=CMATMUL(1,VLDQT,FTDDDUM,I,J)
!
          GCURR(293+(I-1)*3+J)=CMATMUL(1,VLUQT,FMUDUM,I,J)
          GCURR(302+(I-1)*3+J)=CMATMUL(1,VLDQT,FMDDUM,I,J)
          GCURR(352+(I-1)*3+J)=CMATMUL(0,VLUQ,MQMDUM,I,J)
          GCURR(370+(I-1)*3+J)=CMATMUL(0,VRU,MUMDUM,I,J)
          GCURR(379+(I-1)*3+J)=CMATMUL(0,VRD,MDMDUM,I,J)
          GCURR(323+(I-1)*3+J)=CMATMUL(1,VLUQT,AUMDUM,I,J)
          GCURR(332+(I-1)*3+J)=CMATMUL(1,VLDQT,ADMDUM,I,J)
!
          GCURR(399+(I-1)*3+J)=CMATMUL(1,VLUQT,TRIUDUM,I,J)
          GCURR(408+(I-1)*3+J)=CMATMUL(1,VLDQT,TRIDDUM,I,J)
          GCURR(429+(I-1)*3+J)=CMATMUL(1,VLUQT,MTSFUDUM,I,J)
          GCURR(438+(I-1)*3+J)=CMATMUL(1,VLDQT,MTSFDDUM,I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE STROTBACK(GCURR,GMASS,SMASS)
!
!Purpose: Rotate all the matrices from the original current basis.
!
      IMPLICIT NONE
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      DOUBLE COMPLEX GMASS(601),GCURR(601)
!
      DOUBLE COMPLEX FU(3,3),FD(3,3),MQ(3,3),MU(3,3),MD(3,3)
      DOUBLE COMPLEX AU(3,3),AD(3,3),LU(3,3),LD(3,3)
      DOUBLE COMPLEX GTPQ(3,3),GTPU(3,3),GTPD(3,3),GTQ(3,3)
      DOUBLE COMPLEX GTSQ(3,3),GTSU(3,3),GTSD(3,3)
      DOUBLE COMPLEX FTUQ(3,3),FTDQ(3,3),FTUU(3,3),FTDD(3,3)
      DOUBLE COMPLEX FMU(3,3),FMD(3,3),MQM(3,3),MUM(3,3),MDM(3,3)
      DOUBLE COMPLEX AUM(3,3),ADM(3,3)
      DOUBLE COMPLEX TRIU(3,3),TRID(3,3),MTSFU(3,3),MTSFD(3,3)
!
      DOUBLE COMPLEX FUP(3,3),FDP(3,3),MQP(3,3),MUP(3,3),MDP(3,3)
      DOUBLE COMPLEX AUP(3,3),ADP(3,3),LUP(3,3),LDP(3,3)
      DOUBLE COMPLEX GTPQP(3,3),GTPUP(3,3),GTPDP(3,3),GTQP(3,3)
      DOUBLE COMPLEX GTSQP(3,3),GTSUP(3,3),GTSDP(3,3)
      DOUBLE COMPLEX FTUQP(3,3),FTDQP(3,3),FTUUP(3,3),FTDDP(3,3)
      DOUBLE COMPLEX FMUP(3,3),FMDP(3,3),MQMP(3,3),MUMP(3,3),MDMP(3,3)
      DOUBLE COMPLEX AUMP(3,3),ADMP(3,3)
      DOUBLE COMPLEX TRIUP(3,3),TRIDP(3,3),MTSFUP(3,3),MTSFDP(3,3)
!
      DOUBLE COMPLEX FUDUM(3,3),FDDUM(3,3)
      DOUBLE COMPLEX MQDUM(3,3),MUDUM(3,3),MDDUM(3,3)
      DOUBLE COMPLEX AUDUM(3,3),ADDUM(3,3),LUDUM(3,3),LDDUM(3,3)
      DOUBLE COMPLEX GTPQDUM(3,3),GTPUDUM(3,3),GTPDDUM(3,3),GTQDUM(3,3)
      DOUBLE COMPLEX GTSQDUM(3,3),GTSUDUM(3,3),GTSDDUM(3,3)
      DOUBLE COMPLEX FTUQDUM(3,3),FTDQDUM(3,3),FTUUDUM(3,3),FTDDDUM(3,3)
      DOUBLE COMPLEX FMUDUM(3,3),FMDDUM(3,3)
      DOUBLE COMPLEX MQMDUM(3,3),MUMDUM(3,3),MDMDUM(3,3)
      DOUBLE COMPLEX AUMDUM(3,3),ADMDUM(3,3)
      DOUBLE COMPLEX TRIUDUM(3,3),TRIDDUM(3,3),MTSFUDUM(3,3)
      DOUBLE COMPLEX MTSFDDUM(3,3)
!
      DOUBLE COMPLEX VLUQ(3,3),VLDQ(3,3)
      DOUBLE COMPLEX VLUQT(3,3),VRUT(3,3),VLDQT(3,3),VRDT(3,3)
      DOUBLE COMPLEX CMATMUL
      INTEGER I,J,SMASS
!
!Set the SU(2) doublet rotation
!If SMASS=1 we rotate to the quark mass basis
!Otherwise we rotate between different current bases.
!
      DO I=1,3
        DO J=1,3
          IF(SMASS.EQ.1)THEN
            VLUQ(I,J)=VLU(I,J)
            VLDQ(I,J)=VLD(I,J)
          ELSE
            IF(SVLQ.EQ.1)THEN
              VLUQ(I,J)=VLU(I,J)
              VLDQ(I,J)=VLU(I,J)
            ELSE
              VLUQ(I,J)=VLD(I,J)
              VLDQ(I,J)=VLD(I,J)
            END IF
          END IF
        END DO
      END DO
!
      DO I=1,601
        GMASS(I)=GCURR(I)
      END DO
!
      DO I=1,3
        DO J=1,3
          VLUQT(I,J)=VLUQ(J,I)
          VRUT(I,J)=VRU(J,I)
          VLDQT(I,J)=VLDQ(J,I)
          VRDT(I,J)=VRD(J,I)
        END DO
      END DO   
!
      DO I=1,3
        DO J=1,3
          FU(I,J)=GCURR(3+(I-1)*3+J)
          FD(I,J)=GCURR(12+(I-1)*3+J)
          MQ(I,J)=GCURR(62+(I-1)*3+J)
          MU(I,J)=GCURR(80+(I-1)*3+J)
          MD(I,J)=GCURR(89+(I-1)*3+J)
          AU(I,J)=GCURR(33+(I-1)*3+J)
          AD(I,J)=GCURR(42+(I-1)*3+J)
!
          LU(I,J)=GCURR(111+(I-1)*3+J)
          LD(I,J)=GCURR(120+(I-1)*3+J)
          GTPQ(I,J)=GCURR(138+(I-1)*3+J)
          GTPU(I,J)=GCURR(156+(I-1)*3+J)
          GTPD(I,J)=GCURR(165+(I-1)*3+J)
          GTQ(I,J)=GCURR(185+(I-1)*3+J)
          GTSQ(I,J)=GCURR(205+(I-1)*3+J)
          GTSU(I,J)=GCURR(214+(I-1)*3+J)
          GTSD(I,J)=GCURR(223+(I-1)*3+J)
          FTUQ(I,J)=GCURR(232+(I-1)*3+J)
          FTDQ(I,J)=GCURR(241+(I-1)*3+J)
          FTUU(I,J)=GCURR(259+(I-1)*3+J)
          FTDD(I,J)=GCURR(268+(I-1)*3+J)
!
          FMU(I,J)=GCURR(293+(I-1)*3+J)
          FMD(I,J)=GCURR(302+(I-1)*3+J)
          AUM(I,J)=GCURR(323+(I-1)*3+J)
          ADM(I,J)=GCURR(332+(I-1)*3+J)
          MQM(I,J)=GCURR(352+(I-1)*3+J)
          MUM(I,J)=GCURR(370+(I-1)*3+J)
          MDM(I,J)=GCURR(379+(I-1)*3+J)
!
          TRIU(I,J)=GCURR(399+(I-1)*3+J)
          TRID(I,J)=GCURR(408+(I-1)*3+J)
          MTSFU(I,J)=GCURR(429+(I-1)*3+J)
          MTSFD(I,J)=GCURR(438+(I-1)*3+J)
        END DO
      END DO
!
!Now rotate the matrices back
!
      DO I=1,3 
        DO J=1,3
          FUDUM(I,J)=CMATMUL(2,FU,VRUT,I,J)
          FDDUM(I,J)=CMATMUL(2,FD,VRDT,I,J)
          MQDUM(I,J)=CMATMUL(0,MQ,VLUQ,I,J)
          MUDUM(I,J)=CMATMUL(0,MU,VRU,I,J)
          MDDUM(I,J)=CMATMUL(0,MD,VRD,I,J)
          AUDUM(I,J)=CMATMUL(2,AU,VRUT,I,J)
          ADDUM(I,J)=CMATMUL(2,AD,VRDT,I,J)
          LUDUM(I,J)=CMATMUL(2,LU,VRUT,I,J)
          LDDUM(I,J)=CMATMUL(2,LD,VRDT,I,J)
          GTPQDUM(I,J)=CMATMUL(0,GTPQ,VLUQ,I,J)
          GTPUDUM(I,J)=CMATMUL(0,GTPU,VRU,I,J)
          GTPDDUM(I,J)=CMATMUL(0,GTPD,VRD,I,J)
          GTQDUM(I,J)=CMATMUL(0,GTQ,VLUQ,I,J)
          GTSQDUM(I,J)=CMATMUL(0,GTSQ,VLUQ,I,J)
          GTSUDUM(I,J)=CMATMUL(0,GTSU,VRU,I,J)
          GTSDDUM(I,J)=CMATMUL(0,GTSD,VRD,I,J)
          FTUQDUM(I,J)=CMATMUL(2,FTUQ,VRUT,I,J)
          FTDQDUM(I,J)=CMATMUL(2,FTDQ,VRDT,I,J)
          FTUUDUM(I,J)=CMATMUL(2,FTUU,VRUT,I,J)
          FTDDDUM(I,J)=CMATMUL(2,FTDD,VRDT,I,J)
          FMUDUM(I,J)=CMATMUL(2,FMU,VRUT,I,J)
          FMDDUM(I,J)=CMATMUL(2,FMD,VRDT,I,J)
          MQMDUM(I,J)=CMATMUL(0,MQM,VLUQ,I,J)
          MUMDUM(I,J)=CMATMUL(0,MUM,VRU,I,J)
          MDMDUM(I,J)=CMATMUL(0,MDM,VRD,I,J)
          AUMDUM(I,J)=CMATMUL(2,AUM,VRUT,I,J)
          ADMDUM(I,J)=CMATMUL(2,ADM,VRDT,I,J)
          TRIUDUM(I,J)=CMATMUL(2,TRIU,VRUT,I,J)
          TRIDDUM(I,J)=CMATMUL(2,TRID,VRDT,I,J)
          MTSFUDUM(I,J)=CMATMUL(2,MTSFU,VRUT,I,J)
          MTSFDDUM(I,J)=CMATMUL(2,MTSFD,VRDT,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          FUP(I,J)=CMATMUL(0,VLUQT,FUDUM,I,J)
          FDP(I,J)=CMATMUL(0,VLDQT,FDDUM,I,J)
          MQP(I,J)=CMATMUL(1,VLUQ,MQDUM,I,J)
          MUP(I,J)=CMATMUL(1,VRU,MUDUM,I,J)
          MDP(I,J)=CMATMUL(1,VRD,MDDUM,I,J)
          AUP(I,J)=CMATMUL(0,VLUQT,AUDUM,I,J)
          ADP(I,J)=CMATMUL(0,VLDQT,ADDUM,I,J)
          LUP(I,J)=CMATMUL(0,VLUQT,LUDUM,I,J)
          LDP(I,J)=CMATMUL(0,VLDQT,LDDUM,I,J)
          GTPQP(I,J)=CMATMUL(1,VLUQ,GTPQDUM,I,J)
          GTPUP(I,J)=CMATMUL(1,VRU,GTPUDUM,I,J)
          GTPDP(I,J)=CMATMUL(1,VRD,GTPDDUM,I,J)
          GTQP(I,J)=CMATMUL(1,VLUQ,GTQDUM,I,J)
          GTSQP(I,J)=CMATMUL(1,VLUQ,GTSQDUM,I,J)
          GTSUP(I,J)=CMATMUL(1,VRU,GTSUDUM,I,J)
          GTSDP(I,J)=CMATMUL(1,VRD,GTSDDUM,I,J)
          FTUQP(I,J)=CMATMUL(0,VLUQT,FTUQDUM,I,J)
          FTDQP(I,J)=CMATMUL(0,VLDQT,FTDQDUM,I,J)
          FTUUP(I,J)=CMATMUL(0,VLUQT,FTUUDUM,I,J)
          FTDDP(I,J)=CMATMUL(0,VLDQT,FTDDDUM,I,J)
          FMUP(I,J)=CMATMUL(0,VLUQT,FMUDUM,I,J)
          FMDP(I,J)=CMATMUL(0,VLDQT,FMDDUM,I,J)
          MQMP(I,J)=CMATMUL(1,VLUQ,MQMDUM,I,J)
          MUMP(I,J)=CMATMUL(1,VRU,MUMDUM,I,J)
          MDMP(I,J)=CMATMUL(1,VRD,MDMDUM,I,J)
          AUMP(I,J)=CMATMUL(0,VLUQT,AUMDUM,I,J)
          ADMP(I,J)=CMATMUL(0,VLDQT,ADMDUM,I,J)
          TRIUP(I,J)=CMATMUL(0,VLUQT,TRIUDUM,I,J)
          TRIDP(I,J)=CMATMUL(0,VLDQT,TRIDDUM,I,J)
          MTSFUP(I,J)=CMATMUL(0,VLUQT,MTSFUDUM,I,J)
          MTSFDP(I,J)=CMATMUL(0,VLDQT,MTSFDDUM,I,J)
        END DO
      END DO
!
      DO I=1,3
        DO J=1,3
          GMASS(3+(I-1)*3+J)=FUP(I,J)
          GMASS(12+(I-1)*3+J)=FDP(I,J)
          GMASS(62+(I-1)*3+J)=MQP(I,J)
          GMASS(80+(I-1)*3+J)=MUP(I,J)
          GMASS(89+(I-1)*3+J)=MDP(I,J)
          GMASS(33+(I-1)*3+J)=AUP(I,J)
          GMASS(42+(I-1)*3+J)=ADP(I,J)
!
          GMASS(111+(I-1)*3+J)=LUP(I,J)
          GMASS(120+(I-1)*3+J)=LDP(I,J)
!
          GMASS(138+(I-1)*3+J)=GTPQP(I,J)
          GMASS(156+(I-1)*3+J)=GTPUP(I,J)
          GMASS(165+(I-1)*3+J)=GTPDP(I,J)
          GMASS(185+(I-1)*3+J)=GTQP(I,J)
          GMASS(205+(I-1)*3+J)=GTSQP(I,J)
          GMASS(214+(I-1)*3+J)=GTSUP(I,J)
          GMASS(223+(I-1)*3+J)=GTSDP(I,J)
          GMASS(232+(I-1)*3+J)=FTUQP(I,J)
          GMASS(241+(I-1)*3+J)=FTDQP(I,J)
          GMASS(259+(I-1)*3+J)=FTUUP(I,J)
          GMASS(268+(I-1)*3+J)=FTDDP(I,J)
!
          GMASS(293+(I-1)*3+J)=FMUP(I,J)
          GMASS(302+(I-1)*3+J)=FMDP(I,J)
          GMASS(352+(I-1)*3+J)=MQMP(I,J)
          GMASS(370+(I-1)*3+J)=MUMP(I,J)
          GMASS(379+(I-1)*3+J)=MDMP(I,J)
          GMASS(323+(I-1)*3+J)=AUMP(I,J)
          GMASS(332+(I-1)*3+J)=ADMP(I,J)
!
          GMASS(399+(I-1)*3+J)=TRIUP(I,J)
          GMASS(408+(I-1)*3+J)=TRIDP(I,J)
          GMASS(429+(I-1)*3+J)=MTSFUP(I,J)
          GMASS(438+(I-1)*3+J)=MTSFDP(I,J)
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE SYMMTEST(IN,OUT)
!
!Purpose: Test if real matrix IN(3,3) is hermitian.
!         Return 1 if not, 0 otherwise.
!
      IMPLICIT NONE
!
      DOUBLE PRECISION IN(3,3)
      INTEGER OUT
!
      OUT=0
      IF(IN(1,2).NE.IN(2,1))OUT=1
      IF(IN(1,3).NE.IN(3,1))OUT=1
      IF(IN(2,3).NE.IN(3,2))OUT=1
!
      RETURN
      END
!
      FUNCTION TCSFMUL(TH,A)
!
!Purpose: To trace the matrix A but taking into account
!         restrictions on the multiple from thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A input as (3x3) complex
!
      IMPLICIT NONE
!
      DOUBLE COMPLEX A(3,3),TCSFMUL
      INTEGER LOOP,TH(3),I,J
!
      TCSFMUL=(0.D0,0.D0)
      DO LOOP=1,3
        TCSFMUL=TCSFMUL+A(LOOP,LOOP)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
!
      FUNCTION TSFMUL(TH,A)
!
!Purpose: To trace the matrix A but taking into account
!         restrictions on the multiple from thetas.
!
!Expects: TH input as (3) integer (1 or 0)
!         A input as (3x3) complex
!
      IMPLICIT NONE
!
      DOUBLE PRECISION A(3,3),TSFMUL
      INTEGER LOOP,TH(3),I,J
!
      TSFMUL=0.D0
      DO LOOP=1,3
        TSFMUL=TSFMUL+A(LOOP,LOOP)*DBLE(TH(LOOP))
      END DO
!
      RETURN
      END
!
      SUBROUTINE UOUTCOUP(G215,Q)
!
!Purpose: To output the up-type Yukawas, in either
!         the mass or current basis
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE COMPLEX GTMP(215),G215(215)
      DOUBLE PRECISION Q,COSB,SINB
      INTEGER I,SWROT
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!This switch can be used to choose to write the rotated values.
!
      SWROT=1
!
!Remember the old gs
!
      DO I=1,215
        GTMP(I)=G215(I)
      END DO
!
!Perform the rotation to the mass basis if nec.
!
      IF(SWROT.EQ.1)CALL ROTBACK215(G215)
!
!Make some adjustments to facilitate comparisons
!
      DO I=1,215
        IF(ABS(G215(I)).LT.1.D-30)THEN
          G215(I)=(1.D-30,1.D-30)
        END IF
      END DO
!
!Print the set of terms out.
!
      IF(Q.GE.QNH-EPS)THEN
        WRITE(35,52)Q,ABS(G215(4)),ABS(G215(5)),ABS(G215(6)),
     $                ABS(G215(7)),ABS(G215(8)),ABS(G215(9)),
     $                ABS(G215(10)),ABS(G215(11)),ABS(G215(12))
      ELSE
        WRITE(35,52)Q,ABS(G215(34))/SINB,ABS(G215(35))/SINB,
     $   ABS(G215(36))/SINB,ABS(G215(37))/SINB,ABS(G215(38))/SINB,
     $   ABS(G215(39))/SINB,ABS(G215(40))/SINB,ABS(G215(41))/SINB,
     $   ABS(G215(42))/SINB
      END IF
!
!Return the gs to their original values
!
      DO I=1,215
        G215(I)=GTMP(I)
      END DO
!
   52 FORMAT(SP,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10
     $,1X,D17.10,1X,D17.10,1X,D17.10,1X,D17.10)
!
      RETURN
      END
!
      SUBROUTINE UPMHCOND2
!
!Purpose: Apply the matching conditions at m_H when running up
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      DOUBLE PRECISION SINB,COSB
      INTEGER I,J
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
      VEVMH=G(428)
      G(111)=DCMPLX(DSQRT((DBLE(G(428)))**2/(1.D0+TANB**2)))
      G(110)=DCMPLX(DSQRT((DBLE(G(428)))**2*TANB**2
     $        /(1.D0+TANB**2)))
!
      DO I=4,12
        G(I)=G(I+108)/SINB
      END DO
      DO I=13,30
        G(I)=G(I+108)/COSB
      END DO
!
!Refine the EWSB parameters if m_H>m_SUSY
!
      CSHSH=G(427)
      IF(RGEMS.LT.QNH)THEN
        G(62)=COSB**2*CSHSH-2.D0*SINB*COSB*CSHLH+SINB**2*CLHLH
        G(61)=SINB**2*CSHSH+2.D0*SINB*COSB*CSHLH+COSB**2*CLHLH
        G(109)=SINB*COSB*(-CSHSH+CLHLH)+(SINB**2-COSB**2)*CSHLH
      END IF
!
      RETURN
      END
!
      SUBROUTINE UPMHCOND(GCON)
!
!Purpose: Apply the matching conditions at m_H when running up
!
      IMPLICIT NONE
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      DOUBLE COMPLEX GCON(215)
      DOUBLE PRECISION SINB,COSB
      INTEGER I
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!This uses the input value of tanb - i.e. tanb at m_H.
!
      VEVMH=GCON(61)
      GCON(33)=DCMPLX(DSQRT((DBLE(VEVMH))**2/(1.D0+TANB**2)))
      GCON(32)=TANB*GCON(33)
!
      DO I=4,12
        GCON(I)=GCON(I+30)/SINB
      END DO
      DO I=13,30
        GCON(I)=GCON(I+30)/COSB
      END DO
!
      RETURN
      END
!
      SUBROUTINE UPMZMHIGH2(NSTEP)
!
!Purpose: To run back up to M_GUT as part of the iterative running.
!         Although different in name to distinguish from UPMZMHIGH,
!         this subroutine does almost exactly the same thing but
!         starts with saved values of GSM at M_T
!
      IMPLICIT NONE
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/WKYUK/LAMTMT,LAMBMZ,LAMTAMZ
      DOUBLE PRECISION LAMTMT,LAMBMZ,LAMTAMZ
      SAVE/WKYUK/
!
      COMMON/ATMZ/G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      DOUBLE COMPLEX G1MZ,G2MZ,G3MZ,VSMMZ,LAMBDAMZ,LAMTMZ
      SAVE/ATMZ/
!
      COMMON/SMRGE/SMRGEMH,SMQSTEP,NU,SMDR2LP
      DOUBLE PRECISION SMRGEMH,SMQSTEP
      INTEGER NU,SMDR2LP
      SAVE/SMRGE/
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/LOOPS/SSQSTEP,SW2LP
      DOUBLE PRECISION SSQSTEP
      INTEGER SW2LP
      SAVE/LOOPS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      DOUBLE COMPLEX W(1803),GSM(32),WSM(96),CID(3,3)
      DOUBLE PRECISION DW(1803),DGSM(32),DWSM(96)
      DOUBLE PRECISION PI,TMT,TZ,TTH(20),THIGH,DT,T,A1I,A2I,Q
      INTEGER NSTEP,NSTEPSM,I,J,II,LOOPNSTEP,NSTEPMT,BELOWMS
      EXTERNAL CRGE601,DRGE601,CSMRGEDR,DSMRGEDR
!
      DATA CID(1,1)/(1.D0,0.D0)/,CID(1,2)/(0.D0,0.D0)/
     $    ,CID(1,3)/(0.D0,0.D0)/
      DATA CID(2,1)/(0.D0,0.D0)/,CID(2,2)/(1.D0,0.D0)/
     $    ,CID(2,3)/(0.D0,0.D0)/
      DATA CID(3,1)/(0.D0,0.D0)/,CID(3,2)/(0.D0,0.D0)/
     $    ,CID(3,3)/(1.D0,0.D0)/
!
!Files for outcoup
!
!      OPEN(71,FILE='out/uau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(72,FILE='out/uad.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(73,FILE='out/uae.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(74,FILE='out/umhumu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(75,FILE='out/umhdmu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(76,FILE='out/umq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(77,FILE='out/umup.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(78,FILE='out/umd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(79,FILE='out/uml.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(80,FILE='out/ume.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(81,FILE='out/umu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(82,FILE='out/ub.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(83,FILE='out/umtsfu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(84,FILE='out/umtsfd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(85,FILE='out/umtsfe.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(86,FILE='out/utriu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(87,FILE='out/utrid.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(88,FILE='out/utrie.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(89,FILE='out/umhud.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(90,FILE='out/uaum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(91,FILE='out/uadm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(92,FILE='out/uaem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(93,FILE='out/umhumum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(94,FILE='out/umhdmum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(95,FILE='out/umqm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(96,FILE='out/umupm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(97,FILE='out/umdm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(98,FILE='out/umlm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(99,FILE='out/umem.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(70,FILE='out/umum.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(69,FILE='out/ubm.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(68,FILE='out/ug1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(67,FILE='out/ug2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(66,FILE='out/ug3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(65,FILE='out/uytau.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(64,FILE='out/uyb.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(63,FILE='out/uyu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(62,FILE='out/um1.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(61,FILE='out/um2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(60,FILE='out/um3.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(59,FILE='out/ugtpq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(58,FILE='out/ugtpu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(57,FILE='out/ugtq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(56,FILE='out/uftuq.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(55,FILE='out/dvu.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(54,FILE='out/dvd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!      OPEN(53,FILE='out/dmuflavcomp.dat',STATUS='UNKNOWN'
!     $ ,FORM='FORMATTED')
!      OPEN(52,FILE='out/dm12comp.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!
      DO I=1,32
        GSM(I)=(0.D0,0.D0)
        DGSM(I)=0.D0
      END DO
!
      BELOWMS=1
      PI=4.D0*DATAN(1.D0)
!
!Set NSTEP to a larger number than the downward running so that I
!do not go past the unification point with MHIGH(UP)>MHIGH(DOWN)
!Fixed so that spacing (DeltaT) between points is constant.
!
      IF(UNI.EQ.1)THEN
        IF(NLTMT.NE.20)THEN
          NSTEP=NSTEP*LOG(QTHSORT(20)/1.D19)/LOG(QTHSORT(20)/MHIGH)
        ELSE
          NSTEP=NSTEP*LOG(MT/1.D19)/LOG(MT/MHIGH)
        END IF
      END IF
!
!The next line sets the upper scale so high that the couplings will
!unify before the limit is reached
!
      IF(UNI.EQ.1)MHIGH=1.D19
      TZ=LOG(MZ/MHIGH)
      TMT=LOG(MT/MHIGH)
      DO I=1,20
        TTH(I)=LOG(QTHSORT(I)/MHIGH)
      END DO
      THIGH=0.D0
!
!Run from M_Z to M_t
!
      NSTEPSM=50
      DT=(TMT-TZ)/FLOAT(NSTEPSM)
      NU=3 !top decouples at M_Z
      GSM(1)=G1MZ
      GSM(2)=G2MZ
      GSM(3)=G3MZ
      GSM(4)=DCMPLX(MWEAK(1))/VSMMZ
      GSM(8)=DCMPLX(MWEAK(2))/VSMMZ
      GSM(12)=LAMTMZ
      GSM(13)=DCMPLX(MWEAK(3))/VSMMZ
      GSM(17)=DCMPLX(MWEAK(4))/VSMMZ
      GSM(21)=DCMPLX(LAMBMZ)
      GSM(22)=DCMPLX(MWEAK(5))/VSMMZ
      GSM(26)=DCMPLX(MWEAK(6))/VSMMZ
      GSM(30)=DCMPLX(LAMTAMZ)
      GSM(31)=LAMBDAMZ
      GSM(32)=VSMMZ
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          DGSM(I)=DBLE(GSM(I))
        END DO
      END IF
!
      DO II=1,NSTEPSM
        T=TZ+(TMT-TZ)*FLOAT(II-1)/FLOAT(NSTEPSM)
        SMQSTEP=MHIGH*EXP(T)
        IF(II.EQ.1)EPS=ABS(SMQSTEP*(EXP(DT)-1)/(6.D0*PI))
        IF(COMP.EQ.0)THEN
          CALL DRKSTP(32,DT,T,DGSM,DSMRGEDR,DWSM)
        ELSE
          CALL CRKSTP(32,DT,T,GSM,CSMRGEDR,WSM)
        END IF
        EPS=-ABS(SMQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
      IF(COMP.EQ.0)THEN
        DO I=1,32
          GSM(I)=DCMPLX(DGSM(I))
        END DO
      END IF
!
!Now introduce the top yukawa and rotate
!
      GSM(12)=DCMPLX(LAMTMT)
      CALL ROTATESM(GSM)
!
!Place SM couplings at M_T into G
!
      DO I=4,30
        IF(I.LT.7)G(I-3)=GSM(I-3)
        G(I+108)=GSM(I)
      END DO
      G(428)=GSM(32)
      G(429)=GSM(31)
!
!Now run the rest of the way past the thresholds with LAMBDA_t
!
      IF(NLTMT.GE.LOCMH)CALL UPMHCOND2
!
!First threshold
!
      IF(NLTMT.NE.20)THEN
        IF(NLTMT.NE.0)THEN
          DT=(TTH(NLTMT+1)-TMT)/FLOAT(NSTEPTHRESH(NLTMT))
          LOOPNSTEP=NSTEPTHRESH(NLTMT)
        ELSE
          NSTEPMT=INT(ABS(DLOG(MT/QTHSORT(1)))*NSTEP/DLOG(MHIGH/MT)*25)
          DT=(TTH(NLTMT+1)-TMT)/FLOAT(NSTEPMT)
          LOOPNSTEP=NSTEPMT
        END IF
!
        DO II=1,LOOPNSTEP
          T=TMT+(TTH(NLTMT+1)-TMT)*FLOAT(II-1)/FLOAT(LOOPNSTEP)
          SSQSTEP=MHIGH*EXP(T)
!
          IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
            BELOWMS=0
            CALL ROTBACK(1)
            G(12)=G(12)/(1.D0-DBLE(RTISA))
            G(21)=G(21)/(1.D0-DBLE(RBISA))
            G(30)=G(30)/(1.D0-DBLE(RLISA))
            IF(RGEMS.LE.QNH)THEN
              G(120)=G(120)/(1.D0-DBLE(RTISA))
              G(129)=G(129)/(1.D0-DBLE(RBISA))
              G(138)=G(138)/(1.D0-DBLE(RLISA))
            END IF
            CALL ROTATE(1)
          END IF
!
          IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
!
!Check the sfermion thresholds
!
          CALL CHINT(SSQSTEP,G)
!          CALL OUTCOUP(MHIGH*EXP(T),1)
!
          IF(COMP.EQ.0)THEN
            DO I=1,601
              DG(I)=DBLE(G(I))
            END DO
            CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
            DO I=1,601
              G(I)=DCMPLX(DG(I))
            END DO
          ELSE
            CALL CRKSTP(601,DT,T,G,CRGE601,W)
          END IF
          EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        END DO
!
!Matching conditions at M_H
!
        IF(NLTMT+1.EQ.LOCMH)CALL UPMHCOND2
!
        DO I=NLTMT+2,20
!
!Don't run between degenerate thresholds
!
          IF(NSTEPTHRESH(I-1).EQ.0)GOTO 20
!
          DT=(TTH(I)-TTH(I-1))/FLOAT(NSTEPTHRESH(I-1))
!          WRITE(*,*)'RUNNING BETWEEN',QTHSORT(I-1),QTHSORT(I)
!
          DO II=1,NSTEPTHRESH(I-1)
            T=TTH(I-1)+(TTH(I)-TTH(I-1))
     $                         *FLOAT(II-1)/FLOAT(NSTEPTHRESH(I-1))
            SSQSTEP=MHIGH*EXP(T)
!
            IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
              BELOWMS=0
              CALL ROTBACK(1)
              G(12)=G(12)/(1.D0-DBLE(RTISA))
              G(21)=G(21)/(1.D0-DBLE(RBISA))
              G(30)=G(30)/(1.D0-DBLE(RLISA))
              IF(RGEMS.LE.QNH)THEN
                G(120)=G(120)/(1.D0-DBLE(RTISA))
                G(129)=G(129)/(1.D0-DBLE(RBISA))
                G(138)=G(138)/(1.D0-DBLE(RLISA))
              END IF
              CALL ROTATE(1)
            END IF
!
            IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
            CALL CHINT(SSQSTEP,G)
!            CALL OUTCOUP(MHIGH*EXP(T),1)
            IF(COMP.EQ.0)THEN
              DO J=1,601
                DG(J)=DBLE(G(J))
              END DO
              CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
              DO J=1,601
                G(J)=DCMPLX(DG(J))
              END DO
            ELSE
              CALL CRKSTP(601,DT,T,G,CRGE601,W)
            END IF
            EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
          END DO
!
   20     IF(LOCMH.EQ.I)CALL UPMHCOND2
        END DO
      END IF
!
!Set the g-tilde terms equal to ( the usual g x identity )
!This prevents numerical accidents which may cause spurious
!results.
!
      DO I=1,3
        DO J=1,3
          G(138+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(147+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(156+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(165+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(174+(I-1)*3+J)=G(1)*DSQRT(3.D0/5.D0)*CID(I,J)
          G(185+(I-1)*3+J)=G(2)*CID(I,J)
          G(194+(I-1)*3+J)=G(2)*CID(I,J)
          G(205+(I-1)*3+J)=G(3)*CID(I,J)
          G(214+(I-1)*3+J)=G(3)*CID(I,J)
          G(223+(I-1)*3+J)=G(3)*CID(I,J)
          G(429+(I-1)*3+J)=G(429+(I-1)*3+J)*DSQRT(
     $                     ABS(G(61)+G(62)-G(351)-G(352))
     $                     /(2.D0*ABS(G(108)**2)))
          G(438+(I-1)*3+J)=G(438+(I-1)*3+J)*DSQRT(
     $                     ABS(G(61)+G(62)-G(351)-G(352))
     $                     /(2.D0*ABS(G(108)**2)))
          G(447+(I-1)*3+J)=G(447+(I-1)*3+J)*DSQRT(
     $                     ABS(G(61)+G(62)-G(351)-G(352))
     $                     /(2.D0*ABS(G(108)**2)))
        END DO
      END DO
      IF(ABS(G(61)+G(62)-G(351)-G(352)).GT.0.D0)THEN
        G(108)=G(108)*DSQRT(ABS(G(61)+G(62)-G(351)-G(352))
     $               /(2.D0*ABS(G(108)**2)))
      ELSE
!
!Fix mu to Mz if its new squared value would have been negative
!
        G(108)=G(108)*MZ/ABS(G(108))
      END IF
!
!Continue running to m_high. I need to be careful about the case that
!all thresholds are below m_t.
!
      IF(NLTMT.NE.20)THEN
        DT=(THIGH-TTH(20))/FLOAT(NSTEP)
      ELSE
        DT=(THIGH-TMT)/FLOAT(NSTEP)
      END IF
!
      DO II=1,NSTEP
        IF(NLTMT.NE.20)THEN
          T=TTH(20)+(THIGH-TTH(20))*FLOAT(II-1)/FLOAT(NSTEP)
        ELSE
          T=TMT+(THIGH-TMT)*FLOAT(II-1)/FLOAT(NSTEP)
        END IF
        SSQSTEP=MHIGH*EXP(T)
        Q=MHIGH*EXP(T)
        IF(COMP.EQ.0)THEN
          A1I=4.D0*PI/DG(1)**2
          A2I=4.D0*PI/DG(2)**2
        ELSE
          A1I=4.D0*PI/DBLE(G(1)**2)
          A2I=4.D0*PI/DBLE(G(2)**2)
        END IF
        IF(A1I.LT.A2I.AND.UNI.EQ.1)THEN
          MHIGH=Q
          GO TO 30
        END IF
!
        IF(BELOWMS.EQ.1.AND.MHIGH*EXP(T+DT).GE.RGEMS)THEN
          BELOWMS=0
          CALL ROTBACK(1)
          G(12)=G(12)/(1.D0-DBLE(RTISA))
          G(21)=G(21)/(1.D0-DBLE(RBISA))
          G(30)=G(30)/(1.D0-DBLE(RLISA))
          IF(RGEMS.LE.QNH)THEN
            G(120)=G(120)/(1.D0-DBLE(RTISA))
            G(129)=G(129)/(1.D0-DBLE(RBISA))
            G(138)=G(138)/(1.D0-DBLE(RLISA))
          END IF
          CALL ROTATE(1)
        END IF
!
        IF(II.EQ.1)EPS=ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
        CALL CHINT(SSQSTEP,G)
!        CALL OUTCOUP(MHIGH*EXP(T),1)
        IF(COMP.EQ.0)THEN
          DO I=1,601
            DG(I)=DBLE(G(I))
          END DO
          CALL DRKSTP(601,DT,T,DG,DRGE601,DW)
          DO I=1,601
            G(I)=DCMPLX(DG(I))
          END DO
        ELSE
          CALL CRKSTP(601,DT,T,G,CRGE601,W)
        END IF
        EPS=-ABS(SSQSTEP*(EXP(DT)-1)/(6.D0*PI))
      END DO
!
!NB: If this message is written, the programme will probably
!    not succeed. It may be better to put a stop statement here.
!
      IF(UNI.EQ.1)WRITE(*,*)'ERROR: UNIFICATION NOT FOUND', 
     $' IN ITERATIVE SECTION'
      MHIGH=EXP(T)*MHIGH
 30   CONTINUE
!
!      CLOSE(52)
!      CLOSE(53)
!      CLOSE(54)
!      CLOSE(55)
!      CLOSE(56)
!      CLOSE(57)
!      CLOSE(58)
!      CLOSE(59)
!      CLOSE(60)
!      CLOSE(61)
!      CLOSE(62)
!      CLOSE(63)
!      CLOSE(64)
!      CLOSE(65)
!      CLOSE(66)
!      CLOSE(67)
!      CLOSE(68)
!      CLOSE(69)
!      CLOSE(70)
!      CLOSE(71)
!      CLOSE(72)
!      CLOSE(73)
!      CLOSE(74)
!      CLOSE(75)
!      CLOSE(76)
!      CLOSE(77)
!      CLOSE(78)
!      CLOSE(79)
!      CLOSE(80)
!      CLOSE(81)
!      CLOSE(82)
!      CLOSE(83)
!      CLOSE(84)
!      CLOSE(85)
!      CLOSE(86)
!      CLOSE(87)
!      CLOSE(88)
!      CLOSE(89)
!      CLOSE(90)
!      CLOSE(91)
!      CLOSE(92)
!      CLOSE(93)
!      CLOSE(94)
!      CLOSE(95)
!      CLOSE(96)
!      CLOSE(97)
!      CLOSE(98)
!      CLOSE(99)
!
      RETURN
      END
!
      SUBROUTINE UPSQM(G,Q,OUT)
!
!Purpose: Construct the 6x6 mass matrix for up-type squarks.
!         The basis is u_l,c_l,t_l,u_r,c_r,t_r.
!
      IMPLICIT NONE
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE COMPLEX G(601),OUT(6,6)
!
      DOUBLE COMPLEX AU(3,3),MQ(3,3),MU(3,3),MTSFU(3,3)
      DOUBLE COMPLEX FUQ(3,3),FUUR(3,3)
      DOUBLE COMPLEX GLP,G2L
      DOUBLE COMPLEX AUT(3,3),MTSFUT(3,3),FUQT(3,3),FUURT(3,3)
      DOUBLE COMPLEX AUS(3,3),MTSFUS(3,3),FUQS(3,3),FUURS(3,3)
      DOUBLE COMPLEX TRIU(3,3),SFUQ(3,3),SFUUR(3,3)
      DOUBLE COMPLEX TRIUT(3,3),SFUQT(3,3),SFUURT(3,3)
      DOUBLE COMPLEX TRIUS(3,3),SFUQS(3,3),SFUURS(3,3)
      DOUBLE COMPLEX ID(3,3),CMATMUL
      DOUBLE PRECISION TANB,SINB,COSB,Q,VU,VD,VSM,VUSQ,VDSQ,VSMSQ
      INTEGER I,J
!
      DATA ID(1,1)/(1.D0,0.D0)/,ID(1,2)/(0.D0,0.D0)/
     $     ,ID(1,3)/(0.D0,0.D0)/
      DATA ID(2,1)/(0.D0,0.D0)/,ID(2,2)/(1.D0,0.D0)/
     $     ,ID(2,3)/(0.D0,0.D0)/
      DATA ID(3,1)/(0.D0,0.D0)/,ID(3,2)/(0.D0,0.D0)/
     $     ,ID(3,3)/(1.D0,0.D0)/
!
      DO I=1,6
        DO J=1,6
          OUT(I,J)=(0.D0,0.D0)
        END DO
      END DO
!
!Insert the soft mass matrices in the quark mass basis
!
      DO I=1,3
        DO J=1,3
          MQ(I,J)=MQQMASS(I,J)
          MU(I,J)=MUQMASS(I,J)
        END DO
      END DO
!
!Calculate cos (beta) and sin (beta)
!
      VU=DBLE(G(110))
      VD=DBLE(G(111))
      TANB=VU/VD
      COSB=DSQRT(1.D0/(1.D0+TANB**2))
      SINB=TANB*COSB
      VUSQ=VU**2
      VDSQ=VD**2
!
      VSM=DBLE(G(428))
      VSMSQ=VSM**2
!
!Convert the other entries in G(601) to the matrices needed
!
      DO I=1,3
        DO J=1,3
          AU(I,J)=G(33+(I-1)*3+J)
          MTSFU(I,J)=G(429+(I-1)*3+J)
          FUQ(I,J)=G(3+(I-1)*3+J) !Quartics are not run independently
          FUUR(I,J)=G(3+(I-1)*3+J)
!
          TRIU(I,J)=G(399+(I-1)*3+J)
          SFUQ(I,J)=G(111+(I-1)*3+J)
          SFUUR(I,J)=G(111+(I-1)*3+J)
        END DO
      END DO
!
!Quartics are not run independently at this time.
!Set other quartics to be their non-tilde counterparts.
!This can be fixed if the quartic running is introduced.
!*****NB: IF THE QUARTIC RUNNING IS INTRODUCED CARE MUST
!         BE TAKEN WITH G(541) AND G(542) SINCE THEY ARE
!         DEFINED DIFFERENTLY IN THE COMPLEX VERSION AS
!         OPPOSED TO THE REAL VERSION. SEE THE NOTE AT
!         THE BEGINNING OF drge601.f
!
      GLP=DSQRT(3.D0/5.D0)*G(1)
      G2L=G(2)
      G(541)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*DSQRT(3.D0/5.D0)*G(1)
      G(542)=SQRT(DCMPLX(COSB**2)-DCMPLX(SINB**2))*G(2)
!
      DO I=1,3
        DO J=1,3
          AUT(I,J)=AU(J,I)
          MTSFUT(I,J)=MTSFU(J,I)
          FUQT(I,J)=FUQ(J,I)
          FUURT(I,J)=FUUR(J,I)
          AUS(I,J)=CONJG(AU(I,J))
          MTSFUS(I,J)=CONJG(MTSFU(I,J))
          FUQS(I,J)=CONJG(FUQ(I,J))
          FUURS(I,J)=CONJG(FUUR(I,J))
!
          TRIUT(I,J)=TRIU(J,I)
          SFUQT(I,J)=SFUQ(J,I)
          SFUURT(I,J)=SFUUR(J,I)
          TRIUS(I,J)=CONJG(TRIU(I,J))
          SFUQS(I,J)=CONJG(SFUQ(I,J))
          SFUURS(I,J)=CONJG(SFUUR(I,J))
        END DO
      END DO
!
!Split matrix into 3x3 blocks.
!First, the top left block
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J)=VUSQ*CMATMUL(0,FUURS,FUURT,I,J)+MQ(I,J)
     $               +(GLP**2/12.D0-G2L**2/4.D0)*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I,J)=VSMSQ*CMATMUL(0,SFUURS,SFUURT,I,J)+MQ(I,J)
     $               -VSMSQ*(G(541)**2/12.D0-G(542)**2/4.D0)*ID(I,J)
          END IF
        END DO
      END DO
!
!Next the bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J)=VD*MTSFUT(I,J)-VU*AUT(I,J)
          ELSE
            OUT(I+3,J)=-(VSM*TRIUT(I,J))
          END IF
        END DO
      END DO
!
!Top right is the dagger of bottom left
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I,J+3)=VD*MTSFUS(I,J)-VU*AUS(I,J)
          ELSE
            OUT(I,J+3)=-(VSM*TRIUS(I,J))
          END IF
        END DO
      END DO
!
!Finally bottom right
!
      DO I=1,3
        DO J=1,3
          IF(Q.GT.QNH)THEN
            OUT(I+3,J+3)=VUSQ*CMATMUL(0,FUQT,FUQS,I,J)+MU(I,J)
     $                   -GLP**2/3.D0*(VUSQ-VDSQ)*ID(I,J)
          ELSE
            OUT(I+3,J+3)=VSMSQ*CMATMUL(0,SFUQT,SFUQS,I,J)+MU(I,J)
     $                   +VSMSQ*G(541)**2/3.D0*ID(I,J)
          END IF
        END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE USMMA(G,Q,SVLQ)
!
!Purpose: Find the up squark mass matrix, sort the eigenvectors and
!         eigenvalues and then (if required) write them to a file
!         called sqm2u.dat
!
      IMPLICIT NONE
!
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,STRADD,SQM2U
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
      DOUBLE PRECISION SUM,Q
      DOUBLE COMPLEX G(601),CUSQM(6,6),TEMP(7)
      DOUBLE COMPLEX EVERTMP(6,6),CWORK1(99),CWORK2(12)
      INTEGER I,J,K,CIERR,SVLQ
!
      DO I=1,3
        DO J=1,3
          T1EVE(I,J)=(0.D0,0.D0)
        END DO
        T1EVA(I)=(0.D0,0.D0)
      END DO
      COSTHT=(0.D0,0.D0)
      SINTHT=(0.D0,0.D0)
      MST1=0.D0
      MST2=0.D0
!
      CALL UPSQM(G,Q,USQM)
      DO I=1,6
        DO J=1,6
          CUSQM(I,J)=USQM(I,J)
        END DO
      END DO
      SQM2U=STRADD(FNRGE,'.sqm2u')
      OPEN(25,FILE=SQM2U,STATUS='UNKNOWN',FORM='FORMATTED')
!
      WRITE(25,*)
      WRITE(25,*)'THE 6X6 UP-TYPE MATRIX AT THE SCALE: ',Q
      WRITE(25,*)'IN THE CHOSEN BASIS, IS:'
!
      WRITE(25,*)
      WRITE(25,15)'u_l','c_l','t_l','u_r','c_r','t_r'
      WRITE(25,*)
!
      WRITE(25,11)'u_l',USQM(1,1),USQM(1,2),USQM(1,3),
     $            USQM(1,4),USQM(1,5),USQM(1,6)
      WRITE(25,11)'c_l',USQM(2,1),USQM(2,2),USQM(2,3),
     $            USQM(2,4),USQM(2,5),USQM(2,6)
      WRITE(25,11)'t_l',USQM(3,1),USQM(3,2),USQM(3,3),
     $            USQM(3,4),USQM(3,5),USQM(3,6)
      WRITE(25,11)'u_r',USQM(4,1),USQM(4,2),USQM(4,3),
     $            USQM(4,4),USQM(4,5),USQM(4,6)
      WRITE(25,11)'c_r',USQM(5,1),USQM(5,2),USQM(5,3),
     $            USQM(5,4),USQM(5,5),USQM(5,6)
      WRITE(25,11)'t_r',USQM(6,1),USQM(6,2),USQM(6,3),
     $            USQM(6,4),USQM(6,5),USQM(6,6)
!
      WRITE(25,*)
!
      IF(SVLQ.EQ.1)THEN
!
!Diagonalise the matrix using the LAPACK routine ZGEEV.
!
        CALL ZGEEV('V','V',6,CUSQM,6,T1EVA,T1EVE,6,EVERTMP,6,CWORK1,99
     $             ,CWORK2,CIERR)
!
        DO J=1,6
          SUM=0.D0
          DO I=1,6
            SUM=SUM+DBLE(T1EVE(I,J)*CONJG(T1EVE(I,J)))
          END DO
          DO I=1,6
            T1EVE(I,J)=T1EVE(I,J)/DSQRT(SUM)
          END DO
        END DO
!
!Find the eigenvectors with the most amount of stop
!
        DO I=1,5
          DO J=I+1,6
            IF(DSQRT(DBLE(T1EVE(3,I)*CONJG(T1EVE(3,I)))+DBLE(T1EVE(6,I)
     $              *CONJG(T1EVE(6,I)))).LT.
     $         DSQRT(DBLE(T1EVE(3,J)*CONJG(T1EVE(3,J)))+DBLE(T1EVE(6,J)
     $              *CONJG(T1EVE(6,J)))))THEN
              DO K=1,6
                TEMP(K)=T1EVE(K,I)
                T1EVE(K,I)=T1EVE(K,J)
                T1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=T1EVA(I)
              T1EVA(I)=T1EVA(J)
              T1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Now find next two with most amount of scalar charm
!
        DO I=3,5
          DO J=I+1,6
            IF(DSQRT(DBLE(T1EVE(2,I)*CONJG(T1EVE(2,I)))+DBLE(T1EVE(5,I)
     $              *CONJG(T1EVE(5,I)))).LT.
     $         DSQRT(DBLE(T1EVE(2,J)*CONJG(T1EVE(2,J)))+DBLE(T1EVE(5,J)
     $              *CONJG(T1EVE(5,J)))))THEN
              DO K=1,6
                TEMP(K)=T1EVE(K,I)
                T1EVE(K,I)=T1EVE(K,J)
                T1EVE(K,J)=TEMP(K)
              END DO
              TEMP(7)=T1EVA(I)
              T1EVA(I)=T1EVA(J)
              T1EVA(J)=TEMP(7)
            END IF
          END DO
        END DO
!
!Sort the first two eigenvectors by mass
!
        IF(ABS(T1EVA(1)).GT.ABS(T1EVA(2)))THEN
          DO I=1,6
            TEMP(I)=T1EVE(I,1)
            T1EVE(I,1)=T1EVE(I,2)
            T1EVE(I,2)=TEMP(I)
          END DO
          TEMP(7)=T1EVA(1)
          T1EVA(1)=T1EVA(2)
          T1EVA(2)=TEMP(7)
        END IF
!
!Sort the next two by mass
!
        IF(ABS(T1EVA(3)).GT.ABS(T1EVA(4)))THEN
          DO I=1,6
            TEMP(I)=T1EVE(I,3)
            T1EVE(I,3)=T1EVE(I,4)
            T1EVE(I,4)=TEMP(I)
          END DO
          TEMP(7)=T1EVA(3)
          T1EVA(3)=T1EVA(4)
          T1EVA(4)=TEMP(7)
        END IF
!
!And the final two
!
        IF(ABS(T1EVA(5)).GT.ABS(T1EVA(6)))THEN
          DO I=1,6
            TEMP(I)=T1EVE(I,5)
            T1EVE(I,5)=T1EVE(I,6)
            T1EVE(I,6)=TEMP(I)
          END DO
          TEMP(7)=T1EVA(5)
          T1EVA(5)=T1EVA(6)
          T1EVA(6)=TEMP(7)
        END IF
!
!Note about sin{\theta_t}: If I define sin{\theta_t} the same way as
!BT, then -sin{\theta_t}=conjg(t1eve(6,1))
!
        COSTHT=CONJG(T1EVE(3,1))
        SINTHT=-CONJG(T1EVE(6,1))
        MST1=DSQRT(ABS(T1EVA(1)))
        MST2=DSQRT(ABS(T1EVA(2)))
!
        WRITE(25,*)'EIGENVECTORS ARE:'
        WRITE(25,*)
        WRITE(25,15)'t_1','t_2','c_1','c_2','u_1','u_2'
        WRITE(25,*)
!     
        WRITE(25,11)'u_l',T1EVE(1,1),T1EVE(1,2),T1EVE(1,3),
     $              T1EVE(1,4),T1EVE(1,5),T1EVE(1,6)
        WRITE(25,11)'c_l',T1EVE(2,1),T1EVE(2,2),T1EVE(2,3),
     $              T1EVE(2,4),T1EVE(2,5),T1EVE(2,6)
        WRITE(25,11)'t_l',T1EVE(3,1),T1EVE(3,2),T1EVE(3,3),
     $              T1EVE(3,4),T1EVE(3,5),T1EVE(3,6)
        WRITE(25,11)'u_r',T1EVE(4,1),T1EVE(4,2),T1EVE(4,3),
     $              T1EVE(4,4),T1EVE(4,5),T1EVE(4,6)
        WRITE(25,11)'c_r',T1EVE(5,1),T1EVE(5,2),T1EVE(5,3),
     $              T1EVE(5,4),T1EVE(5,5),T1EVE(5,6)
        WRITE(25,11)'t_r',T1EVE(6,1),T1EVE(6,2),T1EVE(6,3),
     $              T1EVE(6,4),T1EVE(6,5),T1EVE(6,6)
!
        WRITE(25,*)
        WRITE(25,*)'EIGENVALUES ARE:'
        WRITE(25,16)T1EVA(1),T1EVA(2),T1EVA(3),T1EVA(4),T1EVA(5)
     $             ,T1EVA(6)
        WRITE(25,*)
        WRITE(25,*)'MASSES ARE:'
        WRITE(25,17)DSQRT(DSQRT(DBLE(T1EVA(1)*CONJG(T1EVA(1)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(2)*CONJG(T1EVA(2)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(3)*CONJG(T1EVA(3)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(4)*CONJG(T1EVA(4)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(5)*CONJG(T1EVA(5)))))
     $             ,DSQRT(DSQRT(DBLE(T1EVA(6)*CONJG(T1EVA(6)))))
      END IF
!
!
 11   FORMAT(SP,1P,A,3X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 15   FORMAT(16X,A,21X,A,21X,A,21X,A,21X,A,21X,A)
 16   FORMAT(SP,1P,6X,'(',D10.3,',',D10.3,')',1X,'(',D10.3,',',
     $       D10.3,')',1X,'(',D10.3,',',D10.3,')',1X,'(',
     $       D10.3,',',D10.3,')',1X,'(',D10.3,',',D10.3,')'
     $      ,1X,'(',D10.3,',',D10.3,')')
 17   FORMAT(10X,F11.4,11X,F11.4,11X,F11.4,11X,F11.4,11X
     $      ,F11.4,11X,F11.4)
      WRITE(25,*)
      CLOSE(25)
!
      RETURN
      END
!
      SUBROUTINE VGEN(T12,T13,T23,D13,MAT)
!
!Purpose: To generate an arbitrary unitary matrix
!         Is passed the values of alpha (T12), beta (T13),
!         gamma (T23) and delta (D13)
!
      IMPLICIT NONE
!
      INTEGER I,J
      DOUBLE PRECISION T12,C12,S12,T13,C13,S13,T23,C23,S23,D13
      DOUBLE COMPLEX MAT(3,3),MAT1(3,3),MAT2(3,3),MAT3(3,3)
     $              ,TMPMAT(3,3),CMATMUL
!
!Find the values of sin and cos of the angles
!
      C12=DCOS(T12)
      S12=DSIN(T12)
      C13=DCOS(T13)
      S13=DSIN(T13)
      C23=DCOS(T23)
      S23=DSIN(T23)
!
!These are the individual matrices
!
      MAT1(1,1)=(1.D0,0.D0)
      MAT1(1,2)=(0.D0,0.D0)
      MAT1(1,3)=(0.D0,0.D0)
      MAT1(2,1)=(0.D0,0.D0)
      MAT1(2,2)=DCMPLX(C23)
      MAT1(2,3)=DCMPLX(S23)
      MAT1(3,1)=(0.D0,0.D0)
      MAT1(3,2)=-DCMPLX(S23)
      MAT1(3,3)=DCMPLX(C23)
      MAT2(1,1)=DCMPLX(C13)
      MAT2(1,2)=(0.D0,0.D0)
      MAT2(1,3)=S13*DCMPLX(DCOS(D13),DSIN(D13))
      MAT2(2,1)=(0.D0,0.D0)
      MAT2(2,2)=(1.D0,0.D0)
      MAT2(2,3)=(0.D0,0.D0)
      MAT2(3,1)=-S13*DCMPLX(DCOS(D13),-DSIN(D13))
      MAT2(3,2)=(0.D0,0.D0)
      MAT2(3,3)=DCMPLX(C13)
      MAT3(1,1)=DCMPLX(C12)
      MAT3(1,2)=DCMPLX(S12)
      MAT3(1,3)=(0.D0,0.D0)
      MAT3(2,1)=-DCMPLX(S12)
      MAT3(2,2)=DCMPLX(C12)
      MAT3(2,3)=(0.D0,0.D0)
      MAT3(3,1)=(0.D0,0.D0)
      MAT3(3,2)=(0.D0,0.D0)
      MAT3(3,3)=(1.D0,0.D0)
!
      DO I=1,3
        DO J=1,3
          TMPMAT(I,J)=CMATMUL(0,MAT1,MAT2,I,J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          MAT(I,J)=CMATMUL(0,TMPMAT,MAT3,I,J)
        END DO
      END DO
!
      RETURN
      END
