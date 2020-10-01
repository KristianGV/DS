#include "PILOT.inc"
        REAL FUNCTION SSXINT(XL,F,XR)
C-----------------------------------------------------------------------
C          Integrate F over (XL,XR) using adaptive Gaussian quadrature.
C          TOLABS = 1e-35: absolute error for convergence.
C          TOLREL = 5e-5:  relative error for convergence.
C          TOLBIN = 1e-3:  relative bin size limit. Set contribution to
C                          zero if bin falls below this.
C-----------------------------------------------------------------------
#ifdef IMPNONE_X
      IMPLICIT NONE
#endif
#include "sslun.inc"
        EXTERNAL F
        INTEGER NMAX
        REAL TOLABS,TOLREL,TOLBIN,XMIN,XLIMS(100)
        REAL R(93),W(93)
        INTEGER PTR(4),NORD(4)
        INTEGER ICOUNT,IBAD
        REAL XL,XR,F
        REAL AA,BB,TVAL,VAL,TOL
        INTEGER NLIMS,I,J
C
        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860,.6612093865,.9324695142,
     1 .1252334085,.3678314990,.5873179543,.7699026742,.9041172563,
     1 .9815606342,.0640568929,.1911188675,.3150426797,.4337935076,
     1 .5454214714,.6480936519,.7401241916,.8200019860,.8864155270,
     1 .9382745520,.9747285560,.9951872200,.0323801710,.0970046992,
     1 .1612223561,.2247637903,.2873624873,.3487558863,.4086864820,
     1 .4669029048,.5231609747,.5772247261,.6288673968,.6778723796,
     1 .7240341309,.7671590325,.8070662040,.8435882616,.8765720203,
     1 .9058791367,.9313866907,.9529877032,.9705915925,.9841245837,
     1 .9935301723,.9987710073,.0162767488,.0488129851,.0812974955,
     1 .1136958501,.1459737146,.1780968824,.2100313105,.2417431561,
     1 .2731988126,.3043649444,.3352085229,.3656968614,.3957976498,
     1 .4254789884,.4547094222,.4834579739,.5116941772,.5393881083,
     1 .5665104186,.5930323648,.6189258401,.6441634037,.6687183100,
     1 .6925645366,.7156768123,.7380306437,.7596023411,.7803690438,
     1 .8003087441,.8194003107,.8376235112,.8549590334,.8713885059,
     1 .8868945174,.9014606353,.9150714231,.9277124567,.9393703398,
     1 .9500327178,.9596882914,.9683268285,.9759391746,.9825172636,
     1 .9880541263,.9925439003,.9959818430,.9983643759,.9996895039/
        DATA W/.4679139346,.3607615730,.1713244924,
     1 .2491470458,.2334925365,.2031674267,.1600783285,.1069393260,
     1 .0471753364,.1279381953,.1258374563,.1216704729,.1155056681,
     1 .1074442701,.0976186521,.0861901615,.0733464814,.0592985849,
     1 .0442774388,.0285313886,.0123412298,.0647376968,.0644661644,
     1 .0639242386,.0631141923,.0620394232,.0607044392,.0591148397,
     1 .0572772921,.0551995037,.0528901894,.0503590356,.0476166585,
     1 .0446745609,.0415450829,.0382413511,.0347772226,.0311672278,
     1 .0274265097,.0235707608,.0196161605,.0155793157,.0114772346,
     1 .0073275539,.0031533461,.0325506145,.0325161187,.0324471637,
     1 .0323438226,.0322062048,.0320344562,.0318287589,.0315893308,
     1 .0313164256,.0310103326,.0306713761,.0302999154,.0298963441,
     1 .0294610900,.0289946142,.0284974111,.0279700076,.0274129627,
     1 .0268268667,.0262123407,.0255700360,.0249006332,.0242048418,
     1 .0234833991,.0227370697,.0219666444,.0211729399,.0203567972,
     1 .0195190811,.0186606796,.0177825023,.0168854799,.0159705629,
     1 .0150387210,.0140909418,.0131282296,.0121516047,.0111621020,
     1 .0101607705,.0091486712,.0081268769,.0070964708,.0060585455,
     1 .0050142027,.0039645543,.0029107318,.0018539608,.0007967921/
C
      DATA TOLABS,TOLREL,TOLBIN,NMAX/1.E-35,5.E-5,1E-3,100/
C
      SSXINT=0
      NLIMS=2
      XLIMS(1)=XL
      XLIMS(2)=XR
      ICOUNT = 0
      IBAD=0
      XMIN=TOLBIN*ABS(XR-XL)
C
10    AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2
      BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2
      TVAL=0
      DO 15 I=1,3
15    TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
      TVAL=TVAL*AA
      DO 25 J=1,4
        VAL=0
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
          ICOUNT = ICOUNT + 1
          IF (ICOUNT .GT. 1E5) THEN
            WRITE(LOUT,*) 'ERROR IN SSXINT: SET SSXINT = 0'
            SSXINT=0
            RETURN
          ENDIF
20      VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
          SSXINT=SSXINT+VAL
          NLIMS=NLIMS-2
          IF (NLIMS.NE.0) GO TO 10
          GO TO 999
        ELSEIF(ABS(AA).LT.XMIN.AND.J.EQ.4) THEN
C           Bin is too small -- set VAL = 0. -- FEP
          IBAD=IBAD+1
          NLIMS=NLIMS-2
          IF (NLIMS.NE.0) GO TO 10
          GO TO 999
        ENDIF
25    TVAL=VAL
      IF (NMAX.EQ.2) THEN
        SSXINT=VAL
        GO TO 999
      END IF
      IF (NLIMS.GT.(NMAX-2)) THEN
        WRITE(LOUT,50) SSXINT,NMAX,BB-AA,BB+AA
50      FORMAT (' ERROR IN SSXINT, SSXINT,NMAX,XL,XR=',G15.7,I5,2G15.7)
        RETURN
      END IF
      XLIMS(NLIMS+1)=BB
      XLIMS(NLIMS+2)=BB+AA
      XLIMS(NLIMS)=BB
      NLIMS=NLIMS+2
      GO TO 10
C
999   CONTINUE
C999   IF(IBAD.GT.0) THEN
c        WRITE(LOUT,*) 'WARNING IN SSXINT: BAD CONVERGENCE FOR ',IBAD,
c     $  ' INTERVALS.'
C        WRITE(LOUT,*) ' '
C      ENDIF
      RETURN
      END
