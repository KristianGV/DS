c Function wlgama to calculate lograithm of complex gamma function in double precision.
c Taken from CERN program library.
c
c Added by torsten.bringmann@fys.uio.no (as well as all functions it depends on)


*
* $Id: clogam64.F,v 1.1.1.1 1996/04/01 15:01:55 mclareni Exp $
*
* $Log: clogam64.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:55  mclareni
* Mathlib gen
*
*
      FUNCTION WLGAMA(Z)
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* defc64.inc
*
      COMPLEX*16
     +     WLGAMA
     +    ,WLOGAM
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* defc64.inc
*
      COMPLEX*16
C    +     Z,W,U,V,H,P,R,GCONJG,GCMPLX
     +     Z,  U,V,H,P,R
      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'CLGAMA/WLGAMA')
      DIMENSION C(10)
      PARAMETER (Z1 = 1, HF = Z1/2)
*
* $Id: gcmpfun.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: gcmpfun.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* gcmpfun.inc
*
*
* $Id: def64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: def64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* def64.inc
*
      DOUBLE PRECISION
     +      GREAL,GIMAG,XARG,YARG
*
* $Id: defc64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* defc64.inc
*
      COMPLEX*16
     +      ZARG,GCONJG,GCMPLX
      GREAL( ZARG)=DREAL( ZARG)
      GIMAG( ZARG)=DIMAG( ZARG)
      GCONJG(ZARG)=DCONJG(ZARG)
      GCMPLX(XARG,YARG)=DCMPLX(XARG,YARG)
      DATA PI /3.14159 26535 89793 24D+0/
      DATA C1 /9.18938 53320 46727 42D-1/
      DATA C2 /1.14472 98858 49400 17D+0/
      DATA C( 1) / 8.33333 33333 33333 33D-2/
      DATA C( 2) /-2.77777 77777 77777 78D-3/
      DATA C( 3) / 7.93650 79365 07936 51D-4/
      DATA C( 4) /-5.95238 09523 80952 38D-4/
      DATA C( 5) / 8.41750 84175 08417 51D-4/
      DATA C( 6) /-1.91752 69175 26917 53D-3/
      DATA C( 7) / 6.41025 64102 56410 26D-3/
      DATA C( 8) /-2.95506 53594 77124 18D-2/
      DATA C( 9) / 1.79644 37236 88305 73D-1/
      DATA C(10) /-1.39243 22169 05901 12D+0/
C     GREAL(U)=DREAL(U)
C     GIMAG(U)=DIMAG(U)
C     GCONJG(U)=DCONJG(U)
C     GCMPLX(X,Y)=DCMPLX(X,Y)
      ENTRY WLOGAM(Z)
      X=Z
      Y=GIMAG(Z)
      IF(Y .EQ. 0 .AND. -ABS(X) .EQ. INT(X)) THEN
       H=0
       WRITE(ERRTXT,101) X
       CALL MTLPRT(NAME,'C306.1',ERRTXT)
      ELSE
       YA=ABS(Y)
       U=GCMPLX(X,YA)
       IF(X .LT. 0) U=1-U
       H=0
       UR=U
       IF(UR .LT. 7) THEN
        UI=GIMAG(U)
        A=ATAN2(UI,UR)
        H=U
        DO 1 I = 1,6-INT(UR)
        UR=UR+1
        U=GCMPLX(UR,UI)
        H=H*U
    1   A=A+ATAN2(UI,UR)
        H=GCMPLX(HF*LOG(GREAL(H)**2+GIMAG(H)**2),A)
        U=U+1
       ENDIF
       R=1/U**2
       P=R*C(10)
       DO 2 I = 9,2,-1
    2  P=R*(C(I)+P)
       H=C1+(U-HF)*LOG(U)-U+(C(1)+P)/U-H
       IF(X .LT. 0) THEN
        UR=INT(X)-1
        UI=PI*(X-UR)
        X=PI*YA
        T=EXP(-X-X)
        A=SIN(UI)
        T=X+HF*LOG(T*A**2+(HF*(1-T))**2)
        A=ATAN2(COS(UI)*TANH(X),A)-UR*PI
        H=C2-GCMPLX(T,A)-H
       ENDIF
       IF(Y .LT. 0) H=GCONJG(H)
      ENDIF
      WLGAMA=H
      RETURN
  101 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',1P,E15.1)
      END
