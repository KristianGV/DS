C
C
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION H2(X,Y)
C
      REAL*8 X,Y
C
      H2=X*DLOG(X)/(1.D0-X)/(X-Y)+Y*DLOG(Y)/(1.D0-Y)/(Y-X)
      RETURN
      END
