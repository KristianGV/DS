************************************************************************
***  dark matter density profiles written in the form:
***
***    rho(r) = rhos * dsdmdprof(r/rs,alpha)
***
***  namely dsdmdprof is the function g(x), e.g., in Eq. (10) of 
***  astro-ph/0207125.
***      
***  NOTE: version valid for a Einasto profile only      
***      
***  input:
***    x [1]
***    alpha [1]: Einasto index      
************************************************************************
      real*8 function dsdmdprof_ein(x,alpha)
      implicit none
      real*8 x,alpha
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprof_ein=dexp(-2.d0/alpha*(x**alpha-1.d0))
      end


************************************************************************
***  first derivative with respect to x of the function dsdmdprof(x) as 
***  implemented in the definition of the density profile:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***      
***  NOTE: version valid for a Einasto profile only      
***      
***  input:
***    x [1]
***    alpha [1]: Einasto index      
************************************************************************
      real*8 function dsdmdprofdx_ein(x,alpha)
      implicit none
      real*8 x,alpha
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprofdx_ein=-2.d0*x**(alpha-1.d0)
     &       *dexp(-2.d0/alpha*(x**alpha-1.d0))
      end


************************************************************************
***  second derivative with respect to x of the function dsdmdprof(x) as 
***  implemented in the definition of the density profile:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***      
***  NOTE: version valid for a Einasto profile only      
***      
***  input:
***    x [1]
***    alpha [1]: Einasto index      
************************************************************************
      real*8 function dsdmd2profdx2_ein(x,alpha)
      implicit none
      real*8 x,alpha
      if(x.lt.1.d-16) x=1.d-16
      dsdmd2profdx2_ein=(-2.d0*(alpha-1.d0)*x**(alpha-2.d0)
     &       +4.d0**x**(2.d0*alpha-2.d0))
     &       *dexp(-2.d0/alpha*(x**alpha-1.d0))
      end
      


************************************************************************
***  volume integral I1 of the dark matter density profile:
***
***    dsI1(x) = \int_0^x dxp xp**2 dsdmprof(xp)
***
***  namely, e.g., Eq. (22) of astro-ph/0207125 with n=1.
***      
***  NOTE: version valid for a Einasto profile only      
***      
***  input:
***    x [1]
***    alpha [1]: Einasto index      
************************************************************************
      real*8 function dsI1_ein(x,alpha)
      implicit none
      real*8 x,alpha,gamma,gammln,gammp,c
      if(x.lt.1.d-16) x=1.d-16
      c=2.d0/alpha
      gamma=dexp(gammln(3.d0/alpha))
      dsI1_ein=1.d0/alpha/c**(3.d0/alpha)*dexp(c)*gamma
     &    *gammp(3.d0/alpha,c*x**alpha)
      end
      
      
************************************************************************
*** auxiliary functions for dsI1_ein, based on NUMERICAL RECIPES: 
************************************************************************
      double precision FUNCTION gammln(xx)
      double precision xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+dlog(stp*ser/x)
      return
      END
ccc
ccc
      double precision FUNCTION gammp(a,x)
      double precision a,x
CU    USES gcf,gser
      double precision gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0) then
        write(*,*) 'DS: bad arguments in gammp, x, a = ',x,a
        stop
      endif
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      endif
      return
      END
ccc
ccc
      double precision FUNCTION gammq(a,x)
      double precision a,x
CU    USES gcf,gser
      double precision gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0) then
        write(*,*) 'DS: bad arguments in gammq, x, a = ',x,a
        stop
      endif
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
ccc
ccc      
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      double precision a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=300,EPS=3.d-16,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      double precision an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i*1.d0-a)
        b=b+2.d0
        d=an*d+b
        if(dabs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(dabs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(dabs(del-1.d0).lt.EPS)goto 1
 11   continue
      write(*,*) 'DS: in gcf, a too large or ITMAX too small'
      write(*,*) 'DS: a, ITMAX = ',a,ITMAX
      stop
 1    gammcf=dexp(-x+a*dlog(x)-gln)*h
      return
      END
ccc
ccc
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      double precision a,gamser,gln,x,EPS
      PARAMETER (ITMAX=300,EPS=3.d-16)
CU    USES gammln
      INTEGER n
      double precision ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0) then
          write(*,*) 'DS: in gser x < 0 = ',x
          stop
        endif
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*EPS)goto 1
 11   continue
      write(*,*) 'DS: in gser, a too large or ITMAX too small'
      write(*,*) 'DS: a, ITMAX = ',a,ITMAX
      stop
 1    gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      END
      
