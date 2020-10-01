!Open-interval Romberg integration for power-law power spectra
!Adapted from Numerical Recipes.
!Author: Pat Scott (p.scott@imperial.ac.uk)
!Date: Some time in 2010/11.
!Added to DarkSUSY: Oct 13 2014


      double precision function qromo(func,n,alpha,kratio,a,b,choose)
      implicit none
      integer JMAX,JMAXP,K,KM
      double precision n,alpha,kratio,a,b,func,ss,EPS
      external func,choose
      parameter (EPS=1.d-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      integer j
      double precision dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,n,alpha,kratio,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) then
            qromo=ss
            return
          endif
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.d0
11    continue
      write(*,*) 'too many steps in qromo'
      stop ! TB, removed pause
      end

      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      integer n,NMAX
      double precision dy,x,y,xa(n),ya(n)
      parameter (NMAX=10)
      integer i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
            write(*,*) 'failure in polint'
            stop ! removed pause
          endif 
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
