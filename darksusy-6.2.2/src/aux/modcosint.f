!Returns modified cosine integral MCi(x) = 1/x^2 * (Ci(x) - lnx - EULER), along with 1/x^2 * Si(x)
!Adapted from Numerical Recipes.
!For x < 0 the routine returns MCi(-x) and the user must supply the -i pi/x^2 him/herself
!si returns a very large value for x=0, instead of crashing or outputting Inf
!Author: Pat Scott (p.scott@imperial.ac.uk)
!Date: Some time in 2010/11.
!Added to DarkSUSY: Oct 11 2014

      subroutine ModCosInt(x,ci,si)
      implicit none
      double precision, intent(IN) :: x
      integer MAXIT,i,k
      double precision ci,si,EPS,PIBY2,FPMIN,TMIN, a,err,fact,sign,sum,sumc,sums,t,term,absc,EULER
      parameter (EPS=6.d-8,MAXIT=100,PIBY2=1.5707963d0,EULER=0.57721566490153286060d0,FPMIN=1.d-300,TMIN=2.d0)
      complex (KIND=8) h,b,c,d,del
      logical odd
      absc(h)=abs(real(h))+abs(aimag(h))
      t=abs(x)
      if(t.eq.0.)then
        si=1.d0/FPMIN
        ci=-0.25d0
        return
      endif
      if(t.gt.TMIN)then
        b=cmplx(1.,t)
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=2,MAXIT
          a=-(i-1.d0)**2
          b=b+2.d0
          d=1.d0/(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(absc(del-1.d0).lt.EPS)goto 1
11      continue
        pause 'cf failed in cisi'
1       continue
        h=cmplx(cos(t),-sin(t))*h
        ci=(-real(h)-log(t)-EULER)/(t*t)
        si=(PIBY2+aimag(h))/(t*t)
      else
        if(t.lt.sqrt(FPMIN))then
          sumc=-0.25d0
          sums=1.d0/t
        else
          sum= 1.d0/t
          sums=1.d0/t
          sumc=-0.25d0
          sign=-1.d0
          fact=0.5d0
          odd=.true.
          do 12 k=3,MAXIT
            fact=fact*t/k
            term=fact/k
            sum=sum+sign*term
            err=term/abs(sum)
            if(odd)then
              sign=-sign
              sums=sum
              sum=sumc
            else
              sumc=sum
              sum=sums
            endif
            if(err.lt.EPS)goto 2
            odd=.not.odd
12        continue
          pause 'maxits exceeded in ModCosInt'
        endif
2       si=sums
        ci=sumc
      endif
      if(x.lt.0.)si=-si
      return

      end subroutine ModCosInt

