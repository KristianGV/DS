
c_______________________________________________________________________
c  Function dsfistat.
c
c  input:
c    P=p_Y/T  -  Momentum of particle Y / temp
c    X=x_Y/T
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-14
c=======================================================================
      real*8 function dsfistat(P,X,x1,x2,eta1,eta2,stat)
      implicit none
      integer stat
      real*8 P,X,x1,x2,eta1,eta2,E,E1,E2,pcm,dsfipcm,C1,C2,ex1,ex2,ex3
      real*8 num,denom

      pcm=dsfipcm(X,x1,x2)
      E=sqrt(P*P+X*X)
      E1=sqrt(pcm*pcm+x1*x1)
      E2=sqrt(pcm*pcm+x2*x2)
      C1=E/X
      C2=P/X
      ex1=eta1*exp(-C1*E1)
      ex2=eta2*exp(-C1*E2)
      ex3=exp(-C2*pcm)

      num=1+log((1-ex1*ex3)*(1-ex2*ex3)/(1-ex1/ex3)/(1-ex2/ex3))
     &/(2*pcm*C2)
      denom=1-eta1*eta2*exp(-E)

c     REMEBER TO CHECK THAT EXP AND LOG IS WELL BEHAVED!!!
      if(stat.eq.0 .or. (eta1.eq.0 .and. eta2.eq.0)) then
            dsfistat=1
      elseif(pcm*C2 .lt. 1E-5*(x1+x2)*C1) then
            dsfistat=(1+ex1/(1-ex1)+ex2/(1-ex2))/denom
      else
            dsfistat=num/denom
      endif
      return
      end 