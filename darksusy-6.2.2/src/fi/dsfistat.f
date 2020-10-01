
c_______________________________________________________________________
c  Function dsfistat calculates function that contains all the 
c     statistical mechanical information used in G_Y->a,b
c
c  input:
c     P=p_Y/T  -  Momentum of mediator / temp
c     X=m_Y/T  -  mass of mediator / temp
c     xi=mi/T  -  mass of out particle i/ temp
c     etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c     pluss (minus) for bosons (fermions)
c     stat - determines if quantum corrections is to be considered 
c     (stat=0 for Maxwell-Boltzmann statistics) (stat\neq 0 for FD/BE)
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-14
c=======================================================================
      real*8 function dsfistat(P,X,x1,x2,eta1,eta2)
      implicit none

      include 'dsficom.h'

      real*8 P,X,x1,x2,eta1,eta2,E,E1,E2,pcm,dsfipcm,C1,C2,ex1,ex2,ex3
      real*8 num,denom
c     pcm=p_cm/T
      pcm=dsfipcm(X,x1,x2)
      E=sqrt(P*P+X*X)
      E1=sqrt(pcm*pcm+x1*x1)
      E2=sqrt(pcm*pcm+x2*x2)
      C1=E/X
      C2=P/X
      ex1=eta1*exp(-C1*E1)
      ex2=eta2*exp(-C1*E2)
      ex3=exp(-C2*pcm)

      if(P.gt.1E3) then
            num=1+log(1/((1-(E1/pcm)*(1-E2/pcm))))
      else
            num=1+log((1-ex1*ex3)*(1-ex2*ex3)/(1-ex1/ex3)/(1-ex2/ex3))
     &/(2*pcm*C2)
      end if
      denom=1-eta1*eta2*exp(-E)

c     REMEBER TO CHECK THAT EXP AND LOG IS WELL BEHAVED!!!
c     When stat=0 <=> eta1=eta2=0 => S=1
      if(stat.eq.0 .or. (eta1.eq.0 .and. eta2.eq.0)) then
            dsfistat=1
c     If p_cm*c2<<(x1+x2)*C1, we get 0/0. Therefore L'Hôpitals rule is used
      elseif(pcm*C2 .lt. 1E-5*(x1+x2)*C1) then
            dsfistat=(1+ex1/(1-ex1)+ex2/(1-ex2))/denom
      else
            dsfistat=num/denom
      endif
      return
      end 