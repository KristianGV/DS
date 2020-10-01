
c_______________________________________________________________________
c  Function dsk1 calculates function K^~_1 used in 1->2 decay to FIMP.
c   it reduces to usual Bessel function of second kind of degree one
c   if eta_i=0
c
c   Input:
c     xi = mi/T
c     etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c     pluss (minus) for bosons (fermions)
c     stat - determines if quantum corrections is to be considered 
c        (stat=0 for Maxwell-Boltzmann statistics) (stat\neq 0 for FD/BE)
c     k1 - output value
c     
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      subroutine dsfik1(x1,x2,x3,eta1,eta2,eta3,k1)   
      implicit none

      include 'dsficom.h'

      real*8, intent(in) :: x1, x2, x3, eta1, eta2, eta3
      real*8, intent(out) :: k1
      real*8 a,b, sum, eps, dsfik1int, dsbessek1
      real*8,external :: dsfik1int_simp


            
      a=1E-30; b=1
      eps = 1E-10

      x1_k1=x1;x2_k1=x2;x3_k1=x3;eta1_k1=eta1;eta2_k1=eta2
      eta3_k1=eta3
c     If stat=0 we use MB <=> K^~_1 reduces to Bessel
      if(stat.eq.0) then
            k1=dsbessek1(x1)/exp(x1)
      else
            call dgadap(a,b,dsfik1int_simp,eps,sum)
            k1=x1*sum
      endif
      end subroutine dsfik1

c     Integrand used to calculate K^~_1
      real*8 function dsfik1int(x,x1,x2,x3,eta1,eta2,eta3)
      implicit none


      real*8  x,x1,x2,x3,eta1,eta2,eta3,z, J, int, c1, dsfistat,a
      z=1/x
      c1=sqrt(z*z-1)
      J=c1*exp(-x1*z)/(1-eta1*exp(-x1*z))
      dsfik1int=J*dsfistat(x1*c1,x1,x2,x3,eta2,eta3)*z*z
      return
      end 

c     Putting variables into common blocks to get function of one variable
c     used to integrate.
      real*8 function dsfik1int_simp(x)
      implicit none

      include 'dsficom.h'

      real*8 x,xtemp,dsfik1int

      dsfik1int_simp=dsfik1int(x,x1_k1,x2_k1,x3_k1,eta1_k1,
     &eta2_k1,eta3_k1)
      return 
      end
      




