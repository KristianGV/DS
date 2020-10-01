
c_______________________________________________________________________
c  Function dsk1.
c
c  
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      subroutine dsfik1(x1,x2,x3,eta1,eta2,eta3,stat,k1)   
      implicit none

      integer, intent(in) :: stat
      real*8, intent(in) :: x1, x2, x3, eta1, eta2, eta3
      real*8, intent(out) :: k1
      real*8 a,b, dsf_int, eps, dsfik1int, dsbessek1
      real*8,external :: dsfik1int_simp
      
      a=1E-9; b=1
      eps = 0.000001
      if(stat.eq.0) then
            k1=dsbessek1(x1)/exp(x1)
      else
            k1=dsf_int(dsfik1int_simp,a,b,eps)
      endif
      end subroutine dsfik1


      real*8 function dsfik1int(z,x1,x2,x3,eta1,eta2,eta3,stat)
      implicit none
      real*8  x1,x2,x3,eta1,eta2,eta3,z, J, int, c1, dsfistat,a
      integer stat

      c1=sqrt(1/z/z-1)
      J=c1*exp(-x1/z)/(1-eta1*exp(-x1/z))
      int=J*dsfistat(x1*c1,x1,x2,x3,eta2,eta3,stat)/z/z
      dsfik1int=int
      return
      end 

      real*8 function dsfik1int_simp(z)
      implicit none
      real*8 z,x1,x2,x3,eta1,eta2,eta3, dsfik1int
      integer stat

      common /dsfik1var/ x1,x2,x3,eta1,eta2,eta3,stat
      dsfik1int_simp=dsfik1int(z,x1,x2,x3,eta1,eta2,eta3,stat)
      return 
      end
      




