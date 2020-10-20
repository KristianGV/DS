c_______________________________________________________________________
c  Function dsfidecab calculates the freeze-in decay abundance.
c
c  Input:
c     TR - reheating temperature 
c     M - mediator mass
c     w - mediator decay width
c     g - mediator internal degrees of freedom
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfidecab(Tmin,TR,w,M,g,eta)
      implicit none

      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 TR,w,g,a,b,Tmin,eps,dsf_int2,sum,M,eta,abserr,resabs,resasc 
      real*8, external :: dsfidecint_simp

      M_dec=M; eta_dec=eta
      a=log(Tmin);b=log(TR);eps = 1E-5
      ! a=Tmin;b=TR;eps = 1E-5
      call dqk21(dsfidecint_simp,a,b,sum,abserr,resabs,resasc)  

      ! call dgadap(a,b,dsfidecint_simp,eps,sum)
c      sum=dsf_int2(dsfidecint_simp,a,b,eps)
      dsfidecab=w*g*M_dec*M_dec*sum/(pi*pi)

      return
      end