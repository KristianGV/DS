c_______________________________________________________________________
c  Function dsfidecab.
c
c  
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfidecab(TR,M,w,g)
      implicit none

      include 'dsmpconst.h'

      real*8 TR,w,g,M,a,b,Tmin,eps,dsf_int
      real*8, external :: dsfidecint_simp
      Tmin=M/100
      a=log(Tmin);b=log(TR);eps = 0.00001

      dsfidecab=w*g*M*M*dsf_int(dsfidecint_simp,a,b,eps)/(pi*pi)

      return
      end