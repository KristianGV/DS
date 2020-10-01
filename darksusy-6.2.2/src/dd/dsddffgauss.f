      subroutine dsddffgauss(q,a,z,ff,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Calculates F^2(q) 
c     
c     F(q) is the gaussian form factor
c
c     Input: a: Mass number of Element
c            z: Atomic Number of Element
c            q: Momentum transfer in GeV
c
c     Output: ff: Form Factor squared (normalized to F^2(0) = 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-------------------------------------------------------

      implicit none
      include 'dsmpconst.h'
      real*8 q,ff
      integer a,z,ierr
      real*8 r,x
      
      ierr=0

      r = 0.89d0*exp(log(dble(a))/3.d0)+0.3d0
      x = q*r*fermiGeV
      ff = exp(-x*x/3.d0)
      
      return
      end
