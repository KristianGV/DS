      subroutine dsddffgould(q,a,z,ff,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Calculates F^2(q) 
c     
c     F(q) is the old obsolete (incorrect) form factor from Gould
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
      real*8 mni,e0,dele
      
      ierr=0

      mni = (z*m_p_amu+(a-z)*m_n_amu)*atomicmassunit
      e0=3./2./mni*0.038938     ! gives units of gev; see gould (a8)
      e0=e0/(0.91*mni**(0.33333)+0.3)**2
      dele=q**2/(2.d0*mni)
      ff=exp(-dele/e0)
      
      return
      end
