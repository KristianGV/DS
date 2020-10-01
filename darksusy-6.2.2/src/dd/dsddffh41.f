      subroutine dsddffh41(q,a,z,ff,ierr)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Calculates F^2(q) 
c     
c     F(q) is Helm type form factor computed using
c      as in DarkSUSY 4.1 (with corrected expansion of f, PG 20080216)
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
      real*8 mni,r,s,r1tmp,r1,x,y,f,x2
      
      ierr=0
      
      if (a.eq.1) then
         ff = 1.d0
      else
         mni = (z*m_p_amu+(a-z)*m_n_amu)*atomicmassunit
         r = (0.91d0*exp(log(mni)/3.d0)+0.3d0)
         s = 1.d0
         r1tmp = r*r-5.d0*s*s
         if (r1tmp.gt.0.d0) then
            r1 = sqrt(r1tmp)
         else
            r1 = r
         endif 
         x = dabs(q*r1*fermiGeV)
         y = q*s*fermiGeV
         if (x.gt.5.d-8) then
            f = (3.d0*(sin(x)-x*cos(x))/x**3)
         else
            x2 = x*x
            f = 1.d0+x2*(-0.1d0+x2*(1.d0/280.d0+x2*(-1.d0/15120.d0+
     &           x2/1330560.d0)))
         endif
         ff = f**2*exp(-y**2)
      endif
      
      return
      end
