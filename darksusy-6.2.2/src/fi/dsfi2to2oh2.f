c_______________________________________________________________________
c  Function dsfi2to2oh2 calculates relic abundance from 2->2.
c     
c    input:
c       TR - reheating temperature
c       Tmin - minimum temperature
c       mi - mass for particle i
c       etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c          pluss (minus) for bosons (fermions)    
c       gi - internal degrees of freedom for particle i
c       c12 - constant = 1/2 if particle 1= particle 2, else c12=1.
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfi2to2oh2(TR,Tmin,m1,m2,eta1,eta2,etaX,g1,g2,c12)
      implicit none

      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 TR,Tmin,s0,GeVg,m,dsfi2to2ab,eta1,eta2,etaX,g1,g2,c12,m1,m2


c     Using that entropy today is s0=2.8912E9 m^-3 = 2.8912E3 cm^-3
      s0=2.8912E3
c     Using that 1 GeV = 1.78266192E-24 g -> GeV/g=1.78266192E-24
c     Since rho_conh2 is iven in units g cm^-3, we have s0/rho_conh2 is in units g^-1
      GeVg=1.78266192E-24
c     Using that s^-1 = 6.582119E-25 GeV
      dsfi2to2oh2=mdm*dsfi2to2ab(TR,Tmin,m1,m2,eta1,eta2,etaX,g1,g2,c12)
     &*s0/rho_conh2*GeVg
      return 
      end