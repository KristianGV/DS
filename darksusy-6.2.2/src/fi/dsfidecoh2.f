c_______________________________________________________________________
c  Function dsfidecoh2 calculates relic density for 1->2
c     
c     Input:
c       TR - Reheating temperature in GeV
c       w - Decay width Y->x,x in GeV
c       M - mediator (Y) mass
c       mdm - dm mass
c       g - internal degrees of freedom to mediator
c       etai=\pm exp^(mu_i/T)  - parameter in phase-space distribution
c       pluss (minus) for bosons (fermions)
c       stat - determines if quantum corrections is to be considered 
c       (stat=0 for Maxwell-Boltzmann statistics) (stat\neq 0 for FD/BE)
c  
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfidecoh2(Tmin,TR,w,M,g,eta)
      implicit none

      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 TR,TRgev,w,g,dsfidecab,s0,Gevg,M,eta,Tmin
c     Using that entropy today is s0=2.8912E9 m^-3 = 2.8912E3 cm^-3
      s0=2.8912E3
c     Using that 1 GeV = 1.78266192E-24 g
c     Since rho_conh2 is given in units g cm^-3, we have s0/rho_conh2 is in units g^-1
      Gevg=1.78266192E-24

      dsfidecoh2= mdm*dsfidecab(Tmin,TR,w,M,g,eta)*s0/rho_conh2*Gevg

      return 
      end