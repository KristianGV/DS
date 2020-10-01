c_______________________________________________________________________
c  Function dsfidecoh2.
c     
c     Input:
c       TR - Reheating temperature in GeV
c       w - Decay width Y->x,x in s^-1
c       M - mediator (Y) mass
c       g - internal degrees of freedom to mediator
c       eta - etas*e^mu/T
c       stat - 0 for no statistics.
c  
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfidecoh2(TR,M,w,g,eta,stat)
      implicit none

      include 'dsmpconst.h'

      real*8 TR,TRgev,M,w,wgev,g,eta,dsfidecab,s0,gtoGeV
      integer stat
c     Using that entropy today is s0=2.8912E9 m^-3 = 2.8912E3 cm^-3
      s0=2.8912E3
c     Using that 1 GeV = 1.78266192E-24 g
c     Since rho_conh2 is iven in units g cm^-3, we have s0/rho_conh2 is in units g^-1
      gtoGeV=1.78266192E-24
c     Using that s^-1 = 6.582119E-25 GeV
      wgev=6.582119d-25 *w
c     Using that 1 GeV = 1.160451812E13 Kelvin
C      TRgev=TR/1.160451812D13
      dsfidecoh2=m*dsfidecab(TR,M,wgev,g,eta,stat)*s0/rho_conh2*gtoGeV

      return 
      end