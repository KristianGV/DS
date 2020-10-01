***********************************************************************
*** dskdboltz implements the Boltzmann equation for the evolution of
*** the WIMP temperature as dy/dx = rhs, expressed in the variables
*** (see arXiv:1603.04884)
***
***    y = m0 * T_DM * s**(-2/3)
***    x = m0 / T
***
*** NB: src_models/xxx/mh/dskdparticles MUST be called before using this!
***    (typically done via a call to dskdtkd)
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
*** updates: 2013-06-11 (removed model-dependence)
***          2017-05-11 (changed definition of y)
***          2018-05-28 (allow for Tdark != Tphoton) 
***********************************************************************

      subroutine dskdboltz(x,y,rhs)
      implicit none
      include 'dsmpconst.h'
      include 'dskdcom.h'

      real*8 x,y,rhs
      real*8 xd, dsrdxi
      real*8  sqrtstar,heff,yeq

c... functions
      logical dsisnan
      real*8 dskdgammarate


      rhs=0d0
      if (dsisnan(x).or.(x.le.0.d0)) return

      call dsrddof(m0/x,sqrtstar,heff)
      xd = x/dsrdxi(x)
c      yeq = x / (heff*2.d0*pi**2/45.d0)**(2.d0/3.d0)
      yeq = x**2/xd / (heff*2.d0*pi**2/45.d0)**(2.d0/3.d0)

c... implement (A6) in 1603.04884 (which takes into account erratum to hep-ph/0612238)
      rhs =(sqrtstar/heff)*dskdgammarate(x)
     &      /(sqrt(4.*pi**3/45.d0)/mpl*m0**2/x**2) ! this is H, up to geff^0.5
     &      *(yeq - y)/x  

      return
      end





        

