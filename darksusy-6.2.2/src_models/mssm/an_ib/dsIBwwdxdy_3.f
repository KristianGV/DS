*****************************************************************************
*** auxiliary routine called by dsIBwwdxdy
*** author: Torsten Bringmann, 2007-02-16
*****************************************************************************

      real*8 function dsIBwwdxdy_3(x,y,m0,mw,mc1,mc2)
      implicit none
      real*8 x,y,m0,mw,mc1,mc2

      dsIBwwdxdy_3 = 
     -   (-64*m0**3*mc2*(mw**8*(-6*mw**2 + mc1**2*(7 + x)) + 
     -      256*m0**10*(-1 + x)*(x - y)**2*y**2 + 
     -      m0**2*mw**6*(4*mc1**2*(11*x**2 + x*(6 - 4*y) - 12*y) + 
     -         mw**2*(7 - 31*x - 24*x**2 + 32*y)) - 
     -      64*m0**8*(x - y)*y*
     -       (-4*mc1**2*(-1 + x)*(x - y)*y + 
     -         mw**2*(3*x**2 + 4*y*(1 + 2*y) - 2*x*(1 + 6*y))) - 
     -      4*m0**4*mw**4*(4*mc1**2*
     -          (-2*x**3 + 2*x*(1 - 3*y)*y - 2*y**2 + x**2*(1 + 8*y))+ 
     -         mw**2*(8*x**3 + 4*(3 - 4*y)*y - 6*x*(1 + 2*y) + 
     -            x**2*(-3 + 16*y))) - 
     -      16*m0**6*mw**2*(4*mc1**2*y*
     -          (3*x**3 - 4*y**2 + 2*x*y*(3 + 2*y) - x**2*(2 + 7*y))+ 
     -         mw**2*(2*x*(1 - 19*y)*y + x**3*(-2 + 8*y) + 
     -            2*y**2*(-1 + 16*y) + x**2*(1 + 8*y - 8*y**2)))))/
     -  (mw**4*(-2*mc1**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc2**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (mw**2 + 4*m0**2*(x - y))*(mw**2 - 4*m0**2*y)*
     -    (-2*mc1**2 + mw**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mw**2 + m0**2*(-2 + 4*y)))
      return
      end   ! dsIBwwdxdy_3

