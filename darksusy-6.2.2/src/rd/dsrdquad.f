      real*8 function dsrdquad(x0,y0,x1,y1,x2,y2,x,reslin,err)
**********************************************************************
*** Function dsrdquad performs a quadratic interpolation between the
*** points (x0,y0), (x1,y1) and (x2,y2). It will return value from
*** the interpolation at point x. x is assumed to be x0 <= x <= x2.
*** The points are assumed to be ordered in x, x0<x1<x2.
*** Upon return, err is the relative error compared to a linear      
*** interpolation. The setup is done with Newton divided differences.
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2018-10-27      
**********************************************************************      

      implicit none

      real*8 x0,y0,x1,y1,x2,y2,x,err
      real*8 f0,f1a,f1b,f2
      real*8 res,reslin

      f0=y0
      f1a=(y1-y0)/(x1-x0)
      f1b=(y2-y1)/(x2-x1)
      f2=(f1b-f1a)/(x2-x0)
      res=f0+f1a*(x-x0)+f2*(x-x0)*(x-x1) ! quad interpolation

c..Compare with linear
      if (x.le.x1) then
         reslin=y0+(y1-y0)/(x1-x0)*(x-x0)
      else
         reslin=y1+(y2-y1)/(x2-x1)*(x-x1)
      endif
      err=(res-reslin)/reslin
      dsrdquad=res
      
      return
      end
      
