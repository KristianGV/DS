
*************************
      real*8 function dssem_epfunc(r)
      implicit none
      include 'dsmpconst.h'

      real*8 r,dssem_earthmass

      dssem_epfunc=dssem_earthmass(r)*1.d-3*GNewton/max(r,100.0d0)**2

      return
      end
