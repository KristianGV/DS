
*************************
      real*8 function dssem_edfunc(r)
      implicit none
      include 'dsmpconst.h'

      real*8 r,dssem_earthdens

      dssem_edfunc=dssem_earthdens(r)*1000.0d0*4.0d0*pi*r**2

      return
      end
