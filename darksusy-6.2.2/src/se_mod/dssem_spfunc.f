
*************************
      real*8 function dssem_spfunc(r)
      implicit none
      include 'dsmpconst.h'
      
      real*8 r,dssem_sunmass

      dssem_spfunc=dssem_sunmass(r)*1.d-3*GNewton/max(r,100.0d0)**2
c      write(*,*) 'dssem_spfunc: ',r,dssem_spfunc

      return
      end
