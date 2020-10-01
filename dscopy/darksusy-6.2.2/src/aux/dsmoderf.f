      real*8 function dsmoderf(x)
c =================================================================
c error function 
c modified by l. bergstrom 98-09-15
c modified by p. gondolo 2000-07-19
c see test output below
c used for damour-krauss calculations
c =================================================================
      implicit none
      include 'dsmpconst.h'
      real*8 erf,x
      dsmoderf=dsqrt(pi)/2.d0*erf(x)
      return
      end
