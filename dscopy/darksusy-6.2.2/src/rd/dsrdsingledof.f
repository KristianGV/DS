
c_______________________________________________________________________
c  Function dsrdsingledof returns the energy density of a thermally
c  distributed single particle degree of freedom, normalized to the
c  energy density of a single bosonic degree of freedom in the highly
c  relativistic limit. At high tempertures, the return value is thus
c  1.0 (0.875) for bosons (fermions).
c
c  input:
c    x    - particle mass divided by temperature
c    stat - bosonic (stat=1) or fermionic (stat=2) d.o.f.
c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no), 2018-08-30
c=======================================================================
      real*8 function dsrdsingledof(x,stat)
      implicit none
      include 'dsmpconst.h'
      real*8 x
      integer stat
      real*8 tmp, dsbessek0, dsbessek1
      integer i

      dsrdsingledof = 0.0d0
      if (x.gt.30d0) return
      tmp = 0.0d0

      if (stat.eq.1) then   ! Boson
        if (x.lt.0.1) then
          tmp = 1.0d0
        else
          do i=1,4
            tmp = tmp + 
     &            (3*x**2/(1.*i**2)*dsbessek0(i*x)
     &             + (6*x/(1.*i**3)+x**3/(1.*i)) * dsbessek1(i*x) )
     &                *exp(-i*x)   
          enddo
          tmp = tmp*15/pi**4/0.9967  ! last factor ensures correct 
                                     ! normalization for x->0
        endif      
      elseif (stat.eq.2) then  ! Fermion
        if (x.lt.0.1) then
          tmp = 0.875d0
        else
          do i=1,4
            tmp = tmp + (-1.)**(i+1)*
     &            (3*x**2/(1.*i**2)*dsbessek0(i*x)
     &             + (6*x/(1.*i**3)+x**3/(1.*i)) * dsbessek1(i*x) )
     &                *exp(-i*x)                      
          enddo
          tmp = tmp*15/pi**4*1.00116 ! last factor ensures correct 
                                     ! normalization for x->0
        endif
      else
        tmp = 0.0d0
        write(*,*) 'Warning: dsrdsingledof called with illegal option stat = ',stat
        stop
      endif

      dsrdsingledof = tmp

      return
      end


      
