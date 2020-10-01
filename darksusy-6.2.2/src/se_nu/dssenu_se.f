       real*8 function dssenu_se(x,vbar)
c----------------------------------------------------------------------
c      dssenu_se
c      x: mx/m(i)
c      vbar: three-dimensional velocity dispersion of WIMPs in the halo
c      lars bergstrom 1995-12-12
c----------------------------------------------------------------------
       implicit none
       real*8 x,vbar,vesc,a
c...vesc is an average escape velocity that makes this approximation as
c...good as possible
       vesc=13.2d0 
       a=1.5d0*x/(x-1.d0)**2*(vesc**2/vbar**2)  ! Eq. (9-22) in jkg
       dssenu_se=(a**1.5d0/(1+a**1.5d0))**(1.d0/1.5d0) ! Eq. (9-21) in jkg
c       write(*,*) 'dssenu_se: ',x,vbar,dssenu_se
       return
       end
