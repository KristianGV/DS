      subroutine dsddvearth(t,vearth)
c
c speed of the earth relative to the galaxy in km/s at 
c time t in days from 12:00 UT Dec 31, 1999
c
c formulas from Green, astro-ph/0304446, originally by Lewin and Smith
c
      implicit none
      include 'dshmcom.h'
      real*8 t,vearth,vea(3),vlsr(3),vsun(3),vtot(3)
      real*8 ve,lmeansun,g,lambda,ellipticity,lambda0,
     &     betax,betay,betaz,lambdax,lambday,lambdaz,velambda
      integer modulation
      common /ddmodul/ modulation
      vlsr(1) = 0.d0
      vlsr(2) = v_sun ! 220.d0
      vlsr(3) = 0.d0
      vsun(1) = 10.0d0
      vsun(2) = 5.2d0
      vsun(3) = 7.2d0
      if (modulation.eq.0) then
         vea(1) = 0.d0
         vea(2) = 0.d0
         vea(3) = 0.d0
      else
         ve = v_earth
         lmeansun = 280.460d0+0.9856003d0*t
         g = 357.528d0+0.9856003d0*t
         lambda = lmeansun+1.915d0*dsin(g)+0.020d0*dsin(2.d0*g)
         ellipticity = 0.016722d0
         lambda0 = 13.0d0
         betax = -5.5303d0
         betay = 59.575
         betaz = 29.812
         lambdax = 266.141
         lambday = -13.3485
         lambdaz = 179.3212
         velambda = ve*(1.d0-ellipticity*dsin(lambda-lambda0))
         vea(1) = velambda*dcos(betax)*dsin(lambda-lambdax)
         vea(2) = velambda*dcos(betay)*dsin(lambda-lambday)
         vea(3) = velambda*dcos(betaz)*dsin(lambda-lambdaz)
      endif
      vtot(1) = vlsr(1)+vsun(1)+vea(1)
      vtot(2) = vlsr(2)+vsun(2)+vea(2)
      vtot(3) = vlsr(3)+vsun(3)+vea(3)
      vearth = sqrt(vtot(1)**2+vtot(2)**2+vtot(3)**2)
      return
      end


