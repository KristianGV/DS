
***********************************************************************
*** l.b. and j.e. 1999-04-06
*** auxiliary function for r-integration
*** input: radius in centimeters
*** output: integrand in cm^-1 s^-1
***********************************************************************

      real*8 function dssenu_ceint(r,foveru)
      implicit none
      real*8 mx,sigsi,mu,muplus,muminus,ma,rx
      real*8 v,res,dssenu_ceint2,r,dssem_earthvesc,
     &  umin,umax,foveru,sigsd
      integer imass,vtype
      external foveru
      common/seint/mx,ma,sigsi,sigsd,rx,vtype
      common/seint2/imass
      external dssenu_ceint2
      include 'dshmcom.h'
      include 'dsmpconst.h'

      rx=r
      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0
c...begin velocity integration
c...
c...determine integration limits
      v=dssem_earthvesc(r/100.0d0)   ! escape velocity at r

c...For a general velocity distribution
      umin=0.d0                  ! lower velocity limit
      umax=sqrt(mu/muminus**2)*v ! upper velocity limit
c...Note, we should only allow such velocities that we scatter to
c...velocities lower than the escape velocity. This means that we have to
c...make sure that (mu/muplus^2 > u^2 / (u^2 + v^2 )). This is the same as
c...the condition umax above.

c      write(*,*) 'umin=',umin,'  umax=',umax
      if (umin.lt.umax) then
        call dssenu_hiprecint2(dssenu_ceint2,foveru,umin,umax,res)
        res=res*1.0d5    ! km/s -> cm/s from u integration
        dssenu_ceint=res*r*r*4.*pi  ! assume spherical earth
      else
        dssenu_ceint=0.0d0
      endif

c      if (umax.gt.30.0d0) then
c        write(*,*) 'umax=',umax,'  res=',res
c      endif

      return
      end
