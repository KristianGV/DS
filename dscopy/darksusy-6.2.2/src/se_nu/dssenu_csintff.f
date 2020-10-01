***********************************************************************
*** l.b. and j.e. 1999-04-06
*** auxiliary function for r-integration
*** input: radius in centimeters
*** output: integrand in cm^-1 s^-1
*** Adapted for the Sun by J. Edsjo, 2003-11-26
*** Modified to allow for numerical form factor integration by
*** J. Edsjo, 2010-05-28
***********************************************************************

      real*8 function dssenu_csintff(r,foveru)
      implicit none
      real*8 mx,sigsi,mu,muplus,muminus,ma,rx
      real*8 v,res,dssenu_csintff2,r,max,dssem_sunvesc,
     &  umin,umax,foveru,sigsd,vp
      integer vtype
      external foveru
      common/seint/mx,ma,sigsi,sigsd,rx,vtype
      external dssenu_csintff2
      include 'dssecom.h'
      include 'dsmpconst.h'

      rx=r
      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0
c...begin velocity integration
c...
c...determine integration limits
      v=dssem_sunvesc(r/100.0d0)   ! escape velocity at r

c...Now use a different maximal velocity after the scatter to allow
c...for cutting away WIMPs that reach too far out (to Jupiter).
      vp=sqrt(v**2-veout**2)

c...general velocity limits
      umin=0.d0                  ! lower velocity limit
c      umax=sqrt(mu/muminus**2)*v ! upper velocity limit
c...Correct umax for a differnt choice of maximum velocity
c...veout is typically set to 0 for standard treatment or to
c...the escape veolocity at Jupiter for a more conservative estimate
c...
c...The 10000.d0 maximum is put to avoid numerical problems when mx and
c...ma matches. In principle, that should be set to the maximum velocity
c...for which f(u) is non-zero, but 10000.d0 is a safe choice here.
      if (abs(muminus).gt.1.d-10) then
        umax=sqrt(max(vp**2*mu/muminus**2-veout**2,0.d0))
        umax=min(umax,10000.d0)
      else
        umax=10000.d0
      endif
c...Note, we should only allow such velocities that we scatter to
c...velocities lower than the escape velocity. This means that we have to
c...make sure that (mu/muplus^2 > u^2 / (u^2 + v^2 )). This is the same as
c...the condition umax above.

      if (umin.lt.umax) then
        call dssenu_hiprecint2(dssenu_csintff2,foveru,umin,umax,res)
        res=res*1.0d5    ! km/s -> cm/s from u integration
        dssenu_csintff=res*r*r*4.*pi  ! assume spherical sun
      else
        dssenu_csintff=0.0d0
      endif

      return
      end
