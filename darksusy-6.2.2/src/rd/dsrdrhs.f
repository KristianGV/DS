      subroutine dsrdrhs(x,wrate,lambda,yeq,nfcn)
c_______________________________________________________________________
c  adimensional annihilation rate lambda in the boltzmann equation
c    y' = -lambda (y**2-yeq**2) and equilibrium dm density in units
c    of the entropy density.
c  input:
c    x - mass/ PHOTON temperature (real)
c    wrate - invariant annihilation rate (real)
c  output:
c    lambda - adimensional parameter in the evolution equation (real)
c    yeq - equilibrium number/entropy densities (real)
c    nfcn - number of calls to wrate (integer)
c  common:
c    'dsrdcom.h' - included common blocks
c  uses qrkck or dgadap.
c  called by dsrdeqn.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994-1996
c  modified: joakim edsjo (edsjo@fysik.su.se) 98-04-28
c  modified 98-04-28: y_eq corrected by a factor of 2 (je)
c  modified 18-05-28: allowed for Tdark != Tphoton (see dsrdxi for description)
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 dsrdthav
      real*8 x,wrate,lambda,yeq,heff,
     &  k2,dsbessek2
      integer nfcn,i
      external wrate
      real*8 xd,dsrdxi
      real*8 t,grav,wav,sqrtgstar,eqcnst
c     grav=mplanck*sqrt(pi/45), eqcnst=45/(4*pi^4)
      parameter (grav=3.2262726d18) 
      parameter (eqcnst=0.1154923d0)

c--------------------------------------------------------------------dof
      t=mco(1)/x
      call dsrddof(t,sqrtgstar,heff)
c      write (*,'(1x,a,1x,e14.8,1x,i3,1x,i3,2(1x,e14.8)))') 
c     & 'dsrdrhs: t,klo,khi,sqrtgstar,heff',
c     & t,klo,khi,sqrtgstar,heff

c... make sure to distinguish between Tphoton and T dark
c... the latter enters in Yeq and thermal average, but NOT in (the rest of) lambda
      xd = x/dsrdxi(x)

c-------------- thermally averaged cross section times relative velocity

      wav=dsrdthav(xd,wrate)

c...This piece of code is just for testing      
c      wav=dsrdthav(27.d0,wrate)
c      do k=1,nrd
c         write(47,*) k,ppp(k),rdwx(k),
c     &        dsrdbreit_wigner(rgev(1),rwid(1),resalpha(1),mco(1),
c     &           ppp(k))*resnorm(1)
c      enddo
c      do k=0,1000
c         pa=14.6+dble(k)/1000.d0*0.5d0
c         write(53,*) k,pa,dsrdwx(pa),
c     &        dsrdbreit_wigner(rgev(1),rwid(1),resalpha(1),
c     &        mco(1),pa)*resnorm(1),dsanwx(pa)
c      enddo
c      
c      stop     ! testing code until here

         
      if (rderr.ne.0) return

c---------------------------------------------------------------- lambda

      lambda=grav*sqrtgstar*wav*mco(1)/(x**2)
      lambda=max(lambda,1.d-100) ! JE Correction for wrate=0 for low p

c--------------------------------------------------- equilibrium density

      yeq=0.d0

      do i=1,nco
        k2=exp(-xd*mco(i)/mco(1))*dsbessek2(xd*mco(i)/mco(1))
        yeq=yeq+(mco(i)/mco(1))**2*mdof(i)*k2
      enddo
c      yeq=eqcnst*yeq*x**2/heff  
      yeq=eqcnst*yeq/xd*x**3/heff  ! NOte that s is still definied wrt to standard T!
      
      return
      end
