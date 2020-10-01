c.. Program som beraknar massan i vintergatan innanfor en sfar med
c.. radien dmin kpc runt solen (med antaganden att profilen ar NFW
c.. fran Battaglia et al, MNRAS 370:1055 2006 -- se dslocal.h).
c..
c.. (Mass of DM in sphere around Sun assuming NFW profile as per
c.. Battaglia+)
c..
c.. Input: dmin (kpc)                 radius of sphere
c.. Output: dsucmh_localmass (M_sun)  DM mass within that sphere
c..
c.. Author: Sofia Sivertsson August/September 2009
c.. Modified: Pat Scott (pat@fysik.su.se) 100529, 110127, 141011

      real*8 function dsucmh_localmass(dmin)
      implicit none
      include 'dslocal.h'
      include 'dsmpconst.h'
      external integrand
      real*8 dmin                ! sfärens radius omkring solen i kpc
      real*8 dsf_int

      if (dmin .lt. (r_virial+r_sun_gal)/(1.d3*pc)) then
        dmin_internal = dmin
        dsucmh_localmass=2.d0*pi*dsf_int(integrand,0.d0,pi,1.d-4)/Msun
      else
        dsucmh_localmass=M_200/Msun
      endif

      return

      end

*************************************************************************
*** integrand ***********************************************************
*************************************************************************

      real*8 function integrand(theta)
      implicit none
      include 'dslocal.h'
      real*8 theta
      external inreintegrand
      real*8 dsf_int2

      theeta=theta
      integrand=dsf_int2(inreintegrand,0.d0,dmin_internal*1.d3*pc,1.d-5)

      return
      end


*************************************************************************
*** inreintegrand (densityfunktion) *************************************
*************************************************************************

      real*8 function inreintegrand(r_sun)
      implicit none
      include 'dslocal.h'
      real*8 density
      real*8 r_sun ! distance from Sun.
      real*8 r_gal ! distance from Galactic centre
      real*8 r_s

      r_s = r_s3**(1.d0/3.d0)
      r_gal=sqrt(r_sun_gal**2+r_sun**2-2.d0*r_sun_gal*r_sun*cos(theeta))

      if (r_gal .lt. r_virial) then
        density=rhobg*del/(r_gal/r_s*(1.d0+r_gal/r_s)**2)
      else
        !Make sure to zero the density outside the virial radius
        density=0.d0
      endif

      !localdens=rhobg*del/(r_sun_gal/r_s*(1d0+r_sun_gal/r_s)**2)
      inreintegrand=density*r_sun**2*sin(theeta)

      return
      end

