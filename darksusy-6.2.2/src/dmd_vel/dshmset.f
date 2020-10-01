      subroutine dshmset(c)
****************************************************************
*** subroutine dshmset:                                      ***
*** initialize the density profile and or the small clump    ***
*** probability distribution                                 ***
*** type of halo:                                            ***
***   hclumpy=1 smooth, hclumpy=2 clumpy                     ***
***                                                          ***
*** a few sample cases are given; specified values of the    ***
*** local halo density 'rho0' and of the length scale        ***
*** parameter 'a' should be considered just indicative       ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
*** small modif: paolo gondolo 00-07-19                      ***
*** mod: 03-11-19 je, 04-01-13 pu, 09-05-07 ps, 09-08-08 ps  ***
****************************************************************
      implicit none
      include 'dshmcom.h'
      include 'dsmpconst.h'
      character*(*) c
      logical firstdd
      data firstdd/.true./
      save firstdd

c...Set default flags
      udfload=.true.  ! load udffile (if chosen by veldf='user') on next
                      ! call to dshmudf.f
      udfearthload=.true. ! load udfearthfile (if chosen by veldfearth='user')
                          ! on next call to dshmudfearth.f

      isodfload=.true.   ! load isodf file on next call to dshmisotrnum.f

c...Three-dimensional velocity dispersion of the WIMPs in the halo
c...Note that in a simple isothermal sphere, the circular speed
c...(approx. solar speed v_sun) is sqrt(2/3)*vd_3d
      vd_3d = 270.0d0     ! WIMP 3D velocity dispersion, km/s
      v_sun = 220.0d0     ! circular velocity at the Sun location, km/s
      v_obs = v_sun       ! Velocity of observer
      vgalesc = 600.d0          ! galactic escape speed in km/s
c...Observer speed w.r.t. the halo, i.e. including Sun + Earth speed,
c...yearly average,
      vobs = 264.d0       ! observer speed w.r.t. the halo in km/s
      v_earth=29.78d0    ! Earth speed in the solar system in km/s

c a few options for the density profile:

c ... navarro-frenk-white profile, smooth profile
      if (c.eq.'nfwsm'.or.c.eq.'default') then
        hclumpy=1
           ! still needed in jfactor routines 
        r_0=8.d0                ! sun galactocentric distance (kpc)
           ! still needed in jfactor routines 
        rho0=0.3d0              ! local halo density (gev cm^-3)
           ! still needed in jfactor routines
c        haloshape='spherical' !here you choose between spherical and axisymm.
c        halotype='albega' !this picks the alpha-beta-gamma profile, set by:
c        alphah=1.d0        
c        betah=3.d0
c        gammah=1.d0
c        ah=20.d0             ! length scale (kpc)
c        Rref=r_0
c        rhoref=rho0         ! normalizing the profile to the local halo density
c        rhcut=1.d-5         ! cut radius (kpc) 
c        haloid='nfwsm'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

      else
         write (*,*) 'dshmset: unrecognized option ',c
         stop
      endif

      return
      end
