**********************************************************************
*** dshmucmhsetprofile: store profile name for UCMH
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Adapted from Pat Scott's implementation with dshmset
**********************************************************************
      subroutine dshmucmhsetprofile(c)
      implicit none
      character*(*) c
      include 'dsucmh.h'

c......UCMH spherical infall profile
      if (c(1:4).eq.'UCMH') then
c        r_0=4.d0            ! distance to UCMH (kpc)
c        rho0=0.3d0          ! local halo density (gev cm^-3)
c        haloshape='spherical' !here you choose between spherical and axisymm.
c        halotype=c          !this picks the UCMH profile, set by:
        ucmhprofile=c(5:)   ! anything except 'plain' needs external numerical routines
        Mh_ucmh=0.d0        ! UCMH mass (solar masses) - must be set model by model
        Rh_ucmh=0.d0        ! UCMH outer radius (kpc) - must be set model by model
        rhcut=0.d0          ! core radius (kpc) - must be set model by model
      else
        write(*,*) 'ERROR in dshmucmhsetprofile: unknown profile ',c
        stop
      endif

      return
      end
