!Set up the density profile for UCMHs, returning core and outer radii.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: Oct 2014
!
!Input:   profile         chosen density profile for UCMHs
!                         General:
!                          "plain": uncontracted, analytical profile with mass as per mucmh
!                         For UCMHs from phase transition <trans> in {"EW", "QCD", "ee"}: 
!                          <trans>_01: F=0.001, core radius 0.1%
!                          <trans>_11: F=0.01, core radius 0.1%
!                          <trans>_unc: uncontracted, numerical
!         z_i             redshift of formation of the the UCMH seed
!         z_c             latest allowed redshift of collapse in order to form UCMHs 
!         zstop           redshift at which UCMHs cease to grow 
!         mchi            WIMP mass (GeV)
!         sigchi          WIMP ann. cross-section (cm^3 s^-1)
!
!In/Out:  mucmh           UCMH mass (M_Sun); input if profile = "plain", output otherwise (as UCMHs are implied to have 
!                          originated in some specific phase transition in this case).
!
!Output:  rucmh           UCMH radius (kpc)
!         rmin            radius below which the radial infall assumption is violated (kpc)
!         rcut            radius at which WIMP self-annihilation cuts off the density profile (kpc)

      subroutine dsucmh_init_profile(profile, z_i, z_c, zstop, mchi, sigchi, mucmh, rucmh, rmin, rcut)

      implicit none

      include 'dsucmh.h'
      include 'dslocal.h'
      include 'dsmpconst.h'

      character (len=6), intent(IN) :: profile
      double precision, intent(IN) :: z_i, z_c, zstop, mchi, sigchi
      double precision, intent(INOUT) :: mucmh
      double precision, intent(OUT) :: rucmh, rmin, rcut


      double precision, parameter :: racc = 1.d-5 !Accuracy parameter for calculation of rmin
      integer, parameter :: mmax = 1000 !max number of iterations allowed trying to find r_cut
      double precision :: zstart, dsucmh_numerical_dens, testval, testdens, topval, bottomval, dsageatz
      integer :: m
      logical :: uncontracted 
      
      !Set zstart, the redshift of creation of the UCMH seed (e.g. string loop), to the lesser of zeq and z_i
      zstart = merge(zeq,z_i,zeq.lt.z_i)

      !Parse the requested profile
      uncontracted = (trim(profile) .eq. 'plain' .or. profile(len(trim(profile))-2:) .eq. 'unc')
      if (profile(1:2) .eq. 'EW') then
        mucmh = 5.6d-19 * (zeq + 1.d0) / (zstop + 1.d0)  ! set by mass scale of perturbations from EW transition)
      else if (profile(1:2) .eq. 'ee') then 
        mucmh = 0.33d0  * (zeq + 1.d0) / (zstop + 1.d0)  ! set by mass scale of perturbations from e+e- annihilation epoch
      else if (profile(1:3) .eq. 'QCD') then
        mucmh = 1.1d-9  * (zeq + 1.d0) / (zstop + 1.d0)  ! set by mass scale of perturbations from QCD transition
      else
        if (trim(profile) .ne. 'plain') then
          write(*,*) 'Unrecognised profile in dsucmh_init_profile.'
          stop
        endif  
      endif

      ! Choose UCMH halo profile
      call dshmucmhsetprofile('UCMH'//trim(profile))

      rucmh = 0.019d0 * mucmh**(0.3333d0) / (1.d0 + zstop)            !Rh in kpc
      rmin = 5.15d-3*(z_c+1.d0)**(-2.43)*((zstop+1.d0)*mucmh)**0.27d0 !radius in kpc at which the radial infall approximation is violated
      !Older; only valid for z_stop = 10
      !rmin = 9.9d-3 * (z_c+1.d0)**(-2.43) * mucmh**0.27d0            !radius in kpc at which the radial infall approximation is violated
      !Oldest; more conservative treatment
      !sigma_DM = 14.d0 * ((1.d0 + z_c)*1.d-3)**(-0.5d0) 
     &!         * mucmh**0.28d0                                       !sigma_DM in cm/s
      !rmin = sigma_DM**(8./7.) * (rucmh * cmperpc * 1.d3)**(11./7.)  !radius in kpc at which the radial 
     &!      / (GNewton * mucmh * SolarMass)**(4./7.) / cmperpc * 1.d-3 ! infall approximation is violated
      if (uncontracted) then
        maxdens = mchi/sigchi/(t_0-t_eq*((zeq+1.d0)/(zstart+1.d0))**1.5d0)
      else
        maxdens = mchi/sigchi/(t_0-dsageatz(zstop))
      endif

      
      if (trim(profile) .eq. 'plain') then ! Use a non-contracted analytical profile

        rcut = 3.d0 * sigchi / cmperpc**3 * 1.d-9 * mucmh * dmfrac * GeVperSolarMass * (t_0-t_eq*((zeq+1.d0)/(zstart+1.d0))**1.5d0)
     &   / (16.d0 * pi * mchi * rucmh**0.75d0)
        rcut = rcut ** (4.d0/9.d0)

      else ! Use adiabatically-contracted numerical profiles

        topval = rucmh

        if (maxdens .lt. dsucmh_numerical_dens(topval)) then !Whole UCMH needs to be truncated
                
          rcut = rucmh

        else 
                
          bottomval = 0.d0
                
          if (maxdens .gt. dsucmh_numerical_dens(bottomval)) then !No truncation (unusual...)
                
            write(*,*) 'Warning: UCMH is not being truncated at all - this is suspicious...'
            rcut = 0.d0
                
          else !Find the radius at which the UCMH needs to be truncated

            do m = 1, mmax
              testval = (topval+bottomval)/2.d0
              testdens = dsucmh_numerical_dens(testval)
              if (testdens .gt. maxdens) then
                bottomval = testval
              else
                topval = testval
              endif
              if (abs(testdens - maxdens) .lt. abs(racc*maxdens)) exit
            enddo         

            if (m .eq. mmax + 1) then
              write(*,*) 'rcut not found after ',mmax,' steps.  Exiting...'
              stop
            endif

            rcut = testval

          endif

        endif

      endif
      
      !Set the internal representation of the profile correctly
      call dshmucmhset(mucmh, rucmh, max(rcut, rmin), 4.d0)

      end subroutine dsucmh_init_profile

