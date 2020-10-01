!Theoretical prediction of cosmological fraction of UCMHs from power-law Gaussian power spectra.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
! 
!Input:   n               spectral index of primordial spectrum of perturbations (dimensionless)
!         mucmh           present-day UCMH mass (M_Sun)
!         z_c             latest allowed redshift of collapse in order to form UCMHs (dimensionless) 
!         zstop           redshift at which UCMHs cease to grow (dimensionless)
!         rough           Rely on Green-Liddle approximation for mass variance (only to be used for back-comparison)
!
!Output:  f (return val)  cosmological fraction of UCMHs today (dimensionless)
!         k               wavenumber of perturbation (Mpc^-1; 0.5d3/lambda)

      double precision function dsucmh_f_powerlaw(n,mucmh,z_i,z_c,zstop,rough,k)

      implicit none

      include 'dsmpconst.h'

      logical, intent(IN) :: rough
      double precision, intent(IN) :: n, mucmh, z_i, z_c, zstop
      double precision, intent(OUT) :: k
      double precision, parameter :: third = 1.d0/3.d0, twothird = 2.d0/3.d0      
      double precision :: mHor, deltam, deltamin, dsucmh_dmin, dsucmh_beta_powerlaw_running, zstart

      !Set zstart, the redshift at which linear accretion begins, to the lesser of zeq and z_i
      zstart = merge(zeq,z_i,zeq.lt.z_i)

      !Work out mass scale of perturbations (assumed constant from Horizon entry until equality)      
      deltam = mucmh * (zstop + 1.d0)/(zstart + 1.d0)  !Note that mucmh is DM+baryons, but deltam is just DM as the baryons do not collapse with the DM

      !Work out length scale of perturbation
      k = 1.d3*(1.3d2/deltam*omegacdmh2/0.112d0)**third

      !Calculate the minimum density contrast required to form a UCMH
      deltamin = dsucmh_dmin(k,z_c)

      !Work out the horizon mass for the scale in question, if needed for the rough beta calculation.
      mHor = merge(Mhzeq**third * (deltam / dmfrac)**twothird, 0.d0, rough) !Using Eqs 2 & 3 in Scott & Siversston (arXiv v5)

      !Calculate the UCMH fraction resulting from the specified spectral index     
      dsucmh_f_powerlaw = dsucmh_beta_powerlaw_running(n,0.d0,k,deltamin,third,mHor,rough) * dmfrac*(zstart + 1.d0)/(zstop + 1.d0)

      end function dsucmh_f_powerlaw


