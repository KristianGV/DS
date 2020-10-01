!Theoretical prediction of cosmological fraction of UCMHs from running power-law Gaussian power spectra with a step.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
! 
!Input:   n               spectral index of primordial spectrum of perturbations (dimensionless)
!         alpha           slope of running spectral index (dimensionless) 
!         p               step size in primordial spectrum of density perturbations (dimensionless)
!         k_s             wavenumber at which spectrum has a step discontinuity (Mpc^-1)
!         mucmh           present-day UCMH mass (M_Sun)
!         z_c             latest allowed redshift of collapse in order to form UCMHs (dimensionless) 
!         zstop           redshift at which UCMHs cease to grow (dimensionless)
!
!Output:  f (return val)  cosmological fraction of UCMHs today (dimensionless)
!         k               wavenumber of perturbation (Mpc^-1; 0.5d3/lambda)

      double precision function dsucmh_f_powerlaw_step(n,alpha,p,k_s,mucmh,z_i,z_c,zstop,k)

      implicit none

      include 'dsmpconst.h'

      double precision, intent(IN) :: n, alpha, p, k_s, mucmh, z_i, z_c, zstop
      double precision, intent(OUT) :: k
      double precision, parameter :: third = 1.d0/3.d0      
      double precision :: deltam, deltamin, dsucmh_dmin, dsucmh_beta_powerlaw_step, zstart

      !Set zstart, the redshift at which linear accretion begins, to the lesser of zeq and z_i
      zstart = merge(zeq,z_i,zeq.lt.z_i)

      !Work out mass scale of perturbations (assumed constant from Horizon entry until equality)      
      deltam = mucmh * (zstop + 1.d0)/(zstart + 1.d0)  !Note that mucmh is DM+baryons, but deltam is just DM as the baryons do not collapse with the DM

      !Work out length scale of perturbation
      k = 1.d3*(1.3d2/deltam*omegacdmh2/0.112d0)**third

      !Calculate the minimum density contrast required to form a UCMH
      deltamin = dsucmh_dmin(k,z_c)

      !Calculate the UCMH fraction resulting from the specified spectral index and slope
      dsucmh_f_powerlaw_step = dsucmh_beta_powerlaw_step(n,alpha,p,k_s,k,deltamin,third)*dmfrac*(zstart + 1.d0)/(zstop + 1.d0)

      end function dsucmh_f_powerlaw_step

