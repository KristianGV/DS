!Theoretical prediction of cosmological fraction of UCMHs from hierarchical non-Gaussianity.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2012
! 
!Input:   log10P_R        log_10 of the Gaussian amplitude of primordial curvature perturbations at scale k (dimensionless)
!         M_3             Dimensionless skewness of the perturbation amplitude probability distribution (ie. degree of non-Gaussianity)
!         mucmh           present-day UCMH mass (M_Sun)
!         z_c             latest allowed redshift of collapse in order to form UCMHs (dimensionless) 
!         zstop           redshift at which UCMHs cease to grow (dimensionless)
!
!Output:  f (return val)  cosmological fraction of UCMHs today (dimensionless)
!         k               wavenumber of perturbation (Mpc^-1; 0.5d3/lambda)
!         log10P_delta    log_10 of the amplitude of primordial density perturbations at scale k (dimensionless)

      double precision function dsucmh_f_ng_hierarchical(log10P_R, M_3, mucmh, z_i, z_c, zstop, k, log10P_delta)

      implicit none

      include 'dsmpconst.h'

      double precision, intent(IN) :: log10P_R, M_3, mucmh, z_i, z_c, zstop
      double precision, intent(OUT) :: k, log10P_delta
      double precision, parameter :: third = 1.d0/3.d0, log10pt191 = log10(0.191d0), log10pt191timespt81 = log10(0.81d0*0.191d0)  
      double precision :: deltam, deltamin, dsucmh_dmin, dsucmh_betaH16, log10correction, nu, dsucmh_mvar_generalised, zstart

      !Set zstart, the redshift at which linear accretion begins, to the lesser of zeq and z_i
      zstart = merge(zeq,z_i,zeq.lt.z_i)

      !Work out mass scale of perturbations (assumed constant from Horizon entry until equality)      
      deltam = mucmh * (zstop + 1.d0)/(zstart + 1.d0)  !Note that mucmh is DM+baryons, but deltam is just DM as the baryons do not collapse with the DM

      !Work out length scale of perturbation
      k = 1.d3*(1.3d2/deltam*omegacdmh2/0.112d0)**third

      !Calculate the minimum density contrast required to form a UCMH
      deltamin = dsucmh_dmin(k,z_c)

      !Work out the mass variance and convert it to nu parameter
      nu = deltamin/sqrt(dsucmh_mvar_generalised(log10P_R))

      !Calculate the UCMH fraction resulting from the specified Gaussian + non-Gaussian power     
      dsucmh_f_ng_hierarchical = dsucmh_betaH16(M_3,nu) * dmfrac * (zstart + 1.d0)/(zstop + 1.d0)

      log10correction = merge(log10pt191,log10pt191timespt81,k.gt.keq)
      log10P_delta = log10P_R + log10correction

      end function dsucmh_f_ng_hierarchical

