!DarkSUSY example program to calculate things about UCMHs
!
!Authors: Pat Scott
!          pscott@imperial.ac.uk
!         Madeleine Anthonisen
!          maddyanthonisen@hotmail.com
!Date: Oct 2014 - April 2015

      program dstest_ucmh

      implicit none

      include 'dsmpconst.h'

      double precision :: mchi, sigchi,  bfs(29)
      double precision :: sensitivity, CL_ptsrc, skyfrac, ethreshold
      double precision :: zstop, z_c, z_i, dist, reqCL
      character(len=6) :: profile
      logical          :: quiet
      double precision :: n, alpha, k_s, log10P_R    

      double precision :: R, mucmh, k, flux, Eflux, rucmh, rcut, rmin
      double precision :: fmax(4), fmax_best, fneq1_15, f_dummy, zeroin
      double precision :: nmax_constant, nmax_running, pmax, log10P_R_max, log10P_delta_max, log10M3_h_max, log10M3_f_max, gmu_best 
      double precision :: dsucmh_fmax_gamma_galptsrc, dsucmh_fmax_gamma_exgalptsrc, dsucmh_fmax_gamma_diffuse, dsucmh_fmax_reion
      double precision :: dsucmh_f_powerlaw, dsucmh_f_powerlaw_running, dsucmh_f_powerlaw_step, dsucmh_f_generalised
      double precision :: diff_fmax_powerlaw, diff_fmax_powerlaw_running, diff_fmax_powerlaw_step, diff_fmax_generalised
      double precision :: diff_fmax_positive_ng_hierarchical, diff_fmax_positive_ng_feeder
      double precision :: dsucmh_f_cs, mdummy, rcore, f_cs, f_diff, Gmu,diff_fmax_cosmic_strings
      double precision :: alpha_cs, beta, gamma, N_cs, K_cs, sigma
      external diff_fmax_powerlaw, diff_fmax_powerlaw_running, diff_fmax_powerlaw_step, diff_fmax_generalised
      external diff_fmax_positive_ng_hierarchical, diff_fmax_positive_ng_feeder, diff_fmax_cosmic_strings
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_gaussps/alpha,k_s,n
      common /variables_for_ng/log10P_R
      common /variables_for_cs/R,mchi,sigchi,bfs,sensitivity,CL_ptsrc,
     & skyfrac,ethreshold,dist,reqCL,alpha_cs,beta,gamma,N_cs,K_cs,sigma,quiet,profile

      integer :: i
      integer, parameter :: imax = 20
      double precision, parameter :: R_cs_min = 1.d-14, R_cs_max=0.14   !Cosmic string radii to cover; kpc
      double precision, parameter :: mucmh_min = 4.d-4, mucmh_max=1.d2 !UCMH masses to cover; solar masses (constraints can cover 1.4d-9 to 1.d12)
      integer :: pdgann

      !WIMP model
      mchi = 10.d0                  !WIMP mass; GeV
      sigchi = 3.d-26               !Annihlation cross-section; cm^3/s
c      bfs(25) = 1.d0                !bbar
      pdgann=5                  ! b b-bar
      

      !Detector    
      sensitivity = 4.d-9           !Fermi 1yr point source limits @5sigma 
      CL_ptsrc = 1.d0-5.73303d-7    !Fermi 1yr point source limits @5sigma
      skyfrac = 1.0d0               !Fraction of the sky probed (all-sky survey mode)
      ethreshold = 0.1d0            !Threshold energy for calculating integrated gamma-ray fluxes; GeV

      !UCMH structure and analysis details
      zstop = 10.d0                 !Redshift of cessation of accretion of additional matter by UCMHs
      z_c = 1000.d0                 !(Latest) collapse redshift of UCMHs
      z_i = zeq                     !Initial redshift of UCMH seed creation (if z_i > zeq, this is just set to zeq)
      dist = 0.1d0                  !Example UCMH distance from Earth for flux calculation; kpc
      reqCL = 0.954d0               !Requested confidence level of limits
      profile = 'plain'             !'plain' or something like ['EW','QCD','ee']_['01','11','unc']
      quiet = .true.                !Suppress warnings about the limit on f being >1 and therefore invalid.

      !Spectral Model
      alpha = 1.d-2                 !Slope of spectral index assumed for running spectra and step spectra
      n = 0.968d0                   !Spectral index of perturbations assumed for step spectra
      k_s = 1.d-7                   !Assumed scale at which step takes place for step spectra
      log10P_R = -7.1d0             !log10(Assumed amplitude of Gaussian part of curvature perturbation for non-Gaussian spectra)

      !Cosmic string parameters
      alpha_cs = 0.05d0             !Ratio of loop length to horizon scale at formation (dimensionless)
      beta = 2.0d0*pi               !Ratio of loop length to loop radius (dimensionless)
      gamma = 10.d0*pi              !Cosmic string loop decay constant (dimensionless)
      N_cs = 40.d0                  !Normalisation constant from cosmic string simulations (number of loops created per horizon volume) 
      K_cs = 1000.d0                !Number of loop radii a loop is permitted to travel before UCMH formation is deemed impossible
      sigma = 0.3d0                 !Root of the variance of the cosmic string loop velocity distribution (fraction of c; value from arXiv:1107.2751)

      !Start up DarkSUSY
      call dsinit
      write(*,*)

!Set up the WIMP model being used.  Load a SUSY model instead as per standard DarkSUSY examples if that's what you prefer.
      call dsgivemodel_generic_wimp(mchi,.true.,sigchi,pdgann,1.d-7)
c      call dsansetup_wimp(mchi, sigchi, bfs)

      !If you're interested in UCMHs from phase transitions, initialise the numerical profiles for them.
      call dsucmh_initdens_ptrans
      
      write(*,*) 'Calculating UCMH limits on various cosmological power spectra.'
      write(*,'(A82,A60)') " M_0(M_Sun)    k(Mpc^-1)     f_max         n_max(PL)     n_max(run.PL) p_max(step)",
     & "P_R_max       P_delta_max   M_3_max(hrcl) M_3_max(feeder)"
      !Loop over a few UCMH masses, just for example purposes.  For UCMHs from cosmic strings, looping over R is more appropriate (see below).
      do i=1,imax        
        mucmh = 10.d0**(log10(mucmh_min) + dble(i-1)/dble(imax-1)*log10(mucmh_max/mucmh_min))

        !1. Initialise the UCMH profile, retrieving its inner and outer radii
        call dsucmh_init_profile(profile, z_i, z_c, zstop, mchi, sigchi, mucmh, rucmh, rmin, rcut)

        !2. Get the gamma-ray fluxes coming from such a UCMH, if it is situated at a distance dist from Earth
        call dsucmh_flux(mucmh, rucmh, max(rcut,rmin), dist, ethreshold, flux, Eflux)

        !3a. Get the limit from Galactic point source searches on the fraction of DM in UCMHs of mass mucmh, at CL reqCL
        fmax(1) = dsucmh_fmax_gamma_galptsrc(flux, dist, skyfrac, sensitivity, CL_ptsrc, reqCL, quiet)
        !3b. Get the limit from extra-galactic point source searches on the fraction of DM in UCMHs of mass mucmh, at CL reqCL
        fmax(2) = dsucmh_fmax_gamma_exgalptsrc(flux, dist, skyfrac, sensitivity, CL_ptsrc, reqCL, quiet)
        !3c. et the limit from diffuse emission from UCMHs in the Milky Way on the fraction of DM in UCMHs of mass mucmh, at CL reqCL
        fmax(3) = dsucmh_fmax_gamma_diffuse(Eflux, dist, skyfrac, reqCL, quiet)
        !3d. Get the approximate limit from reionisation on the fraction of DM in UCMHs of mass mucmh (does not require step 2)
        fmax(4) = dsucmh_fmax_reion(mchi, z_i, zstop, quiet)
        !3e. Choose the best limit
        fmax_best = minval(fmax)

        !4a. Find the corresponding limit on a constant spectral index
        nmax_constant = zeroin(0.5d0,1.8d0,diff_fmax_powerlaw,1.d-5)  
        !4b. Find the corresponding limit on a running spectral index
        nmax_running = zeroin(0.5d0,1.8d0,diff_fmax_powerlaw_running,1.d-5)        
        !4c. Find the corresponding limit on a step in the power spectrum
        pmax = zeroin(1.d-1,1.d2,diff_fmax_powerlaw_step,1.d-5)
        !4d. Find the corresponding limit on the generalised curvature power and amplitude of density perturbations
        log10P_R_max = zeroin(0.d0,-12.d0,diff_fmax_generalised,1.d-5)
        f_dummy = dsucmh_f_generalised(log10P_R_max,mucmh,z_i,z_c,zstop,k,log10P_delta_max) !(This is just log10P_R_max-->log10P_delta_max)
        !4e. Find the corresponding limit on the degree of positive non-Gaussianity, assuming hierarchical scaling
        log10M3_h_max = zeroin(-5.d0,-0.4d0,diff_fmax_positive_ng_hierarchical,1.d-5)
        !4f. Find the corresponding limit on the degree of positive non-Gaussianity, assuming feeder scaling
        log10M3_f_max = zeroin(-5.d0,-0.4d0,diff_fmax_positive_ng_feeder,1.d-5)

        !5. Compute the UCMH fraction given a power law spectrum of perturbations with e.g. index n = 1.15 (only requires step 1)
        fneq1_15 = dsucmh_f_powerlaw(1.15d0,mucmh,z_i,z_c,zstop,.false.,k)

        write(*,'(10E14.7)') mucmh, k, fmax_best, nmax_constant, nmax_running, pmax, 
     &                       10.d0**log10P_R_max, 10.d0**log10P_delta_max, 10.d0**log10M3_h_max, 10.d0**log10M3_f_max

      enddo

      write(*,*)
      write(*,*) " Calculating UCMH limits on cosmic strings."      
      write(*,*) " R (kpc)       M_0 (M_Sun)      Gmu_max"      
      !Loop over R
      do i=1,imax        
        R = 10.d0**(log10(R_cs_min) + dble(i-1)/dble(imax-1)*log10(R_cs_max/R_cs_min))
        !Root-finding method to constrain G\mu for each R 
        gmu_best = zeroin(1.d-11,1.d-5,diff_fmax_cosmic_strings,1.d-12)       
        write(*,'(3E14.7)') R, mucmh, gmu_best
      enddo

      end program dstest_ucmh


      double precision function diff_fmax_powerlaw(n)
      implicit none 
      double precision, intent(IN) :: n
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, alpha, k_s, n_ignore, dsucmh_f_powerlaw
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_gaussps/alpha,k_s,n_ignore
      diff_fmax_powerlaw = dsucmh_f_powerlaw(n,mucmh,z_i,z_c,zstop,.false.,k) - fmax_best
      end function diff_fmax_powerlaw


      double precision function diff_fmax_powerlaw_running(n)
      implicit none 
      double precision, intent(IN) :: n
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, alpha, k_s, n_ignore, dsucmh_f_powerlaw_running
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_gaussps/alpha,k_s,n_ignore
      diff_fmax_powerlaw_running = dsucmh_f_powerlaw_running(n,alpha,mucmh,z_i,z_c,zstop,k) - fmax_best
      end function diff_fmax_powerlaw_running


      double precision function diff_fmax_powerlaw_step(p)
      implicit none 
      double precision, intent(IN) :: p
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, alpha, k_s, n, dsucmh_f_powerlaw_step
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_gaussps/alpha,k_s,n
      diff_fmax_powerlaw_step = dsucmh_f_powerlaw_step(n,alpha,p,k_s,mucmh,z_i,z_c,zstop,k) - fmax_best
      end function diff_fmax_powerlaw_step


      double precision function diff_fmax_generalised(log10P_R)
      implicit none 
      double precision, intent(IN) :: log10P_R
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, log10P_delta, dsucmh_f_generalised
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      diff_fmax_generalised = dsucmh_f_generalised(log10P_R,mucmh,z_i,z_c,zstop,k,log10P_delta) - fmax_best
      end function diff_fmax_generalised


      double precision function diff_fmax_positive_ng_hierarchical(log10M_3)
      implicit none 
      double precision, intent(IN) :: log10M_3
      double precision, parameter :: ln10 = log(10.d0)
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, log10P_R, log10P_delta, dsucmh_f_ng_hierarchical
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_ng/log10P_R
      diff_fmax_positive_ng_hierarchical = dsucmh_f_ng_hierarchical(log10P_R,exp(ln10*log10M_3),mucmh,z_i,z_c,zstop,k,log10P_delta) 
     & - fmax_best
      end function diff_fmax_positive_ng_hierarchical


      double precision function diff_fmax_negative_ng_hierarchical(log10M_3)
      implicit none 
      double precision, intent(IN) :: log10M_3
      double precision, parameter :: ln10 = log(10.d0)
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, log10P_R, log10P_delta, dsucmh_f_ng_hierarchical
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_ng/log10P_R
      diff_fmax_negative_ng_hierarchical = dsucmh_f_ng_hierarchical(log10P_R,-exp(ln10*log10M_3),mucmh,z_i,z_c,zstop,k,log10P_delta) 
     & - fmax_best
      end function diff_fmax_negative_ng_hierarchical


      double precision function diff_fmax_positive_ng_feeder(log10M_3)
      implicit none 
      double precision, intent(IN) :: log10M_3
      double precision, parameter :: ln10 = log(10.d0)
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, log10P_R, log10P_delta, dsucmh_f_ng_feeder
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_ng/log10P_R
      diff_fmax_positive_ng_feeder = dsucmh_f_ng_feeder(log10P_R,exp(ln10*log10M_3),mucmh,z_i,z_c,zstop,k,log10P_delta) 
     & - fmax_best
      end function diff_fmax_positive_ng_feeder


      double precision function diff_fmax_negative_ng_feeder(log10M_3)
      implicit none 
      double precision, intent(IN) :: log10M_3
      double precision, parameter :: ln10 = log(10.d0)
      double precision :: fmax_best, mucmh, z_i, z_c, zstop, k, log10P_R, log10P_delta, dsucmh_f_ng_feeder
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_ng/log10P_R
      diff_fmax_negative_ng_feeder = dsucmh_f_ng_feeder(log10P_R,-exp(ln10*log10M_3),mucmh,z_i,z_c,zstop,k,log10P_delta) 
     & - fmax_best
      end function diff_fmax_negative_ng_feeder


      !Calculate the difference between the observed limit on the differential fraction of DM mass in UCMHs, and the predictions of cosmic strings
      double precision function diff_fmax_cosmic_strings(Gmu)
      implicit none
      double precision, intent(IN) :: Gmu
      double precision :: mchi, sigchi,  bfs(29)
      double precision :: sensitivity, CL_ptsrc, skyfrac, ethreshold
      double precision :: zstop, z_i, z_c, dist, reqCL, alpha_cs, beta, gamma, N_cs, K_cs, sigma
      character(len=6) :: profile
      logical          :: quiet
      double precision :: R, mucmh, k, flux, Eflux, rucmh, rcut, rmin, fmax(4), fmax_best, rcore, df_cs, df_max
      double precision :: dsucmh_fmax_gamma_galptsrc, dsucmh_fmax_gamma_exgalptsrc, dsucmh_fmax_gamma_diffuse, dsucmh_fmax_reion
      double precision :: dsucmh_f_cs
      common /variables_ucmh/z_i,z_c,zstop,fmax_best,mucmh,k
      common /variables_for_cs/R,mchi,sigchi,bfs,sensitivity,CL_ptsrc,skyfrac,ethreshold,dist,reqCL,
     &                          alpha_cs,beta,gamma,N_cs,K_cs,sigma,quiet,profile
      integer :: regime
      double precision:: ronteq, gmu_c1, gmu_c2, gmu_c3, gmu_c4,gmu_c5, gmu_c6, gmu_c7, gmu_c8, gmu_c9
     
        !Find the UCMH mass (and collapse regime) for the chosen G\mu, R and zstop
      call dsucmh_mass_cs(R, Gmu, z_c, zstop, alpha_cs, beta, gamma,
     &     mucmh, ronteq, regime, gmu_c1, gmu_c2, gmu_c3, gmu_c4, gmu_c5,
     &     gmu_c6, gmu_c7, gmu_c8, gmu_c9)
        if (regime .gt. 5) stop 'This example is not valid for loops formed during matter domination, as z_i is assumed >= z_eq.'       

        !Initialise the UCMH profile, retrieving its inner and outer radii
        call dsucmh_init_profile(profile, z_i, z_c, zstop, mchi, sigchi, mucmh, rucmh, rmin, rcut)
        !Choose the larger of R, rcut and rmin as the actual core radius
        rcore = max(R,max(rcut,rmin))
        !Get the gamma-ray fluxes coming from such a UCMH, if it is situated at a distance dist from Earth.
        call dsucmh_flux(mucmh, rucmh, rcore, dist, ethreshold, flux, Eflux)
        
        !Get the limit from Galactic point source searches on the fraction of DM in UCMHs of mass mucmh, at CL reqCL
        fmax(1) = dsucmh_fmax_gamma_galptsrc(flux, dist, skyfrac, sensitivity, CL_ptsrc, reqCL, quiet)
        !Get the limit from extra-galactic point source searches on the fraction of DM in UCMHs of mass mucmh, at CL reqCL
        fmax(2) = dsucmh_fmax_gamma_exgalptsrc(flux, dist, skyfrac, sensitivity, CL_ptsrc, reqCL, quiet)
        !Get the limit from diffuse emission from UCMHs in the Milky Way on the fraction of DM in UCMHs of mass mucmh, at CL reqCL
        fmax(3) = dsucmh_fmax_gamma_diffuse(Eflux, dist, skyfrac, reqCL, quiet)
        !Get the approximate limit from reionisation on the fraction of DM in UCMHs of mass mucmh
        fmax(4) = dsucmh_fmax_reion(mchi, z_i, zstop, quiet)
        !Choose the best limit 
        fmax_best = minval(fmax)
        !Turn this limit into a differential df/dM
        df_max = fmax_best / mucmh

        !Compute the UCMH differential fraction for the chosen G\mu and R from cosmic string loops
        df_cs = dsucmh_f_cs(R, Gmu, z_c, zstop, alpha_cs, beta, gamma, N_cs, K_cs, sigma)
        if (regime .eq. 0) df_cs = 0.d0
               
        !Compute the difference between the observational limits and the theoretical predictions of cosmic strings
        diff_fmax_cosmic_strings = df_cs - df_max
  
      end function diff_fmax_cosmic_strings

