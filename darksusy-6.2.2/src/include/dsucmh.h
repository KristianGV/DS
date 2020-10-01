*         -*- mode: fortran -*-
! UCMH common block.
! Pat Scott (p.scott@imperial.ac.uk)
! Oct 11 2014

      integer, parameter :: PBHfilelen = 2590              !Length of PBH_beta_mass_limits.dat
      double precision :: PBH_masses(PBHfilelen)
      double precision :: PBH_betamaxes(PBHfilelen)
      double precision :: twoucmhsigma_sq,Rh_ucmh,Mh_ucmh,r_min,maxdens
      character(len=20) :: ucmhprofile
      common/ucmhparcom/PBH_masses,PBH_betamaxes,twoucmhsigma_sq,Rh_ucmh,Mh_ucmh,r_min,maxdens,ucmhprofile

      double precision :: rhcut,r_0
      common /ucmhhm/rhcut,r_0
      save /ucmhparcom/,/ucmhhm/
      
