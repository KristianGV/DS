!Returns the expected mass of a UCMH produced in a phase transition.
!Note that the result is DM+baryons, but the internal deltam is just DM as the baryons do not collapse with the DM
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
!
!Inputs: transition  EW, QCD or ee phase transition
!        zstop       redshift at which UCMHs cease to grow (dimensionless)
!Output: Mass of UCMH (M_sun)
      double precision function dsucmh_mass_ptrans(transition, zstop)    
      character(len=*), intent(IN) :: transition
      double precision, intent(IN) :: zstop

      if (trim(transition) .eq. "EW") then
        deltam=5.6d-19      ! mass scale of primordial perturbation (solar masses) (EW)
      else if (trim(transition) .eq. "QCD") then
        deltam=1.1d-9       ! mass scale of primordial perturbation (solar masses) (QCD) 
      else if (trim(transition) .eq. "ee") then
        deltam=0.33d0       ! mass scale of primordial perturbation (solar masses) (ee)
      else 
        write(*,*) 'Unrecognised phase transition in dsucmh_mass_ptrans.'
        stop
      endif

      dsucmh_mass_ptrans = deltam * (zeq + 1.d0) / (zstop + 1.d0) !Mh_UCMH in solar masses

      end function dsucmh_mass_ptrans

