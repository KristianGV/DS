!Limit on cosmological fraction of UCMHs from reionisation.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
! 
!Input:   mchi            WIMP mass (GeV)
!         zstop           redshift at which UCMHs cease to grow (dimensionless)
!         quiet           Suppress warnings
!
!Output:  fmax            maximum cosmological fraction of UCMHs today (dimensionless)

      double precision function dsucmh_fmax_reion(mchi, z_i, zstop, quiet)

      implicit none
      include "dsmpconst.h"

      logical, intent(IN) :: quiet
      double precision, intent(IN) :: mchi, z_i, zstop
      double precision :: massGrowthFactor, zstart

      !Set zstart, the redshift of creation of the UCMH seed (e.g. string loop), to the lesser of zeq and z_i
      zstart = merge(zeq,z_i,zeq.lt.z_i)

      !Calculate the factor by which UCMHs grow in mass from equality to zstop
      massGrowthFactor = (zstart + 1.d0) / (zstop + 1.d0)
 
      !Work out the corresponding maximum fraction of DM contained in UCMHs
      !dsucmh_fmax_reion = massGrowthFactor * 1.d-4*(mchi/100.d0)**1.3d0 !old version (1011.1935v1) 
      dsucmh_fmax_reion = massGrowthFactor * 1.d-2*(mchi/100.d0)         !new version (1011.1935v2)

      !Make sure that fmax is less than or equal to 1 at equality (~where the limit is derived), as UCMHs 
      !must be less than 100% of the dark matter. Otherwise, the limit is invalid for this point, so throw an error.
      if (dsucmh_fmax_reion .gt. massGrowthFactor .and. .not. quiet) then 
        write(*,*) 'Warning: fmax in dsucmh_fmax_reion = ', dsucmh_fmax_reion, '> ', massGrowthFactor,
     &   ', which is the growth factor since equality -- so limits are invalid.'
      endif

      end function dsucmh_fmax_reion

