!Initialisation routine for interpolator in limits on PBH fractions
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: Oct 2014

      subroutine dspbh_init

      implicit none

      include 'dsucmh.h'

      integer :: j
      
      open(unit=2,file="PBH_beta_mass_limits.dat",action='READ')
      read(2,*)
      do j = 1, PBHfileLen
        read(2,*) PBH_masses(j), PBH_betamaxes(j)
      enddo
      close(2)     

      end subroutine
