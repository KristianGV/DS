!Works out the mimimum amplitude perturbation required to form a UCMH
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
! 
!Input:  k          wavenumber corresponding to scale of perturbations (Mpc^-1)
!        z          latest allowed redshift of collapse of UCMHs
!
!Output: dmin       minimum density contrast required to form a UCMH (dimensionless)

      double precision function dsucmh_dmin(k,z)

      implicit none

      include 'dsucmh.h'
      include 'dsmpconst.h'

      double precision, intent(IN) :: k, z
      double precision :: dsucmh_transFunc, dsucmh_mathcalT, dsucmh_aontsq
      external dsucmh_transFunc, dsucmh_mathcalT, dsucmh_aontsq
 
      dsucmh_dmin = 2.d0/9.d0 * (1.5d0*pi)**(2.d0/3.d0) 
      dsucmh_dmin = dsucmh_dmin * dsucmh_transFunc(1.d0) / (k*k*dsucmh_mathcalT(k)) * dsucmh_aontsq(z)

      end function dsucmh_dmin  



