!Compute the square of the ratio of the scalefactor to the age of the Universe
!at a given redshift.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011
!
!Input:   z       redshift    (dimensionless)
!Output:  aontsq  a^2/(ct)^2  (Mpc^-2)

      double precision function dsucmh_aontsq(z)

      implicit none

      include 'dsmpconst.h'

      double precision, intent(IN) :: z
      integer, parameter :: nz = 7
      double precision, parameter ::
     &     zs(nz) = (/50.d0, 100.d0, 130.d0, 150.d0, 200.d0,
     &     500.d0, 1000.d0/)
      double precision, parameter ::
     &     aontsqs(nz) = (/1.78826666d-6, 3.6648633d-6,
     &     4.84663373d-6, 5.65783754d-6, 7.816d-6, 2.274d-5, 5.57d-5/)
      double precision :: weight
      integer :: zlow_index = 0

      !Less accurate at high z.
      !double precision :: ct, dsageatz
      !ct = dsageatz(z) * c_light / mperkpc ! Mpc
      !dsucmh_aontsq = (ct*(z+1.d0))**(-2)

      !More accurate for high z.
      if (z .lt. zs(1) .or. z .gt. zs(nz)) then
        write(*,*) 'z_c = ', z, ' outside permitted range in dsumch_aontsq.'
        stop
      endif
      call dshunt(zs,nz,z,zlow_index)
      if (zlow_index .eq. nz) then
        dsucmh_aontsq = aontsqs(nz)
      else
        weight = (zs(zlow_index+1) - z) / (zs(zlow_index+1) - zs(zlow_index))
        dsucmh_aontsq = exp(weight*log(aontsqs(zlow_index)) + (1.d0 - weight)*log(aontsqs(zlow_index+1)))
      endif

      end function dsucmh_aontsq

