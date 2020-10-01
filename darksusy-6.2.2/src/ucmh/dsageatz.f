!Calculates the approximate age of the Universe at a given redshift. 
!Assumes a flat Universe; accurate to a few percent for z ~< 1000.  Matches Eq 16, astro-ph/0003463.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: April 2015
!
!Input:	  z          redshift
!Output:  dsageatz   age of the Universe at redshift z

      double precision function dsageatz(z)
      
      implicit none
      
      include 'dsmpconst.h'
      
      double precision, intent(IN) :: z
      double precision, parameter :: omega_m = 1.d4 * omegamh2 / H_0**2
      double precision, parameter :: Hubble_time = mperkpc / H_0
      double precision, parameter :: prefactor = 2.d0 * Hubble_time / (3.d0 * (1.d0 - omega_m)**0.5d0)

      dsageatz = prefactor * asinh(((1.d0/omega_m - 1.d0)/(1.d0 + z)**3)**0.5d0)
      
      end function dsageatz
