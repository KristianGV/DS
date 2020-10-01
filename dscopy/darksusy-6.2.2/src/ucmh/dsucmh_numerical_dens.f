!Return the dark matter density a given height above the centre of a UCMH.
!Replace this with a call to your own density function to 
!include other profiles, like adiabatically-contracted ones.
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
!
!Input:  r  height above centre of UCMH (kpc)
!Output: dark matter density (GeV cm^-3)

        double precision function dsucmh_numerical_dens(r)

        implicit none

        double precision :: r, dsucmh_numdens_ptrans
        external dsucmh_numdens_ptrans

        dsucmh_numerical_dens = dsucmh_numdens_ptrans(r)

        end function
