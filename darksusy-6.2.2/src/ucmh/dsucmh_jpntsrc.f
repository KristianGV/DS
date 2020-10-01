!Calculate the J value for point source detection of a single UCMH at distance r_0
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
!
!Inputs: incut  inner radius of UCMH (kpc)
!        inmax  outer radius of UCMH (kpc)
!Output: J factor (cm^-5 GeV^2)

      double precision function dsucmh_jpntsrc(incut,inmax)

      implicit none
      include 'dsucmh.h'
      include 'dsmpconst.h'

      integer :: ptsperprofile, trans_index, prof_index, IER
      double precision :: dsquared, incut, inmax, cutrad, maxrad, TSINTL, Rh_local
      double precision :: radii(3,3,10000), products(3,3,10000), derivs(3,3,10000), sigma(3,3,9999)
      character :: trans*2, prof*3 
      common/dsucmhcom2/radii,products,derivs,sigma
      save/dsucmhcom2/
      common/dsucmhcom3/ptsperprofile
      save/dsucmhcom3/

      cutrad = incut * cmperpc * 1.d3                      !inner radius of ucmh in cm
      maxrad = inmax * cmperpc * 1.d3                      !outer radius of ucmh in cm
      dsquared = r_0*r_0*cmperpc*cmperpc*1.d6              !distance squared in cm

      if (cutrad .gt. maxrad) then

        !Whole halo is truncated - return constant density solution
        dsucmh_jpntsrc = maxdens**2 * maxrad**3 / (3.d0 * dsquared)
        return

      else

        if (ucmhprofile .eq. 'plain') then
          Rh_local = Rh_ucmh * 1.d3 * cmperpc                                  !Rh in cm
          dsucmh_jpntsrc = 3.d0 * (dmfrac * GeVperSolarMass * Mh_ucmh)**2
          dsucmh_jpntsrc = dsucmh_jpntsrc / (128.d0 * pi * pi * ((Rh_local*cutrad)**1.5d0 - cutrad**3))
          dsucmh_jpntsrc = (dsucmh_jpntsrc + maxdens**2 * cutrad**3 / 3.d0) / dsquared
          return
        endif

        trans = ucmhprofile(1:2)
        select case (trans)
        case ('EW')
          trans_index = 1
        case ('QC')
          trans_index = 2
        case ('ee')
          trans_index = 3
        case default
          write(*,*) 'Error in dsucmh_numdens_ptrans: unknown profile ', ucmhprofile
          stop
        end select

        prof = ucmhprofile(len(trim(ucmhprofile))-2:)
        select case (prof)
        case ('_01')
          prof_index = 1
        case ('_11')
          prof_index = 2
        case ('unc')
          prof_index = 3
        case default
          write(*,*) 'Error in dsucmh_numdens_ptrans: unknown profile ', ucmhprofile
          stop
        end select

        dsucmh_jpntsrc = TSINTL(cutrad,maxrad,ptsperprofile,radii(trans_index,prof_index,:),products(trans_index,prof_index,:),
     &    derivs(trans_index,prof_index,:),sigma(trans_index,prof_index,:),IER)
        if (IER .lt. 0) then
          write(*,*) 'Quitting because TSINTL returned error number: ',IER
          stop
        endif

        !convert from cm^3 kg^2 m^-6 to GeV^2 cm^-3
        dsucmh_jpntsrc = dsucmh_jpntsrc * 5.60958912d20 * 5.60958912d20

        !Add in the central core and divide by distance squared
        dsucmh_jpntsrc = (dsucmh_jpntsrc + maxdens**2 * cutrad**3 / 3.d0) / dsquared

      endif

      end function

