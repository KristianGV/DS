!Return the dark matter density a given height above the centre of a UCMH created in a phase transition
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2009
!
!Input:  r  height above centre of UCMH (kpc)
!Output: dark matter density (GeV cm^-3)

      double precision function dsucmh_numdens_ptrans(r)

      implicit none
      include 'dsucmh.h'
      include 'dsmpconst.h'

      double precision :: r, logr, weight, testrad
      integer :: trans_index, prof_index, ptsperprofile, lowerindex, topindex, testindex
      double precision :: logradii(3,3,10000), logdensities(3,3,10000)
      character :: trans*2, prof*3 
      common/dsucmhcom1/logradii,logdensities
      save/dsucmhcom1/
      common/dsucmhcom3/ptsperprofile
      save/dsucmhcom3/

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

        logr = log10(r * 1.d3 * cmperpc)

        if (logr .lt. logradii(trans_index,prof_index,1)) then
          dsucmh_numdens_ptrans = 10.d0**logdensities(trans_index,prof_index,1)
          return
        endif

        !Get the interpolated density using a binary search and linear interpolation
        topindex = ptsperprofile
        lowerindex = 1
        if (logr .gt. logradii(trans_index,prof_index,topindex)) then
          weight = 0.d0
          lowerindex = ptsperprofile
        else
          do
            testindex = (topindex+lowerindex)/2
            testrad = logradii(trans_index,prof_index,testindex)
            if (testrad .gt. logr) then
              topindex = testindex
            else
              lowerindex = testindex
            endif
            if (topindex .eq. lowerindex + 1) exit
          enddo
          weight = (logr - logradii(trans_index,prof_index,lowerindex))/
     &     (logradii(trans_index,prof_index,lowerindex+1) - logradii(trans_index,prof_index,lowerindex))
        endif
     
        dsucmh_numdens_ptrans = 10.d0**((1.d0-weight)*logdensities(trans_index,prof_index,lowerindex) +
     &   weight*logdensities(trans_index,prof_index,lowerindex+1))
             
      end function

