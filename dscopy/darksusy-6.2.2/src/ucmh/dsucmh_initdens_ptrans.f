****************************************************************
*** Numerical density profile initialisation for UCMHs formed***
*** in phase transitions in the early Universe.              ***  
***                                                          ***
*** Author: Pat Scott (pat@fysik.su.se)                      ***      
*** Date: August 2009                                        ***
****************************************************************

      subroutine dsucmh_initdens_ptrans

      implicit none

      integer :: i, j, k, ntransitions, nprofiles, ptsperprofile, IER
      character :: fileprefix*300, filenames_trans(3)*3,
     &  filenames_prof(3)*3, filename*200
      double precision :: logradii(3,3,10000), logdensities(3,3,10000),
     & radii(3,3,10000), densities(3,3,10000)
      double precision :: products(3,3,10000), derivs(3,3,10000),
     &  sigma(3,3,9999), working(19998)
      common/dsucmhcom1/logradii,logdensities
      save/dsucmhcom1/
      common/dsucmhcom2/radii,products,derivs,sigma
      save/dsucmhcom2/
      common/dsucmhcom3/ptsperprofile
      save/dsucmhcom3/
    
      ntransitions = 3
      nprofiles = 3
      filenames_trans = ['EW ', 'QCD', 'ee ']
      filenames_prof = ['01 ', '11 ', 'unc']
      ptsperprofile = 10000
c      fileprefix = dsinstall//'/share/DarkSUSY/ucmh_numdens_ptrans/dens'
      call dsdatafile(fileprefix,'ucmh_numdens_ptrans/dens')
  
      do i = 1, ntransitions
       do j = 1, nprofiles
          filename = trim(fileprefix)//trim(filenames_trans(i))//trim(filenames_prof(j))//'.dat'
          open(unit=2,file=filename,action='READ')
          do k = 1, ptsperprofile
            read(2,*) radii(i,j,ptsperprofile+1-k), densities(i,j,ptsperprofile+1-k)
          enddo
          close(2)
          radii(i,j,:) = radii(i,j,:)*1.d2
          products(i,j,:) = radii(i,j,:)*radii(i,j,:)*densities(i,j,:)*densities(i,j,:)
          logdensities(i,j,:) = log10(densities(i,j,:)*5.60958912d20)
          logradii(i,j,:) = log10(radii(i,j,:))
        
          call TSPSI (ptsperprofile,radii(i,j,:),products(i,j,:),2,0,.false.,.false.,2*(ptsperprofile-1),
     &     working,derivs(i,j,:),sigma(i,j,:),IER) !init interpolator
          if (IER .lt. 0) then
            write(*,*) 'Quitting because TSPSI returned error number: ',IER
            stop
          endif
        enddo
      enddo

      end subroutine





