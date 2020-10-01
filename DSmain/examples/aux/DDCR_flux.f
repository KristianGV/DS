      program DDCR_flux
c
c     This program is an example DarkSUSY main program to calculate the high-
c     energy flux of DM particles that results from ISM cosmic rays scattering 
c     on the standard CDM population in the halo.
c
c     Author: Torsten Bringmann, 2018-07-05 (mod 2019-03-15)
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      include 'dsidtag.h'

      real*8 tstart,tfinish, tmp                                      ! aux
      integer ierr,iwarn,nfc, istatus
      real*8 sigvan, massdm(5), flux(5)
      real*8 sigsi, rholocal
      real*8 TDM, TDMmin, TDMmax
      integer npoints, i, j
      logical he_include

ccc functions
      real*8 dsddDMCRflux
      
ccc output files
      character*80 outfile(2)   
      

      call CPU_TIME(tstart)
      call dsinit


c... whether or not to include the CR He component
      he_include = .true.

      if (he_include) then
        outfile(1)='DMCR_flux.dat'
      else  
        outfile(1)='DMCR_flux_noHe.dat' 
      endif    
       
c... fiducial scattering and annihilation cross sections; 
c... result will not depend on these
      sigsi  = 1.0d-30 ! cm^2
      sigvan = 3.d-26  ! cm^3 s^-1
c... define DM masses to be tabulated (in GeV)
      massdm(1) = 1.0d-3
      massdm(2) = 1.0d-2
      massdm(3) = 1.0d-1
      massdm(4) = 1.0d0
      massdm(5) = 1.0d1
c... range of kinetic energies to tabulate [GeV]
      TDMmin = 1.d-10
      TDMmax = 1.0d1
      npoints = 1000

      open (unit=10,file=outfile(1))
      write(10,*) '# DM masses to be tabulated [GeV]'
      write(10,*) massdm
      write(10,*) 'Tkin [GeV]  | DM flux/cross section [cm^-4 s^-1 GeV^-1] '//
     &            'for DM masses in above order '
      write(10,*)

      do i=1,npoints
        TDM = TDMmin*(TDMmax/TDMmin)**((i-1.)/(npoints-1.))
        do j=1,5            
          call dsgivemodel_generic_wimp(massdm(j),.true., sigvan,5,sigsi*1.0d36)
          if (.not.he_include) then ! exclude He->DM by SD (rather than SI) scattering
            call dsgivemodel_generic_wimp_opt('spin',0.5d0) ! SD couplings require
                                                            ! DM particle with spin!
            call dsgivemodel_generic_wimp_opt('sigsip',0.0d36)
            call dsgivemodel_generic_wimp_opt('sigsin',0.0d36)
            call dsgivemodel_generic_wimp_opt('sigsdp',sigsi*1.0d36)
          endif  
          call dsmodelsetup(ierr,iwarn)
          rholocal=0.3d0 ! GeV/cm^2
          flux(j)=dsddDMCRflux(TDM, rholocal)/sigsi  
        enddo 
        write(10,*) TDM, flux
      enddo

      close(10)

500   call CPU_TIME(tfinish)
      write (*,*)
      write (*,*) '-------------------------------------------------------'
      write (*,*) 'The DarkSUSY example program has finished successfully.'
      write (*,*) 'Total time needed (in seconds): ',tfinish-tstart
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-------------------------------------------------------'
      write (*,*)
      stop
 999  end

      
