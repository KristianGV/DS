      program DDCR_limits
c
c     This program is an example DarkSUSY main program to calclulate various
c     limits on the DM-nucleon scattering cross section that result from
c     the cosmic-ray induced high-velocity flux of DM particles.
c
c     Author: Torsten Bringmann, 2018-11-23
c             mod 2019-10-28 (optimized root finding)
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      include 'dsddcom.h'
      
      real*8 tstart,tfinish                                   ! aux
      integer ierr,iwarn,nfc, istatus, i, npoints
      real*8 sigvan, eps, stepsize, slopesign
      real*8 mdm, mmin, mmax
      real*8 sigsi, sighi, siglo, sigstart, sigtmp
      real*8 count, counthi, countlo, counttmp
      logical sigmamax, overshoot, undershoot


ccc functions
      real*8 dsddDMCRcountrate

ccc output files
      character*80 outfile, experiment
      

c     Here we include the file dsidtag.h which contains information
c     about which particle module is used.
      include 'dsidtag.h'


      call CPU_TIME(tstart)
      call dsinit
      write (*,*) '-------------------------------------------------------'
      write (*,*) 

c... This is how one could change the default treatment of soil absorption, 
c... to fully take into account 'any' Q2-dependence. See header of 
c... dsddTDMattenuation for more details!
c        attenuation_how = 2            

c... choose which experiment (and determine output filename correspondingly)
      experiment = 'Borexino_SD' ! see src/dd/dsddDMCRcountrate.f for options
      sigmamax = .false. ! set to .true. to obtain maximal sigma that can be
                         ! probed (due to soil absorption) 
      
      Deff = 8.02d0 ! kpc; eff distance out to which source term is integrated
                     ! Deff = 0.997 (8.02) kpc corresponds to 1 (10) kpc in real
                     ! distance
      if (sigmamax) then
         outfile='DMCR_limitsMAX'
      else
         outfile='DMCR_limits'      
      endif           
      if (Deff.eq.0.997d0) then
         outfile=outfile(1:index(outfile,' ')-1)//'_1kpc_'//
     &           experiment(1:index(experiment,' ')-1)//'.dat' 
      elseif (Deff.eq.8.02d0) then     
         outfile=outfile(1:index(outfile,' ')-1)//'_10kpc_'//
     &           experiment(1:index(experiment,' ')-1)//'.dat' 
      else
        write(*,*) 'Typo in specifying Deff?'
        write(*,*)
        stop
      endif
           
      open(unit=10,file=outfile)
      write(*,*) 'Determining limits and writing results to file ', outfile
      write(*,*) 'Please wait...'
      write(*,*)

c... scan setup
      mmin     = 1.0d-4  ! minimal and maximal DM masses (GeV)
      mmax     = 2.0d1    
      npoints  = 200
      eps      = 1.0d-4  ! required precision in countrate for limiting cross section
      sigvan   = 3.d-26  ! cm^3 s^-1; results do not depend on this choice
      sigstart = 1.0d-32 ! start value at mmin
      sigsi    = sigstart 
      count    = 1.0d-50 ! init 
    
      do i = 0,npoints-1
        counttmp = count 
        sigtmp = 0.001*sigsi
        sighi = 1.0d-10
        siglo = 1.0d-50
        overshoot = .false. 
        undershoot = .false.
        mdm = mmin*(mmax/mmin)**((1.*i)/(1.*npoints-1.))
        write(*,*) 'mdm = ', mdm
        stepsize = 5.d0 ! initial stepsize for sigsi
100     call dsgivemodel_generic_wimp(mdm,.true., sigvan,5,sigsi*1.0d36)
c... For Borexino, we want for better comparison in any case only SD couplings
c... (-> no CRDM component from scttering on He)
        if (experiment.eq.'Borexino'.or.experiment.eq.'Borexino_SD') then
          call dsgivemodel_generic_wimp_opt('spin',0.5d0) ! SD couplings require
                                                          ! DM particle with spin!
          call dsgivemodel_generic_wimp_opt('sigsip',0.0d36)
          call dsgivemodel_generic_wimp_opt('sigsin',0.0d36)
          call dsgivemodel_generic_wimp_opt('sigsdp',sigsi*1.0d36)
        endif  
        call dsmodelsetup(ierr,iwarn)
        count = dsddDMCRcountrate(experiment)

        if (count.gt.1.0d0) then
          overshoot = .true. 
          sighi = sigsi
          counthi = count
        endif
        slopesign = (count-counttmp)/log(stepsize)! are we on a rising or a falling slope?
        if (count.lt.1.0d0.and.(.not.undershoot).and.((sigmamax.and.slopesign.lt.0.0d0)
     &      .or.((.not.sigmamax).and.slopesign.gt.0.0d0))) then
          undershoot = .true. 
          siglo = sigsi
          countlo = count          
        endif  
        
c... change stepsize adaptively if accuracy goal is not yet met
c        if (abs(count-1.0d0).gt.eps) then
        if ((abs(sighi-siglo)/siglo).gt.eps) then
        
          if (.not.(overshoot.and.undershoot)) then ! don't (yet) enter enter bisection
            if (sigmamax) then
              if ((stepsize.gt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp).or.
     &            (stepsize.lt.1.d0.and.count.gt.1.0d0).or.
     &            (stepsize.lt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp).or.
     &            (stepsize.gt.1.d0.and.count.lt.1.d-6)) 
     &             stepsize=1./stepsize**0.61 ! change search direction and decrease  stepsize
            else
              if ((stepsize.lt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp).or.
     &            (stepsize.gt.1.d0.and.count.gt.1.0d0).or.
     &            (stepsize.gt.1.d0.and.count.lt.1.0d0.and.count.lt.counttmp))
     &             stepsize=1./stepsize**0.71
            endif
            sigsi = sigsi*stepsize
            counttmp = count
            if (abs(1.0-stepsize).lt.eps) goto 200 ! no solution found...              
          else ! bisect
            if (count.gt.1.0d0) then
              sighi = sigsi
              counthi = count
            else
              siglo = sigsi   
              countlo = count         
            endif
            sigsi = sqrt(siglo*sighi)
          endif     
          goto 100 ! re-calculate countrate with new value of sigsi
        endif
      
 200    if (overshoot.and.undershoot) write(10,*) mdm, sigsi
        if (overshoot.and.undershoot) then
          write(*,*) mdm, sigsi
        else
          write(*,*) 'Did not manage to find sigsi for DM mass =', mdm
          write(*,*) '[continue with next mass]'
          sigsi = sigstart       
        endif   
        
 220    continue 
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

      
