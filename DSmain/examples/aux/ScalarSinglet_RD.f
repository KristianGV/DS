      program ScalarSinglet_RD
c
c     This program is an example DarkSUSY main program to calculate the 
c     relic density and kinetic freeze-out temperature for the Scalar 
c     Singlet model.
c
c     Author: Torsten Bringmann, 2018-05-26
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none

      real*8 oh2,xf                                           ! relic density
      real*8 tkd,mcut, mcut2                                  ! kinetic decoupling
      real*8 inputmass, inputlambda                           ! model parameters
      real*8 mmin, mmax, deltaM, oh2goal, doh2goal            ! scan 
      real*8 logdeltam, logdeltamres, stepsize
      integer i, iprint
      real*8 tstart,tfinish                                   ! aux
      integer ierr,iwarn,nfc, istatus

ccc functions
      real*8 dsrdomega, dskdmcut, dskdtkd

ccc output files
      character*80 outfile(2)   
      

c     Here we include the file dsparticles.h which contains common block
c     definitions for particle codes and masses, widths, dof etc.
      include 'dsparticles.h'
      include 'dssm.h'
      include 'dsidtag.h'
      include 'dskdcom.h'


      call CPU_TIME(tstart)
      call dsinit
      write (*,*) '-------------------------------------------------------'
      write (*,*) 

c... fix outpu file names
      outfile(1)='ScalarSinglet_RD_acc.dat'



      inputmass=60.0
      inputlambda=0.000667162


c... enter model parameters and initialize model
      call dsgivemodel_silveira_zee(inputlambda,inputmass)
      call dsmodelsetup(ierr,iwarn)
      write (*, *) 'Scalar Singlet model initialized with mS [GeV], lambda = ',
     &               inputmass, inputlambda           

c... check that there are no issues with the values of the input parameters
      if (ierr.ne.0.or.iwarn.ne.0) then
         write (*,*) 'Errorflags ierr, iwarn = ', ierr, iwarn
         write (*,*) 'Aborting program...'
         goto 500
      endif

      
c... calculate relic density
c      oh2=dsrdomega(0,1,xf,ierr,iwarn,nfc)
      oh2=dsrdomega(0,0,xf,ierr,iwarn,nfc)
      write(*,*)
      write(*,*) '  Relic density: Oh2 = ',oh2,ierr,iwarn
      write(*,*) '  Chemical decoupling at: T_f = ',mass(kdm)/xf,' GeV.'


c...  calculate  kinetic decoupling, and the smallest halo size
      tkd=dskdtkd(1) ! 1=full calculation
      mcut=dskdmcut(mass(kdm),tkd,1) ! 1 = true cut-off
      write(*,*) 
      write(*,*) '  Cutoff mass: M_cut/M_sun = ',mcut 
      write(*,*) '  Kinetic decoupling at: Tkd = ',tkd, ' MeV'
      write(*,*) 



c... scan setup
      mmin = 10.0d0       ! GeV
      mmax = 1000.0d0
      logdeltam = 0.07     ! log(10)-spacing between masses to tabulate
      logdeltamres = 0.005 ! same, in resonance region
      oh2goal = 0.112      ! required relic density
      doh2goal = 1.d-4     ! allowed error in oh2. 
                           ! NB: needs to be rather small for smooth lambda(mS) !
      inputlambda = 0.1    ! start value at mmin


      write(*,*) 
      write(*,50) mmin, mmax
 50   format('Now scanning over all DM masses from ',F6.1,' GeV to ',F6.1,' GeV.')
      write(*,*) 'Calculating lambda and Mcut for correct relic density,'
     &           //' please wait...'


      open (unit=10,file=outfile(1))
      write(10,*) 'mS [GeV]   |   lambda  |   Tkd [MeV] (default)   |  '
     &            //'Mcut (default)   |   Mcut (only light quarks)'
      write(10,*) '-----------+-----------+-------------------------+--'
     &            //'-------------------+----------------------------'
      write(10,*)


      inputmass = mmin  
      iprint = 0          
c... here we start the main loop over all DM masses
 90   continue 

        stepsize = 2.d0 ! inital stepsize for lambda

c... init model and calculate RD
 100    call dsgivemodel_silveira_zee(inputlambda,inputmass)
        call dsmodelsetup(ierr,iwarn)
        if (ierr.ne.0.or.iwarn.ne.0) then
c           write (*, *) 'mS [GeV], lambda = ', inputmass, inputlambda
c           write (*,*) 'Errorflags ierr, iwarn = ', ierr, iwarn
c           write (*,*) 'Aborting program...'
c           goto 500
        endif
c        oh2=dsrdomega(0,1,xf,ierr,iwarn,nfc)
        oh2=dsrdomega(0,0,xf,ierr,iwarn,nfc)

c... change stepsize adaptively if accuracy goal is not yet met
        if (abs(oh2-oh2goal).gt.doh2goal) then
          if (oh2.gt.oh2goal) then
            if (stepsize.lt.1.d0) stepsize=1/stepsize**0.85
          else
            if (stepsize.gt.1.d0) stepsize=1/stepsize**0.85
          endif
          if (abs(1.d0-stepsize).lt.1.d-6) then
            write(*,*) 'Did not manage to find lambda for mS =', inputmass
            write(*,*) '[continue with next mass]'
c            read(*,*)
            inputlambda=0.1
            goto 200
          endif
          inputlambda=inputlambda*stepsize
          goto 100 ! re-calculate RD with new value of lambda
        endif

c        write(*,*) 'TEST : ',inputmass,inputlambda,stepsize,oh2
        iprint = iprint + 1
        if (iprint.gt.20) then
          write(*,*) 'Tabulated lambda up to mS =', inputmass
          iprint=0
        endif  

c... now calculate kinetic decoupling for this value of lambda
        quark_how=2    ! only light quarks above 4*154 MeV
        tkd=dskdtkd(1) 
        mcut2=dskdmcut(mass(kdm),tkd,1)
        quark_how=1    ! default
        tkd=dskdtkd(1) 
        mcut=dskdmcut(mass(kdm),tkd,1)

        write(10,*) inputmass, inputlambda, tkd, mcut, mcut2

 200    continue

        if (inputmass.ge.mmax) goto 300
 
        if (inputmass.lt.35.d0.or.inputmass.gt.130.0d0) then
          inputmass = inputmass*10**logdeltam
        else ! need better resolution
          inputmass = inputmass*10**logdeltamres
        endif

      goto 90 ! loop over DM masses

300   write(*,*)
      write(*,*) 'Done: results written to ', 
     &             outfile(1)(1:index(outfile(1),' ')-1), '.'

500   close(10)
      call CPU_TIME(tfinish)


      write (*,*)
      write (*,*) '-------------------------------------------------------'
      write (*,*) 'The DarkSUSY example program has finished successfully.'
      write (*,*) 'Total time needed (in seconds): ',tfinish-tstart
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-------------------------------------------------------'
      write (*,*)
      stop
 999  end

      
