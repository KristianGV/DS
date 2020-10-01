      program vdSIDM_RD
c
c     This program is an example DarkSUSY main program to calculate the 
c     self-interaction cross section and kinetic decoupling for thermally 
c     produced SIDM. A whole scan takes about 30 minutes to run on a laptop.
c
c     Author: Torsten Bringmann, 2018-05-19
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      include 'dsidtag.h'
      include 'dsmpconst.h'

      real*8 mDMin, mMEDin, gin, gstep, alpha                   ! model parameters
      real*8 oh2, oh2goal, doh2goal, xf                         ! relic density
      real*8 tkd,mcut                                           ! kinetic decoupling
      real*8 vkms, sigT                                         ! self-scattering
      real*8 mDMmin,mDMmax,mmedmin,mmedmax,mDMresmin,mDMresmax  !scan
      real*8 mmedresmin,stepmDM,stepmDMres,stepmmed,stepmmedres
      real*8 deltares, deltares2
      integer nres, npoints, i, scantype
      real*8 tstart,tfinish, tmp                                ! aux
      integer ierr,iwarn,nfc, istatus
      character*80 outfile(2)                                   ! output files
     
ccc functions
      real*8 dsrdomega, dskdmcut, dskdtkd, dssisigtmav

      
      call CPU_TIME(tstart)
      call dsinit

      scantype = 2 ! 1-vector, 2- scalar

c...  output file names
      if (scantype.eq.1) outfile(1)='sigmaT_2D_table.dat'
      if (scantype.eq.2) outfile(1)='sigmaT_2D_table_pwave.dat'


c... settings for observables
      vkms     = 30.0d0   ! 30 is for dwarf galaxy scales
      oh2goal  = 0.112    ! required relic density.
      doh2goal = 3.d-4    ! allowed error in oh2. 
                          ! NB: needs to be rather small for smooth curves!

c... scan setup
      mDMmin      = 0.5d0    ! GeV
      mDMmax      = 2.5d4
      mmedmin     = 0.8d-5
      mmedmax     = 1.0d0
      mDMresmin   = 2.0d1     ! this is where a smaller stepsize is used
      mDMresmax   = 3.0d2     ! for a better sampling of the resonance
      mmedresmin  = 2.5d-3    
c use the following instead for p-wave case
      if (scantype.eq.2) then 
        mDMresmin   = 7.0d0
        mDMresmax   = 2.1d2
        mmedresmin  = 8.0d-4    
      endif
      mmedresmin  = 2.5d-3    
      stepmDM     = 0.14d0    ! in log_10
      stepmmed    = 0.11d0
      stepmDMres  = 0.008d0   ! optimal plot performance for 0.005 -> ~2hrs runtime
      stepmmedres = 0.07d0    ! zoom-in may require slightly better
      gin = 0.01   ! start value at mmin

      write (*,*) '-------------------------------------------------------'
      write(*,*) 
      write(*,60) mDMmin, mDMmax
      write(*,65) mmedmin, mmedmax
 60   format('Now scanning over all DM masses from ',F5.1,' GeV to ',F7.1,' GeV,')
 65   format('and over all mediator masses from ',F8.4,' GeV to ',F5.1,' GeV.')
      write(*,*) 'Calculating (alpha and) sigma_T for correct relic density,'
     &           //' please wait...'


      open (unit=10,file=outfile(1))
      write(10,*) 'mmed [GeV]  | mdm [GeV] | alpha | sigma_T_30 [cm^2/g] | '//
     &            'Mcut [Msun] '
      write(10,*) '-----------+------------+-------+---------------------+-'//
     &            '------------'
      write(10,*)

                        
      mDMin = mDMmin
      npoints = 0
c... here we start the main loop over all points in the (mDM,mmed) plane
 100  continue ! loop over DM masses
      npoints = npoints + 1 ! this is just to control output
      mMEDin = mmedmin 
 110  continue ! loop over mediator masses
      gstep = 2.d0 ! inital stepsize for coupling

c... init model and calculate RD
 200    if (mDMin.le.mmedin) goto 320 ! skip unphysical region
        if (scantype.eq.1) call dsgivemodel_vdSIDM_vector(mDMin, mMEDin, gin) 
        if (scantype.eq.2) call dsgivemodel_vdSIDM_scalar(mDMin, mMEDin, gin) 
        call dsmodelsetup(ierr,iwarn)
        if (ierr.ne.0.or.iwarn.ne.0) then
           if (iwarn.eq.2) then
              write(*,*) 'WARNING: model initialized with very large coupling'
              write(*,*) ' alpha = ', gin**2/4./pi
           else   
             write (*, *) 'mDM,mmed [GeV], g = ', mDMin, mMEDin, gin  
             write (*,*) 'Errorflags ierr, iwarn = ', ierr, iwarn
             write (*,*) 'Aborting program...'
             goto 500
           endif
        endif
        oh2=dsrdomega(0,1,xf,ierr,iwarn,nfc)

c... change stepsize adaptively if accuracy goal is not yet met
        if (abs(oh2-oh2goal).gt.doh2goal) then
          if (oh2.gt.oh2goal) then
            if (gstep.lt.1.d0) gstep=1./gstep**0.65
          else
            if (gstep.gt.1.d0) gstep=1./gstep**0.65
          endif
          if (abs(1.d0-gstep).lt.1.d-5.or.gin.gt.1d1) then
            write(*,*) 'Did not manage to find coupling for mdm, mmed =', 
     &                  mDMin, mMEDin
            write(*,*) '[continue with next mass]'
c            read(*,*)
c            gin=0.1
            goto 300
          endif
          gin=gin*gstep
          goto 200 ! re-calculate RD with new value of gin
        endif

 300    continue

        alpha = gin**2/4./pi
        sigT = dssisigtmav(vkms)
        tkd=dskdtkd(1) ! 1=full calculation
        mcut=dskdmcut(mDMin,tkd,1)
        write(10,*) mMEDin, mDMin, alpha, sigT, mcut

 320    if (mmedin.gt.mmedmax) then
          if (npoints.ge.10) then
            write(*,*) 'DM mass ',mDMin, 'completed.'
            npoints = 0
          endif
          if (mDMin.lt.mDMresmin.or.mDMin.gt.mDMresmax.or.mMEDin.lt.mmedresmin) then
            mDMin  = mDMin*10**stepmDM
          else ! use finer stepsizes for sampling in the resonacne region
            mDMin  = mDMin*10**stepmDMres
          endif
          if (mDMin.gt.mDMmax) goto 400 ! end scan
          goto 100 ! next DM mass
 
        else  ! continue loop over mediator masses, for fixed DM mass      
        
          if (mDMin.lt.mDMresmin.or.mDMin.gt.mDMresmax.or.mMEDin.lt.mmedresmin) then
            mMEDin = mMEDin*10**stepmmed
          else ! use finer stepsizes for sampling in the resonacne region

            nres      = sqrt(6*mDMin*alpha/mMEDin)/pi 
            deltares  = abs(1.-pi**2*nres**2*mMEDin/6./mDMin/alpha)
            deltares2 = abs(1.-pi**2*(nres+1)**2*mMEDin/6./mDMin/alpha)

            if (nres.lt.10.and.deltares.lt.0.2.or.deltares2.lt.0.2.and.
     &          sigT.gt.0.01.and.sigT.lt.1000.0) then
               mMEDin = mMEDin*10**(stepmmedres/10.)        
            else 
               mMEDin = mMEDin*10**stepmmedres        
            endif
          endif
          goto 110
        endif

400   write(*,*)
      write(*,*) 'Done: results written to ', 
     &             outfile(1)(1:index(outfile(1),' ')-1), '.'

500   close(10)
      call CPU_TIME(tfinish)

      write (*,*)
      write (*,*) '-------------------------------------------------------'
      write (*,*) 'This DarkSUSY program has finished successfully.'
      write (*,*) 'Total time needed (in minutes): ',(tfinish-tstart)/60.
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-------------------------------------------------------'
      write (*,*)
      stop
 999  end

      
