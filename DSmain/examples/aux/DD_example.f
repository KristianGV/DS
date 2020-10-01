      program DD_example
c
c     This program is an example DarkSUSY main program to demonstrate 
c     usage of the generic WIMP module for the calculation of direct 
c     detection observables.
c
c     Author: Torsten Bringmann, 2018-11-24
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none      
      real*8 tstart,tfinish                  ! aux
      integer ierr, info, opt
      character charopt
      integer genselfconj,genpdg             ! genWIMP parameters    
      real*8 genmwimp,gensvann,genSI         ! genWIMP parameters   
      integer a(5), z(5), stoich(5)          ! parameters for nuclei 
      real*8 spin
      real*8 sigij(27,27)                    ! cross sections
      real*8 v, e                            ! WIMP-nucleus velocity in km/s
                                             ! recoil energy in keV
      real*8 rsi,rsd                         ! rates

ccc functions
      real*8 dsdmspin

c     Here we include the file dsidtag.h which contains information
c     about which particle module is used.
      include 'dsidtag.h'

      call CPU_TIME(tstart)
      call dsinit
      write (*,*) 
      
c... Let's first read in mass and DD cross sections
 50   write(*,*) 'Enter WIMP mass (GeV): '
      read(5,*) genmwimp
      write(*,*) 'Enter WIMP-nucleon scattering cross section (pb)'
      read(5,*) genSI


c... There are more parameters needed to initialize a generic WIMP,
c... but we will not need them for the calculation of direct detection observables
c... We thus simply set them to some fiducial values
      genselfconj=1
      gensvann=1.0d-26
      genpdg=5
      
c... Now transfer these parameters to the DarkSUSY common blocks      
      call dsgivemodel_generic_wimp(genmwimp,genselfconj,
     &     gensvann,genpdg,genSI)
c... As usual, we have to finalize the initialization of the model
c... by a call to dsmodelsetup
      call dsmodelsetup(ierr,info)

c... As an example, we set up scattering on NaI
      a(1) = 23                 ! Na-23
      z(1) = 11                 ! Na-23
      a(2) = 127                ! I-127
      z(2) = 53                 ! I-127
      stoich(1)=1               ! rel. number of Na nuclei
      stoich(2)=1               ! rel. number of I nuclei
      v = 200.d0                ! WIMP-nucleus velocity in km/s
      e = 10.0                  ! recoil energy in keV

      write(*,*)
100   write (*,*) '-------------------------------------------------------'
      write(*,*) 'Elastic scattering cross sections for'
      write(*,*) '* WIMP-nucleus velocity [km/s] : ',v
      write(*,*) '* Recoil energy [kev] : ',e
      write (*,*) '-------------------------------------------------------'
      call dsddsigma(0.d0,0.d0,a(1),z(1),sigij,ierr)
      write(*,*) 'sigma_SI for Na-23 (pb) = ',1.d36*sigij(1,1)
      write(*,*) 'sigma_SD for Na-23 (pb) = ',1.d36*sigij(4,4)
      call dsddsigma(0.d0,0.d0,a(2),z(2),sigij,ierr)
      write(*,*) 'sigma_SI for I-127 (pb) = ',1.d36*sigij(1,1)
      write(*,*) 'sigma_SD for I-127 (pb) = ',1.d36*sigij(4,4)
      write(*,*)

      write (*,*) '-------------------------------------------------------'
      write(*,*) ' Count-rate for NaI [counts/kg-day-keV] '
      write (*,*) '-------------------------------------------------------'
      call dsdddrde(651.3d0,e,2,a,z,stoich,rsi,rsd,1)
      write(*,*) 'spin-independent = ',rsi
      write(*,*) 'spin-dependent = ',rsd
      write(*,*)      
      write(*,*)      
      
200   write(*,*) 'What do you want to do next ?'
      write(*,*) '  1 = change reference velocity and recoil energy'
      write(*,*) '  2 = change dark matter particle spin'
      write(*,*) '  3 = specify type of cross section entered initially'
      write(*,*) '  [else] = quit'
      read(*,*) opt
      write(*,*)
    
      if (opt.eq.1) then
        write(*,*) 'Enter new reference velocity [km/s] : '
        read(*,*) v
        write(*,*) 'Enter new recoil energy [keV] : '
        read(*,*) e
        goto 100 ! This does note require a new call to dsmodelsetup
      elseif (opt.eq.2) then
        write(*,*) 'Currently the DM spin is set to ', dsdmspin()
        write(*,*) 'What do you want to change this to?' 
        read(*,*) spin
        call dsgivemodel_generic_wimp_opt('spin',spin)
      elseif (opt.eq.3) then
        write(*,*) 'The DM-nucleon cross section is set to sigma = ',genSI
        write(*,*) 'Which of the follwoing options to you prefer ?'
        write(*,*) 'a) sigma = sip = sin [default in generic_wimp]' 
        write(*,*) 'b) sigma = sdp = sdn' 
        write(*,*) 'c) sigma = sip' 
        write(*,*) 'd) sigma = sin' 
        write(*,*) 'e) sigma = sin' 
        write(*,*) 'f) sigma = sin' 
        write(*,*) '[unspecified couplings are set to zero]'
        read(*,*) charopt

c...we first set all cross sections to zero
        call dsgivemodel_generic_wimp_opt('sigsip',0.0d0)
        call dsgivemodel_generic_wimp_opt('sigsin',0.0d0)
        call dsgivemodel_generic_wimp_opt('sigsdp',0.0d0)
        call dsgivemodel_generic_wimp_opt('sigsdn',0.0d0)
c...and then only the chosen cross sections to the specified input value
        if (charopt.eq.'a') then
          call dsgivemodel_generic_wimp_opt('sigsip',genSI)
          call dsgivemodel_generic_wimp_opt('sigsin',genSI)
        elseif (charopt.eq.'b') then
          call dsgivemodel_generic_wimp_opt('sigsdp',genSI)
          call dsgivemodel_generic_wimp_opt('sigsdn',genSI)
        elseif (charopt.eq.'c') then
          call dsgivemodel_generic_wimp_opt('sigsip',genSI)
        elseif (charopt.eq.'d') then
          call dsgivemodel_generic_wimp_opt('sigsin',genSI)
        elseif (charopt.eq.'e') then
          call dsgivemodel_generic_wimp_opt('sigsdp',genSI)
        elseif (charopt.eq.'f') then
          call dsgivemodel_generic_wimp_opt('sigsdn',genSI)
        endif
      else
        goto 500
      endif
   
      call dsmodelsetup(ierr,info)
      if (info.eq.3) write(*,*) 'spin = ',dsdmspin(), 'currently not supported!'
      if (info.eq.10) then
         write(*,*) 'WARNING: Will set spin-dependent couplings to zero'
         write(*,*) 'because the dark matter particles have no spin!'
         write(*,*) '[continue]'
         read(*,*)
      endif
      write(*,*)
      if (ierr.ne.0) goto 200
      goto 100


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

      
