**********************************************************************
*** Subroutine dsinit intializes DarkSUSY. Every DarkSUSY main program
*** should call this routine first thing to get DarkSUSY set up and
*** ready to run
***
***   type : commonly used
***   desc : Initialize DarkSUSY (should always be called)
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: June 11, 2013 
*** (mod Torsten Bringmann, 2015)
**********************************************************************
      subroutine dsinit
      implicit none
      include 'dsio.h'
      include 'dsidtag.h'
      include 'dspbcom.h'
      include 'dssem_sun.h'
      include 'dsmpconst.h'
      include 'dsver.h'

c...Startup
      dsversion=dsver ! move from parameter to variable to allow changes

      write(*,*) 
      write(*,*) 
     &  '*********************************************************'
      write(*,*) 
     &  '*** Welcome to DarkSUSY version'
      write(*,*) '*** ',dsversion
      write(*,*) 
     &  '*********************************************************'
      write(*,*) ' '
      
c... This is reset in dsinit_module. We only include it here to guarantee
c... that the electron mass is always initialized (as it is needed in the 
c... propagation routines 
      m_e = 0.000510999907d0 ! in GeV

c...Initialize specific particle physics module
      call dsinit_module

c...Initialize remaining, particle-independent parts
      call dsrdinit ! initialize RD routines
      call dshmset('default')
      call dskdset('default')
      call dssenu_set('default')
      call dssem_sunset('default')
      call dsanyield_init
      call dsanyield_dbset(0,0.d0) ! default MC and p0bar
      call dsddinit
      call dsdmdinit ! default halo models

c...initialize Sun potential
      sdread=.false. 
      call dssem_sunread 
      
c...Initialization of axisymmetric cosmic-ray routines
      call dscraxiinit
c... + plus a default set of propagation parameters      
      call dscraxiset_default
      
c... program switches
      prtlevel = 0
      idtag = '[no ID tag set]'
      luout = 6  ! unit where messages go

c...load table of nuclides
      call dsreadnuclides

      write(*,*) 'Initialization of DarkSUSY complete.'      
      write(*,*)

      return
      end



