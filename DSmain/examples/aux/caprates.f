**********************************************************************
*** Program to show how capture rates can be calculated.
*** Compared to caprates_ff.f, this program uses the tabulate capture
*** rates. It scans over a set of masses.      
***
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: February, 2018
**********************************************************************

      program caprates

      implicit none

      real*8 mwimp,rho
      real*8 gps,gns,gpa,gna
      integer ng
      parameter (ng=4)
      real*8 gg(ng)
      integer pdgann
      real*8 svann,sigsi,sigsd
      logical selfconj
      real*8 dssenu_capsuntabff
      integer i,n,j
      real*8 capsunsi
      real*8 capsunsd
      

      real*8 rateea,ratesu
      integer istat

c...Initialize DarkSUSY with default settings      
      call dsinit
      call dssenu_set('tabmed')  ! default capture calculation
      call dssem_sunset('default') ! default solar model

c...Set up astrophysics (vel. distribution is currently the default)      
      rho=0.3d0                 ! local halo density in GeV/cm^3
      
c...Set up particle physics model
      selfconj=.true.           ! self-conjugated DM (it own anti-particle)
      svann=3d-26               ! ann cross sect, cm^3/s
      pdgann=5                  ! annihilation to b b-bar (pdg code 5)
      sigsi=1.0d-7              ! SI scattering cross section, pb
      sigsd=1.0d-5              ! SD scattering cross section, pb

      open(unit=47,file='caprates_si.dat',status='unknown',
     &     form='formatted')
      write(47,100) sigsi
      write(47,101)       
 100  format('# For scattering cross section ',E14.7,' pb')
 101  format('#...m_WIMP... C_tot (s^-1)')

      open(unit=48,file='caprates_sd.dat',status='unknown',
     &     form='formatted')
      write(48,100) sigsd
      write(48,101) 

      n=1000
      do i=0,n
         mwimp=10**(dble(i)/dble(n)*4)
         write(*,*) 'WIMP mass: ',mwimp
      
c...Set up generic wimp model with these parameters
         call dsgivemodel_generic_wimp(mwimp,selfconj,svann,pdgann,sigsi)
      
        call sig2g(mwimp,sigsi,sigsi,sigsd,sigsd,gg)
        gps=gg(1)
        gns=gg(2)
        gpa=gg(3)
        gna=gg(4)
         
        capsunsi=dssenu_capsuntabff(mwimp,rho,gps,gns,0.d0,0.d0)
        write(*,*) 'The capture rate is: ',capsunsi,
     &       's^-1.'
        write(47,110) mwimp,capsunsi

        capsunsd=dssenu_capsuntabff(mwimp,rho,0.d0,0.d0,gpa,gna)
        write(*,*) 'The capture rate is: ',capsunsd,
     &       's^-1.'
        write(48,110) mwimp,capsunsd
        
      enddo
      close(47)
      close(48)
      
 110  format(100(1x,E12.6))
      
      end


**********************************************************************

*** converts cross sections (pb) to couplings assuming sign is positive      
      subroutine sig2g(mwimp,sigsip,sigsin,sigsdp,sigsdn,gg)
      implicit none
      include 'dsmpconst.h'

      real*8 mwimp,sigsip,sigsin,sigsdp,sigsdn,gg(4),gps,gns,gpa,gna
      real*8 fkinp,fkinn

      fkinp = 4./pi*(m_p*mwimp/(m_p+mwimp))**2
      fkinn = 4./pi*(m_n*mwimp/(m_n+mwimp))**2


      gps=2.0d0*sqrt(sigsip*1.d-36/fkinp/gev2cm2)
      gns=2.0d0*sqrt(sigsin*1.d-36/fkinn/gev2cm2)

      gpa=sqrt(sigsdp*1.d-36/gev2cm2/fkinp*4.d0/3.d0)
      gna=sqrt(sigsdn*1.d-36/gev2cm2/fkinn*4.d0/3.d0)
      
      gg(1)=gps
      gg(2)=gns
      gg(3)=gpa
      gg(4)=gna

      return
      end
      
