**********************************************************************
*** Program to show how capture rates can be calculated.
*** This program is a bit more advanced than caprates.f as it goes into
*** interal workings of the routines to select exactly which elements to
*** include. It is provided to give an example on how to use routines
*** in more advanced way to get detailed information.      
***
*** This program is intended to be used with the generic_wimp module
*** The program scans a few masses and calculates the total capture rate, and
*** capture rates on individual elements.
*** Compared to the default calculation in dsmain_wimp and dstest, this
*** program does not use tabulated capture rates, instead it always uses
*** the full numerical calculation (with numerical integration over the
*** velocity, radius and momentum transfer and sum over all the chosen
*** elements) to be able to choose different options
*** and choice of elements to include. Hence, the program takes quite a long
*** time to run (typically several hours in the default setup).
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: February, 2018     
**********************************************************************

      program caprates_ff

      implicit none

      real*8 mwimp,rho
      real*8 gps,gns,gpa,gna
      integer ng,scheme
      parameter (ng=4)
      real*8 gg(ng)
      integer pdgann
      real*8 svann,sigsi,sigsd
      logical selfconj
      real*8 dssenu_capsunnumff
      integer i,n,j
      integer zsitot
      parameter (zsitot=9)
      integer zsi(zsitot)
      data zsi/1,2,6,7,8,10,12,26,28/
      real*8 capsunsi(0:zsitot)
      integer zsdtot
      parameter (zsdtot=9)
      integer zsd(zsdtot)
      data zsd/1,2,6,7,8,11,12,13,14/
      real*8 capsunsd(0:zsdtot)
      

      real*8 rateea,ratesu
      integer istat

c...Initialize DarkSUSY with default settings      
      call dsinit
      call dssenu_set('nummed')  ! used with higher level neutrino rate routines
      call dssem_sunset('default') ! default solar model
c...Choose form factors      
      call dsddset('sf_m','best')
      call dsddset('sf_sigma','best') 

c...Set up astrophysics (vel. distribution is currently the default)      
      rho=0.3d0              ! local halo density in GeV/cm^3
      
c...Set up particle physics model
      selfconj=.true.           ! self-conjugated DM (it own anti-particle)
      svann=3d-26               ! ann cross sect, cm^3/s
      pdgann=5                  ! annihilation to b b-bar (pdg code 5)
      sigsi=1.0d-7              ! SI scattering cross section, pb
      sigsd=1.0d-5              ! SD scattering cross section, pb

      open(unit=47,file='caprates_ff_si.dat',status='unknown',
     &     form='formatted')
      write(47,100) sigsi
      write(47,101) (zsi(j),j=1,zsitot)
 100  format('# For scattering cross section ',E14.7,' pb')
 101  format('#...m_WIMP... ....C_tot... C_i for Z_i =',100(1x,I2))

      open(unit=48,file='caprates_ff_sd.dat',status='unknown',
     &     form='formatted')
      write(48,100) sigsd
      write(48,101) (zsd(j),j=1,zsdtot)

      call dssem_sunread        ! Just needed as we call dssenu_selectedelements directly below
      n=10  ! Note, increase this if you want more points in mass
      do i=0,n
         mwimp=10**(dble(i)/dble(n)*4)
         write(*,*) 'WIMP mass: ',mwimp
      
c...Set up generic wimp model with these parameters
         call dsgivemodel_generic_wimp(mwimp,selfconj,svann,pdgann,sigsi)
      
c...Get couplings from particle physics module
c...Alternatively they can be set manually      
c        call dsddgpgn(ng,gg,scheme)
c        if (scheme.eq.1) then     ! for four-fermion couplings
c           gps=gg(1)
c           gns=gg(2)
c           gpa=gg(3)
c           gna=gg(4)
c        else
c           write(*,*) 'ERROR: scheme = ',scheme,
c     &       ' not yet implemented.'
c        endif
        call sig2g(mwimp,sigsi,sigsi,sigsd,sigsd,gg)
        gps=gg(1)
        gns=gg(2)
        gpa=gg(3)
        gna=gg(4)
c...We force a setup with new elements by calling dssenu_set twice
c...Normally this would not be needed, but here we want to work
c...with ther internal settings of the routines to be able to
c...select only one element        
        call dssenu_set('numlo')
        call dssenu_selectelements
        call dssenu_set('nummed')
        call dssenu_selectelements
        capsunsi(0)=dssenu_capsunnumff(mwimp,rho,gps,gns,0.d0,0.d0)
        write(*,*) 'The capture rate is: ',capsunsi(0),
     &       's^-1.'
        do j=1,zsitot
           call dssenu_oneelement(zsi(j),0)
           capsunsi(j)=dssenu_capsunnumff(mwimp,rho,gps,gns,0.d0,0.d0)
        enddo
        write(47,110) mwimp,(capsunsi(j),j=0,zsitot)

c...We force a setup with new elements by calling dssenu_set twice
        call dssenu_set('numlo')
        call dssenu_selectelements
        call dssenu_set('nummed')
        call dssenu_selectelements
        capsunsd(0)=dssenu_capsunnumff(mwimp,rho,0.d0,0.d0,gpa,gna)
        write(*,*) 'The capture rate is: ',capsunsd(0),
     &       's^-1.'
        write(*,*) mwimp,rho,gpa,gna
        do j=1,zsdtot
           call dssenu_oneelement(0,zsd(j))
           capsunsd(j)=dssenu_capsunnumff(mwimp,rho,0.d0,0.d0,gpa,gna)
        enddo
        write(48,110) mwimp,(capsunsd(j),j=0,zsdtot)
        
      enddo
      close(47)
      close(48)
      
 110  format(100(1x,E12.6))
      
      end


**********************************************************************

      subroutine dssenu_oneelement(zsi,zsd)
*** Select just one element for SI and SD for capture rate calculations

      implicit none
      include 'dssem_sun.h'
      include 'dsio.h'
      integer zsi,zsd,l,m,m_1,m_2
      logical inciso

c---spin-independent elements      
      sdelsi=0

      if (zsi.eq.0) goto 200
      do l=zsi,zsi ! only one SI element
         
c...    determine if isotopes should be included
        inciso=.false.
        if (sdmfr(l,1,int(sdn/2)).gt.0.d0) inciso=.true. ! iso information=yes
        if (inciso) then
           m_1=1
           m_2=isomax
        else
           m_1=0
           m_2=0
        endif

        do m=m_1,m_2
           if (sdmfr(l,m,int(sdn/2)).gt.0.d0) then
c              write(*,*) 'AA:  adding isotope ',m
              sdelsi=sdelsi+1
              sdzsi(sdelsi)=l
              sdisosi(sdelsi)=m
              if (prtlevel.ge.2) 
     &          write(*,*) 'Adding SI Z=',l,' i=',m,' A=',sdaa(l,m)
c              write(*,*) 'AA: Adding (Z,I) = ',l,m,' ',
c     &           sdname(l),'-',int(sdaa(l,m)+0.5),
c     &           10**(sdabund(l)-12)*sdaa(l,m)**4
              if (sdelsi.ge.sdelmax) then
                 write(*,*) 'DS ERROR in dssnu_selectelements: ',
     &             'sdelmax too small for SI'
                 stop
              endif
            endif
        enddo
      enddo
      write(*,'(A,I3,A)') 'dssenu_oneelement: ',sdelsi,
     &     ' elements will be included in the Sun',
     &     ' SI capture calculation.'
        
 200  continue
      
c---Spin-dependent elements
      
      sdelsd=0
      if (zsd.eq.0) goto 300
      do l=zsd,zsd ! only one element
         do m=1,isomax
            if (sdsp(l,m).ne.0.d0) then
               if (sdmfr(l,m,int(sdn/2)).gt.0.d0) then
                  sdelsd=sdelsd+1
                  sdzsd(sdelsd)=l
                  sdisosd(sdelsd)=m
                  if (prtlevel.ge.2) 
     &              write(*,*) 'Adding SD Z=',l,' i=',m,' A=',sdaa(l,m)
               endif
               if (sdelsd.ge.sdelmax) then
                  write(*,*) 'DS ERROR in dssenu_selectelements: ',
     &              'sdelmax too small for SD'
                  stop
               endif
            endif
         enddo
      enddo

      write(*,'(A,I3,A)') 'dssenu_selectelements: ',sdelsd,
     &     ' elements will be included in the Sun',
     &     ' SD capture calculation.'

 300  continue
      return
      end

**********

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
      
