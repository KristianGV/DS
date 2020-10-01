***********************************************************************
*** dssenu_capsunnumffi gives the capture rate of WIMPs in the sun
*** given a specified velocity distribution. the integrations over
*** the sun's radius and over the velocity distribution are
*** performed numerically.
*** Compared to dssenu_capsunnumi, the form factors are not assumed to be
*** exponential, and an actual integration over the form factors is
*** performed
*** input: mx [ gev ]
***        gps,gns,gpa,gna - WIMP-nucleon four-fermion couplings
***          (obtained from e.g. dsddgpgn), [GeV^-4]
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
*** Author: Joakim Edsjo
*** Date: 2003-11-26
*** Modified: 2010-05-19 to allow for more elements and arbitrary form factors
***********************************************************************

      real*8 function dssenu_capsunnumffi(mx,gps,gns,gpa,gna,foveru)
      implicit none

      include 'dssem_sun.h'
      include 'dssecap.h'
      include 'dssecom.h'
      include 'dsio.h' ! JE TMP

      real*8 mx,gps,gns,gpa,gna,ma
      real*8 mxx,max,sigsix,rmin,rmax,res,dssenu_csintff,rx,foveru,
     &  sigsdx
      real*8 l2jjpp,l2jjnn,l2jjpn
      integer i,vtype,l,ierr
      external dssenu_csintff,foveru
      common/seint/mxx,max,sigsix,sigsdx,rx,vtype

c...Check if we need to load table
      call dssem_sunread

c...Make sure we have an updated table of elements to include
      call dssenu_selectelements

c...Initialize capture arrays
      do i=1,zmax
         do l=0,isomax
            capsunsi(i,l)=0.d0
            capsunsd(i,l)=0.d0
         enddo
      enddo


c...perform integration as given in Gould, ApJ 521 (1987) 571,
c...but also integrating over q (with form factor) numerically
      mxx=mx
      gpsx=gps
      gnsx=gns
      gpax=gpa
      gnax=gna
      dssenu_capsunnumffi=0.0d0

c...sum over elements
c      write(*,*) 'FFI1: ',mx,gpsx,gnsx,gpax,gnax,sdelsi,sdelsd ! JE TMP
c---Spin-independent
      sctype=1 ! SI
      if (abs(gpsx).gt.0.d0.or.abs(gnsx).gt.0.d0) then
        do i=1,sdelsi
          zx=sdzsi(i)
          ix=sdisosi(i)
          ax=int(sdaa(zx,ix)+0.5d0)
c          write(*,*) 'FFI2: SI: ',i,zx,ix,ax ! JE TMP
          ma=sdma(zx,ix) ! nucleus mass
          max=ma
c...begin radial integration
          rmin=0.d0
          rmax=r_sun*100.0d0  ! sun radius in centimeters
          call dssenu_hiprecint(dssenu_csintff,foveru,rmin,rmax,res)
          capsunsi(zx,ix)=res
          if (ix.gt.0) then ! add to 0 component so that it contains isotope sum
             capsunsi(zx,0)=capsunsi(zx,0)+res
          endif
          dssenu_capsunnumffi=dssenu_capsunnumffi+res
c          write(*,*) 'FFI3 Element i=',i,'  C_SI =',res,'  sum=',dssenu_capsunnumffi
        enddo
      endif

c---Spin-dependent
      sctype=2 ! SD
      if (abs(gpax).gt.0.d0.or.abs(gnax).gt.0.d0) then
        do i=1,sdelsd 
          zx=sdzsd(i)
          ix=sdisosd(i)
          ax=int(sdaa(zx,ix)+0.5d0)
c          write(*,*) 'FFI3: SD: ',i,zx,ix,ax ! JE TMP
          ma=sdma(zx,ix) ! nucleus mass
          max=ma
c...Check that spin-dep form factors exist
c          prtlevel=1 ! JE TMP
          call dsddffsd(1.d-4,ax,zx,l2jjpp,l2jjnn,l2jjpn,ierr)
c          prtlevel=0 ! JE TMP
c          write(*,*) 'XXX SD: Z=',zx,' A=',ax,'ls: ',
c     &        l2jjpp,l2jjnn,l2jjpn ! JE TMP
c... FIXME: More elaborate testing if SD form factor exists
          if (l2jjpp.eq.0.d0.and.l2jjnn.eq.0.d0.and.l2jjpn.eq.0.d0) then
c             write(*,198) 'WARNING in dssenu_capsunnumffi: ',
c     &       'No SD form factor for Z=',zx,' A=',ax ! JE TMP
             res=0.d0
             capsunsd(zx,ix)=-1.d0
          else
c             write(*,199) 'SD form factor found for Z=',zx,' A=',ax ! JE TMP
c...begin radial integration
             rmin=0.d0
             rmax=r_sun*100.0d0  ! sun radius in centimeters
c             prtlevel=1 ! JE TMP
             call dssenu_hiprecint(dssenu_csintff,foveru,rmin,rmax,res)
c             prtlevel=0 ! JE TMP
             capsunsd(zx,ix)=res
c...    add to 0 component so that it contains isotope sum
             if (ix.gt.0) then 
                capsunsd(zx,0)=capsunsd(zx,0)+res
             endif
c             if (res.gt.dssenu_capsunnumffi*0.5d0) then ! JE TMP
c                write(*,*) 'New large SD contribution:'
c                write(*,*) 'SD: z=',zx,' i=',ix,' a=',ax,' ma=',max,
c     &       ' res=',res ! JE TMP
c             endif
             dssenu_capsunnumffi=dssenu_capsunnumffi+res
          endif
c          write(*,*) 'SD: z=',zx,' i=',ix,' a=',ax,' ma=',max,
c     &       ' res=',res ! JE TMP
c          write(*,*) 'Element i=',i,'  C_SD =',res,'  sum=',dssenu_capsunnumi 
        enddo
      endif

c 198   format(A,A,I3,A,I3)
c 199   format(A,I3,A,I3)

      return
      end
