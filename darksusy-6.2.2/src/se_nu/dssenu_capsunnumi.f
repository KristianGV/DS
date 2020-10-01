c  dssenu_capsunnumi.f:
***********************************************************************
*** dssenu_capsunnumi gives the capture rate of WIMPs in the sun
*** given a specified velocity distribution. the integrations over
*** the sun's radius and over the velocity distribution are
*** performed numerically
*** input: mx [ gev ]
***        sigsi [ cm^2 ]
***        sgisd [ cm^2 ]
***        foveru [ cm^-3 (cm s^-1)^-2 ] external function f(u)/u with
***          velocity u [ km s^-1 ] as argument.
*** Author: Joakim Edsjo
*** Date: 2003-11-26
***********************************************************************

      real*8 function dssenu_capsunnumi(mx,sigsi,sigsd,foveru)
      implicit none

      include 'dssem_sun.h'
      include 'dssecom.h'

      real*8 mx,sigsi,sigsd,ma
      real*8 mxx,max,sigsix,rmin,rmax,res,dssenu_csint,rx,foveru,
     &  sigsdx
      integer i,vtype
      external dssenu_csint,foveru
      common/seint/mxx,max,sigsix,sigsdx,rx,vtype

c...Check if we need to load table
      call dssem_sunread

c...Make sure we have an updated table of elements to include
      call dssenu_selectelements

c...perform integration as given in Gould, ApJ 521 (1987) 571.
      mxx=mx
      sigsix=sigsi
      sigsdx=sigsd
      dssenu_capsunnumi=0.0d0

c...sum over elements

c---Spin-independent
      sctype=1 ! SI
      if (sigsi.gt.0.d0) then
        do i=1,sdelsi
          zx=sdzsi(i)
          ix=sdisosi(i)
          ma=sdma(zx,ix) ! nucleus mass
          max=ma
c...begin radial integration
          rmin=0.d0
          rmax=r_sun*100.0d0  ! sun radius in centimeters
          call dssenu_hiprecint(dssenu_csint,foveru,rmin,rmax,res)
          dssenu_capsunnumi=dssenu_capsunnumi+res
c          write(*,*) 'Element i=',i,'  C_SI =',res,'  sum=',dssenu_capsunnumi 
        enddo
      endif

c---Spin-dependent
      sctype=2 ! SD
      if (sigsd.gt.0.d0) then
        do i=1,sdelsd
          zx=sdzsd(i)
          ix=sdisosd(i)
          ma=sdma(zx,ix) ! nucleus mass
          max=ma
c...begin radial integration
          rmin=0.d0
          rmax=r_sun*100.0d0  ! sun radius in centimeters
          call dssenu_hiprecint(dssenu_csint,foveru,rmin,rmax,res)
          dssenu_capsunnumi=dssenu_capsunnumi+res
c          write(*,*) 'Element i=',i,'  C_SD =',res,'  sum=',dssenu_capsunnumi 
        enddo
      endif


      return
      end
