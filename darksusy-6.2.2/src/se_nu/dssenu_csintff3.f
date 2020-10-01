c...
c...auxiliary function for inner integrand
c...This function here is essentially the integrand in Gould (A5),
c...ApJ 321 (1987) 571, but with sigma added and a more general form factor
c...than the exponential one in Gould A5.
c...input: y=Delta-E / (mx*w**2/2)
c...output: integrand in units of cm^2

      real*8 function dssenu_csintff3(y)
      implicit none
      include 'dssem_sun.h'
      include 'dssecom.h'
      include 'dsmpconst.h'

      real*8 y,er
      real*8 mu,mx,muplus,muminus,ma
      real*8 fac1,fac2
      real*8 sigsi,sigsd,r
      complex*16 gg(27,2)
      real*8 sigij(27,27)
      real*8 spinx,dsdmspin
      integer vtype,ierr
      external foveru,dsintcsintff3
      common/seint/mx,ma,sigsi,sigsd,r,vtype

c...Set up g2p calculation      
      gg(1,1)=dcmplx(gpsx,0.d0)
      gg(1,2)=dcmplx(gnsx,0.d0)
      gg(4,1)=dcmplx(gpax,0.d0)
      gg(4,2)=dcmplx(gnax,0.d0)
      spinx=dsdmspin()

c...Setup      
      mu=mx/ma
      muplus=(mu+1.d0)/2.d0
      muminus=(mu-1.d0)/2.d0

      Er=mx*wwx**2*y/(c_light**2*2.d0)*1.d6 ! recoil energy Er=deltae in keV
c      q=sqrt(ma*mx*wwx**2*y/c_light**2)  ! q=sqrt(2*ma*deltae), deltae=mx*ww**2*y/2
      fac1=muplus**2/mu

      if (sctype.eq.1) then     ! SI
         call dsddg2sigma(mx,spinx,wwx,Er,ax,zx,gg,sigij,ierr)
         fac2=sigij(1,1)
      elseif (sctype.eq.2) then ! SD
         call dsddg2sigma(mx,spinx,wwx,Er,ax,zx,gg,sigij,ierr)
         fac2=sigij(4,4)
      else
         write(*,*) 'DS Error in dssenu_csintff3: Invalid sctype: ',sctype
         write(*,*) 'DarkSUSY stopping.'
         stop
      endif

      dssenu_csintff3=fac1*fac2
      return
      end
