c_______________________________________________________________________
c  Function dsfi2to2rhs gives the right hand side of the Boltzmann equation
c   s*dY/dt=RHS. Done by integrating dsfi2to2int over s, and multiplying by
c   temperature dependent factors in the end.
c   
c   input:
c     T - temperature of heath bath
c
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function dsfi2to2rhs(lnT)
      implicit none
      include 'dsmpconst.h'
      include 'dsficom.h'
      include 'dsidtag.h'



      real*8 lnT,T,a,b,eps,geff,sqrtgstar,heff,s_ent,H,HPrime,
     &c,d,sum1,sum2,sum3,dsrdthav,
     &dsbessek2,x,debug1,debug2,debug3,sum,sum_tmp,tmp_upper,
     &tmp_upper1,xf,HsT2
      real*8, external :: dsfi2to2int_simp,dsanwx
      integer narray,ncoann,nrs,nthr,i,j
      parameter (narray=1000)

      real*8 mcoann(narray),dof(narray),tm(narray),
     &rm(narray),rw(narray)

      integer npts2,npts_tmp,iwar,nfc
      real*8 points(narray),points_sorted(narray),tmp_lowest,
     &points_int(narray),abserr,resabs,resasc

      T=exp(lnT)
      ! T=lnT
      call dsrdset('dof','default')
      call dskdgeff(T,geff)
      call dsrddof(T,sqrtgstar,heff)

      s_ent=heff*2*pi*pi/45*T*T*T
      H=sqrt(4*pi**3*geff/45)*T*T/mpl
      HPrime=H*heff/sqrt(geff)/sqrtgstar
      HST2=heff*2*pi*pi/45*T*HPrime

      T_22=T;x1_22=m1_22/T;x2_22=m1_22/T
      x=mdm/T
c     Integration limits
      ! a=1/(100*mdm**2);eps=1.d-5
      ! b=1/(4*mdm**2)
      ! b=1
      if(stat.eq.0) then
c     TRIED TO USE <sv>     
      
            if(x.eq.0) then
              dsfi2to2rhs=0
              write(*,*) "Warning x=0 at T=", T
            else
              dsfi2to2rhs=dsrdthav(x,dsanwx)*
     &(dsbessek2(x)/exp(x))**2/HST2

    !           dsfi2to2rhs=mpl*sqrtgstar*sqrt(pi/45)*(45/4/pi**2)**2
    !  &*(dsbessek2(x)/exp(x))**2*dsrdthav(x,dsanwx)/heff**2/T**3*4*pi**4
    !  &/sqrt(8*pi)                    
            end if


! c------------ MY OWN ROUTINE
!         call dsrdparticles(0,narray,
!      &  selfcon,ncoann,mcoann,dof,nrs,rm,rw,nthr,tm)
!         a=(4*mdm**2)+1.d-30
!         b=1.d20;eps=1.d-5
!         sum=0
! c------------ Makes integration limits break points
!         npts2=2+nthr+(2*nrs)
!         npts_tmp=1
!         if(npts2.eq.2) then
!           points(1)=a
!           points(2)=b
!         else
!           points(1)=a
!         do i=1,nthr
!           do j=1,nrs
!             if((points(npts_tmp).lt.rm(j)**2).and.(rm(j).lt.tm(i))) then
!               npts_tmp=npts_tmp+1
!               points(npts_tmp)=(rm(j)-rw(j))**2
!               npts_tmp=npts_tmp+1
!               points(npts_tmp)=(rm(j)+rw(j))**2
!             end if
!           end do
!           npts_tmp=npts_tmp+1
!           points(npts_tmp)=tm(i)**2
!         end do
!         points(npts2)=b
!         end if
! c------------ sorting points from lowest to largest
!         do i=1,npts2
!           ! finding lowest number in points
!           tmp_lowest=1.d100 
!           do j=1,npts2
!             if(i.eq.1) then
!               if(tmp_lowest.gt.points(j)) then
!                 tmp_lowest=points(j)
!               end if
!             else
!               if((tmp_lowest.gt.points(j)).and.
!      &            (points(j).gt.points_sorted(i-1))) then
!                 tmp_lowest=points(j)
!               end if
!             end if
!           end do
!           points_sorted(i)=tmp_lowest
!         end do
! c ------- Want to integrate over x=1/s insted, need new integration limits
!         do i=1,npts2
!           points_int(i)=1/points_sorted(npts2-(i-1))
!         end do

! c-------- Integrate from xmin to xmax with breakpoints in points_sorted
! c         This means to sum over integral from points_int(i) to points_int(i+1)
!         sum=0
!         do i=1,(npts2-1)
!           call dgadap(log(points_int(i)),log(points_int(i+1)),
!      &dsfi2to2int_simp,eps,sum_tmp)
!     !       call dqk21(dsfi2to2int_simp,points_int(i),
!     !  &points_int(i+1),sum,abserr,resabs,resasc)  

!           sum=sum+sum_tmp
!         end do

!         dsfi2to2rhs=sum/HPrime/s_ent*T

c----------------------------------------------------------------------
      else 
        a=1/(4*mdm**2)
        b=1/((4*mdm**2)+1.d40);eps=1.d-4

        call dgadap(b,a,dsfi2to2int_simp,eps,sum1)
        dsfi2to2rhs=(sum1)/HPrime/s_ent*T
      end if

      return
      end
