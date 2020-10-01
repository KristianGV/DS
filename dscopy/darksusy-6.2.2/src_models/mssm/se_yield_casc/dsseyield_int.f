      real*8 function dsseyield_int(f,a,b,
     &     chi,wwh,m0,m1,m2,
     & e0,eth,thm,fk,fv,seerror)      
c_______________________________________________________________________
c  integrate function f between a and b
c  Note, we give all arguments here instead of having them via common
c  blocks to allow for Open MP parallellization      
c  input
c    integration limits a and b
c  called by dsseyieldfth
c  author: joakim edsjo (edsjo@fysik.su.se) 96-05-16
c  based on paolo gondolos wxint.f routine.
c=======================================================================
      implicit none
c      include 'dsseyieldcom.h'
      include 'dsio.h'
      real*8 f,a,b,tot,eps,st,os,ost,del,sum,x
      real*8 m0,m1,m2,e0,eth,thm
      integer chi,wwh,fk,fv   ! wwh=1 - sun, 2 - earth
      integer jmax,it,l,j,nfcn,jdid,seerror
      external f
c      parameter (a=-1.0,b=1.0,eps=1.0d-4,jmax=20)
      parameter (eps=1.0d-2,jmax=30)  ! je change in eps ps change in jmax
      dsseyield_int=0.d0
      del=b-a
      ost=0.5*del*(f(a,m0,m1,m2,e0,eth,thm,chi,
     &     wwh,fk,fv)+
     &  f(b,m0,m1,m2,e0,eth,thm,chi,
     &     wwh,fk,fv))
      x=0.5*(b+a)
      st=0.5*(ost+del*f(x,m0,m1,m2,e0,eth,thm,chi,
     &   wwh,fk,fv))
      os=(4.0*st-ost)/3.0
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=a+0.5*del
        sum=0.0
        do l=1,it
           sum=sum+f(x,m0,m1,m2,e0,eth,thm,chi,
     &        wwh,fk,fv)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (abs(tot-os).le.eps*abs(os)) then
           dsseyield_int=tot
           return
        endif
     	os=tot
        ost=st
c        type *,'jdid',jdid,' os',os, 'ost',ost
      enddo

!$omp critical (stdout)
      if (prtlevel.ge.2) then
         write(*,*) 'DS WARNING: too many steps in dsseyield_int.'
      endif
!$omp end critical (stdout)

      seerror=1
      dsseyield_int=0.0d0

      end






