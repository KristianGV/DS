      real*8 function dspbbesselzeroj0(s)
************************************************************************
*** s-th zero of J0
************************************************************************
      implicit none
      integer s
ccc
      real*8 store(100000,2)
      common/pbintj0store/store
      real*8 Rstore,pbrhstore
      integer supjstore,supfstore
      common/pbintj0storeck/Rstore,pbrhstore,supjstore,supfstore
ccc
 100  continue
      if(s.gt.supjstore) then
        call dspbset_besselj0store(1,-100.d0)
        goto 100
      endif
      dspbbesselzeroj0=store(s,1)
      return
      end 


      real*8 function dspbbesselfactorj0(s,R)
************************************************************************
*** with nus the s-th zero of J0, this is J0(nus*R/diffRh)/J1(nus)**2
************************************************************************
      implicit none
      include 'dscraxicom.h'
      integer s
      real*8 R
ccc
      real*8 store(100000,2)
      common/pbintj0store/store
      real*8 Rstore,pbrhstore
      integer supjstore,supfstore
      common/pbintj0storeck/Rstore,pbrhstore,supjstore,supfstore
ccc
      if(dabs(R-Rstore).gt.1.d-5*dabs(R+Rstore)) supfstore=0
      if(dabs(pbrhstore-diffRh).gt.1.d-5*dabs(pbrhstore+diffRh)) 
     &   supfstore=0
 100  continue
      if(s.gt.supjstore) then
        call dspbset_besselj0store(1,R)
        goto 100
      endif
 200  continue
      if(s.gt.supfstore) then
        call dspbset_besselj0store(0,R)
        goto 200
      endif
      dspbbesselfactorj0=store(s,2)
      return
      end 



      subroutine dspbset_besselj0store(how,R)
************************************************************************
*** inputs:
***     how: for how=1 add a number=incr of zeros of J0 to store(ii,1) 
***     R: for R>0 add values J0(zero**R/diffRh)/J1(zero)**2 to 
***        store(ii,2)
************************************************************************
      implicit none
      include 'dscraxicom.h'
      integer how,ii,incr,nzeron,sup
      parameter(incr=5000,nzeron=100000)
      real*8 R,jnzero(100000),zero
      real*8 dbesj0,dbesj1
ccc
      real*8 store(100000,2)
      common/pbintj0store/store
      real*8 Rstore,pbrhstore
      integer supjstore,supfstore
      common/pbintj0storeck/Rstore,pbrhstore,supjstore,supfstore
ccc
      if(how.eq.1) then
        sup=supjstore+incr
        if(sup.gt.nzeron) then
          write(*,*) 'DS: in dspbset_besselstore exceeded the dimension'
          write(*,*) 'DS: of the store common block'
          stop
        endif
        call dbzejy(dble(0),sup,1,1.d-16,jnzero)
        do ii=supjstore,sup
          store(ii,1)=jnzero(ii)
        enddo
        supjstore=sup
      endif
      if(R.ge.0.d0) then
        do ii=supfstore,supjstore
          zero=store(ii,1)
          store(ii,2)=dbesj0(zero*R/diffRh)/(dbesj1(zero))**2
        enddo
        supfstore=supjstore
        Rstore=R
        pbrhstore=diffRh
      endif
      return
      end


      subroutine dspbbesselini
************************************************************************
*** intializations for bessel function zeros and bessel transforms
************************************************************************
      implicit none
ccc
      real*8 Rstore,pbrhstore
      integer supjstore,supfstore
      common/pbintj0storeck/Rstore,pbrhstore,supjstore,supfstore
ccc
      supjstore=0
      supfstore=0
      Rstore=-100.d0      
      pbrhstore=-100.d0
      return
      end 



