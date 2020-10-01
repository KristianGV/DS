************************************************************************
*** function dsIB2BRtree returns the branching ratio for the decay rate
*** of possible on-shell states into fermions (at tree level).
***
***   Input: decay - string that describes the decay process (see below)
***             kf - particle code of one of the final state fermions
***                  (integer; see description in dsib2sigmav or dmssm.h)
***
***   Output: partial decay rate at tree level, divided by total width
***           (not necessarily at tree level!) 
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-12-15, update 2105-03-22 (streamlined handling of k2)
*** MG 2015-01-20: added k2 as argument
***                output for ff=k2,k2 or fF=k2,k2+1, resp.
***                sign flip in matrix element for V->ff (see comments)
************************************************************************

      real*8 function dsIB2BRtree(decay,kf)
      implicit none
      include 'dsmssm.h'


c------------------------ variables ------------------------------------

      character*(*) decay
      real*8  result, msq, xim, xip
      integer kf,k1,k2,k3
      
c------------------------ functions ------------------------------------

      real*8 dsIB2msq2gamma, dsIB2msqvector, dsIB2msqscalar

c-----------------------------------------------------------------------

      dsIB2BRtree=0.0d0
      msq=0.0d0
      result=0.0d0
      k2=kf
      if (k2.lt.1.or.k2.gt.12) then
        write(*,*) 'ERROR in dsIB2BRtree: non-existing fermion code kf = ',k2
        return
      endif
      
      
      if (decay.eq.'tWb') then
        if (k2.lt.11) then
          write(*,*) 'ERROR in dsIB2BRtree: tWb cannnot be used '//
     &               'for fermion code kf = ',k2
          return
        endif

        xim=(mass(kw)**2-mass(kb)**2)/mass(kt)**2
        xip=(mass(kw)**2+mass(kb)**2)/mass(kt)**2
        msq=0.5d0*dble(gl(kW,kt,kb))**2*
     &      (mass(kt)**2+mass(kb)**2-2*mass(kw)**2+
     &       (mass(kt)**2-mass(kb)**2)**2/mass(kw)**2)
        result=msq*dsIB2msq2gamma(mass(kt),xip,xim)/2.d0 ! 2.d0 is what top width 
                                                         ! is set to internally.
                                                         ! -> needs FIXING!!!

      elseif (decay.eq.'Hbt') then
        if (k2.lt.11) then
          write(*,*) 'ERROR in dsIB2BRtree: Hbt cannnot be used '//
     &               'for fermion code kf = ',k2
          return
        endif

        xim=(mass(kt)**2-mass(kb)**2)/mass(khc)**2
        xip=(mass(kt)**2+mass(kb)**2)/mass(khc)**2
        msq=dsIB2msqscalar(khc,kt,kb,xip)
        result=msq*dsIB2msq2gamma(mass(khc),xip,xim)/width(khc)

      elseif (decay.eq.'WfF'.or.decay.eq.'HfF') then
        if (k2.gt.10) then
          write(*,*) 'ERROR in dsIB2BRtree: WfF and HfF can not be used for '//
     &               '3rd generation quarks!'
          return
        endif

       
        if (decay.eq.'WfF') k1=kw
        if (decay.eq.'HfF') k1=khc
        k2=kf-mod(kf+1,2)          ! if k2 even, subtract 1 
        k3=k2+1
        if (mass(k1).gt.(mass(k2)+mass(k3))) then
          xim=(mass(k2)**2-mass(k3)**2)/mass(k1)**2
          xip=(mass(k2)**2+mass(k3)**2)/mass(k1)**2
          if (decay.eq.'WfF') msq=dsIB2msqvector(k1,k2,k3,xip,xim)
          if (decay.eq.'HfF') msq=dsIB2msqscalar(k1,k2,k3,xip)
          result=result+msq*dsIB2msq2gamma(mass(k1),xip,xim)/width(k1)
        endif
      
      elseif (decay.eq.'Zff'.or.
     &        decay.eq.'Hff'.or.decay.eq.'hff'.or.decay.eq.'Aff') then

        if (decay.eq.'Zff') k1=kz
        if (decay.eq.'Hff') k1=kh1
        if (decay.eq.'hff') k1=kh2
        if (decay.eq.'Aff') k1=kh3
        if (mass(k1).gt.(2*mass(k2))) then
          xip=2*mass(k2)**2/mass(k1)**2
          if (decay.eq.'Zff') then
            msq=dsIB2msqvector(k1,k2,k2,xip,0.d0)
          else 
            msq=dsIB2msqscalar(k1,k2,k2,xip)
          endif  
          result=result+msq*dsIB2msq2gamma(mass(k1),xip,0.d0)/width(k1)
        endif
        
      else
        write(*,*) 'ERROR in dsIB2BRtree: undefined input decay = ', decay
        return
      endif

      dsIB2BRtree=result

      return
      end




   
c...auxiliary functions used by dsib2BRtree
      real*8 function dsIB2msq2gamma(m,xip,xim)
      implicit none
      include 'dsmpconst.h'
      real*8 m,xip,xim      
      dsIB2msq2gamma=sqrt(1.-2*xip+xim**2)/(16.*pi*m)    
      return
      end
   
      real*8 function dsIB2msqscalar(k1,k2,k3,xip) ! k1 -> k2 k3
      implicit none
      include 'dsmssm.h'
      integer k1,k2,k3
      real*8 xip,xit,res      
      xit=mass(k2)*mass(k3)/mass(k1)**2
      res=(1.-xip)*(abs(gl(k1,k2,k3))**2+abs(gr(k1,k2,k3))**2)
     &     -2.*xit*dble(gr(k1,k2,k3)*conjg(gl(k1,k2,k3))+gl(k1,k2,k3)*conjg(gr(k1,k2,k3))) 
c MG  4 gr gl* -> 2 (gr gl* + gr* gl), color factor
      if(k2.ge.7) res=res*3
      dsIB2msqscalar=mass(k1)**2*res    
      return
      end
   
      real*8 function dsIB2msqvector(k1,k2,k3,xip,xim) ! k1 -> k2 k3
      implicit none
      include 'dsmssm.h'
      integer k1,k2,k3
      real*8 xip,xim,xit,res      
      xit=mass(k2)*mass(k3)/mass(k1)**2
      res=(1.-xip/2.-xim**2/2.)*(abs(gl(k1,k2,k3))**2+abs(gr(k1,k2,k3))**2)*
     &     2./3. +2.*xit*dble(gr(k1,k2,k3)*conjg(gl(k1,k2,k3))+gl(k1,k2,k3)*conjg(gr(k1,k2,k3))) 
c MG  - 4 gr gl -> 2 (gr gl* + gr* gl)   (note sign flip)
      if(k2.ge.7) res=res*3
      dsIB2msqvector=mass(k1)**2*res    
      return
      end
   
