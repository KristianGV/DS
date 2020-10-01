************************************************************************
*** subroutine dsib2ffZamps calculates helicity amplitudes for fFW final 
*** states.
***
***   Input: f                    - final state fermion
***                                 (see header of dsib2sigmav)
***
*** Author: Francesca Calore, 2013-04-10
***         Torsten Bringmann, 2014-02-20
************************************************************************
      subroutine dsib2ffZamps(f)
      implicit none
      include 'dsib2com.h'

      integer f, i, j 
      integer ieff
c-----------------------------------------------------------------------

 
c... s channel 
      call dsib2ffZsFSR(f, 1)
      call dsib2ffZsFSR(f, 2)
      do i = 1, 2
         call dsib2ffZsVIB(f, i) ! sum over kh1, kh2
      end do
      do i = 1, 4
         call dsib2ffZsISR(f, i) ! sum over kn(i)
      end do
      
      
c... t+u channel 
      if(f.eq.1.or.f.eq.3.or.f.eq.5) then   !sneutrinos in t+u channel
        ieff = IB2sfgen 
      else
        ieff = 2*IB2sfgen
      endif
      
      do i = 1,ieff              ! sum over sfermions
         call dsib2ffZtuFSR(f, 1, i)         
         call dsib2ffZtuFSR(f, 2, i)         
         do j = 1,ieff                           
            call dsib2ffZtuVIB(f, i, j)    
         end do
         do j = 1, 4
            call dsib2ffZtuISR(f, 1, i , j)
            call dsib2ffZtuISR(f, 2, i , j)
         end do
      end do     


      return
      end

