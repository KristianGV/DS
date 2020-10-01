************************************************************************
*** subroutine dsib2ffWamps calculates helicity amplitudes for fFW final 
*** states.
***
***   Input: f                    - final state fermion
***                                 (see header of dsib2sigmav)
***
*** Author: Francesca Calore, 2013-04-10
***         Torsten Bringmann, 2014-02-20
************************************************************************
      subroutine dsib2ffWamps(f)
      implicit none
      include 'dsib2com.h'
      
      integer f, i, j 
      integer ieffsf, ieffsff
c-----------------------------------------------------------------------

c... s channel 
      call dsib2ffWsFSR(f, 1)
      call dsib2ffWsFSR(f, 2)
      call dsib2ffWsVIB(f)  
      do i = 1, 2                     ! sum over kcha(i)
         call dsib2ffWsISRH(f, i)
         call dsib2ffWsISRW(f, 1, i)
         call dsib2ffWsISRW(f, 2, i)
      end do  


c... t+u channel 
      if(f.eq.1.or.f.eq.3.or.f.eq.5) then ! sneutrinos in t+u channel 
         ieffsf = IB2sfgen 
         ieffsff = 2*IB2sfgen
      elseif(f.eq.2.or.f.eq.4.or.f.eq.6) then  
         ieffsf = 2*IB2sfgen 
         ieffsff = IB2sfgen  
      else                                ! squarks
         ieffsf = 2*IB2sfgen
         ieffsff = 2*IB2sfgen
      endif

c TB debug 2019-05-31: corrected(?) call to FSR routines      
c      write(*,*) 'ffW FSR: ', f, ieffsff, ieffsf
      do i = 1, ieffsff
         call dsib2ffWtuFSR(f, 1, i)         
      end do
      do i = 1, ieffsf
         call dsib2ffWtuFSR(f, 2, i)         
      end do
c      do i = 1, 2 ! sum over sfermions (sneutrino correctly added)
c         call dsib2ffWtuFSR(f, 1, i)         
c         call dsib2ffWtuFSR(f, 2, i)         
c      end do
c      write(*,*) 'ffW VIB: ', f, ieffsff, ieffsf
      do i = 1,ieffsf 
         do j = 1,ieffsff                           
            call dsib2ffWtuVIB(f, i, j)    
         end do 
      end do
 
c updated code for tuISR (Francesca Calore 01/2015)
c      write(*,*) 'ffW ISR: ', f, ieffsff, ieffsf
      do j =1,2
         do i = 1, ieffsff
            call dsib2ffWtuISR1(f,  i , j)
         end do
         do i = 1, ieffsf
            call dsib2ffWtuISR2(f,  i , j)
         end do
      end do
c      write(*,*) '(ffW done)'
 
      return
      end
