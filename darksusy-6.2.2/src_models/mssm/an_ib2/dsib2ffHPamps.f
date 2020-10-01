************************************************************************
*** subroutine dsib2ffHPamps calculates helicity amplitudes for fFH final 
*** states with charged scalars.
***
***   Input: f                    - final state fermion
***                                 (see header of dsib2sigmav)
***
*** Author: Francesca Calore, 2014-02-XX
************************************************************************
      subroutine dsib2ffHPamps(f)
      implicit none
      include 'dsib2com.h'
      
      integer f, i, j 
      integer ieffsf, ieffsff
c-----------------------------------------------------------------------

c... s channel 
      call dsib2ffHPsFSRH3(f, 1)
      call dsib2ffHPsFSRH3(f, 2)
      call dsib2ffHPsFSRZ(f, 1)
      call dsib2ffHPsFSRZ(f, 2)
      call dsib2ffHPsVIB(f)  
      do i = 1, 2                     ! sum over kcha(i)
         call dsib2ffHPsISRH(f, i)
         call dsib2ffHPsISRW(f, i)
      end do

c... t+u channel 
      if(f.eq.1.or.f.eq.3.or.f.eq.5) then ! !sneutrinos in t+u channel 
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
c      write(*,*) 'ffHp FSR: ', f, ieffsff, ieffsf
      do i = 1, ieffsff
         call dsib2ffHPtuFSR(f, 1, i)         
      end do
      do i = 1, ieffsf
         call dsib2ffHPtuFSR(f, 2, i)         
      end do
c      do i = 1, 2 ! sum over sfermions (sneutrino correctly added)
c         call dsib2ffHPtuFSR(f, 1, i)         
c         call dsib2ffHPtuFSR(f, 2, i)         
c      end do
c      write(*,*) 'ffHp VIB: ', f, ieffsff, ieffsf
      do i = 1,ieffsf 
         do j = 1,ieffsff                           
            call dsib2ffHPtuVIB(f, i, j)    
         end do 
      end do

c      write(*,*) 'ffHp ISR: ', f, ieffsff, ieffsf
      do j =1,2         
         do i = 1, ieffsff  
            call dsib2ffHPtuISR1(f,  i , j) 
         end do
         do i = 1, ieffsf  
            call dsib2ffHPtuISR2(f,  i , j)
         end do
      end do
c      write(*,*) '(ffHp done)'


      return
      end
