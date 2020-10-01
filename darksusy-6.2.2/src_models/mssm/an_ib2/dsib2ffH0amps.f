************************************************************************
*** subroutine dsib2ffH0amps calculates helicity amplitudes for fFH final 
*** states with neutral, CP-even scalars.
***
***   Input: f                    - final state fermion
***                                 (see header of dsib2sigmav)
***          hi=1,2               - SM Higgs (1) or heavy Higgs (2)
***
*** Author: Francesca Calore, 2014-02-20
************************************************************************
      subroutine dsib2ffH0amps(f,hi)
      implicit none
      include 'dsib2com.h'
     
      integer f, hi, i, j 
      integer ieff
c-----------------------------------------------------------------------

c... s channel 
      call dsib2ffH0sFSR(f, 1, hi)
      call dsib2ffH0sFSR(f, 2, hi)

      call dsib2ffH0sVIB(f, hi)  

      do i = 1, 4
         call dsib2ffH0sISR(f, i, hi) ! sum over kn(i)
      end do


c... t+u channel 
      if(f.eq.1.or.f.eq.3.or.f.eq.5) then   !sneutrinos in t+u channel
        ieff = IB2sfgen  
      else
        ieff = 2*IB2sfgen 
      endif
      
      do i = 1,ieff              ! sum over sfermions
         call dsib2ffH0tuFSR(f, 1, i, hi)         
         call dsib2ffH0tuFSR(f, 2, i, hi)         
         do j = 1,ieff                           
            call dsib2ffH0tuVIB(f, i, j, hi)    
         end do
         do j = 1, 4
            call dsib2ffH0tuISR(f, 1, i , j, hi)
            call dsib2ffH0tuISR(f, 2, i , j, hi)
         end do
      end do     

      return
      end
