***********************************************************************
*** Routine dsib2set sets default settings for electroweak internal 
*** bremsstrahlung (IB2) calculations. Must be called before any of the
*** IB2 routines can be used (typically done in dsinit)!
*** 
***  c - character string specifying choice to be made 
***  possible values for c:
***
***  'default' - include all ffB channels (where B is any Boson)
***  'off'     - no electroweak IB contribution
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-03-06
***********************************************************************

      subroutine dsib2set(c)
      implicit none
      include 'dsib2com.h'
      include 'dsmssm.h'
      character*(*) c
      integer i
c-----------------------------------------------------------------------

      IB2sfgen = 1 ! only take into account dominant sfermion flavour
                   ! in propagators (IB2sfgen = 3 for all)

c...all channels
      if (c.eq.'standard') then
          do i=1,612
            IB2flag(i)=1
          enddo 
            ib2svhow=1      ! first integrate over E1, then EV to obtain sv
            ib2dnhow='full' ! use full SUSY expression for dN/dE 
                            ! distributions of particles in 3-body final state
            IB2acc=0.05     ! accuracy goal for sv integration.
            intres=4.d0     ! integration around resonance (in units of m*\Gamma)
            

c...FSR option
      elseif (c.eq.'FSR') then
            ! not yet implemented  
            ib2dnhow='FSR'  ! only include model-independent FSR for dN/dE 
                            ! distributions of particles in 3-body final state
          
c...Off option
      elseif (c.eq.'off'.or.c.eq.'default') then
         do i=1,612
            IB2flag(i)=0
         enddo

      else
         write(*,*) 'ERROR in dsIB2set -- unknown option: ',c
         write(*,*) 'Stopping...'
         stop
      endif

      return
      end





        

