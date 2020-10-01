**********************************************************************
*** dsrdinit initializes relic density routines and sets defaults
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: June 19, 2019      
**********************************************************************

      subroutine  dsrdinit
      implicit none
      include 'dsrdcom.h'

      call dsrdcom              ! make sure we set defaults from data file

      call dsrdset('dof','default') ! set default d.o.f. tables

      rdinit=1234 ! show that we have called dsrdinit

      return
      
      end

