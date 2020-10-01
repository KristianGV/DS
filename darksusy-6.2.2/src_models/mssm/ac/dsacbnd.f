*******************************************************************************
***  subroutine dsacbnd checks collider bounds                              ***
***                                                                         ***
***                                                                         ***
*******************************************************************************
      subroutine dsacbnd(warning)
      implicit none
      include 'dsaccom.h'
      integer warning
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dsacbnd')
     
      
      if (aclabel.eq.'default') then
         call dsacbnd13(warning)
      elseif (aclabel.eq.'preDS6') then
         call dsacbnd12(warning)
      elseif (aclabel.eq.'pdg2011b') then
         call dsacbnd11(warning)
      else
         write(*,*) 'Error in dsacbnd: invalid option: ',aclabel
         stop
      endif
      return
      end
