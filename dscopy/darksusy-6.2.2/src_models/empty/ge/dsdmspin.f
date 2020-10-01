*******************************************************************************
*** Function dsdmspin returns the dark matter spin (in hbar)                ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : dark matter spin                                                ***
***                                                                         ***
*** author: Paolo Gondolo 2016-11-20                                        ***
*******************************************************************************

      real*8 function dsdmspin()
      implicit none

      include 'dsempty.h'

      dsdmspin=spin(kdm)

      return
      end
