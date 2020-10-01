*******************************************************************************
*** Function dsdmspin returns the WIMP spin (in hbar)                       ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : WIMP spin                                                       ***
***                                                                         ***
*** author: Paolo Gondolo 2016-11-20                                        ***
*******************************************************************************

      real*8 function dsdmspin()
      implicit none

      include 'dssilveira_zee.h'

      dsdmspin=spin(kdm)

      return
      end
