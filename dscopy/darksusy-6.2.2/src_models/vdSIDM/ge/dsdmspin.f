*******************************************************************************
*** Function dsdmspin returns the dark matter spin (in hbar)                ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : dark matter spin                                                ***
***                                                                         ***
*** author: Torsten Bringmann 2018-02-18                                    ***
*******************************************************************************

      real*8 function dsdmspin()
      implicit none

      include 'dsvdSIDM.h'

      dsdmspin=spin(kdm)

      return
      end
