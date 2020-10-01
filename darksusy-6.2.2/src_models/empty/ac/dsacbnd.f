*******************************************************************************
***  subroutine dsacbnd checks collider bounds                              ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : Check collider bounds                                           ***
***                                                                         ***
*** author: Torsten.Bringmann@fys.uio.no                                    ***
*** date  : 2015-06-11                                                      ***
*******************************************************************************
      subroutine dsacbnd(excl)
      implicit none
      integer excl

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsacbnd')
  
c... empty model is not excluded
      excl=0
      return
      end
