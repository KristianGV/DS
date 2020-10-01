*******************************************************************************
*** Function dssisigtm provides the momentum-transfer cross section per DM  ***
*** mass, in conventional units of cm**2/g.                                 ***
***                                                                         ***
***  type : interface                                                       ***
***  desc : momentum-transfer cross section                                 ***
***                                                                         ***
***  Input:                                                                 ***
***    vkms - relative velocity of scattering DM particles (in km/s)        ***
***                                                                         ***
***  This interfacte function is required by dssisigtmav in src/.           ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-05-18                                                         ***
*******************************************************************************
      real*8 function dssisigtm(vkms)
      implicit none

      real*8 vkms

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dssisigtm')

      
      return
      end

