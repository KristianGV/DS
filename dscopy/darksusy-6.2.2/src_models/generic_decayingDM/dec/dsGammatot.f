**********************************************************************
*** function dsGammatot returns the *total* decay rate / width of the
*** DM particle.
***                                                                         
***  type : interface                                                       
***                                                                         
*** Units of returned width: 1/s
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2016-02-06
**********************************************************************

      real*8 function dsGammatot()
      implicit none
      include 'dsgeneric_decayingDM.h'

     
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_decayingDM','dssigmav0tot')

c... This common block value is set in dsgivemodel_decayingDM
      dsGammatot=Gammatot
    
      return
      end
