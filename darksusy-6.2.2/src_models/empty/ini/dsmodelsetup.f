*******************************************************************************
***  subroutine dsmodelsetup sets up a particle model                       ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : Sets up a new particle physics model in module                  ***
***                                                                         ***
*** author: Torsten.Bringmann@fys.uio.no                                    ***
*** date 2014-05-09                                                         ***
*******************************************************************************
      subroutine dsmodelsetup(unphys,warning)
      implicit none
      include 'dsempty.h'

      integer unphys,warning
c-----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsmodelsetup')

      call dsnewidnumber ! This should always be the first call in dsmodelsetup,
                         ! and makes sure the model has a new unique ID number.
                         ! Several routines in src_models/ and in src/ 
                         ! *require* this to be set anew for every new model.



      unphys=0
      warning=0

      kdm=1
      mass(kdm)=1.d2

      end


