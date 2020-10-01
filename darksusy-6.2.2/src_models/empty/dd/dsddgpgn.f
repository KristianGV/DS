*******************************************************************************
***  subroutine dsddgpgn returns DM nucleon four-fermion couplings          ***
***                                                                         ***
***  type : interface                                                       ***
***         [NB: this routine will eventually be phased out as *interface*  ***
***              routine. Currently, neutrino telescope routines are the    ***
***              only routines in /src that access it. ]                    ***
***                                                                         ***
***  desc : Four fermion couplings.
***                                                                         ***
***   output:                                                               ***
***     gg complex*16 array of couplings                                     ***
***          first index is operator index (1:ddng), with ddng in dsddcom.h ***
***          second index is proton/neutron (1:2)                           ***
***     units: GeV^-2                                                       ***
***                                                                         ***
***  author: torsten.bringmann@fys.uio.no 18-11-2016                        ***
***    2016-11-20 Paolo Gondolo abandoned scheme                            ***
*******************************************************************************
      subroutine dsddgpgn(gg,ierr) 

      implicit none
      include 'dsddcom.h'
      complex*16 gg(ddng,2)
      integer ierr,i

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsddgpgn')

      ierr=0
      do i=1,ddng
        gg(i,1)=dcmplx(0.d0,0.d0)
        gg(i,2)=dcmplx(0.d0,0.d0)
      enddo

      return
      end
