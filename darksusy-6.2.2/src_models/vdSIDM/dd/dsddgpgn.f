*******************************************************************************
***  subroutine dsddgpgn returns DM nucleon four-fermion couplings          ***
***                                                                         ***
***  type : interface                                                       ***
***         [NB: this routine will eventually be phased out as *interface*  ***
***              routine. Currently, neutrino telescope routines are the    ***
***              only routines in /src that access it. ]                    ***
***                                                                         ***
***   output:                                                               ***
***     gg complex*16 array of couplings                                     ***
***          first index is operator index (1:ddng), with ddng in dsddcom.h ***
***          second index is proton/neutron (1:2)                           ***
***     units: GeV^-2                                                       ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date  : 2018-02-18                                                      ***
*******************************************************************************
      subroutine dsddgpgn(gg,ierr) 
      implicit none
      include 'dsddcom.h'
      
      complex*16 gg(ddng,2)
      integer i, ierr


c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dsddgpgn')

c... currently, the coupling of dark sector to SM particles is assumed to be
c... negligible in the vdSIDM module, and direct detection signals are therefore
c... set to zero. See the generic_WIMP module for a simple way of adding 
c... non-zero CR source terms.

      ierr=0
      do i=1,ddng
        gg(i,1)=dcmplx(0.d0,0.d0)
        gg(i,2)=dcmplx(0.d0,0.d0)
      enddo

      return
      end
