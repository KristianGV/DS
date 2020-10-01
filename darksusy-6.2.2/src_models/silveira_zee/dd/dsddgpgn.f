      subroutine dsddgpgn(gg,ierr)
*******************************************************************************
***  subroutine dsddgpgn returns DM nucleon four-fermion couplings          ***
***                                                                         ***
***  type : interface                                                       ***
***         [NB: this routine will eventually be phased out as *interface*  ***
***              routine. Currently, neutrino telescope routines are the    ***
***              only routines in /src that access it. ]                    ***
***                                                                         ***
***   output:                                                               ***
***     gg     : complex*16 (ddng,2) : four-fermion couplings in GeV^-2       ***
***                                                                         ***
***  author: paolo gondolo (paolo.gondolo@utah.edu) 2016                    ***
*******************************************************************************
      implicit none
      include 'dssilveira_zee.h'
      include 'dsddcom.h'
      complex*16 gg(ddng,2)
      integer i, ierr
      real*8 fn,gn,mn, dsmwimp

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dsanwx')

      ierr=0
      do i=1,ddng
        gg(i,1)=dcmplx(0.d0,0.d0)
        gg(i,2)=dcmplx(0.d0,0.d0)
      enddo

      fn=0.30d0
      mn=0.5d0*(m_n+m_p)
      gn=lambda*fn*mn/dsmwimp()/2.d0/mass(khsm)**2
      gg(1,1)=dcmplx(gn,0.d0)
      gg(1,2)=dcmplx(gn,0.d0)
      gg(4,1)=dcmplx(0.d0,0.d0)
      gg(4,2)=dcmplx(0.d0,0.d0)

      return
      end
