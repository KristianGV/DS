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
***  author: paolo gondolo (paolo@physics.utah.edu) 2004                    ***
***    2008-02-16 paolo gondolo double spin-dependent squark contribution   ***
***    2008-02-16 paolo gondolo undo running quark masses in couplings      ***
***    2016-05-19 Paolo Gondolo, Joakim Edsjo, scheme setup                 ***
***    2016-11-20 Paolo Gondolo abandoned scheme                            ***
***    2018-10-24 Torsten Bringmann, nucleon x-sections now in common block ***
*******************************************************************************
      subroutine dsddgpgn(gg,ierr) 

      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsddcom.h'
      include 'dsmpconst.h'

      complex*16 gg(ddng,2)
      integer ierr,i

      real*8 mwimp,dsmwimp
      real*8 fkinp,fkinn,gps,gns,gpa,gna

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dsddgpgn')

      ierr=0
      do i=1,ddng
        gg(i,1)=dcmplx(0.d0,0.d0)
        gg(i,2)=dcmplx(0.d0,0.d0)
      enddo

      mwimp=dsmwimp()

      fkinp = 4./pi*(m_p*mwimp/(m_p+mwimp))**2
      fkinn = 4./pi*(m_n*mwimp/(m_n+mwimp))**2

c...NOTE: By providing cross sections, we have lost the information
c...on the sign of the g's. We will here assume that they are positive and
c...of the same sign, but this does not need to be the case. 
c...For a more general use of the DD routines, replace this function 
c...and provide the g's directly in any other way.

      gps=2.0d0*sqrt(sigsip/fkinp/gev2cm2)
      gns=2.0d0*sqrt(sigsin/fkinn/gev2cm2)

      gpa=sqrt(sigsdp/gev2cm2/fkinp*4.d0/3.d0)
      gna=sqrt(sigsdn/gev2cm2/fkinn*4.d0/3.d0)

      gg(1,1)=dcmplx(gps,0.d0)
      gg(1,2)=dcmplx(gns,0.d0)
      gg(4,1)=dcmplx(gpa,0.d0)
      gg(4,2)=dcmplx(gna,0.d0)

      return
      end
