*******************************************************************************
*** Function dsanwx provides the  WIMP self-annihilation invariant rate.    ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  This function REPLACES the function with the same name in the          ***
***  generic_WIMP module, to take into account threshold effects when       ***
***  calculating threshold effects.                                         ***
***                                                                         ***
***  See the two ways of making generic_wimp_oh2, as defined in the makefile***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2017-05-11                                                         ***
*******************************************************************************
      real*8 function dsanwx(p)
      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsmpconst.h'

      real*8 p
      real*8 s, vmoeller, dsmwimp, pdv, dsmass, dswidth

      real*8 dsanthreshold

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dsanwx')


      dsanwx=0.0d0
      
      s = 4.0d0*(dsmwimp()**2+p**2)

c... CMS energy must be large enough to produce final states
c... This takes into account that *one* of the final states may be virtual (with mass 0) 
c... NB: This implies that we have to change the definition of thresholds consistently 
c... also in dsrdparticles, and hence replace that function as well -> see below
c      write(*,*) 'AAA ',p,dsmwimp(),dsmass(svch),dsmass(svch)**2,s ! JE TMP
      if (dsmass(svch)**2.gt.s) return 
      
      vmoeller = 4.*p/sqrt(s) ! Moeller velocity in the CMS

      pdv=sqrt(p**2+dsmwimp()**2)/2.0d0 ! p/vmoeller
      
      dsanwx= 4.0d0*pdv*sqrt(s)*(sva + svb*vmoeller**2)/gev2cm3s


c... Everything below this line is new compared to src_models/generic_wimp/an/dsanwx.f, 
c... and implements the threshold correction.

      if (dswidth(svch).gt.0.0d0) then !assume point-like interation -> n=0
         dsanwx = dsanwx
     &       * dsanthreshold(s,dsmass(svch),
     &       dsmass(svch),dswidth(svch),0)
c         write(*,*) 'BBB: ',p,sqrt(s),dsanwx,dsanthreshold(s,dsmass(svch),
c     &        dsmass(svch),dswidth(svch),0),dswidth(svch),
c     &        2.d0*dsmass(svch)-8.d0*dswidth(svch),
c     &        6.d0*dsmass(svch),dsmwimp(),dsmass(svch) ! JE TMP
      endif

            
      return
      end
      

*******************************************************************************
*******************************************************************************
c... -> need to change this to consistently define thresholds      
      subroutine dsrdparticles(option,nsize,ncoann,mcoann,dof,nrs,rm,rw,nthr,tm)
      implicit none
      include 'dsgeneric_wimp.h'
      integer option

      real*8 dsmwimp, dsmass, dswidth

      integer ncoann,nrs,nthr,nsize
      real*8 mcoann(*),dof(*),tm(*),
     &     rm(*),rw(*)
      integer kcoann(nsize)

c----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dsrdparticles')


c...Add DM particle to annihilation list
      ncoann=1
      mcoann(1)=mass(kdm)
      dof(1)=kdof(kdm)
      kcoann(1)=kdm

      nrs=0
      nthr=0
      
c... for a generic WIMP, a threshold can occur if the DM particle
c... is lighter than the final state particle
      if (dsmwimp().lt.dsmass(svch)) then
        nthr=3
        tm(1)=1.0d0*dsmass(svch) ! This is the value for 3-body final states 
        tm(2)=2.0d0*dsmass(svch)-16.0d0*dswidth(svch) ! empirical additional threshold
                                                     ! to help 3-body case
        tm(3)=2.0d0*dsmass(svch) ! This is the value for 2-body final states
      endif


      return

      end

      
