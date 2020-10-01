*******************************************************************************
***  subroutine dsmodelsetup sets up a particle model                       ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2016-05-09                                                         ***
*******************************************************************************
      subroutine dsmodelsetup(istatus,info)
      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsio.h'

      integer istatus,info
c-----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dsmodelsetup')

c... This should always be the first call in dsmodelsetup, and makes sure the 
c... model has a new unique ID number. Several routines in src_models/ and 
c... in src/ *require* this to be set anew for every new model.
      call dsnewidnumber
      
      
      istatus=0
      info=0

      if (mass(kdm).le.0.d0) then
         istatus=-1
         info=1
      endif
      if (sva.lt.0.d0) then !.or.sva+svb.lt.0) then
         istatus=-1
         info=2
      endif
      if ((spin(kdm).ne.0.0d0).and.(spin(kdm).ne.0.5d0)) then ! Spin 1 not yet implemented!
         istatus=-1
         info=3
      endif 
      if (sigsdp.lt.0.d0.or.sigsdn.lt.0.d0.or.sigsip.lt.0.d0.or.sigsin.lt.0.d0) 
     & then
         istatus=-1
         info=4
      endif

      if (sigsdp.ne.0.0d0.or.sigsdn.ne.0.0d0) then
        if (spin(kdm).eq.0.0d0) then
          if (prtlevel.gt.0) then
            write(*,*) 'WARNING: Will set spin-dependent couplings to zero'
            write(*,*) 'because the dark matter particles have no spin!'
          endif
          sigsdp = 0.0d0
          sigsdn = 0.0d0
          info=10
        endif
      endif


      end


