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
      include 'dsvdSIDM.h'
      include 'dsmpconst.h'

      integer istatus,info
c-----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dsmodelsetup')

      call dsnewidnumber ! This should always be the first call in dsmodelsetup,
                         ! and makes sure the model has a new unique ID number.
                         ! Several routines in src_models/ and in src/ 
                         ! *require* this to be set anew for every new model.

      istatus=0
      info=0

      if (mass(kdm).le.0.d0) then
         istatus=-1
         info=1
      endif

c... check that dimensionless coulings are not too large
      if (gDM.gt.sqrt(4*pi).or.gDR.gt.sqrt(4*pi)) then 
         istatus=-1
         info=2
      endif

      width(kdr) = 0.0d0 !

c... DM fermion, vector mediator
      if (DMtype.eq.1.and.Mediatortype.eq.2) then

        width(kmed) = gDR**2*mass(kmed)/12./pi
        if (mass(kdm).lt.(mass(kmed)/2.)) ! add decay width to DM
     &  width(kmed) = width(kmed) +   
     &     gDM**2*(1. + 2*mass(kdm)**2/mass(kmed)**2)*4./3. ! = |M|^2/mass(kmed)**2
     &     *sqrt(1.-4*mass(kdm)**2/mass(kmed)**2)/16./pi    ! phase space factor        

c... DM fermion, scalar mediator
      elseif (DMtype.eq.1.and.Mediatortype.eq.3) then

        width(kmed) = gDR**2*mass(kmed)/16./pi
        if (mass(kdm).lt.(mass(kmed)/2.)) ! add decay width to DM
     &  width(kmed) = width(kmed) +   
     &     gDM**2*(1. - 4*mass(kdm)**2/mass(kmed)**2)*2. ! = |M|^2/mass(kmed)**2
     &     *sqrt(1.-4*mass(kdm)**2/mass(kmed)**2)/16./pi ! phase space factor        
   
      endif

      end


