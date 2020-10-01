*******************************************************************************
***  subroutine dsmodelsetup sets up a particle model                       ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2016-02-06                                                         ***
*******************************************************************************
      subroutine dsmodelsetup(unphys,warning)
      implicit none
      include 'dsgeneric_decayingDM.h'

      integer unphys,warning
      
      real*8 dsmass,dsmwimp
      integer i,j,k
c-----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_decayingDM','dsmodelsetup')

      call dsnewidnumber ! This should always be the first call in dsmodelsetup,
                         ! and makes sure the model has a new unique ID number.
                         ! Several routines in src_models/ and in src/ 
                         ! *require* this to be set anew for every new model.

      unphys=0
      warning=0
      
c... check that the DM mass is positive
      if (dsmwimp().le.0.0d0) then
        unphys=-1
        warning=111
      endif

c... check that lifetime is positive
      if (Gammatot.le.0.0d0) then
        unphys=-1
        warning=warning+222
      endif

c... remove kinematically disallowed channels from the list of all available channels
      i=1
      do 10 while (i.le.numdecch2b)
        if (dsmwimp().lt.(dsmass(dec_2body(i,1))+dsmass(dec_2body(i,2)))) then
          if (decBR(i).ne.0.0d0) warning = 1 ! The model was set up with a non-zero decay
                                             ! rate to a kinematically disallowed channel.
                                             ! This channel will be treated as invisible
          do j=i,numdecch2b-1
            decBR(j) = decBR(j+1)
            do k=1,6
              dec_2body(j,k) = dec_2body(j+1,k) 
            enddo
          enddo
          numdecch2b = numdecch2b - 1
        else
          i=i+1
        endif
  10  continue


c... same, for channels that result in monochromatic CRs
      i=1
      do 20 while (i.le.numyieldch_line)      
        if (dsmwimp().lt.(dsmass(yieldchannels_line(i,1))+
     &                   dsmass(yieldchannels_line(i,2)))) then
          do j=i,numyieldch_line-1
              yieldchannels_line(j,1) = yieldchannels_line(j,1) 
              yieldchannels_line(j,2) = yieldchannels_line(j,2) 
          enddo
          numyieldch_line = numyieldch_line - 1
        else
          i=i+1
        endif
  20  continue


      
      return
      end

