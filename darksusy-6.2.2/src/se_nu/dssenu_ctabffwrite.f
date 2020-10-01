      subroutine dssenu_ctabffwrite(wh,i,file)

***********************************************************************
*** Writes out tabulated capture rates (apart from couplings)
*** Input: wh = 'su' or 'ea' for sun or earth
***        i = table number to write
***        file = file to write to
*** Author: Joakim Edsjo
*** Date: 2015-06-12
***********************************************************************

      implicit none
      include 'dssecap.h'
      include 'dshmcom.h'
      include 'dsver.h'

      integer i,l,j

      character*2 wh

      character*200 file


c...Now write the data file
     
      write(*,*) 'dssenu_ctabffwrite: Opening file ',file
      open(unit=13,file=file,
     &  form='formatted',status='unknown')

      if (wh.eq.'su'.or.wh.eq.'SU') then
        write(13,101) wh,veldf,ncff,dsversion
      else
        write(*,*) 'DS ERROR in dssenu_ctabffwrite:'
        write(*,*) 'Earth not implemented. Stopping.'
        stop
c        write(13,101) wh,veldfearth,dsversion
      endif
 101  format('#',1x,A2,1x,A12,1x,I5,1x,A)

      do l=0,ncff
        if (wh.eq.'su'.or.wh.eq.'SU') then
          write(13,'(6(1x,E14.8))') (ctabffsu(l,i,j),j=1,6)
        else
c          write(13,'(6x,E14.8)') (ctabffea(l,i,j),j=1,6)
        endif
      enddo
      
      close(13)

      return

      end


      

