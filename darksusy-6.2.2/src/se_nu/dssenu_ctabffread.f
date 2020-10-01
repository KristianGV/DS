      subroutine dssenu_ctabffread(wh,i,file)

***********************************************************************
*** Reads in tabulated capture rates (apart from couplings)
*** Input: wh = 'su' or 'ea' for sun or earth (Earth not yet implemented)
***        i = table number to read
***        file = file name to read
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2015-06-12
***********************************************************************

      implicit none
      include 'dssecap.h'
      include 'dshmcom.h'
      include 'dsio.h'

      integer i,l,j,ncfffile

      character*2 wh,whfile
      character*12 vdf,vdfc

      character*200 file
      logical bad

c...Now read in the data file

c...Make sure 'num' and 'numc' is treated consistenly     
      vdfc=veldf
      if (vdfc.eq.'numc') vdfc='num'

      if (prtlevel.ge.2) write(*,*) 'dssenu_ctabffread: Opening file ',file
      open(unit=13,file=file,
     &  form='formatted',status='old',err=200)


      read(13,101) whfile,vdf,ncfffile
 101  format(1x,1x,A2,1x,A12,1x,I5,1x,A)


c...Check consistency of files
      bad=.false.
      if (wh.eq.'su'.or.wh.eq.'SU') then
        if (whfile.ne.'su'.and.whfile.ne.'SU') bad=.true.
        if (vdf.ne.vdfc) bad=.true.
        if (ncff.ne.ncfffile) bad=.true.
        if (bad) then
          write(*,*) 'ERROR in dssenu_ctabffread: ',
     &      'File type mismatch.'
          write(*,*) 'Tried to read wh=',wh,' and vdfc=',
     &      vdfc,' and ncff=',ncff
          write(*,*) 'but found wh=',whfile,' and vdfc=',vdfc,
     &      ' and ncff=',ncfffile,
     &      ' in file. Stopping.'
          stop
        endif
      else                      ! if (wh.eq.'ea'.or.wh.eq.'EA') then
        write(*,*) 'DS ERROR in dssenu_ctabffread:'
        write(*,*) '  Earth not implemented yet. Stopping.'
        stop
c        if (whfile.ne.'ea'.and.whfile.ne.'EA') bad=.true.
c        if (vdf.ne.veldfearth) bad=.true.
c        if (bad) then
c          write(*,*) 'ERROR in dssenu_ctabread: ',
c     &      'File type mismatch.'
c          write(*,*) 'Tried to read wh=',wh,' and veldfearth=',
c     &      veldfearth,'.'
c          write(*,*) 'but found wh=',whfile,' and veldfearth=',vdf,
c     &      'in file. Stopping.'
c          stop
c        endif
      endif

      do l=0,ncff
        if (wh.eq.'su'.or.wh.eq.'SU') then
          read(13,*) (ctabffsu(l,i,j),j=1,6)
        else
c          read(13,*) (ctabffea(l,i,j,j=1,6))
        endif
      enddo

      close(13)
      
      return

 200  continue  ! we get here if file didn't exist

      close(13)

c...Create file
      write(*,*) 'WARNING in dssenu_ctabffread:'
      write(*,*) 'The requested file ',file,' does not exist.'
      write(*,*) 'I will create it for you',
     &  ' (this only needs to be done once),'
      write(*,*)  'but it could take several hours/days.',
     &  ' Be patient please...'
      call dssenu_ctabffcreate(wh,i)
      call dssenu_ctabffwrite(wh,i,file)

      return

      end


      

