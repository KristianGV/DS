      subroutine dsrdreaddof(filename)
c_______________________________________________________________________
c     read table of effective degrees of freedom in the early universe
c     See header of dsrddof for a definition of the quantities provided
c     input:
c       filename - character*(*) - name of file containing table of dof
c     common:
c       'dsrdcom.h' - included common blocks
c  author: paolo gondolo (paolo@physics.utah.edu) 2005
c  mod: 2017-05-12 torsten bringmann (re-added fe = geff)
c  mod: 2018-02-02 torsten bringmann (initialized klo,khi to allow model-independent
c                                     use of dsrddof)
c=======================================================================
      implicit none
      character*(*) filename
      character*300 msg
      character*200 scr
      include 'dsrdcom.h'
      include 'dsio.h'
      integer i
c-----------------------------------------------------------------------
      write (msg,'(a,a)') 'dsrdreaddof: reading ',filename
      if (prtlevel.ge.2) call dswrite(0,0,msg)
      open (unit=99,file=filename)
      read (99,'(A)',err=2000) scr
      read (99,'(A)',err=2000) scr
      read (99,'(A)',err=2000) scr
      i=1
 100  continue
c... added fe (needed e.g. for KD routines)
c... NB: this does not cause problems even if fe(i) is not provided (as in HP tables) 
      read (99,*,end=1000,err=2000) tgev(i),fg(i),fh(i),fe(i)
      i=i+1
      goto 100
 1000 continue
      close(99)
      nf=i-1
      klo=0
      khi=nf
      return
 2000 continue
      close(99)
      write (*,*) 'dsrdreaddof: error while reading ',filename
      end
