*****************************************************************************
***   subroutine dsanreadfiles loads (from disk) the the requested data
***   files. yieldk is
***   the yield type (51,52 or 53 (or 151, 152, 153)) for
***   positron yields, cont. gamma or muon neutrino yields respectively.
***   yieldk is used to check that the provided data file is of the
***   correct type.  if yieldk=51,52 or 53 integrated yields are loaded and
***   if yieldk =151, 152 or 153, differential yields are loaded.
***   Author: Joakim Edsjo
***   edsjo@fysik.su.se date: 96-10-23 (based on dsmuinit.f version
***   3.21) 
***   modified: 98-01-26
***   modified: 09-10-20 pat scott pat@fysik.su.se
***   modified: 04-05-09 Joakim Edsjo, edsjo@fysik.su.se
***   modified: April, 2016, Joakim Edsjo, to read new simulation yields
***     (including dbars)      
*****************************************************************************

      subroutine dsanreadfiles(yieldk)
      implicit none
      include 'dsio.h'
      include 'dsanyieldcom.h'

c------------------------ variables ------------------------------------

      integer i,j,k,l,m,yieldk,yieldkfile,fi,fltype,fl,dbkind
      integer nztmp,ndtmp
      integer chfile
      real*8 mfile
      character*200 filein
      character*255 scr
      character*10 scr2
      character*15 filepref
      character*3 filenr
      character*10 filesuf

      logical first
      data first/.true./
      save first

      if (.not.dsanyieldinitcalled) then
         write(*,*) 'ERROR in dsanreadfiles: dsanyield_init needs to be',
     &     ' called at startup.'
         write(*,*) 'This is usually done by dsinit. Have you',
     &     ' called dsinit?'
         write(*,*) 'Program stopping.'
         stop
      endif

c------------------------------------------------------ load yield tables

      call dsandec(yieldk,fltype,fi)
      if (fltype.le.2) then
         dbkind=1               ! yields
      else
         dbkind=2               ! errors
      endif

c...generate file name
      if (anftype.eq.'a'.or.anftype.eq.'a') then
        filesuf='.dat'
      else
        filesuf='.bin'
      endif
      if (fltype.eq.1.or.fltype.eq.3) then   ! integrated
        filepref='simint'
        write(filenr,'(i3)') mod(yieldk,1000) ! treat error numbers sep.
      else                    ! differential
        filepref='simdiff'
        write(filenr,'(i3)') mod(yieldk,1000) ! treat error numbers sep.
      endif
      if (fltype.le.2) then     ! yields
         call dsdatafile(filein,filepref//filenr//filesuf)
      else                      ! error on yields
         call dsdatafile(filein,filepref//filenr
     &     //'-err'//filesuf)
      endif

c...delete spaces in file name
      fl=200
      do l=1,fl
 40     if (filein(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            filein(m:m)=filein(m+1:m+1)
          enddo
          if (fl.eq.l) goto 50
          goto 40
        endif
      enddo
 50   continue

C========== NEW CODE FROM HERE ==========
      
      if (fltype.eq.1.or.fltype.eq.3) then ! integrated yields
c...integrated yields
c        write(*,*) 'enter file name for integrated yield tables:'
c        read(*,'(a)') filein
c        write(*,*) 'is this file ascii (a) or binary (b)?'
c        read(*,*) filetype

c---------- All yield types but dbar ----------
       if (fi.lt.9.or.fi.gt.13) then ! not dbar
c...clear the table
          do j=1,nch
            do k=1,nmass
              do l=0,zn
                phiint(l,k,j,fi)=0.0d0
              enddo
            enddo
          enddo

          if (prtlevel.gt.1) 
     &         write(*,*) 'loading integrated yield tables from file ',
     &         filein
          if (anftype.eq.'a'.or.anftype.eq.'A') then
            open(unit=13,file=filein,status='old',form='formatted')
            read(13,500) scr2,yieldkfile
            do i=1,11
               read(13,'(a)') scr ! read header lines
            enddo
          else
          open(unit=13,file=filein,status='old',form='unformatted')
          read(13) yieldkfile
          endif
 500      format(1x,a10,1x,i4)

          if (yieldk.ne.yieldkfile) then
            write(*,*) 'DS ERROR in dsanreadfiles: '
            write(*,*) 'the requested yield type is inconsistent with',
     &        ' the provided data file.'
            write(*,*) 'data file: ',filein
            write(*,*) 'requested yieldtype:   ',yieldk
            write(*,*) 'data file''s yieldtype: ',yieldkfile
            write(*,*) 'program stopped.'
            stop
          endif

          yieldtype(fltype,fi)=yieldkfile
          do j=1,nch
            if (prtlevel.gt.2) write(*,*) '      channel number ',j
            do k=1,nmass
              if (k.ge.milow(j)) then
                if (anftype.eq.'a'.or.anftype.eq.'A') then
                  read(13,*) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  read(13,2000) (phiint(l,k,j,fi),l=0,zn-1)
                else
                  read(13) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  read(13) (phiint(l,k,j,fi),l=0,zn-1)
                endif
              endif
            enddo
          enddo
          close(13)

c---------- dbar yields ----------
        else ! dbar
c...clear the table
          do j=1,nch
             do k=1,nmass
                do m=-1,dpn-1
                   do l=0,zndb
                      phiintdb(l,m,k,j,fi,dbkind)=0.0d0
                   enddo
                enddo
             enddo
          enddo

          if (prtlevel.gt.1) 
     &         write(*,*) 'loading integrated yield tables from file ',
     &         filein
          if (anftype.eq.'a'.or.anftype.eq.'A') then
            open(unit=13,file=filein,status='old',form='formatted')
            read(13,510) scr2,yieldkfile
            read(13,510) scr2,dbct
            read(13,501) scr2,dbp0fit(fi)
            read(13,501) scr2,dbp0low(fi)
            read(13,501) scr2,dbp0high(fi)
            read(13,'(a)') scr
            read(scr(19:50),*) nztmp
            read(13,'(a)') scr
            read(scr(19:50),*) ndtmp
c            read(13,502) scr3,nztmp
c            read(13,502) scr3,ndtmp

 510      format(1x,a10,1x,i4)

            if (nztmp.ne.zndb) then
               write(*,*) 'DS Error in dsanreadfiles.',
     &          ' Mismatch in z-bins: ',
     &              nztmp
               stop
            endif
            if (ndtmp.ne.dpn) then
               write(*,*) 'DS Error in dsanreadfiles.',
     &          ' Mismatch in dp-bins: ',
     &              dpn
               stop
            endif
            do i=1,12
               read(13,'(a)') scr ! read header lines
            enddo
          else
            open(unit=13,file=filein,status='old',form='unformatted')
            read(13) yieldkfile
            read(13) dbct
            read(13) dbp0fit(fi)
            read(13) dbp0low(fi)
            read(13) dbp0high(fi)

          endif

 501      format(1x,a10,1x,e13.6)
c 502      format(1x,a17,1x,I4)

          if (yieldk.ne.yieldkfile) then
            write(*,*) 'DS ERROR in dsanreadfiles: '
            write(*,*) 'the requested yield type is inconsistent with',
     &        ' the provided data file.'
            write(*,*) 'data file: ',filein
            write(*,*) 'requested yieldtype:   ',yieldk
            write(*,*) 'data file''s yieldtype: ',yieldkfile
            write(*,*) 'program stopped.'
            stop
          endif

          yieldtype(fltype,fi)=yieldkfile
          do j=1,nch
            if (prtlevel.gt.2) write(*,*) '      channel number ',j
            do k=1,nmass
              if (k.ge.milow(j)) then
                if (anftype.eq.'a'.or.anftype.eq.'A') then
                  read(13,*) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  do m=0,dpn-1
                     read(13,2000)
     &                 (phiintdb(l,m,k,j,fi,dbkind),l=0,zndb-1)
                  enddo
                else
                  read(13) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  do m=0,dpn-1
                    read(13) (phiintdb(l,m,k,j,fi,dbkind),l=0,zndb-1)
                  enddo
                endif
              endif
            enddo
          enddo
          close(13)
        endif


c========== Differential yields ==========
      else   ! differential yields

c...differential yields
c        write(*,*) 'enter file name for differential yield tables:'
c        read(*,'(a)') filein
c        write(*,*) 'is this file ascii (a) or binary (b)?'
c        read(*,*) filetype

c--------- All yield types but dbar
         if (fi.lt.9.or.fi.gt.13) then ! all but dbar
c...clear yield table
          do j=1,nch
            do k=1,nmass
              do l=-1,zn
                phidiff(l,k,j,fi)=0.0d0
              enddo
            enddo
          enddo

          if (prtlevel.gt.1) 
     &         write(*,*) 'loading differential yield tables from file ',
     &         filein
          if (anftype.eq.'a'.or.anftype.eq.'a') then
            open(unit=13,file=filein,status='old',form='formatted')
            read(13,500) scr2,yieldkfile
            do i=1,11
              read(13,'(a)') scr ! header lines
            enddo
          else
            open(unit=13,file=filein,status='old',form='unformatted')
            read(13) yieldkfile
          endif

          if (yieldk.ne.yieldkfile) then
            write(*,*) 'DS ERROR in dsanreadfiles: '
            write(*,*) 'the requested yield type is inconsistent with',
     &        ' the provided data file.'
            write(*,*) 'data file: ',filein
            write(*,*) 'requested yieldtype:   ',yieldk
            write(*,*) 'data file''s yieldtype: ',yieldkfile
            write(*,*) 'program stopped.'
            stop
          endif

          yieldtype(fltype,fi)=yieldkfile
          do j=1,nch
            if (prtlevel.gt.2) write(*,*) '      channel number ',j
            do k=1,nmass
              if (k.ge.milow(j)) then
                if (anftype.eq.'a'.or.anftype.eq.'a') then
                  read(13,*) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  read(13,2000) (phidiff(l,k,j,fi),l=0,zn-1)
                else
                  read(13) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  read(13) (phidiff(l,k,j,fi),l=0,zn-1)
                endif
                do l=0,zn-1 ! correct units of dyield / dz
                  phidiff(l,k,j,fi)=phidiff(l,k,j,fi)/
     &              (dz(l))
                enddo
                phidiff(-1,k,j,fi)=phidiff(0,k,j,fi)
                phidiff(zn,k,j,fi)=phidiff(zn-1,k,j,fi)
              endif
            enddo
          enddo
          close(13)

c---------- dbar yields ----------
        else ! dbar
c...clear yield table
          do j=1,nch
            do k=1,nmass
              do m=-1,dpn-1
                do l=-1,zndb
                  phidiffdb(l,m,k,j,fi,dbkind)=0.0d0
                enddo
              enddo
            enddo
          enddo

          if (prtlevel.gt.1) 
     &         write(*,*) 'loading differential yield tables from file ',
     &         filein
          if (anftype.eq.'a'.or.anftype.eq.'a') then
            open(unit=13,file=filein,status='old',form='formatted')
            read(13,510) scr2,yieldkfile
            read(13,510) scr2,dbct
            read(13,501) scr2,dbp0fit(fi)
            read(13,501) scr2,dbp0low(fi)
            read(13,501) scr2,dbp0high(fi)
            read(13,'(a)') scr
            read(scr(19:50),*) nztmp
            read(13,'(a)') scr
            read(scr(19:50),*) ndtmp
c            read(13,502) scr3,nztmp
c            read(13,502) scr3,ndtmp
            if (nztmp.ne.zndb) then
               write(*,*) 'DS Error in dsanreadfiles.',
     &          ' Mismatch in z-bins: ',
     &              nztmp
               stop
            endif
            if (ndtmp.ne.dpn) then
               write(*,*) 'DS Error in dsanreadfiles.',
     &          ' Mismatch in dp-bins: ',
     &              dpn
               stop
            endif
            do i=1,12
              read(13,'(a)') scr ! header lines
            enddo
          else
            open(unit=13,file=filein,status='old',form='unformatted')
            read(13) yieldkfile
            read(13) dbct
            read(13) dbp0fit(fi)
            read(13) dbp0low(fi)
            read(13) dbp0high(fi)
          endif

          if (yieldk.ne.yieldkfile) then
            write(*,*) 'DS ERROR in dsanreadfiles: '
            write(*,*) 'the requested yield type is inconsistent with',
     &        ' the provided data file.'
            write(*,*) 'data file: ',filein
            write(*,*) 'requested yieldtype:   ',yieldk
            write(*,*) 'data file''s yieldtype: ',yieldkfile
            write(*,*) 'program stopped.'
            stop
          endif

          yieldtype(fltype,fi)=yieldkfile
          do j=1,nch
            if (prtlevel.gt.2) write(*,*) '      channel number ',j
            do k=1,nmass
              if (k.ge.milow(j)) then
                if (anftype.eq.'a'.or.anftype.eq.'a') then
                  read(13,*) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  do m=0,dpn-1
                     read(13,2000)
     &                 (phidiffdb(l,m,k,j,fi,dbkind),l=0,zndb-1)
                  enddo
                else
                  read(13) chfile,mfile
                  if (chfile.ne.j) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected channel ',j,
     &   ' but found channel ',chfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  if (abs(mfile-mi(k))/mi(k).gt.0.01d0) then
                     write(*,*) 
     &   'DS ERROR in dsanreadfiles. Inconsistent haloann file.'
                     write(*,*) 'Expected mass ',mi(k),
     &   ' but found mass ',mfile
                     write(*,*) 'Stopping.'
                     stop
                  endif
                  do m=0,dpn-1
                    read(13) (phidiffdb(l,m,k,j,fi,dbkind),l=0,zndb-1)
                  enddo
                endif
                do m=0,dpn-1
                  do l=0,zndb-1 ! correct units of dyield / dz
                     phidiffdb(l,m,k,j,fi,dbkind)=
     &               phidiffdb(l,m,k,j,fi,dbkind)/
     &                (dbdz(l))
                  enddo
                  phidiffdb(-1,m,k,j,fi,dbkind)=
     &              phidiffdb(0,m,k,j,fi,dbkind)
                  phidiffdb(zndb,m,k,j,fi,dbkind)=
     &              phidiffdb(zndb-1,m,k,j,fi,dbkind)
                enddo
              endif
            enddo
          enddo
          close(13)
        endif

      endif

      if (prtlevel.gt.1) 
     &     write(*,*) 'loading of halo yield tables finished.'

 2000 format(1000(1x,e12.6))

      return

      end
