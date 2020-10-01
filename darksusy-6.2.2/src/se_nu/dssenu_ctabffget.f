      real*8 function dssenu_ctabffget(wh,ct,mx)

***********************************************************************
*** Interpolates in capture rate tables and returns the capture
*** rate (apart from the coupling).
*** This routine is for the full form factor integration.      
*** Input: wh ('su' or 'ea' for sun or earth)
***        ct, copuling type:
***           1 = gps**2      
***           2 = gns**2      
***           3 = gps*gns      
***           4 = gpa**2      
***           5 = gna**2      
***           6 = gpa*gna      
***        mx WIMP mass in GeV
*** Hidden input: velocity distribution model as given in
*** veldf (for the Sun) and veldfearth (for the Earth)
*** and possible escape velocity corrections (Jupiter effects)
*** Author: Joakim Edsjo
*** Date: 2015-06-12
***********************************************************************

      implicit none

      include 'dssecap.h'
      include 'dshmcom.h'
      include 'dssecom.h'

      integer i,ct,mxi,index
      character*2 wh
      character*200 file
      character*10 vec
      character*1 acc
      real*8 mx,mxpl,tmp,tmp1,tmp2
      logical newfile

c...Determine which table is needed based on veldf or veldfearth
c...Generate a file name
      call dsdatafile(file,'ctabff-')
      write(vec,'(F10.2)') veout
      write(acc,'(I1)') sesunacc
      if (wh.eq.'su'.or.wh.eq.'SU') then
        if (veldf.eq.'num'.or.veldf.eq.'numc') then
           call
     &    dscharadd(file,'su-num-'//acc//'-'//vec//'-'//haloid//'.dat')
        else
          call dscharadd(file,'su-'//acc//'-'//vec//'-'//veldf//'.dat')
        endif
      else
         write(*,*) 'DS ERROR in dssenu_ctabffget:'
         write(*,*)
     &    'Earth tables with full FF integration not implemented yet.'
         write(*,*) 'Program stopping.'
         stop
c        call dscharadd(file,'ea-'//veldfearth//'.dat')
      endif

      newfile=.true.
      if (wh.eq.'su'.or.wh.eq.'SU') then
        do i=1,nsuffloaded
          if (file.eq.filesuff(i)) then
            index=i
            newfile=.false.
            goto 10
          endif
        enddo
        nsuffloaded=nsuffloaded+1  ! one more file loaded
        if (nsuffloaded.gt.ntsuff) then
          write(*,*)
     &  'DS WARNING in dssenu_ctabffget: Too many files loaded'
          write(*,*) 'DarkSUSY will still work, but your performance ',
     &      'will not be optimal.'
          write(*,*) 'Increase ntsuff in dssecap.h to fix this problem.'
          nsuffloaded=ntsuff
        endif
        index=nsuffloaded
        filesuff(index)=file
      else
c        do i=1,neaffloaded
c          if (file.eq.fileeaff(i)) then
c            index=i
c            newfile=.false.
c            goto 10
c          endif
c        enddo
c        neaffloaded=neaffloaded+1  ! one more file loaded
c        if (neaffloaded.gt.nteaff) then
c          write(*,*) 'WARNING in dssenu_ctabffget: Too many files loaded'
c          write(*,*) 'DarkSUSY will still work, but your performance ',
c     &      'will not be optimal.'
c          write(*,*) 'Increase nteaff in dssecap.h to fix this problem.'
c          neaffloaded=nteaff
c        endif
c        index=neaffloaded
c        fileeaff(index)=file
      endif

 10   continue

c...if newfile=.false., we have already loaded this file with index index
c...if new=.true., we need to reload the file

c...Load tables if needed
      if (newfile) then
        call dssenu_ctabffread(wh,index,file)
      endif

c...Find entry
      if (mx.lt.1.0d0.or.mx.ge.1.0d5) then
        write(*,*) 'DS WARNING from dssenu_ctabffget: ',
     &    'WIMP mass outside allowed range: ',mx
        dssenu_ctabffget=0.0d0
      endif

      tmp=log10(mx)*dble(ncff)/5.0d0
      mxi=int(tmp)
      mxpl=tmp-mxi

      if (wh.eq.'su'.or.wh.eq.'SU') then

c...Use log interpolation on capture rate (if positive, otherwise linear)
c...Log gives smaller interpolation errors, but entry 3 and 6 could be
c...negative. If it is, use linear for that one.
         tmp1=ctabffsu(mxi,index,ct)
         tmp2=ctabffsu(mxi+1,index,ct)
         if (tmp1.lt.0.d0.or.tmp2.lt.0.d0) then
           dssenu_ctabffget=(1.0d0-mxpl)*tmp1
     &        +mxpl*tmp2
         else
           dssenu_ctabffget=exp((1.0d0-mxpl)*log(tmp1)
     &       +mxpl*log(tmp2))
         endif
         
      else
c...Earth not implemented with full FF integration yet
         dssenu_ctabffget=0.d0
      endif

      return

      end


      

