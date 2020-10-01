      real*8 function dssenu_ctabget(wh,st,mx)

***********************************************************************
*** Interpolates in capture rate tables and returns the capture
*** rate (apart from the cross section)
*** Input: wh ('su' or 'ea' for sun or earth)
***        st, spin-type (1:spin-independent, 2:spin-dependent)
***        mx WIMP mass in GeV
*** Hidden input: velocity distribution model as given in
*** veldf (for the Sun) and veldfearth (for the Earth)
*** and possible escape velocity corrections (Jupiter effects)
*** Author: Joakim Edsjo
*** Date: 2003-11-27
***********************************************************************

      implicit none

      include 'dssecap.h'
      include 'dshmcom.h'
      include 'dssecom.h'

      integer i,st,mxi,index
      character*2 wh
      character*200 file
      character*10 vec
      real*8 mx,mxpl,tmp
      logical newfile

c...Determine which table is needed based on veldf or veldfearth
c...Generate a file name
      call dsdatafile(file,'ctab-')
      write(vec,'(F10.2)') veout
      if (wh.eq.'su'.or.wh.eq.'SU') then
        if (veldf.eq.'num'.or.veldf.eq.'numc') then
          call dscharadd(file,'su-num-'//vec//'-'//haloid//'.dat')
        else
          call dscharadd(file,'su-'//vec//'-'//veldf//'.dat')
        endif
      else
        call dscharadd(file,'ea-'//veldfearth//'.dat')
      endif

      newfile=.true.
      if (wh.eq.'su'.or.wh.eq.'SU') then
        do i=1,nsuloaded
          if (file.eq.filesu(i)) then
            index=i
            newfile=.false.
            goto 10
          endif
        enddo
        nsuloaded=nsuloaded+1  ! one more file loaded
        if (nsuloaded.gt.ntsu) then
          write(*,*)
     &   'DS WARNING in dssenu_ctabget: Too many files loaded'
          write(*,*) 'DarkSUSY will still work, but your performance ',
     &      'will not be optimal.'
          write(*,*) 'Increase ntsu in dssecap.h to fix this problem.'
          nsuloaded=ntsu
        endif
        index=nsuloaded
        filesu(index)=file
      else
        do i=1,nealoaded
          if (file.eq.fileea(i)) then
            index=i
            newfile=.false.
            goto 10
          endif
        enddo
        nealoaded=nealoaded+1  ! one more file loaded
        if (nealoaded.gt.ntea) then
          write(*,*) 'WARNING in dssenu_ctabget: Too many files loaded'
          write(*,*) 'DarkSUSY will still work, but you performance ',
     &      'will not be optimal.'
          write(*,*) 'Increase ntea in dssecap.h to fix this problem.'
          nealoaded=ntea
        endif
        index=nealoaded
        fileea(index)=file
      endif

 10   continue

c...if newfile=.false., we have already loaded this file with index index
c...if new=.true., we need to reload the file

c...Load tables if needed
      if (newfile) then
        call dssenu_ctabread(wh,index,file)
      endif

c...Find entry
      if (mx.lt.1.0d0.or.mx.ge.1.0d5) then
        write(*,*) 'WARNING from dssenu_ctabget: ',
     &    'WIMP mass outside allowed range: ',mx
        dssenu_ctabget=0.0d0
      endif

      tmp=log10(mx)*dble(nc)/5.0d0
      mxi=int(tmp)
      mxpl=tmp-mxi

      if (wh.eq.'su'.or.wh.eq.'SU') then

        if (st.eq.1) then ! spin-independent
c          dssenu_ctabget=(1.0d0-mxpl)*ctabsusi(mxi,index)
c     &      +mxpl*ctabsusi(mxi+1,index)
          dssenu_ctabget=exp((1.0d0-mxpl)*log(ctabsusi(mxi,index))
     &      +mxpl*log(ctabsusi(mxi+1,index)))
        else  ! spin-dependent
c          dssenu_ctabget=(1.0d0-mxpl)*ctabsusd(mxi,index)
c     &      +mxpl*ctabsusd(mxi+1,index)
          dssenu_ctabget=exp((1.0d0-mxpl)*log(ctabsusd(mxi,index))
     &      +mxpl*log(ctabsusd(mxi+1,index)))
        endif

      else

        if (st.eq.1) then ! spin-independent
c          dssenu_ctabget=(1.0d0-mxpl)*ctabea(mxi,index)
c     &      +mxpl*ctabea(mxi+1,index)
          dssenu_ctabget=exp((1.0d0-mxpl)*log(ctabea(mxi,index))
     &      +mxpl*log(ctabea(mxi+1,index)))
        else ! spin-dependent
          dssenu_ctabget=0.0d0
        endif

      endif

      return

      end


      

