      subroutine dsdmddriver_num(iwin,nrein,rein,nchin,chin,labin,reout)
c_______________________________________________________________________
c this is the sample driver to initialize, print or link corresponding
c functions in case of a dark matter density profile defined as the
c interpolation of a table radius versus density profile
c
c inputs:
c   - iwin: integer setting the action of the driver; it can be equal
c     to:
c       idmdcreate -> profile creation is completed here (first settings
c         are in the dsdmddriver subroutine)      
c       idmdparprint -> the current halo profile is printed 
c         depending on the input labin, see below      
c       idmddensph -> the driver links to the spherical density profile
c   - labin: external label of the current halo model: you got here in
c       case of labin containing the string 'num'; if it contains also
c       the string 'mw' the model is intepreted as referring to the
c       Milky Way
c   - rein(nrein) real*8 vector for different inputs according to
c      different iwin values       
c   - chin(nchin) character*10 vector for input parameter labels
c output:         
c   - reout: real*8 output for the different function linking as
c     specified by the input value iwin 
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      include 'dsdvcom.h'
      include 'dsmpconst.h'
      include 'dsio.h'
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
      real*8 dsdmdsetlabel
      external dsdmdsetlabel
ccc      
      integer ii,jj,nt,indexout
      real*8 r,radouttr,radintr,rho0
      logical mwlab
      character*10 label
      logical error,dsisnan
      real*8 xtmp,ytmp,xin,yout
ccc
      if(iwin.eq.idmdcreate.or.iwin.eq.idmdparprint) then
        idmdlab(1)=iradintr
        chdmdlab(1)=chradintr
        idmdlab(2)=iradouttr
        chdmdlab(2)=chradouttr
        idmdlab(3)=iobjdist
        chdmdlab(3)=chobjdist
        idmdlab(4)=itabnum
        chdmdlab(4)=chtabnum
        ndmdlabtot=4
      endif
ccc
      if(iwin.eq.idmdcreate) then
c
c load the first four parameters:
c        
        do ii=1,ndmdlabtot
          do jj=1,nchin
            if(chin(jj).eq.chdmdlab(ii)) then
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(idmdlab(ii),ihalocur)=indexout
              goto 100
            endif
          enddo
          write(*,*) 'DS: in dsdmddriver_num missing initialization'
          write(*,*) 'DS: of the halo parameter : ',chdmdlab(ii)
          write(*,*) 'DS: when setting the halo model : ',labin
          write(*,*) 'DS: program stopped '
          stop
 100      continue
        enddo
c
c save table kind:
c
        tabnumtype=int(dmdparvec(ihalotag(itabnum,ihalocur)))
c
c load the table:
c
        nentries=0
ccc
ccc check whether the number of tabulated points is specified:
ccc      
        call dsdmdnumlab(' ',0,label)
        do ii=1,nchin
          if(label.eq.chin(ii)) then
            nentries=int(rein(ii))
            goto 101
          endif
        enddo
ccc
ccc fill in the table:
ccc      
 101    jj=0
 20     jj=jj+1
        call dsdmdnumlab('x',jj,label)
        do ii=1,nchin
          if(label.eq.chin(ii)) then
            xtabn(jj)=rein(ii)
            goto 10
          endif
        enddo
        goto 30
 10     call dsdmdnumlab('y',jj,label)
        do ii=1,nchin
          if(label.eq.chin(ii)) then
            ytabn(jj)=rein(ii)
            goto 20
          endif
        enddo
        write(*,*) 'DS: problem in model creation with dsdmddriver_num:'
        write(*,*) 'DS: a x-y mismatch for entry = ',jj
        write(*,*) 'DS: program stopped'
        stop
 30     nt=jj-1
        if(nentries.gt.0.and.nentries.ne.nt) then
          write(*,*) 'DS: error in model creation in dsdmddriver_num:'
          write(*,*) 'DS: expecting a table with # of entries = ',
     &      nentries
          write(*,*) 'DS: you loaded instead ',nt
          write(*,*) 'DS: program stopped'
          stop
        endif
        nentries=nt
ccc
ccc check for nan or negative (too small) entries in case of log:
ccc      
        do ii=1,nentries
          error=.false. 
          if(dsisnan(xtabn(ii)).or.dsisnan(xtabn(ii))) error=.true.
          if((tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &      .and.(xtabn(ii).lt.1.d-60.or.ytabn(ii).lt.1.d-60))
     &      error=.true.
          if(error) then
            write(*,*) 'DS: error in model creation in ',
     &        'dsdmddriver_num:'
            write(*,*) 'DS: error in entry = ',ii
            write(*,*) 'DS: in log scale x, y =  ',xtabn(ii),ytabn(ii)
            write(*,*) 'DS: program stopped'
            stop
          endif
        enddo
ccc
ccc sort values of x:
ccc      
        do ii=1,nentries
          do jj=nentries-1,ii,-1
            if(xtabn(jj).gt.xtabn(jj+1)) then
              xtmp=xtabn(jj)
              xtabn(jj)=xtabn(jj+1)
              xtabn(jj+1)=xtmp
              ytmp=ytabn(jj)
              ytabn(jj)=ytabn(jj+1)
              ytabn(jj+1)=ytmp
            endif
          enddo
        enddo
ccc
ccc log for log tabulations:
ccc      
        if(tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &    then
          do ii=1,nentries
            xtabn(ii)=dlog(xtabn(ii)) 
            ytabn(ii)=dlog(ytabn(ii))
          enddo
        endif
ccc
ccc compute derivatives for spline tabulations:
ccc      
        if(tabnumtype.eq.tabnumspl.or.tabnumtype.eq.tabnumlogspl) 
     &     call dsspline(xtabn,ytabn,nentries,1.d31,1.d31,ytabn2)      
c
c load whether it is Milky Way or not, in case it is load rho0
c        
        mwlab=.false.
        if(ihalotag(iwhichobject,ihalocur).eq.imwsuff) mwlab=.true. 
        if(mwlab) then
          xin=dmdparvec(ihalotag(iobjdist,ihalocur))
          if(tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &       xin=dlog(xin)
          if(tabnumtype.eq.tabnumspl.or.tabnumtype.eq.tabnumlogspl) then
            call dssplint(xtabn,ytabn,ytabn2,nentries,xin,yout)
          elseif(tabnumtype.eq.tabnumlin.or.tabnumtype.eq.tabnumloglin)
     &      then
            call dslinint(xtabn,ytabn,nentries,xin,yout)
          else
            write(*,*) 'DS: error in dsdmddriver_num:'
            write(*,*) 'DS: tabnumtype = ',tabnumtype
            write(*,*) 'DS: instead of being in the range 1-4'
            write(*,*) 'DS: program stopped'
            stop
          endif
          if(tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &       yout=dexp(yout)
          rho0=yout  ! GeV cm^-3
          call dsdmdsetind(indexout,rho0)
          ihalotag(irho0,ihalocur)=indexout
        else
c 
c if this is not the Milky Way, rho0 is set to zero
c
          rho0=0.d0
          call dsdmdsetind(indexout,rho0)
          ihalotag(irho0,ihalocur)=indexout
        endif
c
c no tabulation for CR green functions:
c        
        halotabtag(ihalocur)='none'
        return
      endif
ccc
      if(iwin.eq.idmdparprint) then
c         if (prtlevel.gt.2) then
            write(*,*)
            write(*,*) 'DS: the halo model with access label: '
     &        ,halotag(ihalocur)
            write(*,*) 'DS: has a corresponding authomatically generated'
            write(*,*) 'DS: label for tabulation files equal to: '
     &        ,halotabtag(ihalocur)
            write(*,*) 'DS: the list of parameters for this model is: '
            do ii=1,ndmdlabtot
               write(*,*) 'DS: par & value : ',chdmdlab(ii)
     &                    ,dmdparvec(ihalotag(idmdlab(ii),ihalocur)) 
            enddo
            write(*,*) 'DS: this model is defined via a tabulation'
            if(tabnumtype.eq.tabnumspl)
     &         write(*,*) 'DS: radius vs density and spline interpolation'
            if(tabnumtype.eq.tabnumlin)
     &         write(*,*) 'DS: radius vs density and linear interpolation'
            if(tabnumtype.eq.tabnumlogspl)
     &         write(*,*) 'DS: log(r) vs log(rho) and spline interpolation'
            if(tabnumtype.eq.tabnumloglin)
     &         write(*,*) 'DS: log(r) vs log(rho) and linear interpolation'
            write(*,*) 'DS: list radius vs rho'       
            do ii=1,nentries
               write(*,*) ii,xtabn(ii),ytabn(ii)
            enddo   
c        endif
        return
      endif
ccc
      if(iwin.eq.idmddensph) then
c        
c link to spherical dark matter density profile with inner density
c cutoff and outer truncation radius:
c
        r=rein(idvrsph)
        radouttr=dmdparvec(ihalotag(iradouttr,ihalocur))
        if(r.gt.radouttr) then
          reout=0.d0
          return
        endif  
        radintr=dmdparvec(ihalotag(iradintr,ihalocur))
        r=max(r,radintr)
ccc    
        if((tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &     .and.r.lt.1.d-60) then
          write(*,*) 'DS: error in dsdmddriver_num: calling a log'
          write(*,*) 'DS: tabulation with radius = ',r
          write(*,*) 'DS: program stopped'
          stop
        endif
        xin=r
        if(tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &     xin=dlog(xin)
        if(tabnumtype.eq.tabnumspl.or.tabnumtype.eq.tabnumlogspl) then
          call dssplint(xtabn,ytabn,ytabn2,nentries,xin,yout)
        elseif(tabnumtype.eq.tabnumlin.or.tabnumtype.eq.tabnumloglin)
     &    then
          call dslinint(xtabn,ytabn,nentries,xin,yout)
        else
          write(*,*) 'DS: error in dsdmddriver_num:'
          write(*,*) 'DS: tabnumtype = ',tabnumtype
          write(*,*) 'DS: instead of being in the range 1-4'
          write(*,*) 'DS: program stopped'
          stop
        endif
        if(tabnumtype.eq.tabnumlogspl.or.tabnumtype.eq.tabnumloglin)
     &    yout=dexp(yout)
        reout=yout  ! GeV cm^-3
        return
      endif
ccc
      if(iwin.eq.idmdmassph) then
c        
c link to spherical dark matter mass profile with outer truncation
c radius:
c
        r=rein(idvrsph)
        radouttr=dmdparvec(ihalotag(iradouttr,ihalocur))
        write(*,*) 'DS: error in dsdmddriver_num:'
        write(*,*) 'DS: iwin = idmdmassph but the mass profile'
        write(*,*) 'DS: has not been initialized so far'
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc you called this function with an iwin which is not set, an error
ccc message is written and the program is stopped      
ccc      
      write(*,*) 'DS: linking to dsdmddriver_num with iwin = ',iwin
      write(*,*) 'DS: out of the defined range of values'
      write(*,*) 'DS: program stopped'
      stop
      end


      
      subroutine dsdmdnumlab(axis,n,label)
      implicit none
      include 'dsdmddrvrcom.h'
      character*1 axis
      character(*) label
      integer n
ccc      
      if(n.gt.ntabdim) then
        write(*,*) 'DS: calling dsdmdnumlab with n = ',n
        write(*,*) 'DS: larger that the maximum value = ',ntabdim
        write(*,*) 'DS: program stopped'
        stop
      endif
      if(n.eq.0) then
        label='ntab'
        return
      endif  
      if(axis.eq.'x') then
        label='x'
        if(ntabdim.eq.999) then
          write(n3label,1003) n
          call dscharadd(label,n3label)
        elseif(ntabdim.eq.9999) then
          write(n4label,1004) n
          call dscharadd(label,n4label)
        elseif(ntabdim.eq.99999) then
          write(n5label,1004) n
          call dscharadd(label,n5label)
        else
          write(*,*) 'DS: calling dsdmdnumlab without proper setting'
          write(*,*) 'DS: of ntabdim = ',ntabdim
          write(*,*) 'DS: rather than 999, 9999 or 99999'
          write(*,*) 'DS: program stopped'
          stop
        endif
      elseif(axis.eq.'y') then  
        label='y'
        if(ntabdim.eq.999) then
          write(n3label,1003) n
          call dscharadd(label,n3label)
        elseif(ntabdim.eq.9999) then
          write(n4label,1004) n
          call dscharadd(label,n4label)
        elseif(ntabdim.eq.99999) then
          write(n5label,1004) n
          call dscharadd(label,n5label)
        else
          write(*,*) 'DS: calling dsdmdnumlab without proper setting'
          write(*,*) 'DS: of ntabdim = ',ntabdim
          write(*,*) 'DS: rather than 999, 9999 or 99999'
          write(*,*) 'DS: program stopped'
          stop
        endif
      else
        write(*,*) 'DS: calling dsdmdnumlab with axis = ',axis
        write(*,*) 'DS: different from x or y'
        write(*,*) 'DS: program stopped'
        stop
      endif  
 1003 format(i3)
 1004 format(i4)
c 1005 format(i5)
      return
      end
        
        
