      real*8 function dsdbtdaxitab(R,z,tpin,powerin,isin,how,labhalo
     &   ,loadlab)
************************************************************************
*** tabulated version of the db "confinement time" dsdbtdaxi
***
*** input variables:
***   R = radial coordinate for the position at which the flux is
***      measured in kpc (cylindrical coordinate system)
***   z = vertical coordinate for the position at which the flux is
***      measured in kpc (cylindrical coordinate system)
***      NB: only z=0 works at the moment
***   tpin = antideuteron kinetic energy per nucleon in GeV
***   powerin = integer selecting annihilations (=2) or decays (=1)
***   isin = propagation model option : isin=1 -> model with delta 
***      function approximation for the gas disc; isin=2 -> two-zone
***      model with finite thickness for the gas disk
***   how = 1 - calculate the time only for requested tpin
***         2 - table is tabulated on first call, and then 
***             interpolated
***         3 - as 2, but also write the table to disk at the 
***             first call
***         4 - read table from disk on first call, and use the 
***             subsequent calls. If the file does not exist, it 
***             will be created (as in 3). (use as default)
***     the tabulation is created with dsdbtdaxi with nprec and prec as
***     set by the corresponding pbprec and pbnprec in common blocks. 
***     you need to load a new table whenever:
***       - the halo model is changed
***       - the propagation model is changed
***       - R, z or isin are changed
***     if pbtdaxiloadck=.true. there is an automatic check on these
***     and tables are eventually reloaded
***   labhalo = halo model access label, within the set of models 
***             defined in the halo model repository
***   loadlab = logical variable: in case it is set to true the halo
***     labhalo is selected within this function, otherwise it is
***     assumed that it has been loaded before linking to this function
***     such as in the function dspbdphidtaxi      
***
*** output: 10^15 s
************************************************************************
      implicit none
      include 'dscraxicom.h' ! for 
      include 'dsver.h'  ! for dsversion
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      include 'dslabelcom.h'
      include 'dsio.h'
      real*8 R,z,tpin
      integer powerin,isin,how
      character(*) labhalo
      logical loadlab
ccc
      integer dhow,ii,noout,dspbcraxinumber
      real*8 tp,result,dsdbtdaxi
      real*8 xmax,xmin,xint,ymax,ytoll,yfrac
      character*200 dbfile,scr
      character*1 islab,powerlab
      integer i,k,kk,npoints
      real*8 diff,sum,parread(100)
      integer powerread,isread
      real*8 dsdmssetlabel,dscraxirzsetlabel,dspbcraxisetlabel
      external dsdmssetlabel,dscraxirzsetlabel,dspbcraxisetlabel
      logical loadlabloc
      real*8 parvec(6)
      logical check,match
ccc
      real*8 xtd(1000),ytd(1000),ytd2(1000),realntd
      integer isstore,powerstore,dmdihalotagdb,nocraxirz,nodbcraxi,iload
      character*10 dmslabelt
      common/dbtdaxitabcom/xtd,ytd,ytd2,realntd,isstore,powerstore,
     &  dmdihalotagdb,nocraxirz,nodbcraxi,iload,dmslabelt
ccc
      dhow=how
      if(dhow.lt.1.or.dhow.gt.4) then
        write(*,*) 'DS: ERROR in dsdbtdaxitab: the requested how = ',
     &     dhow,' is out of range. Stopping program.'
        stop
      endif
ccc
      if(loadlab) then
        call dsdmdselect_halomodel(labhalo)
ccc
ccc check whether you have called this function for a halo model which
ccc can be used for Milky Way rates or not
ccc      
        if(.not.dmdmw) then
          write(*,*) 'DS: call to dspbtdaxitab with the halo label: ',
     &       labhalo    
          write(*,*) 'DS: which has not been initialized as suitable'
          write(*,*) 'DS: for Milky Way rates. Program stopped'
          stop
        endif
      endif
      loadlabloc=.false.
ccc
      if((dhow.eq.3.or.dhow.eq.4).and.dmdtmp) then
ccc
ccc a temporary halo model cannot be used for an input/output on file
ccc      
        dhow=2 ! create a new table but do not save it to a file
      endif
ccc
      if(dhow.eq.1) then
ccc
ccc check whether this dspbtdaxi has already been computed:
ccc      
        check=.true.
        parvec(2)=R
        parvec(3)=z
        parvec(4)=tpin
        parvec(5)=dble(powerin)
        parvec(6)=dble(isin)
        call dsdbtdaxidb(dmdihalotag,parvec,dmdtmp,check,match)
        if(match) then
          dsdbtdaxitab=parvec(1)
          return
        endif
ccc
ccc link directly to dsdbtdaxi
ccc
c        write(*,*) 'computing db : ',dmdihalotag,tpin
        dsdbtdaxitab=dsdbtdaxi(R,z,tpin,powerin,isin,pbprec,pbnprec
     &   ,labhalo,loadlabloc)
ccc
ccc store the value of dsdbtdaxi:
ccc      
        check=.false.
        parvec(1)=dsdbtdaxitab
        call dsdbtdaxidb(dmdihalotag,parvec,dmdtmp,check,match)
        return
      endif
ccc
ccc check the match with the position in the galaxy
ccc
      call dscraxirzck(R,z,noout)
      if(nocraxirz.ne.noout) then
ccc
ccc get rzaxitag
ccc
        call dsgetlabel2(dscraxirzsetlabel,rzaxitag)
        nocraxirz=noout
        iload=0
      endif
ccc
ccc check the match with the source function. NOTE: this works only
ccc if you corretly incrementing dsdmsno at each change in the source
ccc function
ccc
      if(dmdihalotagdb.ne.dmdihalotag) then
        dmslabelt=dmdlabel
        dmdihalotagdb=dmdihalotag
        iload=0
      endif
ccc
ccc check the match with the propagation parameters
ccc
      noout=dspbcraxinumber()      
      if(noout.le.0) then   ! zero or negative value corresponds
                            ! to no initialization of the default model
        call dspbcraxick
        noout=dspbcraxinumber()
      endif
      if(nodbcraxi.ne.noout) then
        call dsgetlabel2(dspbcraxisetlabel,pbcraxitag)
        nodbcraxi=noout
        iload=0
      endif
ccc
ccc check the match with the "power" option
ccc
      if(powerstore.ne.powerin) then
        powerstore=powerin
        iload=0
      endif
ccc
ccc check the match with the "is" option
ccc
      if(isstore.ne.isin) then
        isstore=isin
        iload=0
      endif
ccc
      if((dhow.eq.3.or.dhow.eq.4).and.iload.ne.1) then
        write(powerlab,'(I1)') powerin
        write(islab,'(I1)') isin
ccc generate file name and store variables
        call dsdatafile(dbfile,'dbtdaxitab-')
        call dscharadd(dbfile,rzaxitag)
        call dscharadd(dbfile,'-')
        call dscharadd(dbfile,pbcraxitag)
        call dscharadd(dbfile,powerlab)
        call dscharadd(dbfile,islab)
        call dscharadd(dbfile,'-')
        call dscharadd(dbfile,dmslabelt)
        call dscharadd(dbfile,'.dat')
      endif
ccc
 5    if((dhow.eq.2.or.dhow.eq.3).and.iload.ne.1) then 
ccc
ccc make a new table:
ccc
        write(*,*) 'DS: starting dsdbtdaxi tabulation'  
        ymax=0.d0
        xmin=0.1d0
        xmax=10000.d0
        ii=0
ccc
ccc start with 100 points between xmin and xmax, cutting eventually
ccc if the result is too small:
ccc
        do i=0,100
          tp=dexp(dlog(xmin)+(dlog(xmax)-dlog(xmin))/100.d0*i)
          result=dsdbtdaxi(R,z,tp,powerin,isin,pbprec,pbnprec,labhalo
     &      ,loadlabloc)
          ii=ii+1
          xtd(ii)=dlog(tp)
          ytd(ii)=result
          if(ytd(ii).gt.ymax) ymax=ytd(ii)
          if(ytd(ii).lt.1.d-8) goto 111
        enddo
 111    npoints=ii
c        write(*,*) 'DS: Number of points: ',npoints
 130    continue
ccc
ccc integrate the table with extra points if the spacing in the vertical
ccc coordinate is larger than (ytoll*100) %, in the region were the 
ccc vertical coordinate is at least yfrac*ymax   
ccc
        ytoll=0.2d0
        yfrac=1.d-3
        do k=1,npoints-1
          if(dabs(ytd(k)-ytd(k+1)).gt.ytoll*max(ytd(k),ytd(k+1))
     &       .and.max(ytd(k),ytd(k+1)).gt.yfrac*ymax) then
            xint=(xtd(k)+xtd(k+1))/2.d0
            tp=dexp(xint)
            result=dsdbtdaxi(R,z,tp,powerin,isin,pbprec,pbnprec,labhalo
     &      ,loadlabloc)
            do kk=npoints,k+1,-1
              xtd(kk+1)=xtd(kk)
              ytd(kk+1)=ytd(kk)
            enddo
            xtd(k+1)=xint
            ytd(k+1)=result
            write(*,*) 'DS: adding tp, result = ',tp,result
            npoints=npoints+1  
            if(npoints.le.1000) then
              goto 130 
            else
              write(*,*) 'DS: in dsdbtdaxitab exceeded the maximum'
              write(*,*) 'DS: dimension allowed for vectors in the'
              write(*,*) 'DS: table common which is set equal to 1000'
              write(*,*) 'DS: program stopped'
              stop
            endif
          endif  
        enddo
        write(*,*) 'DS: dsdbtdaxi tabulation is over'  
      elseif(dhow.eq.4.and.iload.ne.1) then  ! read table from file
        if (prtlevel.ge.2) write(*,*) 'DS: Reading db table from file ',dbfile
        ii=1 
        open(unit=13,file=dbfile,status='old',form='formatted',
     &       err=200)
        read(13,'(a)',err=200,end=200) scr  ! read header line
        read(13,'(a)',err=200,end=200) scr  ! read header line
        call dsreadlabel2(dscraxirzsetlabel,rzaxitag)
        read(13,1004,err=200,end=200) (parread(i),i=1,npar),powerread,
     &    isread
        do i=1,npar
          diff=dabs(par(i)-parread(i))
          sum=dabs(par(i)+parread(i))
          if(diff.gt.1.d-7*sum) goto 201 
        enddo
        if(powerstore.ne.powerread) goto 202
        if(isstore.ne.isread) goto 203
        read(13,'(a)',err=200,end=200) scr ! read header line
        call dsreadlabel2(dspbcraxisetlabel,pbcraxitag)
        read(13,1000,err=200,end=200) (parread(i),i=1,npar)
        do i=1,npar
          diff=dabs(par(i)-parread(i))
          sum=dabs(par(i)+parread(i))
          if(diff.gt.1.d-7*sum) goto 204 
        enddo
        read(13,'(a)',err=200,end=200) scr  ! read header line
        read(13,1000,end=11)  xtd(ii),ytd(ii) 
 10     ii=ii+1
        if(ii.gt.1000) then
          write(*,*) 'DS: in dsdbtdaxitab exceeded the maximum'
          write(*,*) 'DS: dimension allowed for vectors in the'
          write(*,*) 'DS: table common which is set equal to 1000'
          write(*,*) 'DS: program stopped'
          stop
        endif  
        read(13,1000,end=11)  xtd(ii),ytd(ii) 
        goto 10
 1000   format(100(1x,e14.8))
 11     close(13)
        npoints=ii-1  ! for later spline setup
        if (prtlevel.ge.3) write(*,*) 'DS: Done.'
      endif
ccc
      if (dhow.eq.3) then ! write table to disk
        write(*,*) 'DS: Writing db table to file ',dbfile
        open(unit=13,file=dbfile,status='unknown',form='formatted')
        write(13,1001) dsversion
        write(13,1003)
        call dsreadlabel2(dscraxirzsetlabel,rzaxitag)
        write(13,1004) (par(i),i=1,npar),powerstore,isstore
        write(13,1005)
        call dsreadlabel2(dspbcraxisetlabel,pbcraxitag)
        write(13,1000) (par(i),i=1,npar)
        write(13,1002) 
 1001   format('# Made with DarkSUSY version ',A)
 1002   format('#','..log(tp).',1x,'...table(i)...')
 1003   format('#','..Rstore......',1x,'..zstore......',
     &         1x,'.pw.',1x,'.is.')
 1004   format(2(1x,e14.8),1x,I4,1x,I4)
 1005   format('#','..set..of..diffusion..parameters..')
        do i=1,npoints
          write(13,1000)  xtd(i),ytd(i) 
        enddo  
        close(13)
        if (prtlevel.ge.3) write(*,*) 'DS: Done.'
      endif
ccc
      if(iload.ne.1) then
ccc now set up splines
        realntd=dble(npoints)
        call dsspline(xtd,ytd,int(realntd),1.d31,1.d31,ytd2)
        iload=1  ! do not reload on next call
      endif
ccc now do the actual table lookup
      if(tpin.gt.dexp(xtd(int(realntd)))
     &  .or.tpin.lt.dexp(xtd(1))) then
        write(*,*) 'DS: dsdbtdaxitab called for tp = ',tpin
        write(*,*) 'DS: out of the tabulated range :',dexp(xtd(1)),
     &    dexp(xtd(int(realntd)))
        write(*,*) 'DS: dsdbtdaxitab set to zero '
        dsdbtdaxitab=0.d0
      else   
        call dssplint(xtd,ytd,ytd2,int(realntd),dlog(tpin),result)
        dsdbtdaxitab=result
      endif
ccc
      return
ccc
ccc we get here if one among: R, z, is or one of the diffusion parmaters
ccc were changed
ccc
 201  close(13)
      write(*,*) 'DS: link to a dbtdaxi table file ',dbfile
      write(*,*) 'DS: with a mismatch in:'
      write(*,*) 'DS: R, z : ',(par(i),i=1,npar),(parread(i),i=1,npar)
      stop      
 202  close(13)
      write(*,*) 'DS: link to a dbtdaxi table file ',dbfile
      write(*,*) 'DS: with a mismatch in:'
      write(*,*) 'DS: is : ',powerstore,powerread
      stop      
 203  close(13)
      write(*,*) 'DS: link to a dbtdaxi table file ',dbfile
      write(*,*) 'DS: with a mismatch in:'
      write(*,*) 'DS: is : ',isstore,isread
      stop      
 204  close(13)
      write(*,*) 'DS: link to a dbtdaxi table file ',dbfile
      write(*,*) 'DS: with a mismatch in diffusion parameters :'
      do i=1,npar
        write(*,*) 'DS: i parameter : ',i,par(i),parread(i)
      enddo
      stop
ccc
ccc we get here if the file does not exist or it is not in the right
ccc format
ccc
 200  close(13)
      write(*,*) 'DS: The requested dbtdaxi table file ',dbfile
      write(*,*) 'DS: does not exist, or it is not in the right format'
      write(*,*) 'DS: it will be (re)created for you.'
      dhow=3
      goto 5
ccc
      end


      subroutine dsdbtdaxidb(dmdindex,parvec,tmp,check,match)
ccc
ccc tmp.eq.true/false -> temporary halo model/halo model in database
ccc check.eq.true/false -> check match/add new result
ccc      
      implicit none
      include 'dslabelcom.h'
      integer dmdindex
      real*8 parvec(6)
      logical tmp,check,match
      integer nparvec
      parameter(nparvec=6)
ccc      
      integer ndatjmax,npartot
      parameter(ndatjmax=100,npartot=20)
ccc      
      integer ntmpdatj,ndatj,dmdihatmpdatj,dmdihadatj(ndatjmax),
     &  nparextra
      real*8 tmpdatj(npartot,ndatjmax),datj(npartot,ndatjmax),
     &  parextralc(npartot) 
      common/dbtdaxidbcom/tmpdatj,datj,parextralc,dmdihatmpdatj,
     &  dmdihadatj,ntmpdatj,ndatj,nparextra
ccc
      integer ii,jj
      real*8 dummy,dspbcraxisetlabel,dif,sum
      logical first
      data first/.true./
      save first
ccc
      if(first) then
	ndatj=0
	first=.false.
        dummy=dspbcraxisetlabel(labset) ! for npar
        nparextra=npar
        if(nparvec+nparextra.gt.npartot) then
          write(*,*) 'DS: error in dspbtdaxidb: size of common blocks'
          write(*,*) 'DS: not sufficient, since nparvec+nparextra = '
     &         ,nparvec+nparextra
          write(*,*) 'DS: while npartot = ',npartot
          write(*,*) 'DS: program stopped'
          stop
        endif  
      endif
ccc
      if(check) then
        dummy=dspbcraxisetlabel(labin) ! for par(1) -- par(npar)
        if(tmp) then
          if(dmdindex.ne.dmdihatmpdatj) then
            ntmpdatj=0
          else
            do ii=1,ntmpdatj
              match=.true.
              do jj=2,nparvec 
                dif=dabs(tmpdatj(jj,ii)-parvec(jj))
                sum=dabs(tmpdatj(jj,ii)+parvec(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo
              do jj=1,nparextra 
                dif=dabs(tmpdatj(jj+nparvec,ii)-par(jj))
                sum=dabs(tmpdatj(jj+nparvec,ii)+par(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo
              if(match) then
                parvec(1)=tmpdatj(1,ii)
                return
              endif
            enddo
          endif
          match=.false.
          return
        else
          do ii=1,ndatj
            if(dmdindex.eq.dmdihadatj(ii)) then 
              match=.true.
              do jj=2,nparvec
                dif=dabs(datj(jj,ii)-parvec(jj))
                sum=dabs(datj(jj,ii)+parvec(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo  
              do jj=1,nparextra 
                dif=dabs(datj(jj+nparvec,ii)-par(jj))
                sum=dabs(datj(jj+nparvec,ii)+par(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo
              if(match) then
                parvec(1)=datj(1,ii)
                return
              endif
            endif
          enddo
          match=.false.
          return
        endif
      else
        if(tmp) then
          ntmpdatj=ntmpdatj+1
          if(ntmpdatj.gt.ndatjmax) ntmpdatj=1
          dmdihatmpdatj=dmdindex
          do jj=1,nparvec 
            tmpdatj(jj,ntmpdatj)=parvec(jj)
          enddo
          do jj=1,nparextra 
            tmpdatj(jj+nparvec,ntmpdatj)=par(jj)
          enddo
        else
          ndatj=ndatj+1
          if(ndatj.gt.ndatjmax) ndatj=1
          dmdihadatj(ndatj)=dmdindex
          do jj=1,nparvec 
            datj(jj,ndatj)=parvec(jj)
          enddo
          do jj=1,nparextra 
            datj(jj+nparvec,ndatj)=par(jj)
          enddo
        endif  
      endif      
      return
      end
        
      
