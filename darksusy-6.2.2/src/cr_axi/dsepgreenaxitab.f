      real*8 function dsepgreenaxitab(R,z,DeltaVin,powerin,how,
     &  labhalo,loadlab)
************************************************************************
*** tabulated version of the green function dsepgreenaxi(R,z,DeltaV).
*** on first call the table is either computed or read from a file 
*** provided by the user. this is handled through the input variable:
***     how = 2 - table is tabulated on first  call, and then 
***               interpolated
***           3 - as 2, but also write the table to disk at the 
***               first call
***           4 - read table from disk on first call, and use the 
***               subsequent calls. If the file does not exist, it 
***               will be created (as in 3). (use as default)
*** you need to load a new table whenever:
***   - the halo model is changed
***   - the vertical height of the diffusion zone or the inner radial 
***     cutoff are changed
***   - R or z are changed
***     if epgraxiloadck=.true. there is an automatic check on these
***     and tables are eventually reloaded
***
*** other inputs:
***     R - radial position of the observer in kpc
***     z - vertical position of the observer in kpc
***     DeltaVin - increment in the variable vvar in kpc^2       
***     4*DeltaVin = lambda**2, with lambda the diffusion length in kpc 
***     powerin = integer selecting annihilations (=2) or decays (=1)
***     labhalo = halo model access label, within the set of models 
***             defined in the halo model repository  
***     loadlab = logical variable: in case it is set to true the halo
***       labhalo is selected within this function, otherwise it is
***       assumed that it has been loaded before linking to this 
***       function such as in the function dsepdphidpaxi or dsepdndpaxi
***
*** output: dimensionless
************************************************************************
      implicit none
      include 'dscraxicom.h' ! for crproptag
      include 'dsver.h'  ! for dsversion
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      include 'dslabelcom.h'
      include 'dsio.h'
ccc
      real*8 R,z,deltaVin
      integer powerin,how
      character(*) labhalo
      logical loadlab
ccc
      integer noout,dhow,ii,dsepcraxinumber
      real*8 deltav,result,dsepgreenaxi
      real*8 xmax,xmin,xint,ymax,ytoll,yfrac
      character*200 epfile,scr
      character*1 powerlab
      integer i,k,kk,npoints
      real*8 diff,sum,parread(100)
      integer powerread
      real*8 dsdmssetlabel,dscraxirzsetlabel,dsepcraxisetlabel
      external dsdmssetlabel,dscraxirzsetlabel,dsepcraxisetlabel
      logical loadlabloc
ccc
      real*8 xdv(1000),ydv(1000),ydv2(1000),realndv
      integer powerstore,dmdihalotagep,nocraxirz,noepcraxi,iload
      character*10 dmslabelt
      common/epgreenaxitabcom/xdv,ydv,ydv2,realndv,powerstore,
     &     dmdihalotagep,nocraxirz,noepcraxi,iload,dmslabelt
      save /epgreenaxitabcom/
ccc
      dhow=how
      if(dhow.lt.2.or.dhow.gt.4) then
        write(*,*) 'DS: ERROR in dsepgreenaxitab: the requested how = ',
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
          write(*,*) 'DS: call to dsepgreenaxitab with the halo label: ',
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
ccc check the match with the source function:
ccc

      if(dmdihalotagep.ne.dmdihalotag) then
        dmslabelt=dmdlabel
        dmdihalotagep=dmdihalotag
        iload=0
      endif  
ccc
ccc check the match with the propagation parameters. NOTE: this works
ccc only if you properly set the integer dsepcraxino via dsepcraxick      
ccc

      noout=dsepcraxinumber()
      if(noout.le.0) then   ! zero or negative value corresponds
                            ! to no initialization of the default model
        call dsepcraxick
        noout=dsepcraxinumber()
      endif
      if(noepcraxi.ne.noout) then
        call dsgetlabel2(dsepcraxisetlabel,epcraxitag)
        noepcraxi=noout
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

      if((dhow.eq.3.or.dhow.eq.4).and.iload.ne.1) then
        write(powerlab,'(I1)') powerin
ccc generate file name
        call dsdatafile(epfile,'epgretab-')
        call dscharadd(epfile,rzaxitag)
        call dscharadd(epfile,'-')
        call dscharadd(epfile,epcraxitag)
        call dscharadd(epfile,powerlab)
        call dscharadd(epfile,'-')
        call dscharadd(epfile,dmslabelt)
        call dscharadd(epfile,'.dat')
      endif
ccc
 5    if((dhow.eq.2.or.dhow.eq.3).and.iload.ne.1) then 
ccc
ccc make a new table:
ccc
        write(*,*) 'DS: starting dsepgreenaxi tabulation'  
        ymax=0.d0
        xmin=0.1d0
        xmax=10000.d0
 110    continue
        ii=0
ccc
ccc start with 100 points between xmin and xmax, eventually restarting
ccc if xmin is too large, and stopping if xmax is too large:
ccc
        do i=0,100
          deltav=dexp(dlog(xmin)+(dlog(xmax)-dlog(xmin))/100.d0*i)
          result=dsepgreenaxi(R,z,DeltaV,powerin,labhalo,loadlabloc)
c          write(*,*) 'deltav, result = ',deltav,result
          if(i.eq.0.and.dabs(1.d0-result).gt.1.d-2) then
            xmin=xmin/10.d0
            goto 110
          endif
          ii=ii+1
          xdv(ii)=dlog(deltav)
          ydv(ii)=result
          if(ydv(ii).gt.ymax) ymax=ydv(ii)
c          write(*,*) 'ii,deltav, result = ',ii,deltav,result
          if(ydv(ii).lt.1.d-8) goto 111
        enddo
 111    npoints=ii
c        write(*,*) 'Number of points: ',npoints
        ytoll=0.1d0
        yfrac=1.d-4
 130    continue
        do k=1,npoints-1
          if(dabs(ydv(k)-ydv(k+1)).gt.ytoll*max(ydv(k),ydv(k+1)).and.
     &        max(ydv(k),ydv(k+1)).gt.yfrac*ymax) then
            xint=(xdv(k)+xdv(k+1))/2.d0
            deltav=dexp(xint)
            result=dsepgreenaxi(R,z,DeltaV,powerin,labhalo,loadlabloc)
            do kk=npoints,k+1,-1
              xdv(kk+1)=xdv(kk)
              ydv(kk+1)=ydv(kk)
            enddo
            xdv(k+1)=xint
            ydv(k+1)=result
c            write(*,*) 'adding deltav, result = ',deltav,result
            npoints=npoints+1  
            if(npoints.le.1000) then
              goto 130 
            else
              write(*,*) 'DS: in dsepgreenaxitab exceeded the maximum'
              write(*,*) 'DS: dimension allowed for vectors in the'
              write(*,*) 'DS: table common which is set equal to 1000'
              write(*,*) 'DS: program stopped'
              stop
            endif
          endif  
        enddo
        write(*,*) 'DS: dsepgreenaxi tabulation is over'  
      elseif(dhow.eq.4.and.iload.ne.1) then  
ccc
ccc read table from file
ccc
        if (prtlevel.ge.2) write(*,*) 'DS: Reading ep table from file ',epfile
        ii=0
        open(unit=13,file=epfile,status='old',form='formatted',
     &       err=200)
        read(13,'(a)',err=200,end=200) scr  ! read header line
        read(13,'(a)',err=200,end=200) scr  ! read header line
        call dsreadlabel2(dscraxirzsetlabel,rzaxitag)
        read(13,1004,err=200,end=200) (parread(i),i=1,npar),powerread
        do i=1,npar
          diff=dabs(par(i)-parread(i))
          sum=dabs(par(i)+parread(i))
          if(diff.gt.1.d-7*sum) goto 201 
        enddo
        if(powerstore.ne.powerread) goto 202
        read(13,'(a)',err=200,end=200) scr  ! read header line
        call dsreadlabel2(dsepcraxisetlabel,epcraxitag)
        read(13,1000,err=200,end=200) (parread(i),i=1,npar)
        do i=1,npar
          diff=dabs(par(i)-parread(i))
          sum=dabs(par(i)+parread(i))
          if(diff.gt.1.d-7*sum) goto 204 
        enddo
        read(13,'(a)',err=200,end=200) scr  ! read header line
 10     ii=ii+1
        if(ii.gt.1000) then
          write(*,*) 'DS: in dsepgreenaxitab exceeded the maximum'
          write(*,*) 'DS: dimension allowed for vectors in the'
          write(*,*) 'DS: table common which is set equal to 1000'
          write(*,*) 'DS: program stopped'
          stop
        endif  
        read(13,1000,end=11)  xdv(ii),ydv(ii) 
        goto 10
 1000   format(60(1x,e14.8))
 11     close(13)
        npoints=ii-1  ! for later spline setup
        if (prtlevel.ge.3) write(*,*) 'DS: Done.'
      endif
ccc
      if(dhow.eq.3.and.iload.ne.1) then 
ccc
ccc write table to disk
ccc
        write(*,*) 'DS: Writing ep table to file ',epfile
        open(unit=13,file=epfile,status='unknown',form='formatted')
        write(13,1001) dsversion
        write(13,1003)
        call dsreadlabel2(dscraxirzsetlabel,rzaxitag)
        write(13,1004) (par(i),i=1,npar),powerstore
        write(13,1005)
        call dsreadlabel2(dsepcraxisetlabel,epcraxitag)
        write(13,1000) (par(i),i=1,npar)
        write(13,1002) 
 1001   format('# Made with DarkSUSY version ',A)
 1002   format('#','..log(deltav).',1x,'...table(i)...')
 1003   format('#','..Rstore......',1x,'..zstore......',1x,'.pw.')
 1004   format(2(1x,e14.8),1x,I4)
 1005   format('#','..diffhhstore.',1x,'.diffrcepstore')
        do i=1,npoints
          write(13,1000)  xdv(i),ydv(i) 
        enddo  
        close(13)
        if (prtlevel.ge.3) write(*,*) 'DS: Done.'
      endif
ccc
      if(iload.ne.1) then
ccc now set up splines
        realndv=dble(npoints)
        call dsspline(xdv,ydv,int(realndv),1.d31,1.d31,ydv2)
        iload=1  ! do not reload on next call
      endif
c...Now do the actual table lookup
      if(deltaVin.lt.dexp(xdv(1)).and.deltaVin.gt.-1.d-16) then
        if(dabs(1.d0-ydv(1)).lt.1.d-2) then 
          dsepgreenaxitab=1.d0
        else  
          dsepgreenaxitab=ydv(1)
        endif
      elseif(deltaVin.gt.dexp(xdv(int(realndv))).or.deltaVin.lt.0.d0) 
     &  then
        dsepgreenaxitab=0.d0
      else   
        call dssplint(xdv,ydv,ydv2,int(realndv),dlog(deltaVin),result)
        dsepgreenaxitab=result
      endif
      return
ccc
ccc we get here if one among: R, z diffhh and diffrcep were changed
ccc
 201  close(13)
      write(*,*) 'DS: link to a dsepgreenaxi table file ',epfile
      write(*,*) 'DS: with a mismatch in:'
      write(*,*) 'DS: R, z : ',(par(i),i=1,npar),(parread(i),i=1,npar)
      stop      
 202  close(13)
      write(*,*) 'DS: link to a dsepgreenaxi table file ',epfile
      write(*,*) 'DS: with a mismatch in:'
      write(*,*) 'DS: is : ',powerstore,powerread
      stop
 204  close(13)
      write(*,*) 'DS: link to a dsepgreenaxi table file ',epfile
      write(*,*) 'DS: with a mismatch in diffusion parameters :'
      write(*,*) 'DS: diffhh or diffrcep : ',(par(i),i=1,npar),
     &  (parread(i),i=1,npar)
      stop
ccc
ccc we get here if the file does not exist or it is not in the right
ccc format
ccc
 200  close(13)
      write(*,*) 'DS: The requested dsepgreenaxi table file ',epfile
      write(*,*) 'DS: does not exist, or it is not in the right format'
      write(*,*) 'DS: it will be (re)created if for you.'
      dhow=3
      goto 5
      end

