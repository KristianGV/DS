      subroutine dstausetsph(xxmin,xxmax,nstep,ilosi)
************************************************************************
*** subroutine setting the tabulations of dstaufunsph. it needs to be 
*** called with the correct sequence for the variable nstep, within the
*** functions dslosisph and dsomlosisph 
************************************************************************
      implicit none
      include 'dsdmdintcom.h' ! npairs_norm,taucut,lochow
      include 'dslosidrvrcom.h' ! for ilositauset
      include 'dsdvcom.h' ! for dependent coordinates
      real*8 xxmin,xxmax
      integer nstep,ilosi
ccc
      real*8 xxvec(1000),yyvec(1000),yyvec2(1000)
      integer nintxx,howtau
      common/dstauxxtabcom/xxvec,yyvec,yyvec2,nintxx,howtau
ccc
      real*8 rvec(1000),yvec(1000),yvec2(1000),rminl,rmaxl
      integer nintr,nhow,nrset
      common/dstabfuncom/rvec,yvec,yvec2,rminl,rmaxl,nintr,nhow,nrset
ccc
      real*8 xxmini,xxmaxi,taupartial
      integer nstepi
      common/dstaupartialcom/xxmini,xxmaxi,taupartial,nstepi
ccc
      real*8 tauparvec(3)
      integer nstepcall
      common/dstaupartialcom2/tauparvec,nstepcall
ccc
      real*8 xstore(1000),ystore(1000)
      integer nstore
      common/dstaupartialcom3/xstore,ystore,nstore
ccc
      real*8 tollf,tollx,dsomlosisph2,rbvec(0:1000),ybvec(0:1000),resc,
     &  cut,prec,xadd
      integer how,i,dir,nmaxx
      external dstaupartial1,dsomlosisph2
      integer iwin
      real*8 out
ccc
ccc check: you have to call this function with the sequel nstep = 1,2,3
ccc        or nstep = 4 or nstep = 5
ccc
      if(nstep.eq.1.or.nstep.eq.4.or.nstep.eq.5) then
        nstepcall=0
      else
        if(nstep.ne.nstepcall+1) then
          write(*,*) 'DS: calling dstausetsph with wrong nstep sequel'
          write(*,*) 'DS: nstep, nstepcall = ',nstep,nstepcall
          stop
        endif
      endif
ccc
      if(nstep.eq.1) then
        iwin=ilositauset+ilosi
        dvtri(idvx)=xxmin ! this is conventional for 3 dependent variables
        dvtri(idvy)=xxmax
        dvtri(idvz)=1.001d0*nstep
        call dslosidriver(iwin,ntytri,dvtri,out)
ccc      
        if(locilosi.eq.ilosiabscomp) then
ccc
ccc dstaufunsph computed in dslosidriver, do not load the tabulation 
ccc        
          return
        endif    
        taupartial=0.d0
        tollf=100.d0
        how=2
        tollx=(dlog(xxmax)-dlog(xxmin))/1.d2
        nmaxx=100
        prec=1.d-5
        npairs_norm=1.d0
        npairs_norm=dsomlosisph2(xxmin)
        call dsfun_intvec(dsomlosisph2,dsomlosisph2,xxmin,xxmax,tollf
     &   ,tollx,nmaxx,prec,how,rbvec,ybvec,nintr)
ccc
        rvec(nintr+1)=rbvec(nintr)
        yvec(nintr+1)=taupartial
        resc=npairs_norm*3.0856d21
        resc=resc*taureschoweq1
        do i=nintr,1,-1
          rvec(i)=rbvec(i-1)
          yvec(i)=ybvec(i)*resc+yvec(i+1)
        enddo
        nintr=nintr+1
ccc
ccc eventually add extra points and fill in spline tables
ccc
        tollf=1.1d0
        cut=taucut
        dir=1
        how=2
        call dstautabintsph(dsomlosisph2,prec,resc,cut,tollf,how,dir)
ccc
ccc transfer spline tables
ccc
        nintxx=nintr
        do i=1,nintr
          xxvec(i)=rvec(i)
          yyvec(i)=yvec(i)
          yyvec2(i)=yvec2(i)
c          write(*,*) 'step1 : ',i,dexp(xxvec(i)),dexp(yyvec(i))
        enddo
        tauparvec(2)=dexp(yyvec(1))
        nstepcall=1
        howtau=2
        nstore=nintxx
        do i=1,nstore
          xstore(i)=dexp(xxvec(i))
          ystore(i)=dexp(yyvec(i))
        enddo
        return
      elseif(nstep.eq.2.or.nstep.eq.4) then
        iwin=ilositauset+ilosi
        dvtri(idvx)=xxmin ! this is conventional for 3 dependent variables
        dvtri(idvy)=xxmax
        dvtri(idvz)=1.001d0*nstep
        call dslosidriver(iwin,ntytri,dvtri,out)
ccc      
        if(locilosi.eq.ilosiabscomp) then
ccc
ccc dstaufunsph computed in dslosidriver, do not load the tabulation 
ccc        
          return
        endif    
        taupartial=0.d0
        if(nstep.eq.2) taupartial=tauparvec(2)
        nstepi=nstep
        xxmini=xxmin
        tollf=1.1d0
        tollx=(xxmax-xxmin)/100.d0
        how=1
        npairs_norm=1.d0
        npairs_norm=dsomlosisph2(xxmin)
        call dstabfun(dstaupartial1,xxmin,xxmax,tollf,tollx,how)
        nintxx=nintr
        do i=1,nintr
          xxvec(i)=rvec(i)
          yyvec(i)=yvec(i)
          yyvec2(i)=yvec2(i)
c          write(*,*) 'step2/4 : ',i,xxvec(i),yyvec(i)
        enddo
        if(nstep.eq.2) then
          tauparvec(3)=yyvec(nintr)
          nstepcall=2
        endif
        howtau=1
        return
      elseif(nstep.eq.3) then
        if(locilosi.eq.ilosiabscomp) then
ccc
ccc dstaufunsph computed in dslosidriver, do not load the tabulation
ccc but only update nstep:          
ccc        
          iwin=ilositauset+ilosi
          dvtri(idvx)=xxmin ! this is conventional for 3 dependent variables
          dvtri(idvy)=xxmax
          dvtri(idvz)=1.001d0*nstep
          call dslosidriver(iwin,ntytri,dvtri,out)
          return
        endif    
        taupartial=tauparvec(3)
        nintr=nstore
        rvec(1)=xstore(1)
        yvec(1)=taupartial
        do i=2,nstore
          rvec(i)=xstore(i)
          yvec(i)=yvec(i-1)+ystore(i-1)-ystore(i)
        enddo
        do i=0,20
          xadd=dexp(dlog(xxmin)+(dlog(xxmax)-dlog(xxmin))/20.d0*i)
          if(xadd.gt.rvec(nintr)) then
            nintr=nintr+1
            rvec(nintr)=xadd 
            yvec(nintr)=yvec(nintr-1)
          endif
        enddo
        nintxx=nintr
        do i=1,nintr
          xxvec(i)=dlog(rvec(i))
          yyvec(i)=dlog(yvec(i))
c          write(*,*) 'step3 : ',i,dexp(xxvec(i)),dexp(yyvec(i))
        enddo
ccc
ccc fill in the spline:
ccc
        call dsspline(xxvec,yyvec,nintr,1.d31,1.d31,yyvec2)
        nhow=how
        nrset=123456
        howtau=2
      elseif(nstep.eq.5) then
        iwin=ilositauset+ilosi
        dvtri(idvx)=xxmin ! this is conventional for 3 dependent variables
        dvtri(idvy)=xxmax
        dvtri(idvz)=1.001d0*nstep
        call dslosidriver(iwin,ntytri,dvtri,out)
ccc      
        if(locilosi.eq.ilosiabscomp) then
ccc
ccc dstaufunsph computed in dslosidriver, do not load the tabulation 
ccc        
          return
        endif    
        taupartial=0.d0
        tollf=100.d0
        how=2
        prec=1.d-5
        npairs_norm=1.d0
        npairs_norm=dsomlosisph2(xxmin)
        call dsfun_intvec(dsomlosisph2,dsomlosisph2,xxmin,xxmax,tollf
     &   ,prec,how,rbvec,ybvec,nintr)
ccc
        rvec(1)=rbvec(0)
        yvec(1)=taupartial
        resc=npairs_norm*3.0856d21
        resc=resc*taureschoweq1
        do i=2,nintr+1
          rvec(i)=rbvec(i-1)
          yvec(i)=ybvec(i-1)*resc+yvec(i-1)
        enddo
        nintr=nintr+1
ccc
        tollf=1.1d0
        cut=taucut
        dir=2
        how=2
        call dstautabintsph(dsomlosisph2,prec,resc,cut,tollf,how,dir)
c
        nintxx=nintr
        do i=1,nintr
          xxvec(i)=rvec(i)
          yyvec(i)=yvec(i)
          yyvec2(i)=yvec2(i)
c          write(*,*) 'step5 : ',i,dexp(xxvec(i)),dexp(yyvec(i))
        enddo
        howtau=2
        return
      else
        write(*,*) 'DS: dstausetsph called with wrong nstep = ',nstep
        stop
      endif
      end
ccc
ccc
ccc
      subroutine dstautabintsph(fun,prec,resc,cut,tollf,how,dir)
************************************************************************
*** integrates the tabulations of dstaufunsph with extra points in case
*** two nearby points are such that:
***   max(yvec(i)/yvec(i+1),yvec(i+1)/yvec(i)) > tollf  &
***   max(yvec(i),yvec(i+1)) > cut
*** dstaufunsph comes from the integral of the external function fun, 
*** with precision prec. 
*** other inputs:
***   how = 1: add subintervals in linear scale, 
***   how = 2: add subintervals in log scale
***   dir=1,  rvec(nintr) is the minimum
***   dir=2,  rvec(nintr) is the maximum
************************************************************************
      implicit none
      real*8 fun,prec,resc,cut,tollf,tollr
      integer how,dir
      external fun
ccc
      real*8 rvec(1000),yvec(1000),yvec2(1000),rminl,rmaxl
      integer nintr,nhow,nrset
      common/dstabfuncom/rvec,yvec,yvec2,rminl,rmaxl,nintr,nhow,nrset
ccc
      real*8 rint,yint,diff,miny,maxy,rinf,rsup,yinf,ysup,addres,
     &  trap,result,eps
      integer nmax,i,k,kk,nint
      logical addpoint
ccc
ccc maximum number of intervals, this is set internally since it goes in
ccc the dimension of variables in the common block dstabfuncom 
ccc
      nmax=1000
ccc
ccc check the flag how, tollf, tollr
ccc
      if(.not.(how.eq.1.or.how.eq.2)) then
        write(*,*) 'DS: dstautabintsph called with wrong option how : ',
     &    how
        stop 
      endif
      if(tollf.le.1.d0) then
        write(*,*) 'DS: in dstautabintsph tollf = ',tollf,
     &             ' not properly intialized'
        stop
      endif
      if(how.eq.1) then
        tollr=(rvec(nintr)-rvec(1))/1000.d0
      else ! if(how.eq.2) 
        tollr=(dlog(rvec(nintr))-dlog(rvec(1)))/1000.d0
      endif
      nint=nintr
      if(.not.(dir.eq.1.or.dir.eq.2)) then
        write(*,*) 'DS: dstautabintsph called with wrong option dir : ',
     &    dir
        stop 
      endif
ccc
ccc fill in vectors
ccc
 100  continue
      do k=1,nint-1
        miny=min(dabs(yvec(k)),dabs(yvec(k+1)))
        maxy=max(dabs(yvec(k)),dabs(yvec(k+1)))
        addpoint=.false.
        if(maxy.gt.miny*tollf.and.maxy.gt.cut) addpoint=.true.
        if(how.eq.1) then
          diff=dabs(rvec(k)-rvec(k+1))
        else !if(how.eq.2)
          diff=dabs(dlog(rvec(k))-dlog(rvec(k+1)))
        endif
        if(diff.lt.tollr) addpoint=.false.
        if(addpoint) then
          if(how.eq.1) then
            rint=(rvec(k+1)+rvec(k))/2.d0
          else ! if(how.eq.2) 
            rint=dexp((dlog(rvec(k+1))+dlog(rvec(k)))/2.d0)
          endif
          if(dir.eq.1) then
            rinf=rint
            rsup=rvec(k+1)
            addres=yvec(k+1)
          else ! if(dir.eq.2)
            rinf=rvec(k)
            rsup=rint
            addres=yvec(k)
          endif
          yinf=fun(rinf)
          ysup=fun(rsup)
          trap=0.5d0*(dabs(yinf)+dabs(ysup))*(rsup-rinf)
ccc
ccc you might be trying to integrate a fuction identically 0
ccc
          if(dabs(trap).eq.0.d0) then
            result=0.d0
            goto 30
          endif
          eps=prec*trap/10.d0
c          write(*,*) 
c     &      rvec(k+1),rvec(k),yvec(k+1),yvec(k),rinf,rsup,eps,prec
          call dsfun_int(fun,rinf,rsup,eps,prec,result)
 30       yint=resc*result+addres
          do kk=nint,k+1,-1
            rvec(kk+1)=rvec(kk)
            yvec(kk+1)=yvec(kk)
          enddo
          rvec(k+1)=rint
          yvec(k+1)=yint
c          write(*,*) 'add  interval: xint, yint : ',
c     &                 k+1,rvec(k+1),yvec(k+1)
          nint=nint+1  
          if(nint.le.nmax) then
            goto 100 
          else
            write(*,*) 'DS: in dstautabintsph exceeded the maximum dim'
            write(*,*) 'DS: allowed for vectors in the rvec,yvec block'
            write(*,*) 'DS: which is to nmax : ',nmax
            stop
          endif
        endif  
        continue
      enddo
ccc
ccc take out values lt cut at the extremes of the tabulated interval 
ccc
 300  continue
      if(dabs(yvec(nint)).lt.cut) then
        nint=nint-1
        rmaxl=rvec(nint)
        goto 300
      endif
 400  continue
      if(dabs(yvec(1)).lt.cut) then
        do i=1,nint-1
          rvec(i)=rvec(i+1)
          yvec(i)=yvec(i+1)
        enddo
        nint=nint-1
        rminl=rvec(1)
        goto 400
      endif
ccc
ccc if after taking out the zeros nint.le.2 there is a problem, print out the
ccc the initial points and stop
ccc
      if(nint.le.2) then
        write(*,*) 'DS: too many zeros in dstautabintsph'
        stop
      endif
ccc
ccc if how.eq.2 tabulation in log,log:
ccc
      if(how.eq.2) then
        do i=1,nint
          rvec(i)=dlog(rvec(i))
          yvec(i)=dlog(yvec(i))
        enddo
      endif
ccc
ccc fill in the spline:
ccc
      nhow=how
      nintr=nint
      call dsspline(rvec,yvec,nintr,1.d31,1.d31,yvec2)
      nrset=123456
      return
      end
ccc
ccc
ccc
      real*8 function dstaupartial1(xx)
      implicit none
      real*8 xx
      include 'dsdmdintcom.h' ! npairs_norm,taucut,taureschoweq1
ccc
      real*8 xxmini,xxmaxi,taupartial
      integer nstepi
      common/dstaupartialcom/xxmini,xxmaxi,taupartial,nstepi
ccc
      real*8 xxlocmin,xxlocmax,result,dstaupartial2,resc
      integer ihow
ccc
      if(nstepi.eq.2.or.nstepi.eq.4) then
        ihow=1
        xxlocmax=xx
        xxlocmin=xxmini
        result=dstaupartial2(xxlocmin,xxlocmax,ihow)
      else
        write(*,*) 'DS: dstaupartial1 called with wrong nstepi = ',
     &       nstepi
        stop
      endif
      resc=npairs_norm*3.0856d21
      resc=resc*taureschoweq1
      result=result*resc+taupartial
      if(result.lt.taucut) result=0.d0
      dstaupartial1=result
      return
      end
ccc
ccc
ccc
      real*8 function dstaupartial2(xxmin,xxmax,ihow)
      implicit none
      real*8 xxmin,xxmax
      integer ihow
ccc
      real*8 tollf2,tollab,prec,result,eps,dsomlosisph2
      integer nmax
      external dsomlosisph2
ccc
      if(ihow.eq.1) then
        prec=1.d-4
        eps=(xxmax-xxmin)*prec*dsomlosisph2(xxmax)/10.d0
        call dsfun_int(dsomlosisph2,xxmin,xxmax,eps,prec,result)
      elseif(ihow.eq.2) then
        tollf2=100.d0
        tollab=(dlog(xxmax)-dlog(xxmin))/1000.d0
        nmax=1000
        prec=1.d-5
        call dsfun_intpar(dsomlosisph2,dsomlosisph2,xxmin,xxmax,
     &     tollf2,tollab,nmax,prec,ihow,result)
      else
        write(*,*) 'DS: dstaupartial2 called with wrong ihow = ',ihow
        stop
      endif
      dstaupartial2=result
      return
      end

