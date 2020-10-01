      subroutine dstabfun(fun,rmin,rmax,tollf,tollr,how)
c_______________________________________________________________________
c  set a tabulation for the function fun in the interval (rmin, rmax);
c  intermadiate points r_i are such that 
c    max(fun(r_i),fun(r_i+1)) < tollf * min(fun(r_i),fun(r_i+1))
c  input:
c    function to tabulate fun
c    extrema for tabulation rmin,rmax
c    tollerance on tabulation tollf
c    tollerance on how small the interval (r_i,r_i+1) can be tollr
c    maximum numeber of intervals nmax (this must be .le.1000)
c    how = 1: add subintervals and tabulate fun in linear scale 
c    how = 2: add subintervals and tabulate fun in log scale
c  output:
c    result of integration res
c_______________________________________________________________________
      implicit none
      real*8 fun,rmin,rmax,tollf,tollr
      integer how
      external fun
ccc
      real*8 rvec(1000),yvec(1000),yvec2(1000),rminl,rmaxl
      integer nintr,nhow,nrset
      common/dstabfuncom/rvec,yvec,yvec2,rminl,rmaxl,nintr,nhow,nrset
ccc
      real*8 rveccheck(100),yveccheck(100)
      real*8 rint,yint,diff,ymax,zero,miny,maxy
      integer nmax,i,k,kk,nint,nintcheck
      logical dsisnanorinf,addpoint
c
c maximum number of intervals, this is set internally since it goes in
c the dimension of variables in the common block dstabfuncom 
c
      nmax=1000
c
c check the flag how, tollf, tollr
c
      if(.not.(how.eq.1.or.how.eq.2)) then
        write(*,*) 'DS: dstabfun called with wrong option how : ',how
        stop 
      endif
      if(tollf.le.1.d0) then
        write(*,*) 'DS: in dstabfun tollf = ',tollf,
     &             ' not properly intialized'
        stop
      endif
      if(how.eq.1.and.tollr*1.d3*nmax.lt.(rmax-rmin)) then
        write(*,*) 'DS: in dstabfun tollr, nmax not properly intialized'
        write(*,*) 'DS: rmax-rmin, tollr, nmax = ',rmax-rmin,tollr,nmax 
        stop
      endif
      if(how.eq.2.and.tollr*1.d3*nmax.lt.(dlog(rmax)-dlog(rmin))) then
        write(*,*) 'DS: in dstabfun tollr, nmax not properly intialized'
        write(*,*) 'DS: dlog(rmax)-dlog(rmin), tollr, nmax = ',
     &             dlog(rmax)-dlog(rmin),tollr,nmax 
        stop
      endif     
c
c clean the table:
c
      do i=1,1000
        rvec(i)=0.d0
        yvec(i)=0.d0
      enddo
c
c start with 5 points equally spaced in lin or log scale
c
      rminl=rmin   
      rmaxl=rmax   
      nint=5
      nintcheck=nint
      if(nintcheck.gt.100) then
        write(*,*) 'DS: internal inconsistency in dstabfun'
        stop
      endif
      if(how.eq.1) then
        diff=(rmax-rmin)/dble(nint-1)
      else ! if(how.eq.2) 
        diff=(dlog(rmax)-dlog(rmin))/dble(nint-1)
      endif
      ymax=-1.d40
      do i=1,nint
        if(how.eq.1) then
          rvec(i)=rmin+diff*(i-1)
        else ! if(how.eq.2) 
          rvec(i)=dexp(dlog(rmin)+diff*(i-1))
        endif
        yvec(i)=fun(rvec(i))
        if(dsisnanorinf(yvec(i))) then
          write(*,*) 'DS: dstabfun: fun is = ',yvec(i),
     &               ' at x = ',rvec(i)
          write(*,*) 'DS: within the interval of tabulation: ',rmin,rmax
          stop 
        endif
        if(dabs(yvec(i)).gt.ymax) ymax=dabs(yvec(i)) 
c        write(*,*) 'r,res : ',i,rvec(i),yvec(i)
        rveccheck(i)=rvec(i)
        yveccheck(i)=yvec(i)
      enddo
c
c set here what you are calling 'zero'
c
      zero=1.d-30*ymax
c
c fill in vectors
c
 100  continue
      do k=1,nint-1
        miny=min(dabs(yvec(k)),dabs(yvec(k+1)))
        maxy=max(dabs(yvec(k)),dabs(yvec(k+1)))
        addpoint=.false.
        if(maxy.gt.miny*tollf.and.maxy.gt.zero) addpoint=.true.
        if(how.eq.1) then
          diff=dabs(rvec(k)-rvec(k+1))
        else !if(how.eq.2) 
          diff=dabs(dlog(rvec(k))-dlog(rvec(k+1)))
        endif
        if(diff.lt.tollr) addpoint=.false.
        if(maxy.gt.miny*(1.d0+10.d0*(tollf-1.d0))
     &     .and.diff.gt.0.1d0*tollr) addpoint=.true.
        if(addpoint) then
          if(how.eq.1) then
            rint=(rvec(k+1)+rvec(k))/2.d0
          else ! if(how.eq.2) 
            rint=dexp((dlog(rvec(k+1))+dlog(rvec(k)))/2.d0)
          endif
          yint=fun(rint)
c
c if you hit a singularity skip the point
c
          if(dsisnanorinf(yint)) goto 200
          do kk=nint,k+1,-1
            rvec(kk+1)=rvec(kk)
            yvec(kk+1)=yvec(kk)
          enddo
          rvec(k+1)=rint
          yvec(k+1)=yint
c          write(*,*) 'add  interval: xint, yint : ',
c     &                 k+1,rvec(k+1),yvec(k+1)
          nint=nint+1  
          if(nint.lt.nmax) then
            goto 100 
          else
            write(*,*) 'DS: in dstabfun exceeded the maximum dim'
            write(*,*) 'DS: allowed for vectors in the rvec,yvec block'
            write(*,*) 'DS: which is to nmax = ',nmax
            do kk=1,nint-1
              write(*,*) 'DS : ',kk,rvec(kk),yvec(kk)
            enddo
            stop
          endif
        endif  
 200    continue
      enddo
c
c take out eventual zeros at the extremes of the tabulated interval 
c
 300  continue
      if(dabs(yvec(nint)).lt.zero) then
        nint=nint-1
        rmaxl=rvec(nint)
        goto 300
      endif
 400  continue
      if(dabs(yvec(1)).lt.zero) then
        do i=1,nint-1
          rvec(i)=rvec(i+1)
          yvec(i)=yvec(i+1)
        enddo
        nint=nint-1
        rminl=rvec(1)
        goto 400
      endif
c
c if after taking out the zeros nint.le.2 there is a problem, print out the
c the initial points and stop
c
      if(nint.le.2) then
        write(*,*) 'DS: too many zeros in dstabfun'
        do i=1,nintcheck
          write(*,*) 'DS: #,r,res = ',i,rveccheck(i),yveccheck(i)
        enddo
        stop
      endif
c
c add points very close to the extremes of the tabulated interval 
c
      goto 999
      k=1
      if(how.eq.1) then
        rint=rvec(k)+(rvec(k+1)-rvec(k))/50.d0
      else ! if(how.eq.2) 
        rint=dexp(dlog(rvec(k))+(dlog(rvec(k+1))-dlog(rvec(k)))/50.d0)
      endif
      yint=fun(rint)
      write(*,*) rint,yint 
      do kk=nint,k+1,-1
        rvec(kk+1)=rvec(kk)
        yvec(kk+1)=yvec(kk)
      enddo
      rvec(k+1)=rint
      yvec(k+1)=yint
      nint=nint+1  
c
      k=nint-1
      if(how.eq.1) then
        rint=rvec(k+1)-(rvec(k+1)-rvec(k))/50.d0
      else ! if(how.eq.2) 
        rint=dexp(dlog(rvec(k+1))-(dlog(rvec(k+1))-dlog(rvec(k)))/50.d0)
      endif
      yint=fun(rint)
      write(*,*) rint,yint 
      do kk=nint,k+1,-1
        rvec(kk+1)=rvec(kk)
        yvec(kk+1)=yvec(kk)
      enddo
      rvec(k+1)=rint
      yvec(k+1)=yint
      nint=nint+1
 999  continue
c
c if how.eq.2 tabulation in log,log:
c
      if(how.eq.2) then
        do i=1,nint
          rvec(i)=dlog(rvec(i))
          yvec(i)=dlog(max(zero,yvec(i)))
        enddo
      endif
c
c fill in the spline:
c
      nhow=how
      nintr=nint
      call dsspline(rvec,yvec,nintr,1.d31,1.d31,yvec2)
      nrset=123456
      return
      end


      subroutine dstabfunb(fun,rmin,rmax,tollf,tollr,npmin,how)
c_______________________________________________________________________
c  set a tabulation for the function fun in the interval (rmin, rmax);
c  intermadiate points r_i are such that 
c    max(fun(r_i),fun(r_i+1)) < tollf * min(fun(r_i),fun(r_i+1))
c  input:
c    function to tabulate fun
c    extrema for tabulation rmin,rmax
c    tollerance on tabulation tollf
c    tollerance on how small the interval (r_i,r_i+1) can be tollr
c    maximum numeber of intervals nmax (this must be .le.1000)
c    npmin = minimum number of points
c    how = 1: add subintervals and tabulate fun in linear scale 
c    how = 2: add subintervals and tabulate fun in log scale
c  output:
c    result of integration res
c_______________________________________________________________________
      implicit none
      real*8 fun,rmin,rmax,tollf,tollr
      integer npmin,how
      external fun
ccc
      real*8 rvec(1000),yvec(1000),yvec2(1000),rminl,rmaxl
      integer nintr,nhow,nrset
      common/dstabfuncom/rvec,yvec,yvec2,rminl,rmaxl,nintr,nhow,nrset
ccc
      real*8 rveccheck(100),yveccheck(100)
      real*8 rint,yint,diff,ymax,zero,miny,maxy
      integer nmax,i,k,kk,nint,nintcheck
      logical dsisnanorinf,addpoint
c
c maximum number of intervals, this is set internally since it goes in
c the dimension of variables in the common block dstabfuncom 
c
      nmax=1000
c
c check the flag how, tollf, tollr
c
      if(.not.(how.eq.1.or.how.eq.2)) then
        write(*,*) 'DS: dstabfun called with wrong option how : ',how
        stop 
      endif
      if(tollf.le.1.d0) then
        write(*,*) 'DS: in dstabfun tollf = ',tollf,
     &             ' not properly intialized'
        stop
      endif
      if(how.eq.1.and.tollr*1.d3*nmax.lt.(rmax-rmin)) then
        write(*,*) 'DS: in dstabfun tollr, nmax not properly intialized'
        write(*,*) 'DS: rmax-rmin, tollr, nmax = ',rmax-rmin,tollr,nmax 
        stop
      endif
      if(how.eq.2.and.tollr*1.d3*nmax.lt.(dlog(rmax)-dlog(rmin))) then
        write(*,*) 'DS: in dstabfun tollr, nmax not properly intialized'
        write(*,*) 'DS: dlog(rmax)-dlog(rmin), tollr, nmax = ',
     &             dlog(rmax)-dlog(rmin),tollr,nmax 
        stop
      endif     
c
c clean the table:
c
      do i=1,1000
        rvec(i)=0.d0
        yvec(i)=0.d0
      enddo
c
c start with 5 points equally spaced in lin or log scale
c
      rminl=rmin   
      rmaxl=rmax   
      nint=npmin
      nintcheck=nint
      if(nintcheck.gt.100) then
        write(*,*) 'DS: internal inconsistency in dstabfun'
        stop
      endif
      if(how.eq.1) then
        diff=(rmax-rmin)/dble(nint-1)
      else ! if(how.eq.2) 
        diff=(dlog(rmax)-dlog(rmin))/dble(nint-1)
      endif
      ymax=-1.d40
      do i=1,nint
        if(how.eq.1) then
          rvec(i)=rmin+diff*(i-1)
        else ! if(how.eq.2) 
          rvec(i)=dexp(dlog(rmin)+diff*(i-1))
        endif
        yvec(i)=fun(rvec(i))
        if(dsisnanorinf(yvec(i))) then
          write(*,*) 'DS: dstabfun: fun is = ',yvec(i),
     &               ' at x = ',rvec(i)
          write(*,*) 'DS: within the interval of tabulation: ',rmin,rmax
          stop 
        endif
        if(dabs(yvec(i)).gt.ymax) ymax=dabs(yvec(i)) 
c        write(*,*) 'r,res : ',i,rvec(i),yvec(i)
        rveccheck(i)=rvec(i)
        yveccheck(i)=yvec(i)
      enddo
c
c set here what you are calling 'zero'
c
      zero=1.d-30*ymax
c
c fill in vectors
c
 100  continue
      do k=1,nint-1
        miny=min(dabs(yvec(k)),dabs(yvec(k+1)))
        maxy=max(dabs(yvec(k)),dabs(yvec(k+1)))
        addpoint=.false.
        if(maxy.gt.miny*tollf.and.maxy.gt.zero) addpoint=.true.
        if(how.eq.1) then
          diff=dabs(rvec(k)-rvec(k+1))
        else !if(how.eq.2) 
          diff=dabs(dlog(rvec(k))-dlog(rvec(k+1)))
        endif
        if(diff.lt.tollr) addpoint=.false.
        if(maxy.gt.miny*(1.d0+10.d0*(tollf-1.d0))
     &     .and.diff.gt.0.1d0*tollr) addpoint=.true.
        if(addpoint) then
          if(how.eq.1) then
            rint=(rvec(k+1)+rvec(k))/2.d0
          else ! if(how.eq.2) 
            rint=dexp((dlog(rvec(k+1))+dlog(rvec(k)))/2.d0)
          endif
          yint=fun(rint)
c
c if you hit a singularity skip the point
c
          if(dsisnanorinf(yint)) goto 200
          do kk=nint,k+1,-1
            rvec(kk+1)=rvec(kk)
            yvec(kk+1)=yvec(kk)
          enddo
          rvec(k+1)=rint
          yvec(k+1)=yint
c          write(*,*) 'add  interval: xint, yint : ',
c     &                 k+1,rvec(k+1),yvec(k+1)
          nint=nint+1  
          if(nint.lt.nmax) then
            goto 100 
          else
            write(*,*) 'DS: in dstabfun exceeded the maximum dim'
            write(*,*) 'DS: allowed for vectors in the rvec,yvec block'
            write(*,*) 'DS: which is to nmax = ',nmax
            do kk=1,nint-1
              write(*,*) 'DS : ',kk,rvec(kk),yvec(kk)
            enddo
            stop
          endif
        endif  
 200    continue
      enddo
c
c take out eventual zeros at the extremes of the tabulated interval 
c
 300  continue
      if(dabs(yvec(nint)).lt.zero) then
        nint=nint-1
        rmaxl=rvec(nint)
        goto 300
      endif
 400  continue
      if(dabs(yvec(1)).lt.zero) then
        do i=1,nint-1
          rvec(i)=rvec(i+1)
          yvec(i)=yvec(i+1)
        enddo
        nint=nint-1
        rminl=rvec(1)
        goto 400
      endif
c
c if after taking out the zeros nint.le.2 there is a problem, print out the
c the initial points and stop
c
      if(nint.le.2) then
        write(*,*) 'DS: too many zeros in dstabfun'
        do i=1,nintcheck
          write(*,*) 'DS: #,r,res = ',i,rveccheck(i),yveccheck(i)
        enddo
        stop
      endif
c
c add points very close to the extremes of the tabulated interval 
c
      goto 999
      k=1
      if(how.eq.1) then
        rint=rvec(k)+(rvec(k+1)-rvec(k))/50.d0
      else ! if(how.eq.2) 
        rint=dexp(dlog(rvec(k))+(dlog(rvec(k+1))-dlog(rvec(k)))/50.d0)
      endif
      yint=fun(rint)
      write(*,*) rint,yint 
      do kk=nint,k+1,-1
        rvec(kk+1)=rvec(kk)
        yvec(kk+1)=yvec(kk)
      enddo
      rvec(k+1)=rint
      yvec(k+1)=yint
      nint=nint+1  
c
      k=nint-1
      if(how.eq.1) then
        rint=rvec(k+1)-(rvec(k+1)-rvec(k))/50.d0
      else ! if(how.eq.2) 
        rint=dexp(dlog(rvec(k+1))-(dlog(rvec(k+1))-dlog(rvec(k)))/50.d0)
      endif
      yint=fun(rint)
      write(*,*) rint,yint 
       do kk=nint,k+1,-1
        rvec(kk+1)=rvec(kk)
        yvec(kk+1)=yvec(kk)
      enddo
      rvec(k+1)=rint
      yvec(k+1)=yint
      nint=nint+1
 999  continue
c
c if how.eq.2 tabulation in log,log:
c
      if(how.eq.2) then
        do i=1,nint
          rvec(i)=dlog(rvec(i))
          yvec(i)=dlog(max(zero,yvec(i)))
        enddo
      endif
c
c fill in the spline:
c
      nhow=how
      nintr=nint
      call dsspline(rvec,yvec,nintr,1.d31,1.d31,yvec2)
      nrset=123456
      return
      end



