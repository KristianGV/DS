      subroutine dstabfunl(fun,rmin,rmax,fpmin,fpmax,tollf,tollr,how)
c_______________________________________________________________________
c  set a tabulation for the function fun in the interval (rmin, rmax);
c  a sparse regular grid is set and then intermadiate points r_i are
c  added up to reaching the condition that: 
c
c    abs(1-interpolation(r_i)/fun(r_i)) < tollf      
c
c  input:
c    function to tabulate fun
c    extrema for tabulation rmin,rmax
c    fpmin & fpmax: first derivative of fun computed in rmin,rmax; to
c      use natural splines call the subroutine with values larger than
c      1.d30; for the linear interpolation case these are not used
c    tollerance on tabulation tollf
c    tollerance on how small the interval (r_i,r_i+1) can be tollr
c    maximum numeber of intervals nmax (this must be .le.1000)
c    how = 1: add subintervals and tabulate fun in linear scale with
c             linear interpolation 
c    how = 2: add subintervals and tabulate fun in log scale with
c             linear interpolation
c    how = 101: the same as how=1 but with cubic spline interpolation
c    how = 102: the same as how=2 but with cubic spline interpolation
c  output:
c    tabulation stored in dstabfuncom common block
c_______________________________________________________________________
      implicit none
      real*8 fun,rmin,rmax,fpmin,fpmax,tollf,tollr
      integer how
      external fun
ccc
      real*8 rvec(1000),yvec(1000),yvec2(1000),rminl,rmaxl
      integer nintr,nhow,nrset
      common/dstabfuncom/rvec,yvec,yvec2,rminl,rmaxl,nintr,nhow,nrset
ccc
      real*8 rint,yint,diff,diffy,yextr
      integer nmax,i,kk,nint,ki
      logical dsisnanorinf
ccc
      integer ityp,iar,iaf,nbkp
      integer ilint,isint,ialin,ialog 
      parameter(ilint=1,isint=2,ialin=1,ialog=2)
ccc
ccc check how flag and set intern
ccc      
      if(how.lt.100) then
        ityp=ilint      ! linear interpolation
        if(how.eq.1) then
          iar=ialin
          iaf=ialin
        elseif(how.eq.2) then
          iar=ialog
          iaf=ialog
        else
          write(*,*) 'DS: calling dstabfunl with wrong how = ',how
          stop
        endif  
      else
        ityp=isint      ! cubic spline interpolation
        if(how.eq.101) then
          iar=ialin
          iaf=ialin
        elseif(how.eq.102) then
          iar=ialog
          iaf=ialog
        else
          write(*,*) 'DS: calling dstabfunl with wrong how = ',how
          stop
        endif  
      endif
c
c maximum number of intervals, this is set internally since it goes in
c the dimension of variables in the common block dstabfuncom 
c
      nmax=1000
c
c check whether tollr is set within a reasonable range
c
      if(iar.eq.ialin.and.tollr*1.d3*nmax.lt.(rmax-rmin)) then
        write(*,*) 'DS: in dstabfunl tollr not properly intialized'
        write(*,*) 'DS: rmax-rmin, tollr, nmax = ',rmax-rmin,tollr,nmax 
        stop
      endif
      if(iar.eq.ialog.and.tollr*1.d3*nmax.lt.(dlog(rmax)-dlog(rmin)))
     &   then
        write(*,*) 'DS: in dstabfunl tollr not properly intialized'
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
      nint=10
      nbkp=0
      if(iar.eq.ialin) then
        diff=(rmax-rmin)/dble(nint-1)
      elseif(iar.eq.ialog) then
        diff=(dlog(rmax)-dlog(rmin))/dble(nint-1)
      else
        write(*,*) 'DS: in dstabfunl wrong iar = ',iar
        stop
      endif
      do i=1,nint
        if(iar.eq.ialin) then
          rvec(i)=rmin+diff*(i-1)
        elseif(iar.eq.ialog) then
          rvec(i)=dexp(dlog(rmin)+diff*(i-1))
        endif
        yvec(i)=fun(rvec(i))
        if(dsisnanorinf(yvec(i))) then
          write(*,*) 'DS: dstabfunl: fun is = ',yvec(i),
     &               ' at x = ',rvec(i)
          write(*,*) 'DS: within the interval of tabulation: ',rmin,rmax
          stop 
        endif
c        write(*,*) 'r,res : ',i,rvec(i),yvec(i)
      enddo
c
c go to logs if required:
c      
      if(iaf.eq.ialog) then
        do i=1,nint
          if(yvec(i).le.0.d0) then 
            write(*,*) 'DS: dstabfunl: request to tabulate in log scale'
            write(*,*) 'DS: negative fun: ',yvec(i),' at x = ',rvec(i)
            write(*,*) 'DS: within the interval of tab.: ',rmin,rmax
            stop 
          endif
          yvec(i)=dlog(yvec(i))
        enddo   
      endif
ccc
      if(iar.eq.ialog) then
        do i=1,nint
          rvec(i)=dlog(rvec(i))
        enddo   
      endif
c
c fill in vectors
c
      ki=1
 101  continue
      if(ityp.eq.isint.and.nint.ne.nbkp) then
        call dsspline(rvec,yvec,nint,fpmin,fpmax,yvec2)
        nbkp=nint
      endif  
      if(ki.eq.nint) goto 102
c
c if diff.gt.tollr add one point between ki and ki+1     
c
      diff=dabs(rvec(ki)-rvec(ki+1))
      if(diff.lt.tollr) then
        ki=ki+1
        goto 101
      endif
      rint=rvec(ki)+diff/2.d0
      if(iar.eq.ialin) then
        yint=fun(rint)
      elseif(iar.eq.ialog) then
        yint=fun(dexp(rint))
      endif
      if(dsisnanorinf(yint)) then
        if(iar.eq.ialog) rint=dexp(rint)
        write(*,*) 'DS: dstabfunl: fun is = ',yint,' at x = ',rint
        write(*,*) 'DS: within the interval of tabulation: ',rmin,rmax
        stop 
      endif
      if(iaf.eq.ialog) then
        if(yint.le.0.d0) then 
          if(iar.eq.ialog) rint=dexp(rint)
          write(*,*) 'DS: dstabfunl: request to tabulate in log scale'
          write(*,*) 'DS: negative fun: ',yint,' at x = ',rint
          write(*,*) 'DS: within the interval of tab.: ',rmin,rmax
          stop 
        endif
        yint=dlog(yint)
      endif
      if(ityp.eq.isint) then
        call dssplint(rvec,yvec,yvec2,nint,rint,yextr)
      elseif(ityp.eq.ilint) then
        yextr=(yvec(ki+1)+yvec(ki))/2.d0
      endif  
      diffy=dabs(yint/yextr-1.d0)
ccc      
      do kk=nint,ki+1,-1
        rvec(kk+1)=rvec(kk)
        yvec(kk+1)=yvec(kk)
      enddo
      rvec(ki+1)=rint
      yvec(ki+1)=yint
c      write(*,*) 'adding: ',ki,rvec(ki+1),yvec(ki+1)
      nint=nint+1  
      if(nint.ge.nmax) then
        write(*,*) 'DS: in dstabfunl exceeded the maximum dim'
        write(*,*) 'DS: allowed for vectors in the rvec,yvec block'
        write(*,*) 'DS: which is to nmax = ',nmax
        do kk=1,nint-1
          write(*,*) 'DS : ',kk,rvec(kk),yvec(kk)
        enddo
        stop
      endif      
      if(diffy.lt.tollf) ki=ki+2
      goto 101
c      
 102  continue
      nhow=how
      nintr=nint
c
c even in case the linear interpolation used, load a natural spline:
c      
      if(ityp.eq.ilint) then 
        call dsspline(rvec,yvec,nintr,1.d31,1.d31,yvec2)
      endif  
      nrset=123456
      return
      end
