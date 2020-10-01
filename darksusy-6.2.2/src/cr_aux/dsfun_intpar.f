      subroutine dsfun_intpar(f1,f2,a,b,tollf2,tollab,nmax,prec,how,res)
c_______________________________________________________________________
c  integrate function f1 between a and b breaking the integration 
c  interval into subintervals (ap,bp) such that 
c    max(f2(ap),f2(bp)) < tollf2 * min(f2(ap),f2(bp))
c  input:
c    integration function f1
c    weight function f2
c    integration limits a and b
c    tollerance on weight function f2, tollf2
c    tollerance on how small the interval (ap,bp) can be tollab
c    maximum number of intervals nmax 
c    precision on the sub-integral prec
c    how = 1: add subintervals in linear scale, 
c    how = 2: add subintervals in log scale
c  integration performed with the subroutine dsfun_int
c  output:
c    result of integration res
c_______________________________________________________________________
      implicit none
ccc
      real*8 f1,f2,a,b,tollf2,tollab,prec,res
      integer nmax,how
      external f1,f2
ccc
      real*8 xvec(0:nmax),yvec(0:nmax),miny,maxy,diffx,xint,yint,
     &  partial,xinf,xsup,yinf,ysup,trap,eps,result
      logical dsisnanorinf,addpoint
      integer i,step,nint,k,kk,kmax
c
c check the flag how, tollf2, tollab
c
      if(.not.(how.eq.1.or.how.eq.2)) then
        write(*,*) 'DS: dsfun_intpar called with wrong option how : ',
     &              how
        stop 
      endif
      if(tollf2.le.1.d0) then
        write(*,*) 'DS: in dsfun_intpar tollf2 = ',tollf2,
     &             ' not properly intialized'
        stop 
      endif
      if(how.eq.1.and.tollab*1.d3*nmax.lt.(b-a)) then
      write(*,*) 'DS: in dsfun_intpar tollab, ', 
     &           'nmax not properly intialized'
      write(*,*) 'DS: b-a, tollab, nmax : ',b-a,tollab,nmax 
      stop
      endif
      if(how.eq.2.and.tollab*1.d3*nmax.lt.(dlog(b)-dlog(a))) then
      write(*,*) 'DS: in dsfun_intpar tollab, ', 
     &           'nmax not properly intialized'
      write(*,*) 'DS: dlog(b)-dlog(a), tollab, nmax : ',
     &             dlog(b)-dlog(a),tollab,nmax 
      stop
      endif
c
c start with two intervals
c
      xvec(0)=a
      if(how.eq.1) then
        xvec(1)=0.5d0*(a+b)
      else ! if(how.eq.2) 
        xvec(1)=dexp(0.5d0*(dlog(a)+dlog(b)))
      endif
      xvec(2)=b
      do i=0,2
        yvec(i)=f2(xvec(i))
c        write(*,*) 'starting : ',xvec(i),yvec(i) 
      enddo
      nint=2
c
c f2 may be singular in one or both extremes of integration
c
      step=0
 100  continue
      if(dsisnanorinf(yvec(0))) then
        step=step+1
        if(how.eq.1) then
          xvec(0)=xvec(0)+step*tollab
        else ! if(how.eq.2) 
          xvec(0)=dexp(dlog(xvec(0))+step*tollab)
        endif
        if(step.le.3.and.xvec(0).lt.xvec(1)) then
          yvec(0)=f2(xvec(0))
        else
          write(*,*) 'DS: dsfun_intpar: f2(a) is : ',yvec(0)
          stop 
        endif
        goto 100
      endif
      step=0
 200  continue
      if(dsisnanorinf(yvec(2))) then
        step=step+1
        if(how.eq.1) then
          xvec(2)=xvec(2)-step*tollab
        else !if(how.eq.2) 
          xvec(2)=dexp(dlog(xvec(2))-step*tollab)
        endif
        if(step.le.3.and.xvec(2).gt.xvec(1)) then
          yvec(2)=f2(xvec(2))
        else
          write(*,*) 'DS: dsfun_intpar: f2(b) is : ',yvec(2)
          stop 
        endif
        goto 200
      endif
c
c if f2 is singular in (a+b)/2 or exp((log(a)+log(b))/2) there is 
c   a problem
c
      if(dsisnanorinf(yvec(1))) then
        write(*,*) 'DS: dsfun_intpar: f2 is : ',yvec(1),
     &             ' at x = ',xvec(1)
        write(*,*) 'DS: within the interval of integration: ',a,b
        stop 
      endif
c
c fill in vectors
c
 300  continue
      do k=0,nint-1
        miny=min(dabs(yvec(k)),dabs(yvec(k+1)))
        maxy=max(dabs(yvec(k)),dabs(yvec(k+1)))
        addpoint=.false.
        if(maxy.gt.miny*tollf2) addpoint=.true.
        if(how.eq.1) then
          diffx=dabs(xvec(k)-xvec(k+1))
        else !if(how.eq.2) 
          diffx=dabs(dlog(xvec(k))-dlog(xvec(k+1)))
        endif
        if(diffx.lt.tollab) addpoint=.false.
        if(addpoint) then
          if(how.eq.1) then
            xint=(xvec(k+1)+xvec(k))/2.d0
          else ! if(how.eq.2) 
            xint=dexp((dlog(xvec(k+1))+dlog(xvec(k)))/2.d0)
          endif
          yint=f2(xint)
c
c if you hit a singularity skip the point
c
          if(dsisnanorinf(yint)) goto 400
          do kk=nint,k+1,-1
            xvec(kk+1)=xvec(kk)
            yvec(kk+1)=yvec(kk)
          enddo
          xvec(k+1)=xint
          yvec(k+1)=yint
c          write(*,*) 'add  interval: xint, yint : ',
c     &                 k+1,xvec(k+1),yvec(k+1)
          nint=nint+1  
          if(nint.le.nmax) then
            goto 300 
          else
            write(*,*) 'DS: in dsfun_intpar exceeded the maximum dim'
            write(*,*) 'DS: allowed for vectors in the xvec,yvec block'
            write(*,*) 'DS: which is nmax = ',nmax
            stop
          endif
        endif  
 400    continue
      enddo
c 
c decide on which side start the integration
c
      kmax=-100
      maxy=-1.d10
      do k=0,nint
        if(maxy.lt.dabs(yvec(k))) then
          kmax=k
          maxy=dabs(yvec(k))
        endif
c        write(*,*) k,xvec(k),yvec(k)
      enddo
c
      if(kmax.le.nint-kmax) then
        k=0
        partial=0.d0
 20     xinf=xvec(k)
        yinf=yvec(k)
        if(k.eq.0) xinf=a
 10     xsup=xvec(k+1)
        ysup=yvec(k+1)
        if(k+1.eq.nint) xsup=b
        trap=0.5d0*(dabs(yinf)+dabs(ysup))*(xsup-xinf)
        if(trap.le.prec*1.d2*dabs(partial).and.k+1.lt.nint) then
          k=k+1
          goto 10
        endif
c
c you might be trying to integrate a fuction identically 0
c
        if(dabs(trap).eq.0.d0) then
          result=0.d0
          goto 30
        endif
        eps=prec*trap*1.d-2
        call dsfun_int(f1,xinf,xsup,eps,prec,result)
 30     partial=partial+result
c        write(*,*) k,xinf,xsup,result,partial
        if(k+1.lt.nint) then
          k=k+1
          goto 20
        endif
      else
        k=nint
        partial=0.d0
 21     xsup=xvec(k)
        ysup=yvec(k)
        if(k.eq.nint) xsup=b
 11     xinf=xvec(k-1)
        yinf=yvec(k-1)
        if(k-1.eq.0) xinf=a
        trap=0.5d0*(dabs(yinf)+dabs(ysup))*(xsup-xinf)
        if(trap.le.1.d-3*dabs(partial).and.k-1.gt.0) then
          k=k-1
          goto 11
        endif
c
c you might be trying to integrate a fuction identically 0
c
        if(dabs(trap).eq.0.d0) then
          result=0.d0
          goto 31
        endif
        eps=prec*trap*1.d-2
        call dsfun_int(f1,xinf,xsup,eps,prec,result)
 31     partial=partial+result
c        write(*,*) k,xinf,xsup,result,partial
        if(k-1.gt.0) then
          k=k-1
          goto 21
        endif
      endif
      res=partial
      return
      end



      subroutine dsfun_intparb(f1,f2,a,b,tollf2,tollab,nmax,prec,how
     &  ,res)
c_______________________________________________________________________
c  integrate function f1 between a and b breaking the integration 
c  interval into subintervals (ap,bp) such that 
c    max(f2(ap),f2(bp)) < tollf2 * min(f2(ap),f2(bp))
c  input:
c    integration function f1
c    weight function f2
c    integration limits a and b
c    tollerance on weight function f2, tollf2
c    tollerance on how small the interval (ap,bp) can be tollab
c    maximum number of intervals nmax 
c    precision on the sub-integral prec
c    how = 1: add subintervals in linear scale, 
c    how = 2: add subintervals in log scale
c  integration performed with the subroutine dsfun_intb
c  output:
c    result of integration res
c_______________________________________________________________________
      implicit none
ccc
      real*8 f1,f2,a,b,tollf2,tollab,prec,res
      integer nmax,how
      external f1,f2
ccc
      real*8 xvec(0:nmax),yvec(0:nmax),miny,maxy,diffx,xint,yint,
     &  partial,xinf,xsup,yinf,ysup,trap,eps,result
      logical dsisnanorinf,addpoint
      integer i,step,nint,k,kk,kmax
c
c check the flag how, tollf2, tollab
c
      if(.not.(how.eq.1.or.how.eq.2)) then
        write(*,*) 'DS: dsfun_intparb called with wrong option how : ',
     &              how
        stop 
      endif
      if(tollf2.le.1.d0) then
        write(*,*) 'DS: in dsfun_intparb tollf2 = ',tollf2,
     &             ' not properly intialized'
        stop 
      endif
      if(how.eq.1.and.tollab*1.d3*nmax.lt.(b-a)) then
      write(*,*) 'DS: in dsfun_intparb tollab, ', 
     &           'nmax not properly intialized'
      write(*,*) 'DS: b-a, tollab, nmax : ',b-a,tollab,nmax 
      stop
      endif
      if(how.eq.2.and.tollab*1.d3*nmax.lt.(dlog(b)-dlog(a))) then
      write(*,*) 'DS: in dsfun_intparb tollab, ', 
     &           'nmax not properly intialized'
      write(*,*) 'DS: dlog(b)-dlog(a), tollab, nmax : ',
     &             dlog(b)-dlog(a),tollab,nmax 
      stop
      endif
c
c start with two intervals
c
      xvec(0)=a
      if(how.eq.1) then
        xvec(1)=0.5d0*(a+b)
      else ! if(how.eq.2) 
        xvec(1)=dexp(0.5d0*(dlog(a)+dlog(b)))
      endif
      xvec(2)=b
      do i=0,2
        yvec(i)=f2(xvec(i))
c        write(*,*) 'starting : ',xvec(i),yvec(i) 
      enddo
      nint=2
c
c f2 may be singular in one or both extremes of integration
c
      step=0
 100  continue
      if(dsisnanorinf(yvec(0))) then
        step=step+1
        if(how.eq.1) then
          xvec(0)=xvec(0)+step*tollab
        else ! if(how.eq.2) 
          xvec(0)=dexp(dlog(xvec(0))+step*tollab)
        endif
        if(step.le.3.and.xvec(0).lt.xvec(1)) then
          yvec(0)=f2(xvec(0))
        else
          write(*,*) 'DS: dsfun_intparb: f2(a) is : ',yvec(0)
          stop 
        endif
        goto 100
      endif
      step=0
 200  continue
      if(dsisnanorinf(yvec(2))) then
        step=step+1
        if(how.eq.1) then
          xvec(2)=xvec(2)-step*tollab
        else !if(how.eq.2) 
          xvec(2)=dexp(dlog(xvec(2))-step*tollab)
        endif
        if(step.le.3.and.xvec(2).gt.xvec(1)) then
          yvec(2)=f2(xvec(2))
        else
          write(*,*) 'DS: dsfun_intparb: f2(b) is : ',yvec(2)
          stop 
        endif
        goto 200
      endif
c
c if f2 is singular in (a+b)/2 or exp((log(a)+log(b))/2) there is 
c   a problem
c
      if(dsisnanorinf(yvec(1))) then
        write(*,*) 'DS: dsfun_intparb: f2 is : ',yvec(1),
     &             ' at x = ',xvec(1)
        write(*,*) 'DS: within the interval of integration: ',a,b
        stop 
      endif
c
c fill in vectors
c
 300  continue
      do k=0,nint-1
        miny=min(dabs(yvec(k)),dabs(yvec(k+1)))
        maxy=max(dabs(yvec(k)),dabs(yvec(k+1)))
        addpoint=.false.
        if(maxy.gt.miny*tollf2) addpoint=.true.
        if(how.eq.1) then
          diffx=dabs(xvec(k)-xvec(k+1))
        else !if(how.eq.2) 
          diffx=dabs(dlog(xvec(k))-dlog(xvec(k+1)))
        endif
        if(diffx.lt.tollab) addpoint=.false.
        if(addpoint) then
          if(how.eq.1) then
            xint=(xvec(k+1)+xvec(k))/2.d0
          else ! if(how.eq.2) 
            xint=dexp((dlog(xvec(k+1))+dlog(xvec(k)))/2.d0)
          endif
          yint=f2(xint)
c
c if you hit a singularity skip the point
c
          if(dsisnanorinf(yint)) goto 400
          do kk=nint,k+1,-1
            xvec(kk+1)=xvec(kk)
            yvec(kk+1)=yvec(kk)
          enddo
          xvec(k+1)=xint
          yvec(k+1)=yint
c          write(*,*) 'add  interval: xint, yint : ',
c     &                 k+1,xvec(k+1),yvec(k+1)
          nint=nint+1  
          if(nint.le.nmax) then
            goto 300 
          else
            write(*,*) 'DS: in dsfun_intparb exceeded the maximum dim'
            write(*,*) 'DS: allowed for vectors in the xvec,yvec block'
            write(*,*) 'DS: which is nmax = ',nmax
            stop
          endif
        endif  
 400    continue
      enddo
c 
c decide on which side start the integration
c
      kmax=-100
      maxy=-1.d10
      do k=0,nint
        if(maxy.lt.dabs(yvec(k))) then
          kmax=k
          maxy=dabs(yvec(k))
        endif
c        write(*,*) k,xvec(k),yvec(k)
      enddo
c
      if(kmax.le.nint-kmax) then
        k=0
        partial=0.d0
 20     xinf=xvec(k)
        yinf=yvec(k)
        if(k.eq.0) xinf=a
 10     xsup=xvec(k+1)
        ysup=yvec(k+1)
        if(k+1.eq.nint) xsup=b
        trap=0.5d0*(dabs(yinf)+dabs(ysup))*(xsup-xinf)
        if(trap.le.prec*1.d2*dabs(partial).and.k+1.lt.nint) then
          k=k+1
          goto 10
        endif
c
c you might be trying to integrate a fuction identically 0
c
        if(dabs(trap).eq.0.d0) then
          result=0.d0
          goto 30
        endif
        eps=prec*trap*1.d-2
        call dsfun_intb(f1,xinf,xsup,eps,prec,result)
 30     partial=partial+result
c        write(*,*) 'k, xinf, xsup partial, result',k,xinf,xsup,
c     &     partial,result
        if(k+1.lt.nint) then
          k=k+1
          goto 20
        endif
      else
        k=nint
        partial=0.d0
 21     xsup=xvec(k)
        ysup=yvec(k)
        if(k.eq.nint) xsup=b
 11     xinf=xvec(k-1)
        yinf=yvec(k-1)
        if(k-1.eq.0) xinf=a
        trap=0.5d0*(dabs(yinf)+dabs(ysup))*(xsup-xinf)
        if(trap.le.1.d-3*dabs(partial).and.k-1.gt.0) then
          k=k-1
          goto 11
        endif
c
c you might be trying to integrate a fuction identically 0
c
        if(dabs(trap).eq.0.d0) then
          result=0.d0
          goto 31
        endif
        eps=prec*trap*1.d-2
        call dsfun_intb(f1,xinf,xsup,eps,prec,result)
 31     partial=partial+result
c        write(*,*) 'k, xinf, xsup partial, result',k,xinf,xsup,
c     &     partial,result
        if(k-1.gt.0) then
          k=k-1
          goto 21
        endif
      endif
      res=partial
      return
      end


