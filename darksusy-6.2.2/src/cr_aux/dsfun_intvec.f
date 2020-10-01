      subroutine dsfun_intvec(f1,f2,a,b,tollf2,tollab,nmax,prec,how,
     &  abvec,resvec,nk)
c_______________________________________________________________________
c  integral of function f1 between a and b + res1:
c     int_a^b dx f1 + res1 
c  breaking the integration interval into subintervals abvec such that 
c    max(f2(abvec(i)),f2(abvec(i+1))) 
c       < tollf2 * min(f2(abvec(i)),f2(abvec(i+1)))
c  i.e. in the same way as it is done in dsfun_intpar.
c  input:
c    integration function f1
c    weight function f2
c    integration limits a and b
c    tollerance on weight function f2, tollf2
c    tollerance on how small (abvec(i),abvec(i+1)) can be tollab
c    maximum number of intervals nmax 
c    precision on the sub-integral prec
c    how = 1: add subintervals in linear scale, 
c    how = 2: add subintervals in log scale
c  integration performed with the subroutine dsfun_int
c  output:
c    number of subintervals nk 
c    vector abvec such that:
c      abvec(0) = a
c      abvec(nk) = b
c    vector resvec such that:
c      resvec(0) = 0.d0
c      resvec(i) = integral of f1 in the interval (abvec(i-1),abvec(i)) 
c_______________________________________________________________________
      implicit none
      integer how,nk
      real*8 f2,a,b,tollf2,tollab,prec,abvec(0:1000),resvec(0:1000)
      external f1,f2
c
      real*8 xvec(0:1000),yvec(0:1000),miny,maxy,diffx,xint,yint,
     &  xinf,xsup,yinf,ysup,trap,eps,result
      logical dsisnanorinf,addpoint
      integer nmax,i,step,nint,k,kk
c
c check the flag how, tollf2, tollab
c
      if(.not.(how.eq.1.or.how.eq.2)) then
        write(*,*) 'DS: dsfun_intvec called with wrong option how : ',
     &              how
        stop 
      endif
      if(tollf2.le.1.d0) then
        write(*,*) 'DS: in dsfun_intvec tollf2 = ',tollf2,
     &             ' not properly intialized'
        stop 
      endif
      if(how.eq.1.and.tollab*1.d3*nmax.lt.(b-a)) then
      write(*,*) 'DS: in dsfun_intvec tollab, ', 
     &           'nmax not properly intialized'
      write(*,*) 'DS: b-a, tollab, nmax : ',b-a,tollab,nmax 
      stop
      endif
      if(how.eq.2.and.tollab*1.d3*nmax.lt.(dlog(b)-dlog(a))) then
      write(*,*) 'DS: in dsfun_intvec tollab, ', 
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
          write(*,*) 'DS: dsfun_intvec: f2(a) is : ',yvec(0)
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
          write(*,*) 'DS: dsfun_intvec: f2(b) is : ',yvec(2)
          stop 
        endif
        goto 200
      endif
c
c if f2 is singular in (a+b)/2 or exp((log(a)+log(b))/2) there is 
c   a problem
c
      if(dsisnanorinf(yvec(1))) then
        write(*,*) 'DS: dsfun_intvec: f2 is : ',yvec(1),
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
            write(*,*) 'DS: in dsfun_intvec exceeded the maximum dim'
            write(*,*) 'DS: allowed for vectors in the xvec,yvec block'
            write(*,*) 'DS: which is nmax = ',nmax
            stop
          endif
        endif  
 400    continue
      enddo
ccc
      abvec(0)=a
      resvec(0)=0.d0
      k=0
 20   xinf=xvec(k)
      yinf=yvec(k)
      if(k.eq.0) xinf=a
      xsup=xvec(k+1)
      ysup=yvec(k+1)
      if(k+1.eq.nint) xsup=b
      trap=0.5d0*(dabs(yinf)+dabs(ysup))*(xsup-xinf)
c
c you might be trying to integrate a fuction identically 0
c
      if(dabs(trap).eq.0.d0) then
        result=0.d0
        goto 30
      endif
      eps=prec*trap/10.d0
      call dsfun_int(f1,xinf,xsup,eps,prec,result)
 30   resvec(k+1)=result
      abvec(k+1)=xsup
      if(k+1.lt.nint) then
        k=k+1
        goto 20
      endif
      nk=nint
      return
      end


