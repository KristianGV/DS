      subroutine dsfun_int(f,a,b,eps,prec,res)
c_______________________________________________________________________
c  integrate function f between a and b using dqagse 
c  input:
c    integration function: f
c    integration limits: a and b
c    integration precision: eps (rechecked and reset internally)
c    relative precision needed in worst case: prec
c  output:
c    result of integration res
c_______________________________________________________________________
      implicit none
ccc
      real*8 f,a,b,eps,prec,res
      external f
ccc
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsrel,epsabs,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      logical dsisnanorinf 
ccc
      if(b.le.a.or.dabs(eps).lt.1.d-200) then
        res=0.d0
        return
      endif
ccc
      limit=5000
      epsabs=eps
      epsrel=1.d2*epsabs
 10   call dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval,ier,
     &            alist,blist,rlist,elist,iord,last)
      if((dabs(result).eq.0.d0.and.epsabs.le.0.d0).or.
     &   dsisnanorinf(result)) then
        write(*,*) 'DS: SERIUOS WARNING, DO NOT TRUST RESULTS !!!! '
        write(*,*) 'DS: problem in dsfun_int'
        write(*,*) 'DS: a,b,epsabs,eps,result : ',a,b,epsabs,eps,result
      endif
      if(dabs(result).eq.0.d0) then
        res=result
        return
      endif
      if(dabs(result)*prec.lt.epsabs) then
c        write(*,*) 'reload a, epsabs, epsrel',epsabs, epsrel,
c     &  dabs(result)*prec
        epsabs=1.d-2*dabs(result)*prec
        epsrel=1.d2*epsabs
        goto 10
      endif
      if(dabs(result)*prec.gt.1.d2*epsabs) then
        epsabs=dabs(result)*prec
      endif
      eps=epsabs
      res=result
      return
      end

      logical function dsisnanorinf(a)
c_______________________________________________________________________
c  function checking whether the argument a is +infinity or -infinity
c  or it is not a number
c_______________________________________________________________________
      implicit none
      real*8 a
      logical dsisnan,dsisinf
      dsisnanorinf=.false.
      if(dsisnan(a)) dsisnanorinf=.true.
      if(dsisinf(a)) dsisnanorinf=.true.
      end


      logical function dsisinf(a)
c_______________________________________________________________________
c  function checking whether the argument a is +infinity or -infinity
c_______________________________________________________________________
      implicit none
      real*8 a
      dsisinf=.false.
      if(a.eq.0.d0) return
      if(a.eq.2.d0*a) dsisinf=.true.
      end


      subroutine dsfun_intb(f,a,b,eps,prec,res)
c_______________________________________________________________________
c  integrate function f between a and b using dqagseb
c  input:
c    integration function: f
c    integration limits: a and b
c    integration precision: eps (rechecked and reset internally)
c    relative precision needed in worst case: prec
c  output:
c    result of integration res
c_______________________________________________________________________
      implicit none
ccc
      real*8 f,a,b,eps,prec,res
      external f
ccc
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsrel,epsabs,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      logical dsisnanorinf 
ccc
      if(b.le.a.or.dabs(eps).lt.1.d-200) then
        res=0.d0
        return
      endif
ccc
      limit=5000
      epsabs=eps
      epsrel=1.d2*epsabs
 10   call dqagseb(f,a,b,epsabs,epsrel,limit,result,abserr,neval,ier,
     &            alist,blist,rlist,elist,iord,last)
      if((dabs(result).eq.0.d0.and.epsabs.le.0.d0).or.
     &   dsisnanorinf(result)) then
        write(*,*) 'DS: SERIUOS WARNING, DO NOT TRUST RESULTS !!!! '
        write(*,*) 'DS: problem in dsfun_intb'
        write(*,*) 'DS: a,b,epsabs,eps,result : ',a,b,epsabs,eps,result
c        stop
      endif
      if(dabs(result).eq.0.d0) then
        res=result
        return
      endif
      if(dabs(result)*prec.lt.epsabs) then
c        write(*,*) 'reload b, epsabs, epsrel',epsabs, epsrel,
c     &  dabs(result)*prec
        epsabs=1.d-2*dabs(result)*prec
        epsrel=1.d2*epsabs
        goto 10
      endif
      if(dabs(result)*prec.gt.1.d2*epsabs) then
        epsabs=dabs(result)*prec
      endif
      eps=epsabs
      res=result
      return
      end
