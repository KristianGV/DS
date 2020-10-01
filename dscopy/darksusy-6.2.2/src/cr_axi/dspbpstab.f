      subroutine dspbpstab(tp,R,isin,tollfz,tollfr,tollfr2)
************************************************************************
*** tabulation of pbpoint function needed for axisymmetric source
************************************************************************
      implicit none
      include 'dspbcom.h'
      real*8 tp,R,tollfz,tollfr,tollfr2
      integer isin
ccc
      real*8 rcint,zcint
      common/pbcintrcom/rcint,zcint
ccc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
ccc
ccc common block in the tabulations dstabfun and dstabfunb
ccc
      real*8 rvec(1000),yvec(1000),yvec2(1000),rminl,rmaxl
      integer nintr,nhow,nrset
      common/dstabfuncom/rvec,yvec,yvec2,rminl,rmaxl,nintr,nhow,nrset
ccc
      real*8 rmin,rmax,tollf,tollr
      integer npmin,how
      external dspbpsGCfReq0
ccc
      real*8 zgrid(1000,3),rgrids(1000,1000,3),rgridgc(1000,1000,3),
     &  zlinear
      integer nzgrid,nrgrids(1000),nrgridgc(1000)
      common/dspbpstabcom/zgrid,rgrids,rgridgc,zlinear,nrgrids,
     &  nrgridgc,nzgrid
ccc
      real*8 tptab,Rtab,R0maxtab,z0maxtab
      integer iscall
      common/dspbpsGCitcallcom/tptab,Rtab,R0maxtab,z0maxtab,iscall
ccc
      integer nzadd,iadd,i,k,j
      real*8 zvadd(100),zadd,dspbpsGCfReq0
ccc
      real*8 nsfzrvec(1000),nsfzyvec(1000),nsfzyvec2(1000)
      integer nsfz
      common/Sfztabcom/nsfzrvec,nsfzyvec,nsfzyvec2,nsfz
ccc
      real*8 rtmp(1000),ytmp(1000),ytmp2(1000)
      external dspbpsSfz,dspbpsGCfzittab
ccc
      Rs=R
      tps=tp
      is=isin
c
c first set the vertical steps in the grid: 
c
      tollf=tollfz
      how=1
      rmin=0.d0
      rmax=min(0.99d0*pbhh,1.20d0*zcint)
      tollr=(rmax-rmin)/500.d0
      npmin=30
      call dstabfunb(dspbpsGCfReq0,rmin,rmax,tollf,tollr,npmin,how)
c
c extract the vertical steps from the common block
c
      nzgrid=nintr
      do i=1,nzgrid
        zgrid(i,1)=rvec(i)
        zgrid(i,2)=yvec(i)
c        write(*,*) 'z steps i,z0 : ',i,rvec(i)
      enddo
c
c if you are using pbpointe, there may be some problems in the
c tabulation around pbhg, add a few more points around there: 
c
      if(is.eq.2) then
c
c this is just an educated guess, change it if it does not work!!!
c
        nzadd=7
        zvadd(1)=pbhg-0.01d0
        zvadd(2)=pbhg
        zvadd(3)=pbhg+0.01d0
        zvadd(4)=pbhg+0.02d0
        zvadd(5)=pbhg+0.03d0
        zvadd(6)=pbhg+0.05d0
        zvadd(7)=pbhg+0.08d0
        zlinear=zvadd(6)
c
        if(zvadd(1).lt.zgrid(1,1)) then
          write(*,*) 'DS: in dspbpstab wrong hierarchy between'
          write(*,*) 'DS: zvadd(1) and zgrid(1) : ',zvadd(1),zgrid(1,1)
          stop
        endif
        do k=1,nzadd
          zadd=zvadd(k)
          do i=2,nzgrid
            if(zadd.gt.zgrid(i-1,1).and.zadd.lt.zgrid(i,1)) then
              iadd=i
              goto 10
            endif
          enddo
          write(*,*) 'DS: in dspbpstab error for zadd = ',zadd
          stop
 10       nzgrid=nzgrid+1
          do i=nzgrid,iadd+1,-1
            zgrid(i,1)=zgrid(i-1,1)
            zgrid(i,2)=zgrid(i-1,2)
          enddo
          zgrid(iadd,1)=zadd
          zgrid(iadd,2)=dspbpsGCfReq0(zadd)
        enddo
      endif
c
c for each vertical step generate the radial steps
c
      do i=1,nzgrid
        z0s=zgrid(i,1)
c        write(*,*) 'tab i,z0 : ',i,z0s
        tollf=tollfr
        how=1
        rmin=pbr0-1.30d0*rcint
        if(rmin.lt.0.1d0) then
          write(*,*) 'DS: in dspbpstab too large rcint = ',rcint
          write(*,*) 'DS: too small tabulation minimum  = ',rmin
          stop
        endif
        rmax=pbr0+1.30d0*rcint
        tollr=(rmax-rmin)/500.d0
        call dstabfun(dspbpsSfz,rmin,rmax,tollf,tollr,how)
c
c extract the radial steps from the common block
c
        nrgrids(i)=nintr
        do j=1,nrgrids(i)
          rgrids(i,j,1)=rvec(j)
          rgrids(i,j,2)=yvec(j)
          rgrids(i,j,3)=yvec2(j)
        enddo
c
c tabulate the average over theta of dspbpsGCfRfz versus the radial 
c coordinate in the GC reference frame:
c
c load first the tabulation for dspbpsSfztab
c
        z0s=zgrid(i,1)
        nsfz=nrgrids(i)
        do j=1,nrgrids(i)
          nsfzrvec(j)=rgrids(i,j,1)
          nsfzyvec(j)=rgrids(i,j,2)
          nsfzyvec2(j)=rgrids(i,j,3)
        enddo
        tollf=tollfr2
        how=1
        rmin=0.01d0
        rmax=1.20d0*rcint
        tollr=(rmax-rmin)/500.d0
        call dstabfun(dspbpsGCfzittab,rmin,rmax,tollf,tollr,how)
c
c extract the radial steps from the common block and make a tabulation
c as if negative radii have to be included as well
c
        do j=1,nintr
          rtmp(j+nintr+1)=rvec(j)
          ytmp(j+nintr+1)=yvec(j)
          rtmp(nintr-j+1)=-rvec(j)
          ytmp(nintr-j+1)=yvec(j)
        enddo
        rtmp(nintr+1)=0.d0
        ytmp(nintr+1)=dexp(zgrid(i,2))
        nrgridgc(i)=2*nintr+1
        call dsspline(rtmp,ytmp,nrgridgc(i),1.d30,1.d30,ytmp2)
        do j=1,nrgridgc(i)
          rgridgc(i,j,1)=rtmp(j)
          rgridgc(i,j,2)=ytmp(j)
          rgridgc(i,j,3)=ytmp2(j)
        enddo
      enddo
      R0maxtab=rgridgc(1,nrgridgc(1),1)
      z0maxtab=zgrid(nzgrid,1)
      tptab=tp
      Rtab=R
      iscall=isin
      return
      end


      real*8 function dspbpsStab(L,z0)
************************************************************************
*** tabulated version of dspbpsS(L,z0,tp) with tp=tps
************************************************************************
      implicit none
      real*8 L,z0
ccc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
ccc
      real*8 zgrid(1000,3),rgrids(1000,1000,3),rgridgc(1000,1000,3),
     &  zlinear
      integer nzgrid,nrgrids(1000),nrgridgc(1000)
      common/dspbpstabcom/zgrid,rgrids,rgridgc,zlinear,nrgrids,
     &  nrgridgc,nzgrid
ccc
      integer i,j,ntmp
      real*8 xtmp(1000),ytmp(1000),ytmp2(1000),xxtmp(2000),yytmp(1000),
     &  yytmp2(1000),result
ccc
      do i=1,nzgrid
        xxtmp(i+nzgrid-1)=zgrid(i,1)
        do j=1,nrgrids(i)
          xtmp(j)=rgrids(i,j,1)
          ytmp(j)=rgrids(i,j,2)
          ytmp2(j)=rgrids(i,j,3)
        enddo
        call dssplint(xtmp,ytmp,ytmp,nrgrids(i),L,yytmp(i+nzgrid-1))
        if(i.gt.1) then
          xxtmp(nzgrid-i+1)=-xxtmp(i+nzgrid-1)
          yytmp(nzgrid-i+1)=yytmp(i+nzgrid-1)
        endif
      enddo
      ntmp=2*nzgrid-1
      if(is.eq.2.and.dabs(z0).gt.zlinear) then
        call dslinint(xxtmp,yytmp,ntmp,z0,result)
      else
        call dsspline(xxtmp,yytmp,ntmp,1.d30,1.d30,yytmp2)
        call dssplint(xxtmp,yytmp,yytmp2,ntmp,z0,result)
      endif
      dspbpsStab=result
      return
      end


      real*8 function dspbpsGCitcall(R0,z0,tp,R,isin)
************************************************************************
*** this function is linking to dspbpsGCittab, loading the tabulation
*** in case it is needed
************************************************************************
      implicit none
      real*8 R0,z0,tp,R
      integer isin
      logical tabcall
      real*8 tollfz,tollfr,tollfr2,dspbpsGCittab
ccc
      real*8 tptab,Rtab,R0maxtab,z0maxtab
      integer iscall
      common/dspbpsGCitcallcom/tptab,Rtab,R0maxtab,z0maxtab,iscall
ccc
      tabcall=.false.
      if(R0.gt.R0maxtab) tabcall=.true.
      if(z0.gt.z0maxtab) tabcall=.true.
      if(dabs(tp-tptab).gt.1.d-3*dabs(tp+tptab)) tabcall=.true.
      if(dabs(R-Rtab).gt.1.d-3*dabs(R+Rtab)) tabcall=.true.
      if(isin.ne.iscall) tabcall=.true.
      if(tabcall) then
ccc
ccc this should be good enough, cross check
ccc
        tollfz=1.01d0
        tollfr=1.01d0
        tollfr2=1.005d0
        call dspbpstab(tp,R,isin,tollfz,tollfr,tollfr2)
      endif
      dspbpsGCitcall=dspbpsGCittab(R0,z0)
      return
      end



      real*8 function dspbpsGCittab(R0,z0)
************************************************************************
*** this function is the tabulated version of dspbpsGCit
************************************************************************
      implicit none
      Real*8 R0,z0
ccc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
ccc
      real*8 zgrid(1000,3),rgrids(1000,1000,3),rgridgc(1000,1000,3),
     &  zlinear
      integer nzgrid,nrgrids(1000),nrgridgc(1000)
      common/dspbpstabcom/zgrid,rgrids,rgridgc,zlinear,nrgrids,
     &  nrgridgc,nzgrid
ccc
      integer i,j,ntmp
      real*8 xtmp(1000),ytmp(1000),ytmp2(1000),xxtmp(2000),yytmp(1000),
     &  yytmp2(1000),result
ccc
      do i=1,nzgrid
        xxtmp(i+nzgrid-1)=zgrid(i,1)
        do j=1,nrgridgc(i)
          xtmp(j)=rgridgc(i,j,1)
          ytmp(j)=rgridgc(i,j,2)
          ytmp2(j)=rgridgc(i,j,3)
        enddo
        call dssplint(xtmp,ytmp,ytmp,nrgridgc(i),R0,yytmp(i+nzgrid-1))
        if(i.gt.1) then
          xxtmp(nzgrid-i+1)=-xxtmp(i+nzgrid-1)
          yytmp(nzgrid-i+1)=yytmp(i+nzgrid-1)
        endif
      enddo
      ntmp=2*nzgrid-1
      if(is.eq.2.and.dabs(z0).gt.zlinear) then
        call dslinint(xxtmp,yytmp,ntmp,z0,result)
      else
        call dsspline(xxtmp,yytmp,ntmp,1.d30,1.d30,yytmp2)
        call dssplint(xxtmp,yytmp,yytmp2,ntmp,z0,result)
      endif
      dspbpsGCittab=result
      return
      end


      real*8 function dspbpsSfztab(L)
************************************************************************
*** tabulated version of dspbpsSfz, spline tables are set in 
*** dspbpstab and there is no cross check here that a consistent 
*** tabulation is being used
************************************************************************
      implicit none
      real*8 L,result
ccc
      real*8 nsfzrvec(1000),nsfzyvec(1000),nsfzyvec2(1000)
      integer nsfz
      common/Sfztabcom/nsfzrvec,nsfzyvec,nsfzyvec2,nsfz
ccc
      call dssplint(nsfzrvec,nsfzyvec,nsfzyvec2,nsfz,L,result)
      dspbpsSfztab=result
      return
      end

 
      real*8 function dspbpsGCfRfztab(theta0)
************************************************************************
*** version of dspbpsGCfRfz linking to a tabulation rather than
*** dspbpss, spline tables are set in dspbpstab and there is no
*** cross check here that a consistent tabulation is being used
************************************************************************
      implicit none
      include 'dspbcom.h'
      real*8 theta0
      real*8 L,dspbpsSfztab
ccc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
ccc
      L=dsqrt(Rs**2+R0s**2-2.d0*Rs*R0s*dcos(theta0))
      dspbpsGCfRfztab=dspbpsSfztab(L)
      return
      end


      real*8 function dspbpsGCfzittab(R0)
************************************************************************
*** this function is the angular average of exp(dspbpsGCfRfztab) as
*** a function of the radial coordinate R0 in the cylindrical GC frame
*** spline tables are set in dspbpstab and there is no
*** cross check here that a consistent tabulation is being used
************************************************************************
      implicit none
      include 'dsmpconst.h' 
      real*8 R0
      real*8 low,up,eps,prec,result,dspbpsGCfzittabint
      external dspbpsGCfzittabint
ccc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
ccc
      R0s=R0
      low=0.d0
      up=pi
      prec=1.d-5
      eps=dabs(dspbpsGCfzittabint(up))*(up-low)*prec*1.d-2
      call dsfun_int(dspbpsGCfzittabint,low,up,eps,prec,result)
      dspbpsGCfzittab=result/up
      return
      end
ccc
      real*8 function dspbpsGCfzittabint(theta0)
      implicit none
      real*8 dspbpsGCfRfztab,theta0
      dspbpsGCfzittabint=dexp(dspbpsGCfRfztab(theta0))
      return
      end



      real*8 function dspbpsS(L,z0,tp,isin)
************************************************************************
*** this is dlog of dspbtdpsc assuming the observer, i.e.  
*** the sun is in the origin, i.e. (x,y,z)=(0,0,0), and for the 
*** position of the source specified by L = sqrt(x0**2+y0**2) and z0
***
***     L,z0 - position of the source (kpc)
***     tp - antiproton kinetic energy (gev)
***
*** output: in dlog(kpc^-3 10^15 s)
************************************************************************
      implicit none
      real*8 L,z0,tp
      integer isin
      real*8 x0,dspbtdpsc
      integer is
      common/pboptioncom/is
      is=isin
      x0=L
      dspbpsS=dlog(dspbtdpsc(0.d0,0.d0,0.d0,x0,0.d0,z0,tp,is))
      return
      end


      real*8 function dspbpsSfz(L)
************************************************************************
*** this is dlog of dspbtdpsc assuming the observer, i.e.  
*** the sun is in the origin, i.e. (x,y,z)=(0,0,0), and for the 
*** position of the source specified by L = sqrt(x0**2+y0**2) and
*** z0s=z0
***
***     L,z0s - position of the source (kpc)
***     tps - antiproton kinetic energy (gev)
***
*** output: in dlog(kpc^-3 10^15 s)
************************************************************************
      implicit none
      real*8 L
      real*8 x0,dspbtdpsc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
      x0=L
      dspbpsSfz=dlog(dspbtdpsc(0.d0,0.d0,0.d0,x0,0.d0,z0s,tps,is))
      return
      end


      real*8 function dspbpsGCfReq0(z0)
************************************************************************
*** this is dlog of dspbtdpsc assuming the observer, i.e.  
*** the sun is in the origin, i.e. (x,y,z)=(0,0,0), and for the 
*** position of the source specified by z0 and the radial coordinate
*** R0 in the cylindrical GC coordinate frame set to zero
***
***     z0 - position of the source (kpc)
***     tps - antiproton kinetic energy (gev)
***
*** output: in dlog(kpc^-3 10^15 s)
************************************************************************
      implicit none
      include 'dspbcom.h'
      real*8 z0
      real*8 x0,dspbtdpsc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
      x0=Rs
      dspbpsGCfReq0=dlog(dspbtdpsc(0.d0,0.d0,0.d0,x0,0.d0,z0,tps,is))
      return
      end



      real*8 function dspbpsGCfRfz(theta0)
************************************************************************
*** this is dlog of dspbtdpsc assuming the observer, i.e. the sun is in
*** the origin, i.e. (x,y,z)=(0,0,0), and for the 
*** position of the source specified by z0s and the radial coordinate
*** and angular, R0s and theta in the cylindrical GC coordinate frame
***
***     R0s,z0s,theta0 - position of the source (kpc,kpc,rad)
***     tps - antiproton kinetic energy (gev)
***
*** output: in dlog(kpc^-3 10^15 s)
************************************************************************
      implicit none
      include 'dspbcom.h'
      real*8 theta0
      real*8 x0,dspbtdpsc
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
      x0=dsqrt(Rs**2+R0s**2-2.d0*Rs*R0s*dcos(theta0))
      dspbpsGCfRfz=dlog(dspbtdpsc(0.d0,0.d0,0.d0,x0,0.d0,z0s,tps,is))
      return
      end


      real*8 function dspbpsGCit(R0,z0,tp,isin)
************************************************************************
*** this is the angular average of exp(dspbpsGCfRfz)
***
*** output: in kpc^-3 10^15 s
************************************************************************
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h' 
      real*8 R0,z0,tp
      integer isin
      real*8 tps,Rs,z0s,R0s
      common/dspbpsscom/tps,Rs,z0s,R0s
      integer is
      common/pboptioncom/is
      real*8 low,up,eps,prec,result,dspbpsGCitint
      external dspbpsGCitint
      tps=tp
      z0s=z0
      R0s=R0
      is=isin
      low=0.d0
      up=pi
      prec=1.d-5
      eps=dabs(dspbpsGCitint(up))*(up-low)*prec*1.d-2
      call dsfun_intb(dspbpsGCitint,low,up,eps,prec,result)
      dspbpsGCit=result/up
      return
      end
ccc
      real*8 function dspbpsGCitint(theta0)
      implicit none
      real*8 dspbpsGCfRfz,theta0
      dspbpsGCitint=dexp(dspbpsGCfRfz(theta0))
      return
      end




