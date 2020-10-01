      real*8 function vvar(pp,iv)
************************************************************************
*** input: pp - momentum (GeV)
***        iv = 1 - links to vvarnum
***        iv = 2 - links to vvarana
*** output: v = int_0^{u(p)} d\tilde{u} D(\tilde{u}) (kpc^2)
************************************************************************
      implicit none
      real*8 pp
      integer iv
      real*8 vvarnum,vvarana
ccc
      if(iv.eq.1) then
        vvar=vvarnum(pp) ! kpc^2
      elseif(iv.eq.2) then
        vvar=vvarana(pp) ! kpc^2
      else
        write(*,*) 'DS: in vvar invalid iv = ',iv
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end



      real*8 function vvarana(pp)
************************************************************************
*** input: pp - momentum (GeV)
*** output: v = int_0^{u(p)} d\tilde{u} D(\tilde{u}) (kpc^2)
***         analytical inversion in case of diffusion coefficient set
***         nkdiff = 1 or 4 & betalabel=.false., and momentum loss
***         rate being blossmean * pp**2 
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h'
      real*8 pp ! momentum in GeV
      real*8 uu,uvarana
ccc
      uu=uvarana(pp) ! 10^15 s
      if(nkdiff.eq.1) then
        vvarana=k0halo*(10.d0/blossmean/kdiffrig0)**kdiffdelta
     &                   *uu**(1.d0-kdiffdelta)/(1.d0-kdiffdelta)
      elseif(nkdiff.eq.4) then
        vvarana=k0halo*(uu+(10.d0/blossmean/kdiffrig0)**kdiffdelta
     &                   *uu**(1.d0-kdiffdelta)/(1.d0-kdiffdelta))
                  ! 10^27 cm^2 s^-1 10^15 s = 10^42 cm^2
      else
        write(*,*) 'DS: vvarana called with invalid nkdiff = ',nkdiff
        write(*,*) 'DS: program stopped'
        stop
      endif
      vvarana=vvarana/kpc**2 ! kpc^2
      return
      end



      real*8 function vvarnum(pp)
************************************************************************
*** input: pp - momentum (GeV)
*** output: v = int_0^{u(p)} d\tilde{u} D(\tilde{u}) (kpc^2)
***             assuming fully general momentum loss rate and spatial
***             diffusion coefficient.
*** A numerical integral is needed for the inversion. On the first the
*** routine checks whether the tabulation can be loaded from the file
*** 'vofpnum.dat' 
*** The tabulation needs to be reloaded each time the functions for  
*** the diffusion coefficient or the energy loss rate are changed;
*** their value for pp = 1 GeV is checked in dsepdndpaxi
************************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 pp
ccc
      real*8 xpvarp(50001),ypvarp(50001),ypvarp2(50001)
     & ,xpvarpb(50001),ypvarpb(50001),ypvarpb2(50001)
      common/vofenumcom/xpvarp,ypvarp,ypvarp2,xpvarpb,ypvarpb,ypvarpb2
ccc
      real*8 bloss1i,diff1i,bloss1,diff1,
     &  bloss50i,diff50i,bloss50,diff50,
     &  pvarpoints,vofpvarmin,vofpvarmax
      integer vofeset
      common/vofenumstore/bloss1i,diff1i,bloss1,diff1,bloss50i,diff50i,
     &  bloss50,diff50,pvarpoints,vofpvarmin,vofpvarmax,vofeset
ccc
      integer k,j
      real*8 loclow,locup,par,logpvar,up,low,result,eps,prec
      external vvarnum_int
      character*100 vofefile
      logical getout
ccc
ccc the check on whether you need to reload tabulations has been moved 
ccc up to dsepdndpaxi
ccc
c      call vvarnumsetup
ccc
      if(vofeset.ne.123456) then
        call dsdatafile(vofefile,'vofpnum.dat')
        open(unit=13,file=vofefile,status='old',form='formatted',
     &       err=200)
        read(13,1000,err=200,end=200) bloss1i,diff1i,bloss50i,diff50i
        vofeset=123456
        call vvarnumsetup
        if(vofeset.eq.0) goto 200
        vofeset=0
        read(13,1000,err=200,end=200) pvarpoints,vofpvarmin,vofpvarmax
        do k=1,int(pvarpoints)
          read(13,1000,err=200,end=200) xpvarp(k),ypvarp(k)
     &          ,xpvarpb(k),ypvarpb(k)
          if(k.gt.1.and.(dabs(xpvarp(k)-xpvarp(k-1)).lt.1.d-16
     &       .or.dabs(xpvarpb(k)-xpvarpb(k-1)).lt.1.d-16)) then
            write(*,*) 'DS: the file vofe has two equal x entries : '
            write(*,*) 'DS: k, xp_{k-1} xp_k : ',k,xpvarp(k-1),xpvarp(k)
            write(*,*) 'DS: k, xpb_{k-1} xpb_k : ',k,xpvarpb(k-1)
     &           ,xpvarpb(k)
            stop
          endif
        enddo
 1000   format(60(1x,e18.12))
        close(13)
        call dsspline(xpvarp,ypvarp,int(pvarpoints),1.d31,1.d31,ypvarp2)
        call dsspline(xpvarpb,ypvarpb,int(pvarpoints),1.d31,1.d31,
     &              ypvarpb2)
        vofeset=123456
        goto 100
ccc
 200    close(13)
        write(*,*) 'DS: The vofe file : ',vofefile
        write(*,*) 'DS: does not exist or is not in the required format'
        write(*,*) 'DS: (re)generate it'
        write(*,*) 'DS: vofe tabulation started'
        call vvarnumsetup
        bloss1i=bloss1
        diff1i=diff1
        bloss50i=bloss50
        diff50i=diff50
        pvarpoints=2000
        vofpvarmin=0.001d0
        vofpvarmax=10000.d0
        do k=1,int(pvarpoints)-2
          logpvar=dlog(vofpvarmin)+(dlog(vofpvarmax)-dlog(vofpvarmin))
     &         /(pvarpoints-3.d0)*(k-1)
          low=dexp(logpvar)
          up=vofpvarmax+500.d0
          par=0.d0
          getout=.false.
          do j=1,1000
            loclow=low*10.d0**(j-1)
            locup=low*10.d0**j
            if(locup.gt.up) then 
              locup=up
              getout=.true. 
            endif
            eps=1.d-5
            prec=1.d-5
            call dsfun_int(vvarnum_int,loclow,locup,eps,prec,result)
            result=result ! 10^43 cm^2
            result=result*10.d0/kpc**2 ! kpc^2
            par=par+result
c            write(*,*) k,j,loclow,locup,result,par
            if(getout) goto 111
          enddo
 111      xpvarp(k+1)=logpvar
          ypvarp(k+1)=dlog(par)
c          write (*,*) k,logpvar,ypvarp(k+1),par
        enddo
        write(*,*) 'DS: the vofe tabulation is over'
        xpvarp(1)=xpvarp(2)-(xpvarp(3)-xpvarp(2))/1.d3
        ypvarp(1)=ypvarp(2)
        xpvarp(int(pvarpoints))=xpvarp(int(pvarpoints)-1)
     &    +(xpvarp(int(pvarpoints)-1)
     &      -xpvarp(int(pvarpoints)-2))/1.d3
        ypvarp(int(pvarpoints))=ypvarp(int(pvarpoints)-1)
c        write(*,*) int(pvarpoints),xpvarp(1)
c     &     ,xpvarp(int(pvarpoints))
        call dsspline(xpvarp,ypvarp,int(pvarpoints),1.d31,1.d31,ypvarp2)
        do k=2,int(pvarpoints)-1
          xpvarpb(k)=ypvarp(int(pvarpoints)-k+1)
          ypvarpb(k)=xpvarp(int(pvarpoints)-k+1)
        enddo
        xpvarpb(1)=xpvarpb(2)-(xpvarpb(3)-xpvarpb(2))/1.d3
        ypvarpb(1)=ypvarpb(2)
        xpvarpb(int(pvarpoints))=xpvarpb(int(pvarpoints)-1)
     &    +(xpvarpb(int(pvarpoints)-1)
     &      -xpvarpb(int(pvarpoints)-2))/1.d3
        ypvarpb(int(pvarpoints))=ypvarpb(int(pvarpoints)-1)
        call dsspline(xpvarpb,ypvarpb,int(pvarpoints),1.d31
     &    ,1.d31,ypvarpb2)
        open(unit=13,file=vofefile,status='unknown',form='formatted')
        write(13,1000) bloss1i,diff1i,bloss50i,diff50i
        write(13,1000) pvarpoints,vofpvarmin,vofpvarmax
        do k=1,int(pvarpoints)
          write(13,1000) xpvarp(k),ypvarp(k)
     &      ,xpvarpb(k),ypvarpb(k)
        enddo
        close(13)
        vofeset=123456
        goto 100
      endif
ccc
 100  continue
ccc
ccc actual look up in the table
ccc
      if(dlog(pp).ge.xpvarp(1).and.dlog(pp)
     &   .le.xpvarp(int(pvarpoints))) then
        call dssplint(xpvarp,ypvarp,ypvarp2,int(pvarpoints),dlog(pp),
     &              result)
        vvarnum=dexp(result)
      else
        write(*,*) 'DS: vvarnum called out of the defined range
     &              , pp =',pp
        write(*,*) 'DS: pvarmin, pvarmax = ',dexp(xpvarp(1))
     &       ,dexp(xpvarp(int(pvarpoints)))
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end
ccc
ccc
      real*8 function vvarnum_int(pp)
      implicit none
      real*8 pp,dskdiff,dseppdotm
      vvarnum_int=dskdiff(pp,1,1.d0)/dseppdotm(pp,1) 
           ! 10^27 cm^2 s^-1 * 10^16 GeV^-1 s = 10^43 cm^2 GeV^-1
      return
      end



      subroutine vvarnumsetup
************************************************************************
*** This subroutine checks whether the tabulation of vvarnum is loaded
*** with the currently implemented diffusion coefficient and energy
*** loss rate, it does it loading their values at pp = 1 GeV and 
*** checking against stored coefficients
************************************************************************
      implicit none
      real*8 dseppdotm,dskdiff,diff,sum
ccc
      real*8 bloss1i,diff1i,bloss1,diff1,
     &  bloss50i,diff50i,bloss50,diff50,
     &  pvarpoints,vofpvarmin,vofpvarmax
      integer vofeset
      common/vofenumstore/bloss1i,diff1i,bloss1,diff1,bloss50i,diff50i,
     &  bloss50,diff50,pvarpoints,vofpvarmin,vofpvarmax,vofeset
ccc
      bloss1=dseppdotm(1.d0,1)
      diff=dabs(bloss1i-bloss1)
      sum=dabs(bloss1i+bloss1)
      if(diff.gt.1.d-10*sum) vofeset=0 
      diff1=dskdiff(1.d0,1,0.5d0)
      diff=dabs(diff1i-diff1)
      sum=dabs(diff1i+diff1)
      if(diff.gt.1.d-10*sum) vofeset=0 
      bloss50=dseppdotm(50.d0,1)
      diff=dabs(bloss50i-bloss50)
      sum=dabs(bloss50i+bloss50)
      if(diff.gt.1.d-10*sum) vofeset=0 
      diff50=dskdiff(50.d0,1,0.5d0)
      diff=dabs(diff50i-diff50)
      sum=dabs(diff50i+diff50)
      if(diff.gt.1.d-10*sum) vofeset=0 
ccc
      return
      end


      real*8 function pvarnum(v)
************************************************************************
*** input: v = int_0^{u(p)} d\tilde{u} D(\tilde{u}) (kpc^2)
*** output: p = momentum (GeV)
************************************************************************
      implicit none
ccc
      real*8 xpvarp(50001),ypvarp(50001),ypvarp2(50001)
     & ,xpvarpb(50001),ypvarpb(50001),ypvarpb2(50001)
      common/vofenumcom/xpvarp,ypvarp,ypvarp2,xpvarpb,ypvarpb,ypvarpb2
ccc
      real*8 bloss1i,diff1i,bloss1,diff1,
     &  bloss50i,diff50i,bloss50,diff50,
     &  pvarpoints,vofpvarmin,vofpvarmax
      integer vofeset
      common/vofenumstore/bloss1i,diff1i,bloss1,diff1,bloss50i,diff50i,
     &  bloss50,diff50,pvarpoints,vofpvarmin,vofpvarmax,vofeset
ccc
      real*8 v,result,dummy,vvarnum
ccc
      if(vofeset.ne.123456) then
        dummy=vvarnum(1.d0)
      endif
      if(dlog(v).ge.xpvarpb(1).and.dlog(v)
     &   .le.xpvarpb(int(pvarpoints))) then
        call dssplint(xpvarpb,ypvarpb,ypvarpb2,int(pvarpoints)
     &     ,dlog(v),result)
        pvarnum=dexp(result)
      else
        write(*,*) 'DS: pvarnum called out of the defined range, v =',v
        write(*,*) 'DS: vmin, vmax = ',dexp(xpvarpb(1))
     &       ,dexp(xpvarpb(int(pvarpoints)))
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end


      real*8 function uvar(pp,iv)
************************************************************************
*** input: pp - momentum (GeV)
***        iv = 1 - links to uvarnum
***        iv = 2 - links to uvarana
*** output: u = int_p^\inf d\tilde{p} 1/[-pdot(\tilde{p})] (10^15 s)
************************************************************************
      implicit none
      real*8 pp
      integer iv
      real*8 uvarnum,uvarana
ccc
      if(iv.eq.1) then
        uvar=uvarnum(pp) ! 10^15 s
      elseif(iv.eq.2) then
        uvar=uvarana(pp) ! 10^15 s
      else
        write(*,*) 'DS: in uvar invalid iv = ',iv
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end


      real*8 function uvarana(pp)
************************************************************************
*** input: pp - momentum (GeV)
*** output: u = int_p^\inf d\tilde{p} 1/[-pdot(\tilde{p})] (10^15 s)
***             assuming the momentum loss rate scales with pp^2, i.e.
***             -pdot(pp) = blossmean *pp^2
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 pp ! momentum in GeV
ccc
      uvarana=1.d0/blossmean/pp*10.d0 ! 10^15 s
      return
      end



      real*8 function uvarnum(pp)
************************************************************************
*** input: pp - momentum (GeV)
*** output: u = int_p^\inf d\tilde{p} 1/[-pdot(\tilde{p})]
***             assuming fully general momentum loss rate
*** A numerical integral is needed for the inversion. On the first the
*** routine checks whether the tabulation can be loaded from the file
*** 'uofpnum.dat' 
************************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 pp
ccc
      real*8 xpvarp(50001),ypvarp(50001),ypvarp2(50001)
     & ,xpvarpb(50001),ypvarpb(50001),ypvarpb2(50001)
      common/uofenumcom/xpvarp,ypvarp,ypvarp2,xpvarpb,ypvarpb,ypvarpb2
ccc
      real*8 bloss1i,bloss1,bloss50i,bloss50,pvarpoints,uofpvarmin,
     &  uofpvarmax
      integer uofeset
      common/uofenumstore/bloss1i,bloss1,bloss50i,bloss50,pvarpoints,
     &  uofpvarmin,uofpvarmax,uofeset
ccc
      integer k,j
      real*8 loclow,locup,par,logpvar,up,low,result,eps,prec
      external uvarnum_int
      character*100 uofefile
      logical getout
ccc
ccc the check on whether you need to reload tabulations has been moved 
ccc up to dsepdndpaxi
ccc
c      call uvarnumsetup
ccc
      if(uofeset.ne.123456) then
        call dsdatafile(uofefile,'uofpnum.dat')
        open(unit=13,file=uofefile,status='old',form='formatted',
     &       err=200)
        read(13,1000,err=200,end=200) bloss1i,bloss50i
        uofeset=123456
        call uvarnumsetup
        if(uofeset.eq.0) goto 200
        uofeset=0
        read(13,1000,err=200,end=200) pvarpoints,uofpvarmin,uofpvarmax
        do k=1,int(pvarpoints)
          read(13,1000,err=200,end=200) xpvarp(k),ypvarp(k)
     &      ,xpvarpb(k),ypvarpb(k)
          if(k.gt.1.and.(dabs(xpvarp(k)-xpvarp(k-1)).lt.1.d-16
     &       .or.dabs(xpvarpb(k)-xpvarpb(k-1)).lt.1.d-16)) then
            write(*,*) 'DS: the file uofe has two equal x entries : '
            write(*,*) 'DS: k, xp_{k-1} xp_k : ',k,xpvarp(k-1),xpvarp(k)
            write(*,*) 'DS: k, xpb_{k-1} xpb_k : ',k,xpvarpb(k-1)
     &           ,xpvarpb(k)
            stop
          endif
        enddo
 1000   format(60(1x,e18.12))
        close(13)
        call dsspline(xpvarp,ypvarp,int(pvarpoints),1.d31,1.d31,ypvarp2)
        call dsspline(xpvarpb,ypvarpb,int(pvarpoints),1.d31,1.d31,
     &              ypvarpb2)
        uofeset=123456
        goto 100
ccc
 200    close(13)
        write(*,*) 'DS: The uofe file : ',uofefile
        write(*,*) 'DS: does not exist or is not in the required format'
        write(*,*) 'DS: (re)generate it'
        write(*,*) 'DS: uofe tabulation started'
        call uvarnumsetup
        bloss1i=bloss1
        bloss50i=bloss50
        pvarpoints=2000
        uofpvarmin=0.001d0
        uofpvarmax=10000.d0
        do k=1,int(pvarpoints)-2
          logpvar=dlog(uofpvarmin)+(dlog(uofpvarmax)-dlog(uofpvarmin))
     &         /(pvarpoints-3.d0)*(k-1)
          low=dexp(logpvar)
          up=uofpvarmax+500.d0
          par=0.d0
          getout=.false.
          do j=1,1000
            loclow=low*10.d0**(j-1)
            locup=low*10.d0**j
            if(locup.gt.up) then 
              locup=up
              getout=.true. 
            endif
            eps=1.d-5
            prec=1.d-5
            call dsfun_int(uvarnum_int,loclow,locup,eps,prec,result)
            result=result ! 10^15 s
            par=par+result
c            write(*,*) k,j,loclow,locup,result,par
            if(getout) goto 111
          enddo
 111      xpvarp(k+1)=logpvar
          ypvarp(k+1)=dlog(par)
c          write (*,*) k,logpvar,ypvarp(k+1),par
        enddo
        write(*,*) 'DS: the uofe tabulation is over'
        xpvarp(1)=xpvarp(2)-(xpvarp(3)-xpvarp(2))/1.d3
        ypvarp(1)=ypvarp(2)
        xpvarp(int(pvarpoints))=xpvarp(int(pvarpoints)-1)
     &    +(xpvarp(int(pvarpoints)-1)
     &      -xpvarp(int(pvarpoints)-2))/1.d3
        ypvarp(int(pvarpoints))=ypvarp(int(pvarpoints)-1)
c        write(*,*) int(pvarpoints),xpvarp(1)
c     &     ,xpvarp(int(pvarpoints))
        call dsspline(xpvarp,ypvarp,int(pvarpoints),1.d31,1.d31,ypvarp2)
        do k=2,int(pvarpoints)-1
          xpvarpb(k)=ypvarp(int(pvarpoints)-k+1)
          ypvarpb(k)=xpvarp(int(pvarpoints)-k+1)
        enddo
        xpvarpb(1)=xpvarpb(2)-(xpvarpb(3)-xpvarpb(2))/1.d3
        ypvarpb(1)=ypvarpb(2)
        xpvarpb(int(pvarpoints))=xpvarpb(int(pvarpoints)-1)
     &    +(xpvarpb(int(pvarpoints)-1)
     &      -xpvarpb(int(pvarpoints)-2))/1.d3
        ypvarpb(int(pvarpoints))=ypvarpb(int(pvarpoints)-1)
        call dsspline(xpvarpb,ypvarpb,int(pvarpoints),1.d31
     &    ,1.d31,ypvarpb2)
        open(unit=13,file=uofefile,status='unknown',form='formatted')
        write(13,1000) bloss1i,bloss50i
        write(13,1000) pvarpoints,uofpvarmin,uofpvarmax
        do k=1,int(pvarpoints)
          write(13,1000) xpvarp(k),ypvarp(k)
     &      ,xpvarpb(k),ypvarpb(k)
        enddo
        close(13)
        uofeset=123456
        goto 100
      endif
ccc
 100  continue
ccc
ccc actual look up in the table
ccc
      if(dlog(pp).ge.xpvarp(1).and.dlog(pp)
     &   .le.xpvarp(int(pvarpoints))) then
        call dssplint(xpvarp,ypvarp,ypvarp2,int(pvarpoints),dlog(pp),
     &              result)
        uvarnum=dexp(result)  ! 10^15 s
      else
        write(*,*) 'DS: uvarnum called out of the defined range
     &              , pp =',pp
        write(*,*) 'DS: pvarmin, pvarmax = ',dexp(xpvarp(1))
     &       ,dexp(xpvarp(int(pvarpoints)))
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end
ccc
ccc
      real*8 function uvarnum_int(pp)
      implicit none
      real*8 pp,dseppdotm
      uvarnum_int=1.d0/dseppdotm(pp,1)*10.d0 
           ! 10^15 GeV^-1 s
      return
      end



      real*8 function pofuvar(u,iv)
************************************************************************
*** input: u = int_p^\inf d\tilde{p} 1/[-pdot(\tilde{p})] (10^15 s)
***        iv = 1 - links to pofuvarnum
***        iv = 2 - links to pofuvarana
*** output: pp - momentum (GeV)
************************************************************************
      implicit none
      real*8 u
      integer iv
      real*8 pofuvarnum,pofuvarana
ccc
      if(iv.eq.1) then
        pofuvar=pofuvarnum(u) ! GeV
      elseif(iv.eq.2) then
        pofuvar=pofuvarana(u) ! GeV
      else
        write(*,*) 'DS: in uvar invalid iv = ',iv
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end



      real*8 function pofuvarana(uu)
************************************************************************
*** input: uu - u variable (10^15 s)
*** output: p = uvar^-1(uu) (GeV)
************************************************************************
      implicit none
      include 'dscraxicom.h'
      real*8 uu
ccc
      pofuvarana=1.d0/blossmean/uu*10.d0 ! GeV
      return
      end



      real*8 function pofuvarnum(u)
************************************************************************
*** input: v = int_0^{u(p)} d\tilde{u} D(\tilde{u}) (kpc^2)
*** output: p = momentum (GeV)
************************************************************************
      implicit none
ccc
      real*8 xpvarp(50001),ypvarp(50001),ypvarp2(50001)
     & ,xpvarpb(50001),ypvarpb(50001),ypvarpb2(50001)
      common/uofenumcom/xpvarp,ypvarp,ypvarp2,xpvarpb,ypvarpb,ypvarpb2
ccc
      real*8 bloss1i,bloss1,bloss50i,bloss50,pvarpoints,uofpvarmin,
     &  uofpvarmax
      integer uofeset
      common/uofenumstore/bloss1i,bloss1,bloss50i,bloss50,pvarpoints,
     &  uofpvarmin,uofpvarmax,uofeset
ccc
      real*8 u,result,dummy,uvarnum
ccc
ccc the check on whether you need to reload tabulations has been moved 
ccc up to dsepdndpaxi
ccc
c      call uvarnumsetup
ccc
      if(uofeset.ne.123456) then
        dummy=uvarnum(1.d0)
      endif
      if(dlog(u).ge.xpvarpb(1).and.dlog(u)
     &   .le.xpvarpb(int(pvarpoints))) then
        call dssplint(xpvarpb,ypvarpb,ypvarpb2,int(pvarpoints)
     &     ,dlog(u),result)
        pofuvarnum=dexp(result)
      else
        write(*,*) 'DS: pofuvarnum called out of range, u =',u
        write(*,*) 'DS: umin, umax = ',dexp(xpvarpb(1))
     &       ,dexp(xpvarpb(int(pvarpoints)))
        write(*,*) 'DS: program stopped'
        stop
      endif
      return
      end


      subroutine uvarnumsetup
************************************************************************
*** This subroutine checks whether the tabulation of uvarnum is loaded
*** with the currently implemented diffusion coefficient and energy
*** loss rate, it does it loading their values at pp = 1 GeV and 
*** checking against stored coefficients
************************************************************************
      implicit none
      real*8 dseppdotm,diff,sum
ccc
      real*8 bloss1i,bloss1,bloss50i,bloss50,pvarpoints,uofpvarmin,
     &  uofpvarmax
      integer uofeset
      common/uofenumstore/bloss1i,bloss1,bloss50i,bloss50,pvarpoints,
     &  uofpvarmin,uofpvarmax,uofeset
ccc
      bloss1=dseppdotm(1.d0,1)
      diff=dabs(bloss1i-bloss1)
      sum=dabs(bloss1i+bloss1)
      if(diff.gt.1.d-10*sum) uofeset=0
      bloss50=dseppdotm(50.d0,1)
      diff=dabs(bloss50i-bloss50)
      sum=dabs(bloss50i+bloss50)
      if(diff.gt.1.d-10*sum) uofeset=0 
ccc
      return
      end
