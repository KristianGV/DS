      real*8 function dspbtdps(x,y,z,x0,y0,z0,tp,isin)
************************************************************************
*** antiproton "confinement time per annihilation volume" for a static
*** point source within which WIMP pair annihilations take place.
***
*** inputs:
***   x,y,z - position of the observer (kpc)
***   x0,y0,z0 - position of the source (kpc)
***   tp - antiproton kinetic energy (gev)
***   isin = propagation model option : isin=1 -> model with delta 
***     function approximation for the gas disc; isin=2 -> two-zone
***     model with finite thickness for the gas disk; in both cases 
***     the radial boundary conditions are neglected.
***
*** output: in kpc^-3 10^15 s
***
*** the antiproton flux from wimp pair annihilation is obtained by
*** multiplying the result of this function by the factor:
***   units * 1/(4 pi) * vp * rho2vol * (sigma v)/(2*mwimp**2) * dn/dtp
*** where:
***   vp = vp(tp) = antiproton velocity for given antiproton kinetic
***      energy tp in unit of c 
***   rho2vol = int dr^3_s \rho_s^2(r_s) = integral over the source
***      volume of the wimp density squared in kpc^3 GeV^2 cm^-6
***   sigma v = dssigmav0tot in cm^3 s^-1
***   mwimp = WIMP mass in GeV
***   dn/dtp = dn/dtp(tp) = antiproton spectrum in GeV^-1
***   units = factor taking into account units =
***    1.d15 s sr^-1 2.99792d10 cm s^-1 GeV^2 cm^-6 cm^3 s^-1 GeV^-3
***      = 2.99792d25 cm^-2 s^-1 GeV^-1 sr^-1
*** 
*** for dark matter particle decays, the antiproton flux is obtained 
*** by multiplying the result of this function by the factor:
***   units * 1/(4 pi) * vd * rhovol * dec-rate/mwimp * dN/dtp
*** having defined (as opposed to the previous case):
***   rho2vol = int dr^3_s \rho_s(r_s) = integral over the source
***      volume of the source density in kpc^3 GeV cm^-3
***   dec-rate = dsdecratewimp = dark matter particle decay rate in s^-1
***   mwimp = decaying particle mass in GeV
***   units = factor taking into account units =
***    1.d15 s sr^-1 2.99792d10 cm s^-1 GeV cm^-3 s^-1 GeV^-2
***      = 2.99792d25 cm^-2 s^-1 GeV^-1 sr^-1
***
************************************************************************
      implicit none
      include 'dsmpconst.h'
      real*8 x,y,z,x0,y0,z0,tp
      integer isin
      real*8 dspbtdpsc
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
ccc set the antiproton mass and atomic number
ccc
      mnuc=m_p
      anuc=1
      znuc=1
ccc
      dspbtdps=dspbtdpsc(x,y,z,x0,y0,z0,tp,isin)
      return
      end
ccc
ccc
ccc
      real*8 function dspbtdpsc(x,y,z,x0,y0,z0,tp,isin)
      implicit none
      real*8 x,y,z,x0,y0,z0,tp
      integer isin
      real*8 dspbpoint,dspbpointe,res
ccc
      if(isin.eq.1) then
        res=dspbpoint(x,y,z,x0,y0,z0,tp)
      elseif(isin.eq.2) then
        res=dspbpointe(x,y,z,x0,y0,z0,tp)
      else
        write(*,*) 'DS: wrong is in dspbtdpsc, is = ',isin
        write(*,*) 'DS: program stopped'
        stop
      endif
      dspbtdpsc=res
      return
      end



      real*8 function dspbpoint(x,y,z,x0,y0,z0,tp)
************************************************************************
*** inputs:
***     x,y,z - position of the observer (kpc)
***     x0,y0,z0 - position of the source (kpc)
***     tp - kinetic energy per nucleon (gev)
***
*** output: in kpc^-3 10^15 s
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h'
ccc
      real*8 vbar,pebar,dxy,z0i
      common/dspbpointintcom/vbar,pebar,dxy,z0i
ccc
      real*8 x,y,z,x0,y0,z0,tp
      logical conv
      real*8 pp,rig,ee,dpbar,dskdiff,beta,axsec,dspbsigmavpbar,low,
     &  up,result,result2,dsdbsigmavdbar
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
ccc check whether you are actually outside the diffusion region:
ccc      
      if(dabs(z0).ge.diffhh) then
        dspbpoint=0.d0
        return
      endif
ccc
ccc so far this function has been implemented for z=0 only
ccc
      if(dabs(z).gt.1.d-16) then
        write(*,*) 'DS: dspbpoint has not been implemented for z.ne.0'
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc set tp dependent variables:
ccc
      ee=anuc*tp+mnuc
      pp=dsqrt(dabs(ee**2-mnuc**2))
      rig=pp/dble(znuc)
      beta=pp/ee
      dpbar=dskdiff(rig,1,beta)
      vbar=diffcvel/dpbar    ! 10^5 cm/s * 10^-27 s/cm^-2 = 10^-22 cm^-1
      vbar=vbar*(kpc/10.d0)  ! kpc^-1
      if(anuc.eq.1) then
        axsec=dspbsigmavpbar(ee)  ! mb * 1.d10 cm*s^-1 
      elseif(anuc.eq.2) then
        axsec=dsdbsigmavdbar(ee)  ! mb * 1.d10 cm*s^-1
      else
        write(*,*) 'DS: call to dspbpoint without proper'
        write(*,*) 'DS: intialization of mnuc, anuc : ',mnuc,anuc
        write(*,*) 'DS: program stopped'
        stop
      endif
      pebar=2.d0*diffhg*diffng*axsec/dpbar  
        ! kpc * cm^-3 * 1.d-27 cm^2 * 1.d10 cm*s^-1 * 10^-27 s/cm^-2
        ! 10^-44 * kpc cm^-2
      pebar=pebar*(kpc/10.d0)**2 ! kpc kpc^-2 = kpc^-1
      dxy=dsqrt((x-x0)**2+(y-y0)**2)  ! kpc
      z0i=z0                                 
ccc
      conv=.false.
      low=0.d0
      result=0.d0
      result2=0.d0
      if(dxy.lt.4.d0*z0i) then
        low=0.d0
        up=2.d0/dxy
        call dspbpointexpint(up,result,result2,conv)
        low=2.d0/dxy
      endif
      if(conv) goto 100
      call dspbpointj0int(low,result,result2)
 100  continue
      dspbpoint=result/(2.d0*dpbar)*dexp(-vbar*dabs(z0i)/2.d0)  
                   ! kpc^-1 10^-27 cm^-2 s =kpc^-1 10^15 s (10^21 cm)^-2
      dspbpoint=dspbpoint*kpc**2/(2.d0*pi) ! kpc^-3 10^15 s
      return
      end
ccc
ccc
ccc
      subroutine dspbpointexpint(upin,result,result2,conv)
      implicit none
      real*8 upin,result,result2
      logical conv
ccc
      real*8 vbar,pebar,dxy,z0i
      common/dspbpointintcom/vbar,pebar,dxy,z0i
ccc
      integer ii
      real*8 low,up,eps,prec,res,res2,sum,sum2,sumcheck,dspbpoint_int,
     &  dspbpoint_int2
      external dspbpoint_int,dspbpoint_int2
      logical getout
      conv=.false.
      getout=.false.
      low=0.d0
      sum=0.d0
      sum2=0.d0
      sumcheck=0.d0
      prec=1.d-3
      ii=1
 10   continue
      up=4.d0/z0i*ii
      if(up.ge.upin) then
        up=upin
        getout=.true.
      endif
      eps=dabs(dspbpoint_int(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
      call dsfun_int(dspbpoint_int,low,up,eps,prec,res)
      eps=dabs(dspbpoint_int2(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
      call dsfun_int(dspbpoint_int2,low,up,eps,prec,res2)
      sum=sum+res ! kpc^-1
      sum2=sum2+res2
c      write(*,*) '1: low,up,res,sum : ',low,up,res,sum,res2,sum2
      if(ii.gt.3.and.dabs(sumcheck-sum).lt.1.d-12*abs(sum)) then
        conv=.true.
        getout=.true.
      endif
      sumcheck=sum
      ii=ii+1
      low=up
      if(.not.getout) goto 10
      result=sum
      result2=sum2
      return
      end
ccc
ccc
ccc
      subroutine dspbpointj0int(lowin,result,result2)
      implicit none
      real*8 lowin,result,result2
ccc
      real*8 vbar,pebar,dxy,z0i
      common/dspbpointintcom/vbar,pebar,dxy,z0i
ccc
      integer s,sskip,stake
      real*8 dspbbesselzeroj0
ccc
      real*8 low,up,eps,prec,res,res2,par,par2,sum,sum2,nusk,
     &  dspbpoint_int,dspbpoint_int2
      external dspbpoint_int,dspbpoint_int2
ccc
      low=lowin
      par=result
      par2=1.d0/dsqrt(dxy**2+z0i**2)-result2
      stake=2
      sskip=0
      do s=1,100000
        sskip=sskip+1
        if(sskip.eq.stake) then
          sskip=0
        else
          goto 20
        endif
        nusk=dspbbesselzeroj0(s)
        up=nusk/dxy 
        if(up.le.low) then
          write(*,*) 'DS: in dspbpointj0int up lower than low : ',up,low
          stop
        endif
        prec=1.d-3
        eps=dabs(dspbpoint_int(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
        call dsfun_int(dspbpoint_int,low,up,eps,prec,res)
        eps=dabs(dspbpoint_int2(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
        call dsfun_int(dspbpoint_int2,low,up,eps,prec,res2)
        par=par+res ! kpc^-1
        par2=par2-res2
        sum=par2+par
c        write(*,*) '2: low,up,res,sum : ',low,up,res,par,res2,par2,sum
        if(s.gt.5.and.dabs(sum-sum2).lt.1.d-8*abs(sum)) goto 100
        sum2=sum
        low=up
 20     continue
      enddo
 100  result=sum
      return
      end
ccc
ccc
ccc
      real*8 function dspbpoint_int(kappa)
      implicit none
      include 'dscraxicom.h' ! for diffhh
      real*8 kappa,lambda,partial,ratio,dbesj0,arg
      real*8 vbar,pebar,dxy,z0i
      common/dspbpointintcom/vbar,pebar,dxy,z0i
ccc
      if(vbar.gt.1.d-16) then
        lambda=dsqrt(kappa**2+vbar**2/4.d0)
        if(kappa.gt.1.d0) ratio=dsqrt(1.d0+(vbar/kappa)**2/4.d0)
      else
        lambda=kappa
        if(kappa.gt.1.d0) ratio=1.d0
      endif
      arg=lambda*(diffhh-dabs(z0i))
      if(arg.lt.37.d0) then
        if(dabs(z0i).gt.0.1d0*diffhh) then
          partial=dsinh(lambda*(diffhh-dabs(z0i)))/dsinh(lambda*diffhh)
        else
          partial=dcosh(lambda*dabs(z0i))
     &            -dsinh(lambda*dabs(z0i))/dtanh(lambda*diffhh)
        endif
      else
        partial=dexp(-lambda*dabs(z0i))
      endif
      if(kappa.lt.1.d0) then
        partial=partial
     &          /((pebar+vbar)/2.d0+lambda/dtanh(lambda*diffhh)) ! kpc
        partial=partial*kappa*dbesj0(kappa*dxy) ! dimensionless
      else
        partial=partial
     &          /((pebar+vbar)/(2.d0*kappa)+ratio/dtanh(lambda*diffhh)) ! kpc
        partial=partial*dbesj0(kappa*dxy) ! dimensionless
      endif
      dspbpoint_int=partial
      return
      end
ccc
ccc
ccc
      real*8 function dspbpoint_int2(kappa)
      implicit none
      real*8 kappa,dbesj0
      real*8 vbar,pebar,dxy,z0i
      common/dspbpointintcom/vbar,pebar,dxy,z0i
ccc
      dspbpoint_int2=dexp(-kappa*dabs(z0i))*dbesj0(kappa*dxy)
      return
      end



      real*8 function dspbpointe(x,y,z,x0,y0,z0,tp)
************************************************************************
*** inputs:
***     x,y,z - position of the observer (kpc)
***     x0,y0,z0 - position of the source (kpc)
***     tp - kinetic energy per nucleon (gev)
***
*** output: in kpc^-3 10^15 s
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h'
ccc
      real*8 vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
      common/dspbpointeintcom/vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
ccc
      real*8 x,y,z,x0,y0,z0,tp
      logical conv
      real*8 pp,ee,dpbarh,dpbarg,dskdiff,beta,axsec,dspbsigmavpbar,
     &  low,up,result,result2,rig,dsdbsigmavdbar
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
ccc check whether you are actually outside the diffusion region:
ccc      
      if(dabs(z0).ge.diffhh) then
        dspbpointe=0.d0
        return
      endif
ccc
ccc so far this function has been implemented for z=0 only
ccc
      if(dabs(z).gt.1.d-16) then
        write(*,*) 'DS: dspbpoint has not been implemented for z.ne.0'
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc set tp dependent variables:
ccc
      ee=anuc*tp+mnuc
      pp=dsqrt(dabs(ee**2-mnuc**2))
      rig=pp/dble(znuc)
      beta=pp/ee
      dpbarh=dskdiff(rig,1,beta)
      vbar=diffcvel/dpbarh    ! 10^5 cm/s * 10^-27 s/cm^-2 = 10^-22 cm^-1 
      vbar=vbar*(kpc/10.d0)  ! kpc^-1
      if(anuc.eq.1) then
        axsec=dspbsigmavpbar(ee)  ! mb * 1.d10 cm*s^-1 
      elseif(anuc.eq.2) then
        axsec=dsdbsigmavdbar(ee)  ! mb * 1.d10 cm*s^-1
      else
        write(*,*) 'DS: call to dspbpointe without proper'
        write(*,*) 'DS: intialization of mnuc, anuc : ',mnuc,anuc
        write(*,*) 'DS: program stopped'
        stop
      endif
      pebarh=diffnh*axsec/dpbarh  
        ! cm^-3 * 1.d-27 cm^2 * 1.d10 cm*s^-1 * 10^-27 s/cm^-2
        ! 10^-44 * cm^-2
      pebarh=pebarh*(kpc/10.d0)**2 ! kpc^-2 = kpc^-2
      dpbarg=dskdiff(rig,2,beta)
      pebarg=diffng*axsec/dpbarg  
        ! cm^-3 * 1.d-27 cm^2 * 1.d10 cm*s^-1 * 10^-27 s/cm^-2
        ! 10^-44 * cm^-2
      pebarg=pebarg*(kpc/10.d0)**2 ! kpc^-2 = kpc^-2
      ratiodpbar=dpbarg/dpbarh
      dxy=dsqrt((x-x0)**2+(y-y0)**2)  ! kpc
      z0i=z0                                 
ccc
      conv=.false.
      low=0.d0
      result=0.d0
      result2=0.d0
      if(dxy.lt.4.d0*z0i) then
        low=0.d0
        up=2.d0/dxy
        call dspbpointexpinte(up,result,result2,conv)
        low=2.d0/dxy
      endif
      if(conv) goto 100
      call dspbpointj0inte(low,result,result2)
 100  continue
      if(dabs(z0i).ge.diffhg) 
     &  result=result*dexp(-vbar*(dabs(z0i)-diffhg)/2.d0)  ! kpc^-1
      dspbpointe=result/(2.d0*dpbarh)
                   ! kpc^-1 10^-27 cm^-2 s =kpc^-1 10^15 s (10^21 cm)^-2  
      dspbpointe=dspbpointe*kpc**2/(8.d0*datan(1.d0)) ! kpc^-3 10^15 s
      return
      end
ccc
ccc
ccc
      subroutine dspbpointexpinte(upin,result,result2,conv)
      implicit none
      real*8 upin,result,result2
      logical conv
ccc
      real*8 vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
      common/dspbpointeintcom/vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
ccc
      integer ii
      real*8 low,up,eps,prec,res,res2,sum,sum2,sumcheck,dspbpoint_inte,
     & dspbpoint_int2e
      external dspbpoint_inte,dspbpoint_int2e
      logical getout
      conv=.false.
      getout=.false.
      low=0.d0
      sum=0.d0
      sum2=0.d0
      sumcheck=0.d0
      prec=1.d-3
      ii=1
 10   continue
      up=4.d0/z0i*ii
      if(up.ge.upin) then
        up=upin
        getout=.true.
      endif
      eps=dabs(dspbpoint_inte(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
c      write(*,*) 'first in dspbpointexpinte'
      call dsfun_int(dspbpoint_inte,low,up,eps,prec,res)
c      write(*,*) 'second in dspbpointexpinte'
      eps=dabs(dspbpoint_int2e(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
      call dsfun_int(dspbpoint_int2e,low,up,eps,prec,res2)
      sum=sum+res ! kpc^-1
      sum2=sum2+res2
c      write(*,*) '1: low,up,res,sum : ',low,up,res,sum,res2,sum2
      if(ii.gt.3.and.dabs(sumcheck-sum).lt.1.d-12*abs(sum)) then
        conv=.true.
        getout=.true.
      endif
      sumcheck=sum
      ii=ii+1
      low=up
      if(.not.getout) goto 10
      result=sum
      result2=sum2
      return
      end
ccc
ccc
ccc
      subroutine dspbpointj0inte(lowin,result,result2)
      implicit none
      include 'dscraxicom.h'
      real*8 lowin,result,result2
ccc
      real*8 vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
      common/dspbpointeintcom/vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
ccc
      integer s,sskip,stake
      real*8 dspbbesselzeroj0
ccc
      real*8 low,up,eps,prec,res,res2,par,par2,sum,sum2,nusk,
     &  dspbpoint_inte,dspbpoint_int2e
      external dspbpoint_inte,dspbpoint_int2e
ccc
      low=lowin
      par=result
      if(dabs(z0i).ge.diffhg) then
        par2=2.d0/dsqrt(dxy**2+z0i**2)/(1.d0+ratiodpbar)-result2
      else
        par2=1.d0/dsqrt(dxy**2+z0i**2)/ratiodpbar-result2
      endif
      stake=2
      sskip=0
      do s=1,100000
        sskip=sskip+1
        if(sskip.eq.stake) then
          sskip=0
        else
          goto 20
        endif
        nusk=dspbbesselzeroj0(s)
        up=nusk/dxy 
        if(up.le.low) then
          write(*,*) 'DS: in dspbpointj0int up lower than low : ',up,low
          stop
        endif
        prec=1.d-3
        eps=dabs(dspbpoint_inte(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
c        write(*,*) 'first in dspbpointj0inte, eps: ',eps
        call dsfun_int(dspbpoint_inte,low,up,eps,prec,res)
        eps=dabs(dspbpoint_int2e(low+(up-low)*0.9999d0))
     &     *(up-low)*prec*1.d-2
c        write(*,*) 'second in dspbpointj0inte, eps: ',eps
        call dsfun_int(dspbpoint_int2e,low,up,eps,prec,res2)
        par=par+res ! kpc^-1
        par2=par2-res2
        sum=par2+par
c        write(*,*) '2: low,up,res,sum : ',low,up,res,par,res2,par2,sum
        if(s.gt.5.and.dabs(sum-sum2).lt.1.d-8*abs(sum)) goto 100
        sum2=sum
        low=up
 20     continue
      enddo
 100  result=sum
      return
      end
ccc
ccc
ccc
      real*8 function dspbpoint_inte(kappa)
      implicit none
      include 'dscraxicom.h'
      real*8 kappa,lambdah,lambdag,partial,ratioh,ratiog,ratiolambda,
     &  dbesj0,arg
ccc
      real*8 vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
      common/dspbpointeintcom/vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
ccc
      if(vbar.gt.1.d-16.or.pebarh.gt.1.d-16) then
        lambdah=dsqrt(kappa**2+pebarh+vbar**2/4.d0)
        if(kappa.gt.1.d0) 
     &    ratioh=dsqrt(1.d0+pebarh/kappa**2+(vbar/kappa)**2/4.d0)
      else
        lambdah=kappa
        if(kappa.gt.1.d0) ratioh=1.d0
      endif
ccc
      if(pebarg.gt.1.d-16) then
        lambdag=dsqrt(kappa**2+pebarg)
        if(kappa.gt.1.d0) ratiog=dsqrt(1.d0+pebarg/kappa**2)
      else
        lambdag=kappa
        if(kappa.gt.1.d0) ratiog=1.d0
      endif
      if(kappa.gt.1.d0) then
        ratiolambda=dsqrt(1.d0+pebarh/kappa**2+(vbar/kappa)**2/4.d0)
     &    /dsqrt(1.d0+pebarg/kappa**2)
        partial=1.d0/(vbar/2.d0/kappa
     &                +ratioh/dtanh(lambdah*(diffhh-diffhg))
     &                +ratiodpbar*ratiog*dtanh(lambdag*diffhg))
      else
        ratiolambda=lambdah/lambdag
        partial=kappa/(vbar/2.d0+lambdah/dtanh(lambdah*(diffhh-diffhg))
     &                 +ratiodpbar*lambdag*dtanh(lambdag*diffhg))
      endif
ccc
      if(dabs(z0i).ge.diffhg) then
        arg=lambdah*(diffhh-dabs(z0i))
        if(arg.lt.37.d0) then
          partial=partial*dsinh(arg)
     &             /dsinh(lambdah*(diffhh-diffhg))/dcosh(lambdag*diffhg)
        else
          arg=-lambdah*dabs(z0i)+diffhg*(lambdah-lambdag)
          partial=partial*2.d0*dexp(arg)
        endif
      else
        arg=lambdag*(diffhg-dabs(z0i))
        if(arg.lt.37.d0) then
          partial=partial*(dcosh(arg)/dcosh(lambdag*diffhg)+
     &                     dsinh(arg)/dcosh(lambdag*diffhg)/ratiodpbar*
     &     (vbar/2.d0/lambdag
     &      +ratiolambda/dtanh(lambdah*(diffhh-diffhg))))
        else
          arg=-lambdag*dabs(z0i)
          partial=partial*dexp(arg)*(1.d0+1.d0/ratiodpbar*
     &     (vbar/2.d0/lambdag
     &      +ratiolambda/dtanh(lambdah*(diffhh-diffhg))))
        endif
      endif
      dspbpoint_inte=partial*dbesj0(kappa*dxy) 
      return
      end
ccc
ccc
ccc
      real*8 function dspbpoint_int2e(kappa)
      implicit none
      include 'dscraxicom.h'
      real*8 kappa,dbesj0
ccc
      real*8 vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
      common/dspbpointeintcom/vbar,pebarh,pebarg,ratiodpbar,dxy,z0i
ccc
      if(dabs(z0i).ge.diffhg) then
        dspbpoint_int2e=dexp(-kappa*dabs(z0i))*dbesj0(kappa*dxy)
     &                *2.d0/(1.d0+ratiodpbar)
      else
        dspbpoint_int2e=dexp(-kappa*dabs(z0i))*dbesj0(kappa*dxy)
     &                /ratiodpbar
      endif
      return
      end
