      real*8 function dspbtdps_td(t,tp)
************************************************************************
*** antiproton "confinement time per annihilation volume" for a 
*** non-static point source within which WIMP pair annihilations take 
*** place.
***
*** inputs:
***   it assumes that the observer is located at x,y,z=0       
***     t - time of observation (10^15 s), having assumed t=0 for the
***         time at which the source is the closest to the observer
***     tp - antiproton kinetic energy (gev)
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
***   units * 1/(4 pi) * vp * rhovol * dec-rate/mwimp * dn/dtp
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
      real*8 t,tp
      real*8 dspbtdpsc_td
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
      dspbtdps_td=dspbtdpsc_td(t,tp)
      return
      end
ccc
ccc
ccc
      real*8 function dspbtdpsc_td(t,tp)
      implicit none
      real*8 t,tp
      real*8 tinf,tup,eps,prec,x0,y0,z0,ti,result,par,tinfloc,tuploc,tf
      real*8 deltatlowest,deltatlower,deltat
      integer ii
      logical getout
      parameter (deltatlowest=1.d-7,deltatlower=1.d-6) ! 10^15 s
      external dspbps_tdi_int
ccc
      real*8 tc,tpc
      common/dspbps_tdicom/tc,tpc
ccc
      tc=t
      tpc=tp
      call dsorbitps(0.d0,x0,y0,z0,ti,tf)
      tinf=ti+deltatlowest
      if(t.lt.tinf+deltatlower) then
        write(*,*) 'DS: calling dspbtdpsc_td at time t = ',t,' before'
        write(*,*) 'DS: the source entered the diff. region at t0 =',ti
        dspbtdpsc_td=0.d0
        return
      endif
      eps=1.d-5
      prec=1.d-5
      tup=t
      tup=min(tup,tf-deltatlowest)
ccc
      deltat=0.5d0
      if(tup*tinf.lt.0.d0) then
        par=0.d0
        getout=.false.
        ii=0
 101    tinfloc=tinf+deltat*ii
        tuploc=tinf+deltat*(ii+1)
        if(tuploc.ge.0.d0) then
          tuploc=0.d0
          getout=.true.
        endif
        call dsfun_int(dspbps_tdi_int,tinfloc,tuploc,eps,prec,result)
c        write(*,*) 'a1: ',tinfloc,tuploc,par,result
        par=par+result
        if(.not.getout) then
          ii=ii+1
          goto 101
        endif
        getout=.false.
        ii=0
 102    tinfloc=deltat*ii
        tuploc=deltat*(ii+1)
        if(tuploc.ge.tup) then
          tuploc=tup
          getout=.true.
        endif
        call dsfun_int(dspbps_tdi_int,tinfloc,tuploc,eps,prec,result)
c        write(*,*) 'a2: ',tinfloc,tuploc,par,result
        par=par+result
        if(.not.getout) then
          ii=ii+1
          goto 102
        endif
      else
        par=0.d0
        getout=.false.
        ii=0
 103    tinfloc=tinf+deltat*ii
        tuploc=tinf+deltat*(ii+1)
        if(tuploc.ge.tup) then
          tuploc=tup
          getout=.true.
        endif
        call dsfun_int(dspbps_tdi_int,tinfloc,tuploc,eps,prec,result)
c        write(*,*) 'b: ',tinfloc,tuploc,par,result
        par=par+result
        if(.not.getout) then
          ii=ii+1
          goto 103
        endif
      endif
      dspbtdpsc_td=par
      return
      end
ccc
ccc
      real*8 function dspbps_tdi_int(t0)
      implicit none
      real*8 t0,x0,y0,z0,ti,tf,dspbps_td
ccc
      real*8 tc,tpc
      common/dspbps_tdicom/tc,tpc
ccc
      call dsorbitps(t0,x0,y0,z0,ti,tf)
      dspbps_tdi_int=dspbps_td(0.d0,0.d0,0.d0,tc,x0,y0,z0,t0,tpc)
      return
      end


      real*8 function dspbps_td(x,y,z,t,x0,y0,z0,t0,tp)
************************************************************************
*** inputs:
***     x,y,z,t - position of the observer (kpc) and time of 
***               observation (10^15 s)
***     x0,y0,z0,t0 - position of the source (kpc) at the time t0 when
***                   impulse source was on (10^15 s)
***     tp - antiproton kinetic energy (gev)
***
*** NOTE: this function assumes that mnuc,anuc,znuc in the common block
*** bnucleoncom have been correctly set before calling this function
***
*** output: in kpc^-3
************************************************************************
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h' 
      real*8 x,y,z,t,x0,y0,z0,t0,tp
      real*8 pp,ee,dpbar,dskdiff,beta,rig,axsec,dspbsigmavpbar,
     &  dsdbsigmavdbar,sum,znp,znn,add,sign,dspbsgtd,deltat,
     &  deltat2,sumcheck
      integer n
ccc
      real*8 mnuc
      integer anuc,znuc
      common/pbnucleoncom/mnuc,anuc,znuc
ccc
ccc check whether you are actually outside the diffusion region:
ccc      
      if(dabs(z0).ge.diffhh) then
        write(*,*) 'DS: dspbps_td for source outside the diff. region'
        write(*,*) 'DS: z0 = ',z0,' rather than within ',-diffhh,diffhh
        dspbps_td=0.d0
        return
      endif
ccc
ccc set tp dependent variables:
ccc
      ee=anuc*tp+mnuc
      pp=dsqrt(dabs(ee**2-mnuc**2))
      rig=pp/dble(znuc)
      beta=pp/ee
      dpbar=dskdiff(rig,1,beta)
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
c      10^27 cm^2/s/(10^21 cm)^2 = 1/(10^15 s)
ccc
      sum=dspbsgtd(z,t,z0,t0,dpbar,axsec)
      sumcheck=dabs(sum)
c      write(*,*) 0,sum,sum
      n=1
      sign=-1.d0
 100  znp=sign*z+2.d0*n*diffhh
      znn=sign*z-2.d0*n*diffhh
      add=sign*(dspbsgtd(znp,t,z0,t0,dpbar,axsec)
     &          +dspbsgtd(znn,t,z0,t0,dpbar,axsec))
      sum=sum+add ! 1/sqrt(10^15 s)
c      write(*,*) n,sum,add
      n=n+1
      sign=-sign
      if(n.lt.5.or.dabs(add).gt.1.d-10*abs(sum)) goto 100
      if(dabs(sum).lt.1.d-15*dabs(sumcheck)) then
        write(*,*) 'DS: convergence problem in z sum in dspbps_td'
        write(*,*) 'DS: sqrt(4*D*(t-t0)), diffhh: ',
     &    dsqrt(4.d0*dpbar/kpc**2*(t-t0)),diffhh
        dspbps_td=0.d0
        return
      endif
      deltat=(t-t0)
      deltat2=kpc**2*((x-x0)**2+(y-y0)**2)/(4.d0*dpbar)
      if(deltat.gt.4.d-3*deltat2) then
        dspbps_td=dexp(-deltat2/deltat)/8.d0/pi/deltat/dpbar**1.5d0*sum
c 1/(10^15 s)^(3/2)/(10^27 cm^2/s)^(3/2) = (10^21 cm)^-3
        dspbps_td=dspbps_td*kpc**3  ! kpc^-3
      else
        dspbps_td=0.d0
      endif
      return
      end
ccc
ccc
      real*8 function dspbsgtd(z,t,z0,t0,dpbar,axsec)
ccc output in 1/sqrt(10^15 s)
      implicit none
      include 'dsmpconst.h' 
      include 'dscraxicom.h'
      real*8 z,t,z0,t0,dpbar,axsec,dsderfunc
      real*8 b,q,deltat,deltat2
ccc
      b=1.d-2*kpc*diffhg*diffng*axsec/dsqrt(dpbar)  ! b=P(E)/2/sqrt(D)
c   10^23 cm * cm^-3 * mb *1.d10 cm*s^-1 /sqrt(10^27 cm^2/s)
c   1.d23*1.d-27*1.d10*1.d-6/dsqrt(10^15*s) = 1/dsqrt(10^15*s) 
      q=kpc*(dabs(z)+dabs(z0))/dsqrt(dpbar)   ! q=(|z|+|z_0|)/sqrt(D)
c   10^21 cm/sqrt(10^27 cm^2/s) = sqrt(10^15 s)
      deltat=(t-t0)
      deltat2=kpc**2*(z-z0)**2/(4.d0*dpbar)
      if(deltat.gt.4.d-3*deltat2) then
        dspbsgtd=dexp(-deltat2/deltat)/dsqrt(pi*deltat)
     &             -b*dexp(b*q+b**2*deltat)
     &               *dsderfunc(b*dsqrt(deltat)+q/2.d0/dsqrt(deltat))
      else
        dspbsgtd=0.d0
      endif
      return
      end
ccc
ccc
      real*8 function dsderfunc(x)
      implicit none
      real*8 x,t,z
      z=dabs(x)
      t=1.d0/(1.d0+0.5d0*z)
      dsderfunc=t*dexp(-z*z-1.26551223d0+t*(1.00002368d0
     & +t*(0.37409196d0+t*(0.09678418d0+t*(-0.18628806d0
     & +t*(0.27886807d0+t*(-1.13520398d0+t*(1.48851587d0
     & +t*(-0.82215223d0+t*0.17087277)))))))))
      if(x.lt.0.d0) dsderfunc=2.d0-dsderfunc
      return
      end
