***********************************************************************
*** dskdgammarate returns the momentum transfer rate, as defined in
*** (A6) of arXiv:1603.04884. Input is the *photon* temperature,
***
***    x = mDM / Tphoton
***
*** NB: src_models/xxx/mh/dskdparticles MUST be called before using this!
***    (typically done via a call to dskdtkd)
***
*** author: torsten.bringmann@fys.uio.no, 2016-12-21
*** updates: 2018-05-28 (allow for Tdark != Tphoton) 
***********************************************************************

      real*8 function dskdgammarate(x)
      implicit none

      include 'dsmpconst.h'
      include 'dskdcom.h'

      real*8 x

      integer nsub, i
      real*8  c,subint(20)
      real*8 xd, dsrdxi

c... these parameters are needed by dqagse
      real*8 aint,bint,epsabs,result,abserr
      integer neval,ier, limit
      parameter (limit=40)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last

      logical dsisnan
      real*8 dskdcint
      external dskdcint

      dskdgammarate=0.0d0
      if (dsisnan(x).or.(x.le.0.d0)) return

c... the momentum transfer rate actually depends *exclusively* on the heat bath
c... temperature of the scattering partners, not on that pf the photons
      xd = x/dsrdxi(x)
      Tint=m0/xd

c... integrate over scattering momenta from T/10 to 30 T, which is a *very* good
c... approximation to the full integration [0,infinity],  taking account
c... of resonances by dividing the interval into suitable subintervals

      nsub=1
      aint=Tint/10.
      bint=30*Tint
c      write(*,*)
c      write(*,*) 'Total range : ',aint,bint
c      write(*,*) 'Resonances : ',nKDres,resk
c      write(*,*) '---------------'
      subint(1)=aint
      i=0        ! i indexes the resonance closest to start of present interval
 100  i=i+1
c      write(*,*) 'nsub,i,subint(i): ',nsub,i,subint(i)
      if (i.gt.nKDres.or.(resk(i)/2d0).ge.bint) goto 200 
      if (subint(nsub).ge.2*resk(i)) goto 100
      nsub=nsub+1
      if (subint(nsub-1).le.(resk(i)/2d0)) then
        subint(nsub)=resk(i)/1.9999
        i=i-1    ! the intervall covering this resonance will 
                 ! only be closed in next step
        goto 100
      endif
      if (i.eq.nKDres.or.(resk(i+1)/2d0).gt.2*resk(i)) then
        subint(nsub)=2.0001*resk(i)
      else
        subint(nsub)=(resk(i)+resk(i+1))/2.
      endif
      if (nsub.ge.limit-1.or.subint(nsub).ge.bint) then
        nsub=nsub-1
        goto 200
      else
        goto 100
      endif
c... no more resonances of interest found -- define last subinterval
  200 subint(nsub+1)=bint

      c=0
      do 300 i=1,nsub
        aint=subint(i)
        bint=subint(i+1)
c        write(*,*) 'Now integrating between ',aint,bint
        epsabs=5.d-2*mheps*(bint-aint)*abs(dskdcint(bint)-dskdcint(aint))
        epsabs=epsabs*1.d25 ! integrand is multiplied by 1d25 for better integration properties
        call dqagse(dskdcint,aint,bint,epsabs,1.d-3*mheps,limit,result,
     &         abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        result=result/1.d25 ! correct for rescaling factor in integrand
        if (ier.ne.0.or.result.lt.0) 
     &     write(*,*) 'dskdgammarate: problem in integration of ',
     &     aint,bint,result
c     &     'subinterval ',i, '(of ',nsub,') at T = ',1d3*m0/x,' MeV',
c     &    ' (ier = ',ier,' )'
        if (ier.eq.0.and.result.gt.0) c=c+result
  300 continue

      dskdgammarate=c

      return

      end





        

