***********************************************************************
*** The function dskdtkd returns the kinetic coupling temperature [in MeV]
***
***  input: m0  - DM mass (in GeV)
***         how - integer flag
***
***  possible values for 'how':
***  1 - full calculation 
***      [recommended: following Bringmann, New J. Phys. 11, 105027 (2009);
***       updated in Ref.A of 1603.04884]
***  2 - fast calculation 
***      [neglecting SM masses and resonances in the scattering amplitude]
***  3 - "old" order-of-magnitude estimate 
***      [use only for comparison!]
***
***   type : commonly used
***   desc : Kinetic decoupling temperature
***
***
*** author: torsten bringmann (troms@physto.se), 2010-01-10
*** updates: 2013-06-12 (removed model-dependence)
***          2015-06-06 (added non-standard heat bath temperature)
***********************************************************************

      real*8 function dskdtkd(how)
      implicit none

      include 'dskdcom.h'
      include 'dsmpconst.h'
      include 'dsio.h'

      integer how,i

      integer n,nmin, lspecies,conv
      real*8  tmpres,tmpres2,a,cn,cnsum(4),nj,sqrtgstar,geff, heff, gefftmp
      real*8  dsfac,dsgamma
      real*8  xstart,xf,yi,yf, Hrate
      real*8 xi, dsrdxi, s
      real*8  stepguess,stepmin
      real*8 dsmwimp, dskdgammarate

c... these parameters are needed by dqagse
      real*8 aint,bint,epsabs,intres,abserr,dskdcint2
      integer neval,ier, limit
      parameter (limit=20)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last

      external dskdboltz,dskdcint2

      m0=dsmwimp()
      call dskdparticles    ! prepare KD routines

      tmpres=0d0
      dskdtkd=0d0

      if (how.gt.1.and.how.le.3) then !  use fast calculation to get a first estimate

        nmin=6
        do 10 i=1,4
          cnsum(i)=0d0
 10     continue
        do 50 i=1,6      ! add up contributing leptons
          call dskdm2simp(i,cn,n)
          if (cn.gt.1.d-10.and.n.lt.nmin) nmin=n
          
          if (n.le.-2) then
            write(*,*) 'badly defined parameter in dskdm2simp: n=',n
            write(*,*) 'leaving dskdtkd.f...'
            return
          endif
          if (i.le.3) cnsum(1)=cnsum(1)+cn ! neutrino contributions         
          if (i.ge.4) cnsum(i-2)=cnsum(i-3)+cn
  50    continue
        do 60 i=1,nBSMscattlight  ! add up LIGHT BSM scattering partners,
                                  ! treating them like neutrinos (i.e. massless DR)
                                  ! NB: this requires that dskdm2simp *starts*
                                  ! with new massless dof!!!
          call dskdm2simp(100+i,cn,n)
          if (cn.gt.1.d-10.and.n.lt.nmin) nmin=n
          
          if (n.le.-2) then
            write(*,*) 'badly defined parameter in dskdm2simp: n=',n
            write(*,*) 'leaving dskdtkd.f...'
            return
          endif
          cnsum(1)=cnsum(1)+cn       
  60    continue
    
        n=nmin
        nj=(1.-1/2.**(n+3))*dsfac(n+4)
        if (n.le.4) nj=nj*zeta(n+4)

        lspecies=1       ! start by assuming scattering only with neutrinos

 100    if (lspecies.eq.1) geff=3.35d0      ! only photons and neutrinos,
                                            ! reheating taken into account
        if (lspecies.eq.2) geff=10.75d0     ! rel. e+,e-,nu and photons
        if (lspecies.eq.3) geff=14.2d0      !  "   " and muons
        if (lspecies.eq.4) geff=83.4d0      ! all leptons, geff at 1 GeV

        a = sqrt(10/2.**9/pi**9)/sqrt(geff) 
     -      *cnsum(lspecies)*nj*mpl/m0
        tmpres = m0 / dsgamma((n+1d0)/(n+2d0))
     -             / ( (a/(n+2.))**(1/(n+2.)) )

c... check whether resulting Tkd is consistent with assumed number
c... of leptonic scattering partners
        if (lspecies.eq.1.and.tmpres.gt.(1.2*0.511d-3)) then
           lspecies=2
c           tmpres=1.2*0.511d-3 ! guess tkd=1.2*m_e if only electrons contribute
           if (cnsum(lspecies).gt.0.d0) goto 100
        endif
        if (lspecies.eq.2.and.tmpres.gt.(1.2*105.6d-3)) then
           lspecies=3
c           tmpres=1.2*105.6d-3 ! guess tkd=1.2*m_mu if only mus contribute
           if (cnsum(lspecies).gt.0.d0) goto 100
        endif
        if (lspecies.eq.3.and.tmpres.gt.(1.2*1.78)) then
           lspecies=4
c           tmpres=1.2*1.78  ! guess tkd=1.2*m_tau if *only* taus contribute
           if (cnsum(lspecies).gt.0.d0) goto 100
        endif

c... Slightly improve quick estimate for g
  110   tmpres2=tmpres
        call dskdgeff(tmpres2,gefftmp)
        tmpres=tmpres2*(1.+sqrt(gefftmp/geff)**(1./(2.+n)))/2.
        geff=gefftmp
        if (abs(tmpres-tmpres2)/tmpres.gt.mheps/2.) goto 110
      endif


c        if (how.eq.2.and.tmpres.gt.2*tqcdmin) write(*,*)
c     -    'WARNING: probably too high Tkd=',tmpres*1d3,
c     -    ' -- try full calculation'//
c     -    '(how=1 in dskdtkd) instead!'

      if (how.eq.1) then      ! full calculation  

c... to obtain an initial value for the integration of the Boltzmann eq.,
c... check where the scattering rate is no longer *much* larger than the Hubble rate
c..  xstart=m0/tmpres/100.
        xstart = 10. ! conservative initial guess of when we still are in EQ
  150   xf=1.5*xstart
        Hrate = (m0/xstart)**2*sqrt(4.*pi**3/45.d0*geff)/mpl
        if (Hrate.gt.dskdgammarate(xstart)/1.d1) then
           if (Hrate.gt.dskdgammarate(xstart)*1.d2) then ! we have probably never been in EQ...
             tmpres = 1.0d19
             goto 500
           endif       
           if (prtlevel.gt.0) then
             write(*,*) 'WARNING in dskdtkd: Initial scattering rate too small!'
             write(*,*) 'Calculation of Tkd unreliable'
           endif
        endif
        if (Hrate.lt.dskdgammarate(xstart)/1.d4) then
          xstart=xf
          call dskdgeff(m0/xstart,geff)
          goto 150        
        endif
c        call dskdboltz(xf,0.d0,rhs)
c        if (rhs.gt.1/mheps/5.) then
c          xstart=xf
c          goto 150
c        endif
c        call dskddof(1.d3*m0/xstart,sqrtg,sqrtgt)
        yi = xstart / (geff*2.d0*pi**2/45.d0)**(2.d0/3.d0)  ! ~value (geff instead of heff) in equilibrium
         
c... the integration of the Boltzmann Equation will proceed stepwise until
c... the required precision in Tkd is reached
  200   xf=1.2*xstart
        yf=yi    ! starting with value yi, yf will contain
                 ! result of integration on return
        stepguess=xstart*mheps/10.
        stepmin=xstart*1.D-20
        call dskdinty(yf,xstart,xf,mheps/10.,stepguess,stepmin,dskdboltz,
     &                ier)
        if (ier.ne.0) then
           write(*,*) 'dskdtkd: Cannot integrate Boltzmann Eq.!'
           write(*,*) 'Ti,Tistart,m0/yf,m0/yi : ',
     &                1d3*m0/xstart,5d3*tmpres,m0/yf,m0/yi
c          if (xf.lt.1.05*xstart) then
c            write (*,*) 'dskdtkd: Cannot integrate Boltzmann Eq.!'
c            stop
c          endif
c          xf=xf/1.7
c          goto 200
        endif
        if (abs(yf-yi)/yf.ge.mheps/2.) then  ! not yet reached Tkd (y asymptotes a constant)
          xstart=xf
          yi=yf
          goto 200
        endif
        call dsrddof(tmpres2,sqrtgstar,heff)
        s=heff*(2.*pi**2/45.d0) ! entropy up to T**3 factor
        xi=dsrdxi(xf) ! temperature ratio between scattering partners and photons
        tmpres=xi*m0/yf/s**0.66666  ! first estimate for Tkd

c... In case of large dg/dT, or dxi/dT, the above estimate for Tkd needs to be improved:
        conv=0
  300   tmpres2=tmpres
        call dsrddof(tmpres2,sqrtgstar,heff)
        s=heff*(2.*pi**2/45.d0)
        xi=dsrdxi(m0/tmpres2)
        ! do a weighted average over new and old estimate
        tmpres=(2*tmpres2 + xi*m0/yf/s**0.66666)/3.
        conv=conv+1
        if (conv.gt.15) write(*,*) 'Slow convergence in dskdTkd ',
     &                             '(how=1):Tkd = ',tmpres2,tmpres
        if (conv.gt.30) stop
        if (abs(tmpres-tmpres2)/tmpres.gt.mheps/2.) goto 300

      endif   ! how=1 (full calculation)


c... "old" estimate used in the literature (use only for comparison!)
      if (how.eq.3) then
        conv=0
  400   Tint=tmpres
        call dskdgeff(Tint,geff)  
        tmpres2=tmpres
        aint=1d-1
        bint=3d1
        epsabs=mheps*(bint-aint)*abs(dskdcint2(bint)-dskdcint2(aint))/2.
        call dqagse(dskdcint2,aint,bint,epsabs,1d-3,20,intres,
     &         abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        if (ier.ne.0) write(*,*) 'dskdTkd: integration problem for ',
     &     'how=3 at T = ',1d3*Tint,' MeV. Try different method!'
        tmpres=(2*tmpres2+m0*(1.089d3*m0/mpl*sqrt(geff)/intres)**(1./(2.+n)))/3.
        conv=conv+1 
        if (conv.gt.10) write(*,*) 'Slow convergence in dskdTkd ',
     &                             '(how=3):Tkd = ',tmpres2,tmpres
        if (conv.gt.20) then
           write(*,*) 'exiting dskdTkd...'
           return
        endif
        if (abs(tmpres-tmpres2)/tmpres.gt.mheps) goto 400

      endif     ! how=3 (order-of-magnitude estimate from older literature)


      if (how.lt.1.or.how.gt.3) then
         write(*,*) 'ERROR in dskdtkd -- unknown option: how=',how
         write(*,*) 'Stopping...'
         stop
      endif

 500  dskdtkd=1d3*tmpres   ! give result in MeV

      return

      end





        

