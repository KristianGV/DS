      real*8 function dsrdomega(option,fast,xf,ierr,iwar,nfc)

**********************************************************************
*** function dsrdomega calculates omega h^2 from thermal freeze-out of
*** a dark matter particle.
***
***  type : commonly used
***  desc : Calculate the relic density $\Omega h^2$ of a dark matter particle
***
*** input:
***   option - integer parameter passed to dsrdparticles provided by particle 
***            module (typically describes which coannihilations to include)
***   fast =   0 - standard accurate calculation (accuracy better than 1%)
***            1 - faster calculation: sets parameters for when to add
***                extra points less tough to avoid excessively
***                adding extra points, expected accurarcy: 1% or better      
***            2 - faster calculation: compared to fast=1, this option
***                adds less points in Weff tabulation in general and      
***                is more elaborate in deciding when to include points
***                close to thresholds and resonances
***                expected accuracy: around 1%      
***            3 - even more aggressive on trying minimize the number
***                of tabulated points
***                expected accuracy: 5-10%      
***            9 - superfast. This method still makes sure to include
***                resonances and threholds, but does not attempt to sample
***                them very well. Should give an order of magnitude estimate
***                expected accuracy: order of magnitude      
***           10 - quick and dirty method, i.e. expand the annihilation
***                cross section in x (not recommended)
***                expected accuracy: can be orders of magnitude wrong
***                for models with strong resonances or thresholds   
***           20 - new method to tabulate Weff along with 
***                <sv> calculation (rather than independently)
***                Also uses Breit-Wigner fit to resonances if the fit
***                is good enough.
***           99 - uses same settings as fast=0, but instead of
***                tabulating dsanwx, it uses it directly. This is (unless
***                the annihilation cross section is fast to calculate)
***                very inefficient and is included as an option mostly
***                for testing.
***
*** implicit input:
***    dsanwx -    external routine provided by particle physics module
***                to return the invariant annihilation rate W_eff
***                as a function of momentum, p_eff
***    dsrdparticles - routine provided by particle physics module
***                to return a list of coannihilating particles,
***                their properties, resonances and thresholds
***      
*** output:
***   dsrdomega - omega h^2 for the TOTAL dark matter density of a given DM
***               candidate (i.e. the contributions from both the DM particle
***               and its antiparticle in case DM is not self-conjugate)
***   xf        = m_WIMP / T_f where T_f is the photon temperature at freeze-out
***   ierr      = error from dsrdens or dsrdqad
***   iwar      = warning from dsrdens of dsrdqad
***   nfc       - number of function calls to the effective annihilation
***               cross section
***
*** authors: joakim edsjo, paolo gondolo
*** date: 98-03-03
*** modified: 98-03-03
***           99-07-30 pg
***           02-02-27 joakim edsjo: including sfermion coanns
***           06-02-22 paolo gondolo: streamlined inclusion of coanns
***           13-10-04 paolo gondolo: totally decouple from mssm
***           18-02-28 torsten bringmann: check whether DM is self-conjugate
**********************************************************************

      implicit none
      include 'dsio.h'
      include 'dsidtag.h'
      include 'dsrdcom.h'

      integer ierr,iwar,option,nfc,fast,selfconj
      real*8 oh2,xf
      real*8 dsanwx,dsrdwx
      external dsanwx,dsrdwx
      logical dsisnan
      logical firsterr
      data firsterr/.true./
      save firsterr

c...Temporary array for communicating with dsrdparticles
      real*8 tmco(tharsi),tdof(tharsi),trm(tharsi),trw(tharsi),
     &     ttm(tharsi)
      integer tnco,tnrs,tnthr

c----------------------------------------------------------------------

      if (rdinit.ne.1234) then
         write(*,*) 'DS ERROR: dsrdinit (called from dsinit) is',
     &        ' not called prior to dsrdomega.'
      endif

c...Note. This routine will change RD common blocks (depending on the
c...fast flag. We will here store current state and then reset it
c...at the end to avoid inconsistent results when setting common
c...block variables outside
      call dsrdstate('save')
      
      call dsrdparticles(option,tharsi,selfconj,tnco,tmco,tdof,
     &     tnrs,trm,trw,tnthr,ttm)


      if (selfconj.ne.1.and.selfconj.ne.2) then
         write(*,*) 'ERROR in dsrdomega: undefined value selfconj = ',
     &     selfconj
         write(*,*)
     &  'This should have been set to 1 (for self-conjugate DM)'
         write(*,*) 'or 2 (to distinguish form anti-DM) when setting up'
         write(*,*) 'the particle model.'
         write(*,*) 'Used particle module: ',moduletag
         stop
      endif

c...Now we have the masses, degrees of freedom, resonance and thresholds
c...from the particle physics module. Transfer to RD common blocks

      call dsrdstart(tnco,tmco,tdof,tnrs,trm,trw,tnthr,ttm)

      oh2=0.0d0
      xf=0.0d0
      nfc=0

      if (fast.eq.0.or.fast.eq.99) then       ! fast=0 or 99 calculation
            waccd=0.005d0
            dpminr=1.d-4
            dpthr=5d-4
            wdiffr=0.05d0
            wdifft=0.02d0
c TB add: need better accuracy for low masses!
            cosmin=0.9999 ! no reason to go only down to 5 degrees... !?
            if (mco(1).lt.1.d0) then
              xinit = 0.01   ! instead of 2.0
              xfinal = 5.0d4 ! instead of 2.0d2
            endif
   
      ! if (fast.eq.0.or.fast.eq.99) then       ! fast=0 or 99 calculation
      !    waccd=0.005d0
      !    dpminr=1.d-4
      !    dpthr=5d-4
      !    wdiffr=0.05d0
      !    wdifft=0.02d0         
      elseif (fast.eq.1) then       ! fast=1 calculation
c         dwopt=.true.
         waccd=0.05d0
         dpminr=5d-4
         dpthr=2.5d-3
         wdiffr=0.5d0
         wdifft=0.1d0
c         brmin=1.0d-4           ! 1.0d-3 doesn't make a big speed difference
      elseif (fast.eq.2) then       ! fast=2 calculation
c         dwopt=.true.
         waccd=0.05d0
         dpminr=1d-3
         dpthr=2.5d-3
         wdiffr=0.5d0
         wdifft=0.1d0
c         brmin=1.0d-4           ! 1.0d-3 doesn't make a big speed difference
      elseif (fast.eq.3) then       ! fast=3 calculation
c         dwopt=.true.
         waccd=0.05d0
         dpminr=2d-3
         dpthr=2.5d-3
         wdiffr=0.5d0
         wdifft=0.1d0
         nlow=20
         nhigh=10
         npres=3
         nthup=3
         umax=10.d0
c         brmin=1.0d-4           ! 1.0d-3 doesn't make a big speed difference
      elseif (fast.eq.9) then       ! fast=9 calculation
c         dwopt=.true.
         waccd=0.2d0
         dpminr=1.0d-2
         dpthr=2.5d-3
         wdiffr=0.5d0
         wdifft=0.1d0
         nlow=10
         nhigh=3
         npres=1
         dpres=0.85d0
         nthup=2
         umax=5.d0
c         brmin=1.0d-4           ! 1.0d-3 doesn't make a big speed difference
      elseif (fast.eq.20) then
         dpminr=0.002d0
c         dpminr=0.0001d0 ! possible option, not too much time consuming
         if (thavint.ne.1) then
            if (firsterr) then
               write(*,*) 'DS WARNING in dsrdomega: ',
     &              'fast=20 only well-tested with thavint=1,'
               write(*,*) 'You have set thavint = ',thavint,', ',
     &      'changing to thavint=1.'
               write(*,*) 'This warning will only be printed once.'
               firsterr=.false.
            endif
            thavint=1
         endif
      endif

c--------------------------------------initialize relic density routines
      rdprt = prtlevel

      rdtag=idtag

      if (fast.ge.0.and.fast.le.9.or.fast.eq.99) then  ! standard accurate method
c         call dsrdens(dsanwx,ncoann,kpart,mcoann,dof,nres,rm,rw,
c     &     nthr,tm,oh2,tf,ierr,iwar)
         call dsrdens(dsanwx,oh2,xf,fast,ierr,iwar)
         nfc=nr
      elseif (fast.eq.10) then   ! quick and dirty method
         call dsrdqad(dsanwx,mco(1),oh2,ierr)
         iwar=0
         nfc=2
      elseif (fast.eq.20) then   ! new method for tabulation
         call dsrdens(dsrdwx,oh2,xf,fast,ierr,iwar)
      else 
         write(*,*) 'ERROR in dsrdomega: invalid option fast=',fast 
         stop
      endif

c...  the following two lines are added by je for test
      if (prtlevel.eq.-1) then
         call dsrdwres            ! je test
         write (6,*) '============== new point =============='
c         call dsrdwrate(6,6,1)
c         call dsrdwrate(6,6,2)
c         call dsrdwrate(6,6,3)
      endif

c...  the following line is added by pg for test 050528
      if (prtlevel.eq.-111) call dsrdwintrp(dsanwx,50)

      if (selfconj.eq.1) then
         dsrdomega = oh2  ! DM is its own anti-particle
      else 
         dsrdomega = 2.0d0*oh2  ! dsrdens only returns contribution from DM 
                                ! particle; need to add that from anti-particle
      endif

c... added by TB      
      if (dsisnan(oh2).or.(oh2.lt.0)) then
        if (prtlevel.ge.1) then
           write (*,*) 
     &     'WARNING: dsrdomega calculated the relic density to be ', oh2
           write(*,*) 'ierr, iwar, nfc = ',ierr, iwar, nfc
           write(*,*) 'Setting oh2 and xf instead to zero...'
        endif     
        dsrdomega=0.0d0
        xf=0.0d0
      endif

c...Now rest common block variables to what they were before entering
c...this routine      
      call dsrdstate('reset')

      end
