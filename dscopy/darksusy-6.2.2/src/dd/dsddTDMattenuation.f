*******************************************************************************
*** Subroutine dsddTDMattenuation relates average DM kinetic energies       *** 
*** before and after propagating through a dense medium                     ***
***                                                                         ***
*** As a default, this routine does NOT take into account a possible        ***
*** Q-dependence of the scattering cross section. Effectively assuming the  ***
*** cross section at zero momentum transfer, it thus typically over-        ***
*** estimates the stopping power for e.g. light mediators.                  ***
***                                                                         ***
*** This behavior can be changed by setting the common block flag           ***
*** 'attenuation_how' to 2 (in dsddinit, it is instead set to 1 as default).***
*** WARNING: While this now should work for 'arbitrary' Q-dependence, it    ***
*** has only been tested for single mediators that are neither too light    ***
*** (<<MeV) nor too heavy (>>10 GeV). Since the numerical integration is    ***
*** potentially unstable, the behaviour of this routine must be carefully   ***
*** outside this regime and/or for more complicated energy dependendences.  ***
***                                                                         ***
***  Input:                                                                 ***
***    Tin    - initial kinetic energy  [GeV]                               ***
***    Tz     - average kinetic energy after scattering in medium [GeV]     ***
***    depth  - penetration depth (detector location) [cm]                  ***
***    how    - convert from Tin to Tz (how=1)                              ***
***             or from Tz to Tin (how=2)                                   ***
***                                                                         ***
***  Output:                                                                ***
***    Tin  - initial kinetic energy  [GeV]                                 ***
***    Tz    - average kinetic energy after scattering in medium [GeV]      ***
***                                                                         ***
***  NB: For how=1, the input value of Tz is overwritten on output          ***
***      For how=2, the input value of Tin is overwritten on output         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-07-07                                                         ***
*** mod  2019-10-031 added full Q2 dependence; argument zlfree -> depth     ***
*******************************************************************************
      subroutine dsddTDMattenuation(Tin, Tz, depth, how)
      implicit none
      include 'dsddcom.h'

      real*8 Tin, Tz, depth
      integer how, ierr, targetoption
      real*8 dsmwimp, mdm, zlfree, sip,sin,sdp,sdn, sigsi
      save zlfree

      real*8 dsddlfreesimp, dsgetTzratio

ccc to make sure that we only tabulate once per model / constraint
      integer dsidnumber
      integer idold
      real*8 depthold
      data idold/-123456789/
      data depthold/1.d100/
      save idold
      save depthold

      common /sigtargetcom/ targetoption
      save /sigtargetcom/

      mdm = dsmwimp()
      if (how.ne.1.and.how.ne.2) goto 20
      
      if (idold.ne.dsidnumber().or.depth.ne.depthold) then ! new model / constraint
        if (attenuation_how.eq.2) then ! fully take into account Q2-dependence
          call dstabulateTinTz(depth)
        else ! take the Q2->0 limit and assume scattering independent of Q2
          call dsddsigmanucleon(0.0d0,0.0d0,sip,sin,sdp,sdn,ierr)
          if ((targetoption/10).ne.2) then ! assume that DM scatters on n and p
            sigsi = (sip+sin)/2.           ! take *average* for simplified mean free path       
            zlfree = depth/dsddlfreesimp(sigsi, depth, 1) 
          else ! Borexino SD limits -> only scattering on protons 
            zlfree = depth/dsddlfreesimp(sdp, depth, 2) ! NB: using sigsi->sdp assumes 
                                                        ! an unrealistically high
                                                        ! stopping power!
            if (targetoption.eq.21) zlfree = 1.0d-5*zlfree ! In reality, spin-dependent 
                                                           ! scattering should have 
                                                           ! 'MUCH' larger mean free path
          endif
        endif
        idold=dsidnumber() 
        depthold=depth
      endif

c... this typically indicates that exp() has taken too large arguments
c... (but should be captured now)
      if (Tin.ne.Tin.or.Tz.ne.Tz) then
        write(*,*) 'WARNING in dsddTDMattenuation: Tin, Tz = ',Tin, Tz
      endif

      if (how.eq.1) then ! convert Tin to Tz
        if (Tin.lt.0.0d0) goto 10
        if (attenuation_how.eq.2) then 
          Tz = Tin*dsgetTzratio(Tin,1)
        else
          if (zlfree.gt.30.0d0) then ! nothing arrives that far...
             Tz = 1.0d-50
          else
             Tz = 2.*mdm*Tin / ((2.*mdm+Tin)*exp(zlfree) - Tin)
          endif 
        endif          
      elseif (how.eq.2) then ! convert Tz to Tin
        if (Tz.lt.0.0d0) goto 10
        if (attenuation_how.eq.2) then 
          Tin = Tz*dsgetTzratio(Tz,2)        
        else
          if (zlfree.gt.30.0d0) then ! nothing arrives that far...
            Tin = 1.0d25 ! 'infinitely' high initial energy needed
          elseif (Tz.gt.(2.*mdm/(exp(zlfree)-1.0d0))) then
            Tin = 1.0d25 ! 'infinitely' high initial energy needed
          else
            Tin = 2.*mdm*Tz / ((2.*mdm+Tz)/exp(zlfree) - Tz)
          endif
        endif          
      endif      

      return

 10   write(*,*) 'FATAL ERROR in dsddTDMattenuation:'
      write(*,*) 'You have supplied a non-positive input energy!'
      stop

 20   write(*,*) 'FATAL ERROR in dsddTDMattenuation:'
      write(*,*) 'unknown option how = ',how
      stop 
      end


******************************************************************************
*** auxiliary routines to tabulate and read out conversion between Tin and Tz
******************************************************************************

ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dstabulateTinTz(depth)
c... tabulates ratios (T0,Tz/T0) and (Tz,Tz/T0)      
ccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 Tin, depth
      
      real*8 lnT0dm, lnTzdm, zi, zf, eps, stepguess,stepmin
      integer maxdat, ier, i
      real*8 lnTdmmin, lnTdmmax
      parameter (maxdat=1000)
      parameter (lnTdmmin=log(1.d-9), lnTdmmax=log(1.d9))
      real*8 LogTinTzdat(maxdat,2,2)
      common /soilcom/ LogTinTzdat
      save /soilcom/
      external dsdlogTdxrhs
      logical dsisnan

      do i = 1, maxdat
        lnT0dm = lnTdmmin + (i-1.)/(1.*maxdat-1.)*(lnTdmmax-lnTdmmin)
        zi = 1.0d0 ! start at 1cm depth to avoid zero density at surface
        zf = depth
        eps = 3.0d-5
        stepguess = depth/20.
        stepmin = eps*depth
        lnTzdm = lnT0dm 
        call dskdinty(lnTzdm,zi,zf,eps,stepguess,stepmin,dsdlogTdxrhs,ier)
        if (lnTzdm.lt.(lnTdmmin-2.).or.ier.ne.0.or.dsisnan(lnTzdm)) lnTzdm = -50.
        LogTinTzdat(i,1,1) = lnT0dm   ! T0 -> Tz/T0
        LogTinTzdat(i,1,2) = lnTzdm - lnT0dm
        LogTinTzdat(i,2,1) = lnTzdm   ! Tz -> T0/Tz
        LogTinTzdat(i,2,2) = lnT0dm - lnTzdm

        if (lnTzdm.ne.lnTzdm.or.lnT0dm.ne.lnT0dm) then
          write(*,*) 'ERROR in dstabulateTinTz: ', 
     &                i, exp(lnT0dm), exp(lnTzdm), lnT0dm, lnTzdm
          stop  
        endif

      enddo  
      end ! dstabulateTinTz
  
      
ccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function dsgetTzratio(T,how) 
ccccccccccccccccccccccccccccccccccccccccccccc
c interpolate. assuming that soilcom is ordered
      implicit none
      
      real*8 T, lnt
      integer how  ! 1: T=T0->Tz; 2: T=Tz->T0
      
      integer maxdat, k, dk, klo(2), khi(2)
      parameter (maxdat=1000)
      real*8 res, LogTinTzdat(maxdat,2,2)
      common /soilcom/ LogTinTzdat
      save /soilcom/
      logical initialized
      data initialized /.false./
      save initialized, klo, khi

      lnt = log(T)  
      res = 0.0d0
      if (.not.initialized) then
        klo(1) = 1
        klo(2) = 1
        khi(1) = maxdat
        khi(2) = maxdat
        initialized = .true.
      endif

      if (lnt.le.LogTinTzdat(1,how,1)) then
         klo(how) = 1
         res = LogTinTzdat(1,how,2)
      elseif (lnt.ge.LogTinTzdat(maxdat,how,1)) then
         khi(how) = maxdat
         res = LogTinTzdat(maxdat,how,2)
      else
        dk=maxdat/100
 100    if (LogTinTzdat(khi(how),how,1).lt.lnt) then
           khi(how) = khi(how) + dk
           dk = 2*dk
           if (khi(how).lt.maxdat) goto 100
           khi(how)=maxdat
        endif
        dk=maxdat/100
 110    if (LogTinTzdat(klo(how),how,1).gt.lnt) then
           klo(how) = klo(how) - dk
           dk = 2*dk
           if (klo(how).gt.1) goto 110
           klo(how)=1
        endif
 120    if (khi(how)-klo(how).gt.1) then
           k=(khi(how)+klo(how))/2
           if (LogTinTzdat(k,how,1).gt.lnt) then
              khi(how)=k
           else
              klo(how)=k
           endif
           goto 120
        endif
         res = LogTinTzdat(klo(how),how,2)+
     &               (LogTinTzdat(khi(how),how,2)-LogTinTzdat(klo(how),how,2))
     &               *(lnt-LogTinTzdat(klo(how),how,1))/
     &               (LogTinTzdat(khi(how),how,1)-LogTinTzdat(klo(how),how,1))
      endif
      dsgetTzratio = 1.0d0
      dsgetTzratio = exp(res)

      return
      end ! dsgetTzratio
      

ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dsdlogTdxrhs(z,lnTdm,dlogTdx)   
ccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'dsddcom.h'
      include 'dsmpconst.h'

      real*8 z,lnTdm,dlogTdx

      integer i
      real*8 r,n(Nelements), res, tmp, Trmax, mdm, mN(Nelements)
      real*8 lnTra, lnTrb
      real*8 dsmwimp, dsddTrmax, dsddDMCRsigtarget, dssem_earthdenscomp, 
     &       dsf_int, dsTdsigdT
      external dsTdsigdT
      integer targetoption, tosav
      common /sigtargetcom/ targetoption
      save /sigtargetcom/
      real*8 Tdmcom
      common /Tdxrhscom/ Tdmcom
      logical initialized
      data initialized /.false./
      save initialized, n, mN

c... these parameters are needed by dqagse
      real*8 abserr
      integer neval,ier, limit
      parameter (limit=30)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last
      
      if (.not.initialized) then
        do i=1,Nelements
          r = r_earth-5.d2 ! no variation in density down to detector location!  
          n(i) = dssem_earthdenscomp(r,an(i)) ! density in cm**-3
          mN(i) = mNaU(i)*atomicmassunit      ! mass in GeV
        enddo 
        initialized = .true.
      endif

      res=0.0d0
      dlogTdx=0.0d0
      if (lnTdm.lt.-50.or.lnTdm.gt.60) return ! DEBUG: why would the latter occur !?
      Tdmcom = exp(lnTdm)
      tosav = targetoption
      mdm = dsmwimp()
      do i=1, 1!Nelements
         Trmax = dsddTrmax(Tdmcom, mdm, mN(i))  ! maximal recoil energy of DM particle
         if (Trmax.ne.Trmax) then
           write(*,*) 'ERROR in dsdlogTdxrhs/Trmax: ', z, lnTdm,Tdmcom,trmax
         endif
         targetoption = 1000 + i ! overwrite common block value
c... this only holds for a constant cross section, where the integral
c... is performed trivially
c         res = res+0.5*n(i)*dsddDMCRsigtarget(Trmax/1.01,Tdmcom)/1.01*Trmax**2
c... In general, we need to integrate numerically:
         lnTra = log(1.0d-10*Trmax) 
         if (lnTra.lt.-50.) lnTra=-50.
         lnTrb = Log(0.999999*Trmax) 
         if (lnTrb.gt.60.) lnTrb=60.
         tmp=0.0d0
c This integration routine is not good enough...              
c              if (lnTrb.gt.lnTra) tmp = dsf_int(dsTdsigdT,lnTra,lnTrb,1.0d-4)              
         if (lnTrb.gt.lnTra) call
c FIXME: get rid of hard-coded values that refer to required precision!         
     &       dqagse(dsTdsigdT,lnTra,lnTrb,1.0d-30,2.0d-4,limit,tmp,abserr,
     &       neval,ier,alist,blist,rlist,elist,iord,last)
         res = res + n(i)*tmp 

      enddo
      targetoption = tosav ! restore original targetoption
      
      dlogTdx = -res/Tdmcom ! return dlogTdx rather than dTdm/dx 
      
      end ! dsdlogTdxrhs

ccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function dsTdsigdT(lnTr)   
c... auxiliary routine for integration      
ccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 Tdmcom, Tr, lnTr, dsddDMCRsigtarget
      common /Tdxrhscom/ Tdmcom
      dsTdsigdT = 0.0d0
      if (lnTr.lt.-50.0.or.lnTr.gt.60.0) then      
        return
      endif
      Tr = exp(lnTr)
c... additional factor of Tr from log-integration      
      dsTdsigdT = Tr**2*dsddDMCRsigtarget(Tr,Tdmcom)
      if (dsTdsigdT.ne.dsTdsigdT) write(*,*) 'ERROR in dsTdsigdT:', Tr,Tdmcom,dsTdsigdT
      return
      end ! dsTdsigdT
