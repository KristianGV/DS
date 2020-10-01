************************************************************************
*** function dsIB2sigmav returns the 3-body annihilation rate (sigma v), 
*** in the v->0 limit, for the following channels:
*** (calling it with IB2chinput=0, the sum over all these channels is 
*** returned)
***
***   Input: IB2chinput  - 0   sum over all states
***                        1XX for ffZ
***                        2XX for fFW
***                        3XX for ffh (XX like for ffZ)
***                        4XX for ffH (XX like for ffZ)
***                        5XX for ffA (XX like for ffZ)
***                        6XX for fFH+- (XX like for fFW)
***
***                         XX   | f \bar f Z      |   f \bar f W
***                        ------+-----------------+-------------------
***                         01   | nu_e   nu_e   Z |   nu_e   e      W-
***                         02   | e      e      Z |   e      nu_e   W+
***                         03   | nu_mu  nu_mu  Z |   nu_mu  mu     W-
***                         04   | mu     mu     Z |   mu     nu_mu  W+
***                         05   | nu_tau nu_tau Z |   nu_tau tau    W-
***                         06   | tau    tau    Z |   tau    nu_tau W+
***                         07   | u      u      Z |   u      d      W-
***                         08   | d      d      Z |   d      u      W+
***                         09   | c      c      Z |   c      s      W-
***                         10   | s      s      Z |   s      c      W+
***                         11   | t      t      Z |   t      b      W-
***                         12   | b      b      Z |   b      t      W+
***                         (note that this numbering must be consistent 
***                          with dsmssm.h !!!)
***
***          onshell - 0 subtracts the on-shell contributions ("2-body final 
***                      states")in the narrow width approximation [default]
***                      (note that the result can be negative!)
***                    1 returns the result including on-shell contributions
***
***   Output: annihilation rate in units cm**3 s**-1 
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-12-07, update 2105-03-22 (streamlined handling of kf)
***         update 2016-02-05 updated to DS6 conventions
************************************************************************

      real*8 function dsib2sigmav(IB2chinput,onshell)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'
      include 'dsmpconst.h'


c------------------------ functions ------------------------------------

      real*8   dsIB2dsde_aux,dsib2intres, dssigmav0, dsIB2BRtree
      external dsIB2dsde_aux

c------------------------ variables ------------------------------------

      integer IB2chinput,onshell
      integer IB2ch,ier,err, kf
      real*8  result,tmpres,mfinal1, mfinal2, mfinal3
      real*8  Emin,Emax     
c-----------------------------------------------------------------------


      dsib2sigmav=0d0
      result=0d0
      IB2ch=IB2chinput

c...set channel when summing over all channels
  10  if (IB2chinput.eq.0) then
        if (IB2ch.eq.611) goto 500       ! end sum 
        if (IB2ch/100.eq.2.or.IB2ch/100.eq.6) then
          IB2ch=IB2ch+2
        elseif (IB2ch.gt.0) then
          IB2ch=IB2ch+1
        endif
        if (mod(IB2ch,100).eq.13) IB2ch=IB2ch+88
        if (IB2ch.eq.0) IB2ch=101
      endif

      tmpres=0d0

c... set masses and internal channel specifications
      call dsib2chinit(IB2ch,err)
      if (err.eq.1) then
         write(*,*) 'ERROR in dsIB2sigmav: unknown channel', IB2ch
         goto 450
      endif  

c... determine kinematics
      if (ib2svhow.eq.1) then        ! do E1 integration first, then EB
         mfinal1=MB
         mfinal2=Mf
         mfinal3=Mff
         c_pfinal='B'
      elseif (ib2svhow.eq.2) then    ! do EV integration first, then Ef
         mfinal1=Mf
         mfinal2=Mff
         mfinal3=MB
         c_pfinal='f'
      elseif (ib2svhow.eq.3) then    ! do EV integration first, then Eff
         mfinal1=Mff
         mfinal2=Mf
         mfinal3=MB
         c_pfinal='fbar'
      else
         write(*,*) 'ERROR in dsIB2sigmav: unknown ib2svhow = ',ib2svhow
         return     
      endif   
      if (2.d0.le.1.0001*(mfinal1+mfinal2+mfinal3)) goto 450  !kinematically not allowed          
      Emin = mfinal1
      Emax= 1.d0+(mfinal1**2-(mfinal2+mfinal3)**2)/4.d0

         
c... full 3-body cross section
      tmpres=gev2cm3s*dsib2intres(dsIB2dsde_aux,c_pfinal,Emin,Emax,IB2acc,ier)
      if (ier.ne.0) write(*,*) 'dsib2sigmav: WARNING - error in ',
     &                          c_pfinal,' integration: ', ier, c_Btype, c_ftype

      tmpres=tmpres/1.d20  !correct for factor introduced in dsib2Msqaux

c... default: subtract on-shell contributions in NWA approximation
c MG changed call to BR calculation (BR for each fermion species individually)
c MG added factor 0.5 in dssigmav(11)=WH=W+H- + W-H+, because single contribution is required 
c MG interchanged 8<->9 in case IB2ch.ge.101.and.IB2ch.le.112
c MG 15-03-18 interchanged 5<->6 in case IB2ch.ge.501.and.IB2ch.le.512

      kf = mod(IB2ch,100) ! fermion identifier used in dmssm.h

      if (onshell.ne.1) then
        if (IB2ch.ge.201.and.IB2ch.le.210) then  ! WfF other than Wtb
          tmpres=tmpres - dssigmav0(24,-24)*dsIB2BRtree('WfF',kf)
     &                  - 0.5*dssigmav0(24,-37)*dsIB2BRtree('HfF',kf)  
        elseif (IB2ch.eq.211.or.IB2ch.eq.212) then      ! Wtb
          tmpres=tmpres - 0.5*dssigmav0(24,-37)*dsIB2BRtree('Hbt',11)
     &                  - dssigmav0(6,-6)*dsIB2BRtree('tWb',12)
        elseif (IB2ch.ge.101.and.IB2ch.le.112) then  ! Zff
          tmpres=tmpres - 2.*dssigmav0(23,23)*dsIB2BRtree('Zff',kf)
     &                  - dssigmav0(23,25)*dsIB2BRtree('hff',kf)  
     &                  - dssigmav0(23,35)*dsIB2BRtree('Hff',kf)  
        elseif (IB2ch.ge.601.and.IB2ch.le.610) then  ! HfF
          tmpres=tmpres - 0.5*dssigmav0(24,-37)*dsIB2BRtree('WfF',kf)
        elseif (IB2ch.ge.501.and.IB2ch.le.512) then  ! Aff
          tmpres=tmpres - dssigmav0(35,36)*dsIB2BRtree('Hff',kf)
     &                  - dssigmav0(25,36)*dsIB2BRtree('hff',kf)  
        elseif (IB2ch.ge.401.and.IB2ch.le.412) then  ! Hff
          tmpres=tmpres - dssigmav0(35,36)*dsIB2BRtree('Aff',kf)
     &                  - dssigmav0(23,35)*dsIB2BRtree('Zff',kf)  
        elseif (IB2ch.ge.301.and.IB2ch.le.312) then  ! hff
          tmpres=tmpres - dssigmav0(25,36)*dsIB2BRtree('Aff',kf)
     &                  - dssigmav0(23,25)*dsIB2BRtree('Zff',kf)  
        endif     
      endif

      result=result+tmpres
      
c... loop for total sigma; count WfF and HfF states twice     
 450  if (IB2chinput.eq.0) then
        if ((IB2ch/100.eq.2.or.IB2ch/100.eq.6).and.ier.eq.0) 
     &    result=result+tmpres
        goto 10  
      endif

 500  dsib2sigmav=result

      return
      end

