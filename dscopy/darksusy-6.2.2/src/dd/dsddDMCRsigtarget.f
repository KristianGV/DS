*******************************************************************************
*** Function dsddDMCRsigtarget provides the differential scattering cross   ***
*** section of a (potentially relativistic) DM on a target (detector)       ***
*** nucleus, d\sigma_N / dTr, where Tr is the nuclear recoil energy.        ***
***                                                                         ***
***  type : replaceable                                                     ***
***  desc : Differntial CR - DM scattering cross section                    ***
***                                                                         ***
***  Input:                                                                 ***
***    Tr    - kinetic energy of nucleus [GeV]                              ***
***              (recoil energy after scattering)                           ***
***    Tdm     - initial kinetic energy of DM particle [GeV]                ***
***                                                                         ***
***  'Hidden' input (common block flags):                                   ***
***    targetoption - determines which of the implemented concrete          ***
***                   target options should be used. Typically set in       ***
***                   dsddDMCRcountrate                                     ***
***                                                                         ***
***  Output:  d\sigma_N / dTr, units: cm**2/ GeV                            ***
***                                                                         ***
***  NB: This is the default function to be supplied as input to            ***
***      dsddDMCRdgammadt (typically simply by calling dsddDMCRcountrate).  ***
***      A particle module or main program can replace this function, in    ***
***      particular to add additional target options. When doing so MAKE    ***
***      SURE IT IS CONSISTENT with the scattering cross section specified  ***
***      in dsddDMCRsigCR!                                                  ***
***      Examples (commented out) are provided for simple scalar and vector ***
***      mediators.                                                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-03-10                                                         ***
*** mod  2019-10-30 (added light mediator examples)                         ***
*******************************************************************************
      real*8 function dsddDMCRsigtarget(Tr,Tdm)
      implicit none
      include 'dsmpconst.h'
      include 'dsio.h'
      include 'dsddcom.h'

      real*8 Tr, Tdm
      integer targetoption
      common /sigtargetcom/ targetoption
      
      real*8 Q2, Trmax, mdm            ! general kinematic quantities
      real*8 sigma0, sip, sin, sdp, sdn, mN, s   
      integer ierr, iel

c... comment in for use with light mediators
c      real*8 mmed, tmp, s
c      common /mymed/ mmed

      integer dsidnumber     ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sip, sin, sdp, sdn

c... functions
      real*8 dsmwimp, dsddTrmax
      

      dsddDMCRsigtarget=0.0d0
      mdm=dsmwimp()       
      
      iel = mod(targetoption,1000) ! used for attenuation in soil
      if ((targetoption/1000).eq.1.and.iel.gt.0.and.iel.le.Nelements) then
        mN = mNaU(iel)*atomicmassunit ! convert masses to GeV
      elseif ((targetoption/10).eq.1) then  ! 'xenon1t'
        mN = 122.05   ! mXe [GeV]
      elseif ((targetoption/10).eq.2) then  ! 'Borexino'
        mN = m_p      ! proton target
      elseif ((targetoption/10).eq.3) then  ! 'MiniBoone'
        mN = m_p      ! proton target
      else
        if (prtlevel.gt.1) then
          write(*,*) 'warning in dsddDMCRsigtarget:'
          write(*,*) 'unimplemented targetoption = ',targetoption
          write(*,*) 'Setting cross section to zero.'
        endif  
        return
      endif  
      Trmax = dsddTrmax(Tdm, mdm, mN) ! maximal recoil energy of DM particle
      if (Tr.gt.Trmax) return
      Q2 = 2*mN*Tr ! this allows to implement Q-dependent scatterings
      s  = (mN+mdm)**2 + 2*mN*Tdm  ! cms energy squared


c... The expressions below implement the simplest case where scattering
c... is isotropic (i.e. no Q-dependence apart from form factor)

      if (idold.ne.dsidnumber()) then ! new model -> calc q->0 scattering cross  
                                      ! sections from particle module *once*
        call dsddsigmanucleon(0.0d0,0.0d0,sip,sin,sdp,sdn,ierr)
        idold=dsidnumber()
      endif

      if ((targetoption/1000).eq.1) then ! only keep SI in soil (=dominant!)
        sigma0 = (sip+sin)/2. ! approximate DM-nucleon crosss ection
        sigma0 = sigma0*an(iel)**2*(mN*(m_p+mdm)/m_p/(mN+mdm))**2
      elseif ((targetoption/10).eq.1) then     ! 'xenon1t'
        sigma0 = (sip+sin)/2.   ! NB: This is not the actual Xe cross section,
                                ! but allows to compare to reported limits 
                                ! *per nucleon*, including form factors
      else  ! 'Borexino' or 'MiniBoone'
        sigma0 = sip+sdp        ! proton target
        sigma0 = sigma0/(1. + Q2/0.77**2)**4 ! Dipole form factor suppression
      endif
c      tmp = sigma0/Trmax ! ensure correct normalization
      dsddDMCRsigtarget = sigma0/Trmax ! ensure correct normalization

c... correct estimate by using full expression for SCALAR mediators
c      dsddDMCRsigtarget = tmp
c     &                * mmed**4/(mmed**2 + Q2)**2 ! naive estimate for light mediator
c     &                * (Q2+4*mdm**2)*(Q2+4*mN**2) /
c     &                (16.*s*mN**2*mdm**2/(mN+mdm)**2)


c... correct estimate by using full expression for VECTOR mediators
c      dsddDMCRsigtarget = tmp 
c     &                * mmed**4/(mmed**2 + Q2)**2 ! naive estimate for light mediator
c     &                * (Q2**2-2*Q2*s+2*(s-mN**2-mdm**2)**2) /
c     &                (8.*s*mN**2*mdm**2/(mN+mdm)**2)


      return
      end

