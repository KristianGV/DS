*******************************************************************************
*** Function dsddDMCRsigCR provides the differential scattering cross       ***
*** section of cosmic rays on DM, d\sigma_CR / dTkin.                       ***
***                                                                         ***
***  type : replaceable                                                     ***
***  desc : Differntial CR - DM scattering cross section                    ***
***                                                                         ***
***  Input:                                                                 ***
***    Tkin    - kinetic energy of DM particle [GeV]                        ***
***              (recoil energy after scattering)                           ***
***    Tcr     - initial kinetic energy of CR particle [GeV]                ***
***    CRtype  - CR type (1=proton, 2=helium)                               ***
***                                                                         ***
***  Output:  d\sigma_CR / dTkin, units: cm**2/ GeV                         ***
***                                                                         ***
***  NB: This particular function implements ISOTROPIC SCATTERING for an    ***
***      energy-independent cross section (but taking into account a        ***
***      dipole-like formfactor suppression). A particle module or main     ***
***      program can replace it, and the DMCR routines will automatically   ***
***      use the corresponding modified dependence on Tkin [Q^2]            ***
***      and Tcr [s].                                                       ***
***      Examples (commented out) are provided for simple scalar and vector ***
***      mediators.                                                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-03-10                                                         ***
*** mod  2019-10-30 (added light mediator examples)                         ***
*******************************************************************************
      real*8 function dsddDMCRsigCR(Tkin,Tcr,CRtype)
      implicit none
      include 'dsmpconst.h'
      include 'dsio.h'

      real*8 Tkin, Tcr
      integer CRtype
      
      real*8 Q2, s, Tkinmax, mcr, mdm               ! general kinematic quantities
      real*8 lambda, sigma0, sigij(27,27), sp, sihe ! needed for simplified dipole scattering
      integer ierr

c... comment in for use with light mediators
c      real*8 mmed, tmp, s
c      common /mymed/ mmed

      integer dsidnumber ! to make sure cross section is only called once per model
      integer idold
      data idold/-123456789/
      save idold, sp, sihe 

c... functions
      real*8 dsmwimp, dsddTrmax
      

c... This general part just sets the kinematics and should not be touched
c... when replacing this function      
      dsddDMCRsigCR=0.0d0
      mdm=dsmwimp()       
      Q2 = 2*mdm*Tkin
      
      if (CRtype.eq.1) then
        mcr = m_p
      elseif (CRtype.eq.2) then
        mcr = m_he
      else
        if (prtlevel.gt.1) then
          write(*,*) 'warning in dsddDMCRsigCRdTkin:'
          write(*,*) 'unimplemented option CRtype = ',CRtype
          write(*,*) 'Setting cross section to zero'
        endif  
        return
      endif  
      s = (mcr+mdm)**2 + 2*mdm*Tcr       ! cms energy squared
      Tkinmax = dsddTrmax(Tcr, mcr, mdm) ! maximal recoil energy of DM particle
      if (Tkin.ge.Tkinmax) return


c... Now implement energy-independent, isotropic scattering
c... This part can be replaced with an arbitrary s and Q2 dependence

      if (idold.ne.dsidnumber()) then ! new model -> calc sp and sihe *once*
        sp   = 0.0d0
        sihe = 0.0d0
        call dsddsigma(0.0d0,0.0d0,1,1,sigij,ierr)  ! cross section at zero energy and Q2
        if (ierr.eq.0) sp  = sigij(1,1)+sigij(4,4)  ! SI+SD for protons   
        call dsddsigma(0.0d0,0.0d0,4,2,sigij,ierr) 
        if (ierr.eq.0) sihe = sigij(1,1)            ! SI for helium 
        idold=dsidnumber()
      endif

      sigma0 = 0.0d0
      if (CRtype.eq.1) then
        lambda = 0.77d0 ! Dipole suppression for protons 
        sigma0 = sp
      else
        lambda = 0.41d0 ! Dipole suppression for He 
        sigma0 = sihe
      endif      
      sigma0 = sigma0/(1. + Q2/lambda**2)**4 ! Dipole suppression of cross section
c      sigma0 = sigma0*(mcr+mdm)**2/s ! large Ecm suppression -> rarely relevant...

c      tmp = sigma0/Tkinmax ! ensure correct normalization
      dsddDMCRsigCR = sigma0/Tkinmax ! ensure correct normalization

c... correct estimate by using full expression for SCALAR mediators
c      dsddDMCRsigCR = tmp
c     &                * mmed**4/(mmed**2 + Q2)**2 ! naive estimate for light mediator
c     &                * (Q2+4*mdm**2)*(Q2+4*mcr**2) /
c     &                (16.*s*mcr**2*mdm**2/(mcr+mdm)**2)

c... correct estimate by using full expression for VECTOR mediators
c      dsddDMCRsigCR = tmp
c     &                * mmed**4/(mmed**2 + Q2)**2 ! naive estimate for light mediator
c     &                * (Q2**2-2*Q2*s+2*(s-mcr**2-mdm**2)**2) /
c     &                (8.*s*mcr**2*mdm**2/(mcr+mdm)**2)

      return
      end

