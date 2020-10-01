*****************************************************************************
*** function dsandydth is the differential yield dyield/dcostheta in the
*** cm system boosted to the lab system (including proper jacobians if
*** we are dealing with a differential yield.
*** the function should be integrated from -1 to 1, by e.g.
*** the routine gadap.
*** units: (annihilation)**-1
*****************************************************************************

      real*8 function dsandydth(cth)
      implicit none
      include 'dsanyieldmodelcom.h'

c------------------------ variables ------------------------------------

      real*8 cth,yield,jac,elab,yield1,yield2
      real*8 e1cm,ethcm,dtldt0,tcm(2),dtdt(2)
      integer istat
      integer yieldkk, yieldpdg, diff

c------------------------ functions ------------------------------------
      real*8 dsanyield_sim
c-----------------------------------------------------------------------

c...calculate the energy of particle 1 in the cm system
      e1cm=phie1  ! is this reasonable even for t b-bar ?

c... calculate the cm energy for this angle that corresponds to requested
c... lab energy phieth (m_mu, m_ep etc neglected, but not m_p)
      dtldt0=0.d0
      if (phitype.eq.0) then ! mass of particle neglected
        dtldt0=abs(phigap*(1.0d0+phibep*cth))
        ethcm=phieth/dtldt0
      elseif (phitype.eq.1) then ! mass not neglected and elab<gamma m
        if (cth.le.phicthmax) then
          elab=phieth+phimp
          tcm(1)=(elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp    ! kinetic energy in parent cm system
          tcm(1)=max(tcm(1),0.0d0) ! MG bug fix (cf. also below)
          dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &      sqrt(tcm(1)**2+2.0d0*tcm(1)*phimp)))

          tcm(2)=(elab+
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp    ! kinetic energy in parent cm system
          tcm(2)=max(tcm(2),0.0d0) ! MG bug fix (cf. also below)
          dtdt(2)=abs(phigap*(1.0d0+phibep*cth*(tcm(2)+phimp)/
     &      sqrt(tcm(2)**2+2.0d0*tcm(2)*phimp)))
c          write(*,*) 'tcm(1) = ',tcm(1),'  dtdt(1) = ',dtdt(1)
c          write(*,*) 'tcm(2) = ',tcm(2),'  dtdt(2) = ',dtdt(2)
        else
          dsandydth=0.0d0
c          write(*,*) 'cth = ',cth,' > phicthmax'
          return
        endif
      else  ! phitype.eq.2  ! mass not neglected, elab>gamma m
        elab=phieth+phimp
        tcm(1)=(elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp    ! kinetic energy in parent cm system
        tcm(1)=max(tcm(1),0.0d0) ! MG bug fix (cf. also below)
        dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &    sqrt(tcm(1)**2+2.0d0*tcm(1)*phimp)))
      endif

c... map to PDG code input required by dsanyield_sim
      yieldkk=mod(phifk,100)
      diff=phifk/100
      if (diff.ne.0.and.diff.ne.1) then
        write(*,*) 'ERROR in dsandydth: unsupported argument phifk =', phifk
        dsandydth=0.0d0
        return
      endif      
      if (yieldkk.eq.51) then
         yieldpdg = -11 ! positron yields
      elseif (yieldkk.eq.52) then 
         yieldpdg = 22  ! cont. gammas     
      elseif (yieldkk.eq.53) then 
         yieldpdg = 14  ! muon neutrinos     
      elseif (yieldkk.eq.54) then 
         yieldpdg = -2212 ! antiproton yields     
      elseif (yieldkk.eq.61) then 
         yieldpdg = -1000010020 ! anti-deuteron yields
         call dsanyield_dbset(61,-100.d0)
      elseif (yieldkk.eq.59) then 
         yieldpdg = -1000010020 ! anti-deuteron yields
         call dsanyield_dbset(59,-100.d0)
      elseif (yieldkk.eq.71) then 
         yieldpdg = 14  ! neutrino yields (same as 53)    
      elseif (yieldkk.eq.72) then 
         yieldpdg = 130072 ! muon yields at creation     
      elseif (yieldkk.eq.73) then 
         yieldpdg = 130073 ! integrated muon yields in ice     
      else
        write(*,*) 'ERROR in dsandydth: unspoorted argument phifk =', phifk
        dsandydth=0.0d0
        return
      endif


c... get the yield at this energy
      if (phitype.eq.0) then
        yield=dsanyield_sim(e1cm,ethcm,chi2pdg(phich),'0',yieldpdg,diff,istat)
        if (phifk.ge.100) then  ! differential yields
          jac=1.0d0/dtldt0
          yield=yield*jac
        endif
      elseif (phitype.eq.1) then
c       MG bug fix
c        write(*,*) 'dshadydth: cthcrit='
c     &              , -sqrt(1d0-elab**2/phigap**2/phimp**2)
c        write(*,*) 'dshadydth: cth=', cth
        if (cth.lt.(-sqrt(1d0-elab**2/phigap**2/phimp**2))) then
          ! tcm+mp > elab for 0<tcm<tcm(1) or tcm>tcm(2)
c          write(*,*) 'dshadydth: tcm(1)=', tcm(1), 'tcm(2)=', tcm(2)     
          yield1=dsanyield_sim(e1cm,tcm(1),chi2pdg(phich),'0',yieldpdg,diff,istat)
          yield2=dsanyield_sim(e1cm,tcm(2),chi2pdg(phich),'0',yieldpdg,diff,istat)
          if (phifk.ge.100) then  ! differential yields
            jac=1.0d0/dtdt(1)
            yield1=yield1*jac
            jac=1.0d0/dtdt(2)
            yield2=yield2*jac
            yield=yield1+yield2
          else ! integrated yield
            yield=dsanyield_sim(e1cm,1d-10,chi2pdg(phich),'0',yieldpdg,diff,istat)
     &            -yield1+yield2
          endif
        else ! tcm+mp > elab for all tcm>0
          if (phifk.ge.100) then  ! differential yields
            yield=0d0
          else ! integrated yield
            yield=dsanyield_sim(e1cm,1d-10,chi2pdg(phich),'0',yieldpdg,diff,istat)
          endif
        endif  
c        old code
c          yield1=dsanyield_sim(e1cm,tcm(1),chi2pdg(phich),'0',yieldpdg,diff,istat)
c          yield2=dsanyield_sim(e1cm,tcm(2),chi2pdg(phich),'0',yieldpdg,diff,istat)
c          if (phifk.ge.100) then  ! differential yields
c            jac=1.0d0/dtdt(1)
c            yield1=yield1*jac
c            jac=1.0d0/dtdt(2)
c            yield2=yield2*jac
c          endif
c        yield=yield1+yield2
      else  ! phytype.eq.2
        yield=dsanyield_sim(e1cm,tcm(1),chi2pdg(phich),'0',yieldpdg,diff,istat)
        if (phifk.ge.100) then  ! differential yields
          jac=1.0d0/dtdt(1)
          yield=yield*jac
        endif
      endif

c... take care of factors depending on energy in neutrino-nucleus cross
c... section and muon range
      if (phifk.eq.72.or.phifk.eq.172) then
        yield=yield*dtldt0   ! nu-nucleon cross section
      elseif (phifk.eq.73.or.phifk.eq.173) then
        yield=yield*(dtldt0**2)  ! nu-nucleon cross section and mu range
      endif

c      if (phifk.eq.72.or.phifk.eq.73) then
c        write(*,*) 'warning in dsandydth: integrated yields for',
c     &    ' channel ',phifk,' has to be checked!'
c      endif

c...  include a factor of 1/2 to get average of yield,
c...  (1/(int_-1^1 1 dcth)=1/2)
c...  another factor of 1/2 from that we only consider
c...  one particle.

      dsandydth=0.25d0*yield
c      write(*,*) 'cth = ',cth,'  dsandydth = ',0.25d0*yield

      end











