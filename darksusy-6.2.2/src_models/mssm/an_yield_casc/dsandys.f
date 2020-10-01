*****************************************************************************
*** function dsandys is the differential yield dyield/dcostheta in the
*** cm system boosted to the lab system (including proper jacobians if
*** we are dealing with a differential yield. all decay channels of the
*** scalar boson in question are summed.
*** the function should be integrated from -1 to 1, by e.g.
*** the routine gadap.
*** units: (annihilation)**-1
*** author: joakim edsjo (edsjo@physto.se)
*** date: 1998
*** modified: 98-04-15
*****************************************************************************

      real*8 function dsandys(cth)
      implicit none
      include 'dsanyieldmodelcom.h'
      include 'dsanyieldcom.h'

c------------------------ variables ------------------------------------

      real*8 cth,yield,jac,elab,yield1,yield2,mp,ep
      real*8 e1cm,ethcm,dtldt0,tcm(2),dtdt(2)
      integer istat,hdi(11),ch
      integer yieldkk, yieldpdg, diff

      data hdi/21,20,23,22,25,24,26,13,12,17,19/ ! fundamental -> Higgs decay

c------------------------ constants ------------------------------------

      real*8 tmin,pmin
      parameter(tmin=1.0d-5)   ! to avoid numerical roundoff problems
      parameter(pmin=1.0d-5)   ! to avoid numerical roundoff problems

c------------------------ functions ------------------------------------
      real*8 dsanyield_sim
c-----------------------------------------------------------------------

c...calculate the energy of particle 1 in the cm system
      e1cm=phim0/2.0d0 ! JE: pick virtual mass

c... calculate the cm energy for this angle that corresponds to requested
c... lab energy phieth (m_mu, m_ep etc neglected, but not m_p)
      dtldt0=0.d0
      if (phitype.eq.0) then ! mass of particle neglected
        dtldt0=abs(phigap*(1.0d0+phibep*cth))
        ethcm=phieth/dtldt0
      elseif (phitype.eq.1) then ! mass not neglected and elab<gamma m
        if (cth.le.phicthmax) then
          elab=phieth+phimp
          tcm(1)=max((elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp,tmin)    ! kinetic energy in parent cm system
          tcm(2)=max((elab+
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp,tmin)    ! kinetic energy in parent cm system

          if (phifk.ge.100) then ! differential yields
            dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &        sqrt(max(tcm(1)**2+2.0d0*tcm(1)*phimp,pmin**2))))
            dtdt(2)=abs(phigap*(1.0d0+phibep*cth*(tcm(2)+phimp)/
     &        sqrt(max(tcm(2)**2+2.0d0*tcm(2)*phimp,pmin**2))))
          else
            tcm(1)=max(tmin,tcm(1))
            tcm(2)=max(tmin,tcm(2))
          endif
c          write(*,*) 'tcm(1) = ',tcm(1),'  dtdt(1) = ',dtdt(1)
c          write(*,*) 'tcm(2) = ',tcm(2),'  dtdt(2) = ',dtdt(2)
        else
          dsandys=0.0d0
c          write(*,*) 'cth = ',cth,' > phicthmax'
          return
        endif
      else  ! phitype.eq.2  ! mass not neglected, elab>gamma m
        elab=phieth+phimp
        tcm(1)=max((elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp,tmin)    ! kinetic energy in parent cm system
        dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &    sqrt(max(tcm(1)**2+2.0d0*tcm(1)*phimp,pmin**2))))
      endif

c... get the yield at this energy
      yield=0.0d0

      do ch=1,11
        mp=msim(ch)
        ep=e1cm
        if (ans0br(hdi(ch),phihno).gt.0.0d0) then

c          if(ep.lt.0.99d0*mp) then
c            write(*,*) 'error in dsandys: ep < 0.99 mp'
c            write(*,*) '  channel = ',ch
c            write(*,*) '  mp = ',mp
c            write(*,*) '  ep = ',ep
c          else
c            ep=max(ep,mp*1.0001d0)
c          endif


          if (ep.gt.ub(ch)) then
            write(*,*) 'error in dsandys: too high ep = ',ep
          endif

c          if (ep.lt.lb(ch)) goto 150

c... map to PDG code input required by dsanyield_sim
      yieldkk=mod(phifk,100)
      diff=phifk/100
      if (diff.ne.0.and.diff.ne.1) then
        write(*,*) 'ERROR in dsandys: unsuppported argument phifk =', phifk
        dsandys=0.0d0
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
         yieldpdg = -1000010020 ! anti-deuteron yield
         call dsanyield_dbset(61,-100.d0)
      elseif (yieldkk.eq.59) then 
         yieldpdg = -1000010020 ! anti-deuteron yield, old sph coal mod
         call dsanyield_dbset(59,-100.d0)
      elseif (yieldkk.eq.71) then 
         yieldpdg = 14  ! neutrino yields (same as 53)    
      elseif (yieldkk.eq.72) then 
         yieldpdg = 130072 ! muon yields at creation     
      elseif (yieldkk.eq.73) then 
         yieldpdg = 130073 ! integrated muon yields in ice     
      else
        write(*,*) 'ERROR in dsandys: unsupported argument phifk =', phifk
        dsandys=0.0d0
        return
      endif

          if (phitype.eq.0) then
            yield1=ans0br(hdi(ch),phihno)*
     &        dsanyield_sim(e1cm,ethcm,chi2pdg(ch),'0',yieldpdg,diff,istat)
            if (phifk.ge.100) then  ! differential yields
              jac=1.0d0/dtldt0
              yield1=yield1*jac
            endif
            yield=yield+yield1
          elseif (phitype.eq.1) then
            yield1=dsanyield_sim(e1cm,tcm(1),chi2pdg(ch),'0',yieldpdg,diff,istat)
            yield2=dsanyield_sim(e1cm,tcm(2),chi2pdg(ch),'0',yieldpdg,diff,istat)
            if (phifk.ge.100) then  ! differential yields
              jac=1.0d0/dtdt(1)
              yield1=yield1*jac
              jac=1.0d0/dtdt(2)
              yield2=yield2*jac
            endif
            yield=yield+ans0br(hdi(ch),phihno)*(yield1+yield2)
          else  ! phitype.eq.2
            yield1=ans0br(hdi(ch),phihno)
     &        *dsanyield_sim(e1cm,tcm(1),chi2pdg(ch),'0',yieldpdg,diff,istat)
            if (phifk.ge.100) then  ! differential yields
              jac=1.0d0/dtdt(1)
              yield1=yield1*jac
            endif
            yield=yield+yield1
          endif
        goto 150  
  150   continue
        endif
      enddo

c... take care of factors depending on energy in neutrino-nucleus cross
c... section and muon range
      if (phifk.eq.72.or.phifk.eq.172) then
        yield=yield*dtldt0   ! nu-nucleon cross section
      elseif (phifk.eq.73.or.phifk.eq.173) then
        yield=yield*(dtldt0**2)  ! nu-nucleon cross section and mu range
      endif

c      if (phifk.eq.72.or.phifk.eq.73) then
c        write(*,*) 'warning in dsandys: integrated yields for',
c     &    ' channel ',phifk,' has to be checked!'
c      endif

c...  include a factor of 1/2 to get average of yield,
c...  (1/(int_-1^1 1 dcth)=1/2).

      dsandys=0.5d0*yield
c      write(*,*) 'cth = ',cth,'  dsandys = ',0.25d0*yield

      end











