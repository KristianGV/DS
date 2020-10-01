*****************************************************************************
*** function dsanyielddec integrates dsandys over cos theta.
*** the scalar decay channels are summed in dsandys
*** This routine is for the fundamental channels from S1, S2 and S3 decay.
*** units: (annihilation)**-1
*****************************************************************************

      real*8 function dsanyielddec(eh,hno,egev,yieldk,istat)
      implicit none
      include 'dsanyieldmodelcom.h'
      include 'dsidtag.h'
      include 'dsmpconst.h'

c------------------------ variables ------------------------------------

      real*8 eh,egev,m0
      real*8 elab,e1
      integer istat,yieldk,hno

      real*8 sum,dth
      parameter (dth=0.00001d0)  ! safety when integrating in theta

c------------------------ functions ------------------------------------

      external dsandys
      real*8 dsandys,dsanyield_int

c-----------------------------------------------------------------------

c...Allow for virtual particles
      m0=ans0m(hno)
      m0=min(m0,eh*0.9999d0) ! if virtual


c...set common block variables: boost parameters etc.
      phim0=m0  ! this is the rescaled mass if virtual (see above)
      phibep=sqrt(eh**2-m0**2)/eh
      phigap=eh/m0
      phihno=hno
      phifk=yieldk
      phieth=egev
      e1=m0/2.0d0

      if (yieldk.eq.54.or.yieldk.eq.154) then ! antiproton
        phimp=m_p
        elab=phieth+phimp

        if (elab.le.phigap*phimp) then
          phitype=1
          phicthmax=-sqrt(max((phigap**2*phimp**2-elab**2)/
     &      (phimp**2*phigap**2*phibep**2),0.0d0))
          phicthmax=max(-1.0d0,phicthmax-dth)
          phicthmin=max((phieth-e1*phigap)/
     &      (phigap*phibep*sqrt(e1**2-phimp**2)),-1.0d0)
        else
          phitype=2
          phicthmax=1.0d0
          phicthmin=max((phieth-e1*phigap)/
     &      (phigap*phibep*sqrt(e1**2-phimp**2)),-1.0d0)
        endif
      else  ! all others where mass can be neglected
        phimp=0.0d0
        phitype=0
        phicthmax=1.0d0
        phicthmin=max((phieth-e1*phigap)/
     &    (e1*phigap*phibep),-1.0d0)
      endif

      if (phicthmin.ne.-1.0d0) then
        phicthmin=min(phicthmin+dth,1.0d0)
      endif

      if (yieldk.le.100) then ! integrated yields ok to go all the way
        phicthmax=1.0d0
      endif

      if (phicthmin.ge.phicthmax) then
        dsanyielddec=0.0d0
        return
      endif

c      write(*,*) 'phitype = ',phitype
c      write(*,*) 'phicthmin = ',phicthmin
c      write(*,*) 'phicthmax = ',phicthmax

        sum=dsanyield_int(dsandys,phicthmin,phicthmax)

c        eps=1.0d-3
c        call dgadap(phicthmin,phicthmax,dsandys,eps,sum)

c        if (sum.lt.1.0.and.sum.gt.0.0) then
c          write(6,*) 'sum less then 1.0 in gadap (phiith)'
c          eps=0.01*sum
c          eps0=eps
c  100     call gadap(0.0,thu,dsandydth,eps,sum)
c          if ((eps0/sum).gt.0.01) then
c            eps0=eps0/2.0
c            eps=eps0
c          endif
c          if ((eps0*2/sum).gt.0.01) goto 100
c        endif

c        dsanyielddec=sum

        dsanyielddec=sum
c        write(*,*) 'dsanyielddec called with'
c        write(*,*) '    e0 = ',e0
c        write(*,*) '    m0 = ',m0
c        write(*,*) '    egevr = ',egevr
c        write(*,*) '    mp1 = ',mp1
c        write(*,*) '    mp2 = ',mp2
c        write(*,*) '    ch = ',ch
c        write(*,*) '  result is ',sum

      end




