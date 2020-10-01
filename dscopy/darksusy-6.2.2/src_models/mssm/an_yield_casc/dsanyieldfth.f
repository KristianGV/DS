*****************************************************************************
*** function dsanyieldfth integrates dsandydth over cos theta.
*** it is the yield from particle 1 (which decays from m0)
*** that is calculated. particle one corresponds to channel ch.
*** units: (annihilation)**-1
*****************************************************************************

      real*8 function dsanyieldfth(e0,m0,mp1,mp2,emuthr,ch,yieldk,
     &  istat)
      implicit none
      include 'dsanyieldcom.h'
      include 'dsanyieldmodelcom.h'
      include 'dsidtag.h'
      include 'dsio.h'
      include 'dsmpconst.h'

c------------------------ variables ------------------------------------

      real*8 e0,m0,mp1,mp2,emuthr,m00,e00,mp1t,mp2t
      integer ch,istat,yieldk

      real*8 sum,dth,elab
      parameter (dth=0.001d0)
      logical wb,chok

c------------------------ functions ------------------------------------

      external dsandydth
      real*8 dsandydth,dsanyield_int

c-----------------------------------------------------------------------

      dsanyieldfth=0.d0
      wb=.true.
      m00=m0
      e00=e0

c...Fix mass for virtual particles
      if (e00.lt.m00) then
         m00=e00*0.9999d0
      endif

c...Set masses temporarily. If we have virtual decay products, we will
c...temporarily set them to have lower masses
      mp1t=mp1
      mp2t=mp2

      if (m00.lt.(mp1+mp2)) then
         mp1t=mp1*m00/(mp1+mp2)*0.9999d0
         mp2t=mp2*m00/(mp1+mp2)*0.9999d0
      endif

c...take care of slight numerical inaccuracy problems
c...Obsolte, we now extrapolate
c      if (m00.lt.(mp1+mp2)) then
c        if (m00.gt.0.99*(mp1+mp2)) then
c          m00=(mp1+mp2)*1.0001
c          if(e0.lt.m00) e00=m00
c        else
c          write(*,*) 'error in dsanyieldfth: a particle with mass ',m0
c          write(*,*) 'is let to decay to two particles with mass ',
c     &      mp1,' and ',mp2
c          write(*,*) 'which is not energetically allowed.'
c          write(*,*) 'particle 1 has code ',ch
c          write(*,*) 'the yield from this decay is put to 0.0'
c          write(*,*) 'model: ',idtag
c          dsanyieldfth=0.0
c          return
c        endif
c      endif

c...calculate boost parameters of parent for common block
      phibep=sqrt(e00**2-m00**2)/e00
      phigap=e00/m00

c...set common block variables
      phim0=m00
      phie0=e00
      phim1=mp1t
      phim2=mp2t
      phie1=(m00**2+mp1t**2-mp2t**2)/(2.0d0*m00) ! cm energy of 1
      phie2=(m00**2+mp2t**2-mp1t**2)/(2.0d0*m00) ! cm energy of 2
      phieth=emuthr
      phich=ch
      phifk=yieldk

      if (yieldk.eq.54.or.yieldk.eq.154) then ! antiproton
        phimp=m_p
        elab=phieth+phimp

        if (elab.le.phigap*phimp) then
          phitype=1
          phicthmax=-sqrt(max((phigap**2*phimp**2-elab**2)/
     &      (phimp**2*phigap**2*phibep**2),0.0d0))
          phicthmax=max(-1.0d0,phicthmax-dth)
          phicthmin=max((phieth-phie1*phigap)/
     &      (phigap*phibep*sqrt(phie1**2-phimp**2)),-1.0d0)
        else
          phitype=2
          phicthmax=1.0d0
          phicthmin=max((phieth-phie1*phigap)/
     &      (phigap*phibep*sqrt(phie1**2-phimp**2)),-1.0d0)
        endif
      elseif (yieldk.eq.dbflxk.or.yieldk.eq.(dbflxk+100)) then ! antideuteron
        phimp=m_d ! deuteron mass
        elab=phieth+phimp

        if (elab.le.phigap*phimp) then
          phitype=1
          phicthmax=-sqrt(max((phigap**2*phimp**2-elab**2)/
     &      (phimp**2*phigap**2*phibep**2),0.0d0))
          phicthmax=max(-1.0d0,phicthmax-dth)
          phicthmin=max((phieth-phie1*phigap)/
     &      (phigap*phibep*sqrt(phie1**2-phimp**2)),-1.0d0)
        else
          phitype=2
          phicthmax=1.0d0
          phicthmin=max((phieth-phie1*phigap)/
     &      (phigap*phibep*sqrt(phie1**2-phimp**2)),-1.0d0)
        endif
      else  ! all others where mass can be neglected
        phimp=0.0d0
        phitype=0
        phicthmax=1.0d0
        phicthmin=max((phieth-phie1*phigap)/
     &    (phie1*phigap*phibep),-1.0d0)
      endif

      if (phicthmin.ne.-1.0d0) then
        phicthmin=min(phicthmin+dth,1.0d0)
      endif

      if (yieldk.le.100) then ! integrated yields ok to go all the way
        phicthmax=1.0d0
      endif

      if (phicthmin.ge.phicthmax) then
        dsanyieldfth=0.0d0
        return
      endif


c...check if mneu within correct bounds
      chok = .true.
      if ((ch.eq.1.or.ch.eq.2.or.ch.eq.3.or.ch.eq.4
     &  .or.ch.eq.5.or.ch.eq.7.or.ch.eq.10.or.ch.eq.11)
     &   .and.phie1.lt.lb(ch))
     &  chok=.false.
      if ((ch.eq.6.or.ch.eq.8.or.ch.eq.9).and.phie1.lt.(0.99*lb(ch)))
     &  chok=.false.

      if (.not.chok) then
c        if (prtlevel.gt.0) then
c          write(6,*)
c          write(6,5000)
c     +     'warning in dsanyieldfth: a cm energy of ',
c     +      2.0d0*phie1,' gev wants to be used,'
c          write(6,5010) 'while the lower bound for channel ',ch,
c     +      ' is ',2.0d0*lb(ch),' gev.'
c          write(6,*) 'the yield is put to 0 for these too low energies.'
c          write(6,*) 'the results can thus only be trusted as a',
c     +      ' lower bound.'
c        endif
c        wb=.false.
c        dsanyieldfth=0.0d0
        istat=(istat/2)*2+1
      endif

      if (phie1.gt.ub(ch)) then
        if (prtlevel.gt.0) then
          write(6,*)
          write(6,5000) 'warning in dsanyieldfth: a cm energy of ',
     +      2.0d0*phie1,' gev wants to be used,'
          write(6,5010) 'while the upper bound for channel ',ch,
     +      ' is ',2.0d0*ub(ch),' gev.'
          write(6,5020) 'a cm energy of ',2.0d0*ub(ch),' gev is used',
     +      ' instead for these too high energies.'
          write(6,*) 'the results can thus only be trusted as a',
     +      ' lower bound.'
          endif
        istat=(istat/2)*2+1
      endif

      if (wb) then
        if (phicthmin.lt.0.99d0.and.phicthmax.gt.0.99d0) then
          if (dsandydth(0.99d0).eq.0.0d0) then
            sum=dsanyield_int(dsandydth,0.99d0,phicthmax)+
     &        dsanyield_int(dsandydth,phicthmin,0.99d0)
          else
            sum=dsanyield_int(dsandydth,phicthmin,phicthmax)
          endif
        else
          sum=dsanyield_int(dsandydth,phicthmin,phicthmax)
        endif

c        eps=0.01
c        sum=0.0
c        call gadap(0.0,thu,dsandydth,eps,sum)
c
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

        dsanyieldfth=sum  ! JE Corr 080114: 1.d-15 removed

      endif

 5000 format(' ',a,f8.2,a)
 5010 format(' ',a,i2,a,f8.2,a)
 5020 format(' ',a,f8.2,a,a)

      end
