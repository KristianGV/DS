      subroutine random_model()
c
c     To generate MSSM7 model parameters in a random way
c
      implicit none
      include 'dsidtag.h'
      include 'dsmssm.h'
      real*8 dsrndlog,dsrndlin,dsrndsgn
      integer first,n,idum,i
      real*8 msq,atm,abm
      real*8 mumin,mumax,m1min,m1max,m2min,m2max,mamin,mamax,
     &     tbmin,tbmax,msqmin,msqmax,atmmin,atmmax,abmmin,abmmax
      real*8 dsgf2s2thw
      data first/0/, n/0/, idum/0/
      save first,n,idum
c... First time set random number seed
      if (first.eq.0) then
         idum = -3782932
         n=0
         first=1
      endif
c...  Limits of the Mu-M2-MA-tan(beta)-Msq-At-Ab region explored
      mumin = 10.d0
      mumax = 10000.d0
      m2min = 10.d0
      m2max = 10000.d0
      m1min = 10.d0
      m1max = 10000.d0
      mamin = 10.d0
      mamax = 1000.d0
      tbmin = 1.001d0
      tbmax = 60.d0
      msqmin = 50.d0
      msqmax = 1000.d0
      atmmin = -3.d0
      atmmax = 3.d0
      abmmin = -3.d0
      abmmax = 3.d0
c... random point
      n=n+1
      write (idtag,'(a,i8.8)') 'rndm',n
      mu = dsrndsgn(idum) * dsrndlog(idum,abs(mumin),abs(mumax))
      m2 = dsrndsgn(idum) * dsrndlog(idum,abs(m2min),abs(m2max))
      s2wmz=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),3)
      m1 = m2 * 5.0d0/3.0d0 * s2wmz/(1.d0-s2wmz)
      ma = dsrndlog(idum,mamin,mamax)
      tanbe = dsrndlog(idum,tbmin,tbmax)
      msq = dsrndlog(idum,msqmin,msqmax)
      atm = dsrndlin(idum,atmmin,atmmax)
      abm = dsrndlin(idum,abmmin,abmmax)
c... further assumptions
      ! W mass for unitarity of tree-level annihilation amplitudes
      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      ! GUT relation on gaugino mass parameters
      m3 =  m2 * (alph3mz*s2wmz)/alphem    ! alph3->alph3mz 020912 (JE)
      ! Slepton and squark mass matrices at the weak scale
      mass2q(3) = msq**2
      mass2u(3) = msq**2
      mass2d(3) = 2.0d0*mass2q(3)-mass2u(3)
      mass2l(3) = mass2d(3)
      mass2e(3) = mass2d(3)
      asofte(3) = 0.d0
      asoftu(3) = atm*msq
      asoftd(3) = abm*msq
      do i=1,2
         mass2q(i) = mass2d(3)
         mass2u(i) = mass2d(3)
         mass2d(i) = mass2d(3)
         mass2l(i) = mass2d(3)
         mass2e(i) = mass2d(3)
         asofte(i) = 0.d0
         asoftu(i) = 0.d0
         asoftd(i) = 0.d0
      enddo
      return
      end

