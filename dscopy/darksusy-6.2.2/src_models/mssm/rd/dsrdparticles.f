      subroutine dsrdparticles(option,nsize,
     &  selfcon,ncoann,mcoann,dof,nrs,rm,rw,nthr,tm)
**********************************************************************
*** subroutine dsrdparticles returns which particles are included in
***   the calculation of the WIMP relic density (coannihilations,
***   resonances, thresholds)
***                                                                     
***  type : interface                                                   
***
***  desc : Particles included in relic density calculation
***  desc : (coannihilations, resonances and thresholds)
***
*** input:
***   option = 0 - no coann
***             1 - include all relevant coannihilations (charginos, 
***                 neutralinos and sleptons; Delta m/m<1.5)
***             2 - include only coannihilations betweeen charginos
***                 and neutralinos (Delta m/m<1.5)
***             3 - include only coannihilations between sfermions
***                 and the lightest neutralino (Delta m/m<1.5)
***           101 - include all relevant coannihilations (charginos, 
***                 neutralinos and sleptons; Delta m/m<2.1)
***           102 - include only coannihilations betweeen charginos
***                 and neutralinos (Delta m/m<2.1)
***           103 - include only coannihilations between sfermions
***                 and the lightest neutralino (Delta m/m<2.1)
***   nsize = size of arrays, do not fill larger than this      
***
*** output 
***    selfcon  - specifies whether DM is selfconjugate (1) or not (2) 
***    ncoann - number of particles coannihilating
***    mgev  - relic and coannihilating mass in gev
***    dof   - internal degrees of freedom of the particles
***    nrs   - number of resonances to take special care of
***    rm    - mass of resonances in gev
***    rw    - width of resonances in gev
***    nt    - number of thresholds to take special care of
***            do not include coannihilation thresholds (that's automatic)
***    tm    - sqrt(s) of the thresholds in gev
*** This output is returned (used by the RD routines in src/)
*** and also, together with kcoann, put in common block needed by dsandwdcos.
*** author: paolo gondolo
*** derived from dsrdomega      
*** date: 13-10-04
*** Modified: Joakim Edsjo, edsjo@fysik.su.se, 2015-12-08      
*** Modified: Torsten Bringmann 2018-02-28 (added selfcon)
**********************************************************************

      implicit none
      include 'dsmssm.h'
      include 'dsandwcom.h'
      integer option

      integer selfcon, ncoann,nrs,nthr,i,j,nsize
      real*8 mcoann(*),dof(*),tm(*),
     &     rm(*),rw(*)
      real*8 tmp
      integer kcoann(nsize),ktmp
      integer coproc
      real*8 mcofr

c----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dsrdparticles')

      if (nsize.lt.comax) then
         write(*,*) 'DS ERROR in dsparticles:'
         write(*,*) 'Size of array nsize=',nsize,
     &     ' in call to dsparticles.'
         write(*,*) 'Increase nsize (thavsiz). Program stopping.'
         stop
      endif
      
c... so far, the MSSM module only handles neutralino LSPs. This is
c... a self-conjugate particle
      selfcon = 1

      
c...This routine performs the following
c...1) goes through all particles in the MSSM that we want to
c...include in the coannihilation setup, resonances and thresholds
c...2) Return these (for use in dsrdomega in src), and
c...3) passes them, and kcoann, to common blocks needed for dsandwdcos.      
      
c...Check which coannihilation processes to include
c...The coproc switch determines which processes are included
c...Bit Decimal Coannihilation processes included
c...  0       1 neu_1, neu_2, neu_3, neu_4, cha_1, cha_2
c...  1       2 sel_1,sel_2,smu_1,smu_2,stau_1,stau_2
      mcofr=1.0d0
      if (option.eq.0) then      ! no coannihilations
        coproc=0
      elseif (option.eq.1) then  ! include all coannihilations
        coproc=3
        mcofr=1.5d0
      elseif (option.eq.2) then  ! include only charg. and neu coanns
        coproc=1
        mcofr=1.5d0
      elseif (option.eq.3) then  ! include only slepton coanns
        coproc=2
        mcofr=1.5d0
      elseif (option.eq.101) then  ! include all coannihilations
        coproc=3
        mcofr=2.1d0
      elseif (option.eq.102) then  ! include only charg. and neu coanns
        coproc=1
        mcofr=2.1d0
      elseif (option.eq.103) then  ! include only slepton coanns
        coproc=2
        mcofr=2.1d0
      else
        write(*,*) 'ERROR in dsrdomega: invalid option=',option
        stop
      endif

c...Add lightest neutralino to annihilation list
      ncoann=1
      mcoann(1)=mass(kn(1))
      dof(1)=kdof(kn(1))
      kcoann(1)=kn(1)

c...Add other neutralinos and charginos if chosen to be included
      if (and(coproc,1).ne.0) then
         mcoann(2)=mass(kn(2))
         mcoann(3)=mass(kn(3))
         mcoann(4)=mass(kn(4))
         dof(2)=kdof(kn(2))
         dof(3)=kdof(kn(3))
         dof(4)=kdof(kn(4))
         if (mcoann(2)/mcoann(1).le.mcofr) then
            ncoann=2
            kcoann(ncoann)=kn(2)
         endif
         if (mcoann(3)/mcoann(1).le.mcofr) then
            ncoann=3
            kcoann(ncoann)=kn(3)
         endif
         if (mcoann(4)/mcoann(1).le.mcofr) then
            ncoann=4
            kcoann(ncoann)=kn(4)
         endif
         if (mass(kcha(1))/mcoann(1).le.mcofr) then
            ncoann=ncoann+1
            mcoann(ncoann)=mass(kcha(1))
            dof(ncoann)=kdof(kcha(1))
            kcoann(ncoann)=kcha(1)
         endif
         if (mass(kcha(2))/mcoann(1).le.mcofr) then
            ncoann=ncoann+1
            mcoann(ncoann)=mass(kcha(2))
            dof(ncoann)=kdof(kcha(2))
            kcoann(ncoann)=kcha(2)
         endif
      endif

c...Add slepton coannihilations if chosen to be included
      if(and(coproc,2).ne.0) then

c...Sleptons
         do i=1,6
            if (mass(ksl(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Slepton ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksl(i))
               dof(ncoann)=kdof(ksl(i))
               kcoann(ncoann)=ksl(i)
            endif
         enddo

c...Sneutrinos
         do i=1,3
            if (mass(ksnu(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Sneutrino ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksnu(i))
               dof(ncoann)=kdof(ksnu(i))
               kcoann(ncoann)=ksnu(i)
            endif
         enddo
c...Up squarks
         do i=1,6
            if (mass(ksqu(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Up squark ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksqu(i))
               dof(ncoann)=kdof(ksqu(i))
               kcoann(ncoann)=ksqu(i)
            endif
         enddo

c...Down squarks
         do i=1,6
            if (mass(ksqd(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Down squark ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksqd(i))
               dof(ncoann)=kdof(ksqd(i))
               kcoann(ncoann)=ksqd(i)
            endif
         enddo

      endif

c...add resonances for new call
      nrs=0
c...z resonance
      if (mass(kz).gt.mcoann(1)*2.0d0) then
        nrs=nrs+1
        rm(nrs)=mass(kz)
        rw(nrs)=width(kz)
      endif
c...h1 resonance
      if (mass(kh1).gt.mcoann(1)*2.0d0) then
        nrs=nrs+1
        rm(nrs)=mass(kh1)
        rw(nrs)=width(kh1)
      endif
c...h2 resonance
      if (mass(kh2).gt.mcoann(1)*2.0d0) then
        nrs=nrs+1
        rm(nrs)=mass(kh2)
        rw(nrs)=width(kh2)
      endif
c...h3 resonance
      if (mass(kh3).gt.mcoann(1)*2.0d0) then
        nrs=nrs+1
        rm(nrs)=mass(kh3)
        rw(nrs)=width(kh3)
      endif
c...coannihilation-resonances
      if (ncoann.gt.1) then
        if (mass(khc).gt.mcoann(1)*2.0d0) then   ! sufficient condition
          nrs=nrs+1
          rm(nrs)=mass(khc)
          rw(nrs)=width(khc)
        endif
        if (mass(kw).gt.mcoann(1)*2.0d0) then ! sufficient condition
          nrs=nrs+1
          rm(nrs)=mass(kw)
          rw(nrs)=width(kw)
        endif
      endif

c...add thresholds for new call
      nthr=0
c...ww-threshold
      if (mass(kw).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(kw)
      endif
c...zz-threshold
      if (mass(kz).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(kz)
      endif
c...tt-bar-threshold
      if (mass(kt).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(kt)
      endif
c...note that coannihilation thresholds are automatically added in
c...dsrdstart.

c...We have now set up everything needed for the RD routines
c...Now make sure to transfer what we need to the annihilation rate
c...routines in the mssm module as well. These go into common blocks
      nco=ncoann
      do i=1,nco
         mco(i)=mcoann(i)
         mdof(i)=dof(i)
         kco(i)=kcoann(i)
      enddo

c...sort
      if (nco.ge.2) then
        do i=1,nco-1
          do j=nco-1,i,-1
            if (mco(j).gt.mco(j+1)) then
              tmp=mco(j+1)
              mco(j+1)=mco(j)
              mco(j)=tmp
              tmp=mdof(j+1)
              mdof(j+1)=mdof(j)
              mdof(j)=tmp
              ktmp=kco(j+1)
              kco(j+1)=kco(j)
              kco(j)=ktmp
            endif
          enddo
        enddo
      endif
      return

      end
