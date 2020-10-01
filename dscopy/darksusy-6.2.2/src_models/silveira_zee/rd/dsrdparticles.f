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
***   option = 0 - default
***            !=0 - include various subsets of coannihilations 
***   nsize = size of arrays, do not fill larger than this      
***
*** output 
***    selfcon  - specifies whether DM is selfconjugate (1) or not (2) 
***    ncoann   - number of particles coannihilating
***    mcoann   - relic and coannihilating mass in gev
***    dof      - internal degrees of freedom of the particles
***    nrs      - number of resonances to take special care of
***    rm       - mass of resonances in gev
***    rw       - width of resonances in gev
***    nt       - number of thresholds to take special care of
***               do not include coannihilation thresholds (that's automatic)
***    tm       - sqrt(s) of the thresholds in gev
*** This output is used by the RD routines in src/rd/
*** author: paolo gondolo
*** derived from dsrdomega      
*** date: 13-10-04
*** Modified: Joakim Edsjo, edsjo@fysik.su.se, 2015-12-08      
*** Modified: Torsten Bringmann 2016-06-28
*** Modified: Torsten Bringmann 2018-02-28 (added selfcon)
**********************************************************************

      implicit none
      include 'dssilveira_zee.h'
      integer option

      integer selfcon, ncoann,nrs,nthr,nsize
      real*8 mcoann(*),dof(*),tm(*),rm(*),rw(*)
      integer kcoann(nsize)
      integer coproc

c----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dsrdparticles')

c...the scalar singlet is a self-conjugate particle
      selfcon = 1


c...This routine performs the following
c...1) goes through all particles in the model that we want to
c...include in the coannihilation setup, resonances and thresholds
c...2) Return these (for use in dsrdomega in src), and
c...3) passes them, and kcoann, to common blocks needed for dsandwdcos.      
      
c...Check which coannihilation processes to include
c...The coproc switch determines which processes are included
c      if (option.eq.0) then      ! no coannihilations
        coproc=0
c      else
c        write(*,*) 'ERROR in dsrdparticles: invalid option=',option
c        stop
c      endif

c...Add WIMP to annihilation list
      ncoann=1
      mcoann(1)=mass(kdm)
      dof(1)=kdof(kdm)
      kcoann(1)=kdm

c...add resonances for new call
      nrs=0
c...higgs resonance
      if (mass(khsm).gt.mcoann(1)*2.0d0) then
        nrs=nrs+1
        rm(nrs)=mass(khsm)
        rw(nrs)=width(khsm)
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
c...hh-threshold
      if (mass(khsm).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(khsm)
      endif

c...note that coannihilation thresholds are automatically added in
c...dsrdstart.

      return
      end
      
