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
      include 'dsvdSIDM.h'
      integer option

      integer selfcon, ncoann,nrs,nthr,nsize
      real*8 mcoann(*),dof(*),tm(*),
     &     rm(*),rw(*)
      integer kcoann(nsize)

c----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dsrdparticles')


c...Determine whether DM is self-conjugate
      if (DMtype.eq.1) then ! Dirac DM
        selfcon = 2
      else
        selfcon = 1
      endif

c...Add DM particle to annihilation list
      ncoann=1
      mcoann(1)=mass(kdm)
      dof(1)=kdof(kdm)
      kcoann(1)=kdm

      nrs=0
      nthr=0
      
c... a threshold (as well as resonance) can occur if the DM particle
c... is lighter than the mediator particle
      if (mass(kdm).lt.mass(kmed)) then
        nthr=1
        tm(nthr)=2.0d0*mass(kmed)
        nrs=1
        rm(nrs)=mass(kmed)
        rw(nrs)=width(kmed)
      endif

c... add thresholds at very small velocities, just to make sure this
c... is well enough sampled for large Sommerfeld-enhancement
      nthr=nthr+2
      tm(nthr-1) = 2.0d0*mass(kdm)*(1.0d0 + 0.0001**2)
      tm(nthr) = 2.0d0*mass(kdm)*(1.0d0 + 0.001**2)

      return

      end
