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
**** author: Torsten Bringmann 
*** (derived from mssm/dsrdparticles)
*** date: 14-04-06
**********************************************************************

      implicit none
      include 'dsempty.h'
      integer option

      integer selfcon,ncoann,nrs,nthr,nsize
      real*8 mcoann(*),dof(*),tm(*),
     &     rm(*),rw(*)
      integer kcoann(nsize)

c----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsrdparticles')

c...assume DM to be self-conjugate
      selfcon = 1

c...Add DM particle to annihilation list
      ncoann=1
      mcoann(1)=mass(kdm)
      dof(1)=kdof(kdm)
      kcoann(1)=kdm

      nrs=0
      nthr=0


      return

      end
