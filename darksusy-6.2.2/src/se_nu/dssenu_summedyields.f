*****************************************************************************
*** function dssenu_summedyields gives the total yield of muons 
*** (mu- and mu+) above threshold for
*** a given WIMP mass or the differential muon yield for a given
*** energy and a given angle. 
*** kind =  1 for integrated yields (above emu and below theta)
***         2 for differential yileds (in energy and angle)
***         3 for mixed yields, diff. in energy, integrated up to theta      
*** rtype = 1 for neutrino (nu_mu +/or nu_mu-bar) yields
***         2 for mu- +/or mu+ yields at neutrino-nucleon vertex (cont. events)
***         3 for mu- +/or mu+ yields at imaginary plane in detector.
***            (through-going events)
*** ptype = 1 for particles only (nu_mu or mu-)
***         2 for anti-particles only (nu_mu-bar or mu+)
***         3 for summed yields (both particles and anti-particles)
*** wh='su' corresponds to annihilation in the sun and wh='ea' corresponds
*** to annihilation in the earth. if istat=1 upon return,
*** some inaccesible parts the differential muon spectra has been wanted,
*** and the returned yield should then be treated as a lower bound.
*** if istat=2 energetically forbidden annihilation channels have been
*** wanted. if istat=3 both of these things has happened.
*** units: 1.0e-30 m**-2 (annihilation)**-1 for integrated yield.
***        1.0e-30 m**-2 gev**-1 (degree)^-1 (annihilation)**-1 for
***        differential yield.
*** An additionl unit of m**-1 for rtype 1-2.
*** author: joakim edsjo  edsjo@fysik.su.se
*** date: 96-03-19
*** modified: 96-09-03 to include new index order
*** modified: 97-12-03 to include new muyield routines (v3.21)
*** Modified: 08-04-02 to use new wimpsim routines with neutrino
*** oscillations etc.
*** Modified: 11-04-23 to allow rates to be calculated for only
***  particles or only antiparticles (pat scott, patscott@physics.mcgill.ca)
*****************************************************************************

      real*8 function dssenu_summedyields(emu,theta,wh,kind,rtype,
     & ptype,istat)
      implicit none
      include 'dsseyieldcom.h'

c------------------------ functions ------------------------------------

      real*8 dsseyield

c------------------------ variables ------------------------------------

      real*8 emu,yield,theta
      integer istat,kind,rtype,ptype,t1,t2
      character*2 wh

c----------------------------------------------- set-up common variables

c...loop through different channels and calculate the yield above threshold
c...for each channel.
c      call wirate(6,6,1)

      if (rtype.eq.1) then
         t1=3
         t2=4
      elseif (rtype.eq.2) then
         t1=9
         t2=10
      elseif (rtype.eq.3) then
         t1=13
         t2=14
      else
         write(*,*) 'DS ERROR in dssenu_summedyields. Wrong rtype = ',
     &     rtype
         stop
      endif
         
      seistat=0

      yield = 0.d0
      if (ptype.eq.1 .or. ptype.eq.3) then
        yield=dsseyield(emu,theta,wh,kind,t1,istat)
        seistat=or(seistat,istat)
      endif

      if (ptype.eq.2 .or. ptype.eq.3) then
        yield=yield+dsseyield(emu,theta,wh,kind,t2,istat)
        seistat=or(seistat,istat)
      endif
      
      dssenu_summedyields=yield
      istat=seistat

      end


















