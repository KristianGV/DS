      subroutine dssuconst_yukawa
c_______________________________________________________________________
c  useful constants
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa corrected (pg)
c  modified: 081203 uncomment below to use masses run to the
c    electroweak scale (not to 2*m_lsp like
c    dssuconst_yukawa_running).  This requires
c    extra or external code to set the values of the
c    running masses. pat scott
c    WARNING: see comment in header of dssuconst_yukawa_running.
c             (Torsten Bringmann)
c=======================================================================
      implicit none
      include 'dsmssm.h'
      real*8 aux
      integer i

      cosbe = 1.0d0 / sqrt(1.0d0+tanbe*tanbe)
      sinbe = tanbe * cosbe
      cos2be = 2.0*cosbe*cosbe-1.0d0
      sin2be = 2.0*sinbe*cosbe

      aux = g2weak/dsqrt(2.d0)/mass(kw)
      do i=1,3
         yukawa(knu(i)) = aux*mass(knu(i))/sinbe
         yukawa(kl(i)) = aux*mass(kl(i))/cosbe
         yukawa(kqu(i)) = aux*mass(kqu(i))/sinbe
         yukawa(kqd(i)) = aux*mass(kqd(i))/cosbe
c        Uncomment below for running masses
c         yukawa(knu(i)) = aux*runmass(knu(i))/sinbe
c         yukawa(kl(i)) = aux*runmass(kl(i))/cosbe
c         yukawa(kqu(i)) = aux*runmass(kqu(i))/sinbe
c         yukawa(kqd(i)) = aux*runmass(kqd(i))/cosbe
      enddo

      end
