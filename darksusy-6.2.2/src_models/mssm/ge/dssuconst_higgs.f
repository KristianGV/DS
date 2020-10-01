      subroutine dssuconst_higgs
c_______________________________________________________________________
c  useful constants
c  common:
c    'dsmssm.h' - file with susy common blocks
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa corrected (pg)
c=======================================================================
      implicit none
      include 'dsmssm.h'

c...Default Higgs trilinear couplings
      lam1 = 0.25d0*(g2weak**2+gyweak**2)
      lam2 = lam1
      lam3 = 0.25d0*(g2weak**2-gyweak**2)
      lam4 = -0.5d0*g2weak**2
      lam5 = 0.d0
      lam6 = 0.d0
      lam7 = 0.d0

      end
