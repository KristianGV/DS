      subroutine dssmconst_ckm
c_______________________________________________________________________
c  useful constants for CKM mixing
c  common:
c    'dssm.h' - file with SM common blocks
c  author: paolo gondolo 1994-1999
c  modified: 031105 neutrino's yukawa corrected (pg)
c  modified: 13-04-2019 moved to SM part (tb)
c=======================================================================
      implicit none
      include 'dssm.h'
      real*8 s12,s23,s13,c12,c13,c23,d13
      complex*16 ed13

      s12 = ckms12
      s23 = ckms23
      s13 = ckms13
      d13 = ckmdelta
      c12 = sqrt(1.d0-s12**2)
      c13 = sqrt(1.d0-s13**2)
      c23 = sqrt(1.d0-s23**2)
      ed13 = dcmplx(cos(d13),sin(d13))
      ckm(1,1) = c12*c13
      ckm(1,2) = s12*c13
      ckm(1,3) = s13/ed13
      ckm(2,1) = -s12*c23-c12*s23*s13*ed13
      ckm(2,2) = c12*c23-s12*s23*s13*ed13
      ckm(2,3) = s23*c13
      ckm(3,1) = s12*s23-c12*c23*s13*ed13
      ckm(3,2) = -c12*s23-s12*c23*s13*ed13
      ckm(3,3) = c23*c13

      return
      end
