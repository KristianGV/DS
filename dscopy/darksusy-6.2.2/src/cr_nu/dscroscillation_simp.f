      real*8 function dscroscillation_simp(nui,phinusource)
*******************************************************************************
***   function dscroscillation_simp returns the average neutrino flux of a  ***
***   given flavour at the surface of the earth, given the neutrino         ***
***   fluxes of all flavours at sorce, assuming a distant location in the   ***
***   galactic galo. This simplified version implements the expressions for ***
***   long baseline oscillations (from hep-ph/0506165), with fixed          ***
***   oscillation parameters as in 1506.02657.                              *** 
***                                                                         ***
***   input: nui         - return neutrino type (1- nue, 2- numu, 3-tau)    ***
***          phinusource - neutrino PLUS antineutrino fluxes at source      ***
***                        (all three flavours, in the same order)          ***   
***                                                                         ***
***   unit of return value: same as input source fluxes                     ***
***                                                                         ***
***   author: Torsten.Bringmann@fys.uio.no                                  ***
***   date: 2019-11-15                                                      ***
*******************************************************************************
      implicit none
      real*8 phinusource(3)
      integer nui

      integer j
      real*8 res, Pab(3,3)
      data Pab/0.573, 0.277, 0.150,    ! Eq. after (6) in 1506.02657
     &         0.277, 0.348, 0.375,    ! NB: column-first order
     &         0.150, 0.375, 0.475/    

      dscroscillation_simp = 0.0d0
      if (nui.lt.1.or.nui.gt.3) return

      res = 0.0d0
      do j = 1,3
        res = res + Pab(nui,j)*phinusource(j)
      enddo

      dscroscillation_simp = res      
    
      return
      end


