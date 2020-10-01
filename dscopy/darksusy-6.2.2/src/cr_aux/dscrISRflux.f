*******************************************************************************
*** Function dscrISRflux provides the differential local interstellar       ***
*** cosmic ray flux, per kinetic energy, computed from the intensity        ***
*** provided by dscrISRintensity.                                           *** 
***                                                                         ***
***  Input:                                                                 ***
***    Tkin  - kinetic energy of CR particle [GeV]                          ***
***    CRtype - 1 for protons                                               ***
***             2 for helium                                                ***
***                                                                         ***
***                                                                         ***
***  Units of output: 1/ [cm**2 s GeV]                                      ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-22                                                         ***
*******************************************************************************
      real*8 function dscrISRflux(Tkin, CRtype)
      implicit none
      include 'dsio.h'
      include 'dsmpconst.h'

      real*8 Tkin, mcr, rcr, qcr
      real*8 dscrISRintensity
      integer CRtype
      
      dscrISRflux=0.0d0
      
      if (CRtype.eq.1) then ! protons
        mcr = m_p
        qcr = 1.0d0
      elseif (CRtype.eq.2) then ! helium
        mcr = m_he
        qcr = 2.0d0
      else
        if (prtlevel.gt.1) then
          write(*,*) 'WARNING: dscrISRflux called with' 
          write(*,*) 'unrecognized parameter CRtype = ',CRtype 
          write(*,*) 'I will return zero...'
        endif
        return
      endif 

      rcr = sqrt(Tkin**2+2*Tkin*mcr)/qcr ! rigidity -> need to divide by charge
                                         ! TB, corrected 08/12/2019
      dscrISRflux = 1.d-4*4*pi*
     &              (Tkin + mcr)/(qcr**2*rcr)*dscrISRintensity(rcr)

      return
      end


