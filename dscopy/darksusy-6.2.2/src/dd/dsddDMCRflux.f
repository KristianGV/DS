*******************************************************************************
*** Function dsddDMCRflux provides the local differential flux of DM        ***
*** particles that results from cosmic rays impinging on DM in the          ***
*** diffusive halo, see Bringmann & Pospelov (2018).                        ***
***                                                                         ***
***  Input:                                                                 ***
***    Tkin    - kinetic energy of DM particle [GeV]                        ***
***    rhoDeff - local DM density [GeV/cm**3] multiplied by the effective   ***
***              distance [kpc] out to which the DM flux is assumed to      ***
***              originate from:                                            ***
***                Deff = (\int dV \rho_DM \rho_CR /(4pi d^2)               ***
***                       / (\rho_DM \rho_CR)_local                         ***
***                                                                         ***
***  Output: d\Phi/dT, units: 1/ [cm**2 s GeV]                              ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-22                                                         ***
***  mod 2019-03-10 (diff. scattering x-section now in dsddDMCRsigCR)       ***
*******************************************************************************
      real*8 function dsddDMCRflux(Tkin, rhoDeff)
      implicit none
      include 'dsmpconst.h'

      real*8 Tkin, rhoDeff
      real*8 norm, Tpmin, Tpmax, intres, res
      real*8 mdm, Tdm
      integer CRtype
      common/DMCRflux/ Tdm, CRtype

c... functions
      real*8 dsddDMCRflux_aux, dsmwimp
      external dsddDMCRflux_aux
      
      dsddDMCRflux=0.0d0
      
      mdm  = dsmwimp()
      Tdm  = Tkin
      norm = kpc*1.0d21*rhoDeff/mdm

c...  proton scattering
      CRtype = 1
      Tpmin = sqrt((m_p-Tkin/2.)**2 + (m_p+mdm)**2*Tkin/2./mdm) - m_p + Tkin/2.
      Tpmin = log(Tpmin) ! change to log to better sample fast falling integrand
      Tpmax = Tpmin+4.0d0
      if (Tpmax.lt.2.5) Tpmax = 2.5 ! make sure to sample peak of spectrum
      call dgadap(Tpmin,Tpmax,dsddDMCRflux_aux,1.0d-5,intres)   
      res = intres

c... helium scattering      
      CRtype = 2
      Tpmin = sqrt((m_he-Tkin/2.)**2 + (m_he+mdm)**2*Tkin/2./mdm) - m_he + Tkin/2.
      Tpmin = log(Tpmin) 
      Tpmax = Tpmin+4.0d0
      if (Tpmax.lt.3.0) Tpmax = 3.0 ! make sure to sample peak of spectrum
      call dgadap(Tpmin,Tpmax,dsddDMCRflux_aux,1.0d-5,intres)
      res = res + intres

      dsddDMCRflux = res*norm*1.0d-30 ! undo spurious normalization 
                                      ! from integration
      return
      end


*******************************************************************************
*** auxiliary routine just for integration
*******************************************************************************
      real*8 function dsddDMCRflux_aux(Tcr)
      implicit none
      real*8 Tcr, Tdm
      real*8 dscrISRflux, dsddDMCRsigCR 
      integer CRtype
      common/DMCRflux/ Tdm, CRtype

c... add large normalization to improve integral conversion   
      dsddDMCRflux_aux = 1.0d30*dscrISRflux(exp(Tcr), CRtype)
     &                   *exp(Tcr) ! from conversion to log(Tcr) 
     &                   *dsddDMCRsigCR(Tdm,exp(Tcr),CRtype)
     
c      write(*,*) Tcr, dsddDMCRflux_aux

      return
      end
