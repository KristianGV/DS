*******************************************************************************
*** Function dsddDMCRdgammadt returns the differential rate of recoils      ***
*** for DM hitting a nucleus at rest. Here, the DM flux is assumed to be    ***
*** the one that results from CR interactions with DM particles in the      ***
*** halo, see Bringmann & Pospelov (2018).                                  ***
***                                                                         ***
***  Input:                                                                 ***
***    Trnuc    - recoil energy of nucleus [GeV]                            ***
***    mN       - mass of target nucleus [GeV]                              ***
***    dsigNdTr - differential scattering cross section on target nucleus   ***
***               per nuclear recoil energy [cm^2/GeV]. Must be supplied as ***
***               EXTERNAL function of momentum transfer Q**2               ***
***    depth    - penetration depth (detector location) [cm]                ***
***    rhoDeff  - local DM density [GeV/cm**3] multiplied by the effective  ***
***               distance [kpc] out to which the DM flux is assumed to     ***
***               originate from:                                           ***
***               Deff = (\int dV \rho_DM \rho_CR /(4pi d^2)                ***
***                       / (\rho_DM \rho_CR)_local                         ***
***                                                                         ***
***  Output: d\Gamma/dTrnuc, units: [1/(s GeV)/nucleus]                     ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-24                                                         ***
*** mod  2019-03-13 nuclear recoil rate as external function of Q**2        ***
*** mod  2019-10-04 changed argument zlfree -> depth
*******************************************************************************
      real*8 function dsddDMCRdgammadt(Trnuc,mN,dsigNdTr,depth,rhoDeff)
      implicit none
      include 'dsmpconst.h'

      real*8 Trnuc, mN, dsigNdTr, depth, rhoDeff
      external dsigNdTr
      
      real*8 mdm, Tdmmin, Tdmmax, Tzmin, intres
      real*8 Trnuccom, depthcom, rhoDeffcom
      common /dgdtcom/ Trnuccom, depthcom, rhoDeffcom
      save /dgdtcom/

c... functions
      real*8 dsddDMCRdgammadtrec_aux, dsmwimp
      external dsddDMCRdgammadtrec_aux
      
      dsddDMCRdgammadt = 0.0d0
      
      mdm=dsmwimp()
      Trnuccom=Trnuc
      depthcom=depth
      rhoDeffcom=rhoDeff

      if (mdm.gt.Trnuc/2.) then
        Tzmin= (mdm-Trnuc/2.)* 
     &       (-1. + sqrt(1. + 2.*Trnuc/mN *(mdm+mN)**2/(2*mdm-Trnuc)**2) )  
      else
        Tzmin= (mdm-Trnuc/2.)* 
     &       (-1. - sqrt(1. + 2.*Trnuc/mN *(mdm+mN)**2/(2*mdm-Trnuc)**2) )  
      endif
        
      call dsddTDMattenuation(Tdmmin, Tzmin, depth, 2) ! get Tmin from Tzmin      
      if (Tdmmin.gt.1.0d10) return ! no way we reach the detector

c TB debug : bypass attenuation factor
c      Tdmmin = Tzmin
c      write(*,*) ' from Tzmin -> Tmin : Tdmmin, Tzmin, depth = ',Tdmmin, Tzmin, depth
      Tdmmax = 6.0d1*Tdmmin ! if chosen too high, integral is unreliable
      Tdmmin = log(Tdmmin)
      Tdmmax = log(Tdmmax)
      call dgadap(Tdmmin,Tdmmax,dsddDMCRdgammadtrec_aux,3.0d-4,intres)

c      write(*,*)
c      write(*,*) 'mdm, MN, Trnuc = ',mdm, MN, Trnuc
c      write(*,*) 'dsddDMCRdgammadt : Tdmmin,Tdmmax,intres = ', Tdmmin,Tdmmax,intres

      dsddDMCRdgammadt=1.0d-30*intres   ! undo rescaling in integrand

      return
      end


*******************************************************************************
*** auxiliary routine just for integration
*******************************************************************************
      real*8 function dsddDMCRdgammadtrec_aux(lnTdm)
      implicit none
      real*8 Tdm, lnTdm,Tz, mdm
      real*8 dsmwimp, dsddDMCRflux, dsddDMCRsigtarget
      real*8 Trnuccom, depthcom, rhoDeffcom
      common /dgdtcom/ Trnuccom, depthcom, rhoDeffcom

      dsddDMCRdgammadtrec_aux = 0.0d0
      mdm=dsmwimp()
      tdm = exp(lnTdm)
      call dsddTDMattenuation(Tdm, Tz, depthcom, 1) ! get Tz from Tdm
c TB debug : bypass attenuation factor
c      Tz = Tdm

c      write(*,*) ' from Tmin -> Tzmin : Tdm, Tz, depthcom = ',Tdm, Tz, depthcom

c... The DM flux needs to be evaluated at the 'top of the atmosphere',
c... but the scattering takes place with the DM energy available at the 
c... detector location!
      dsddDMCRdgammadtrec_aux = tdm*dsddDMCRflux(Tdm, rhoDeffcom)
     &                          * dsddDMCRsigtarget(Trnuccom,Tz) 
     &                          *1.0d30 ! rescale for better numerical result

      return
      end
