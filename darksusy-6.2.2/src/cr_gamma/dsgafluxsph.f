      real*8 function dsgafluxsph(egev,diff,xi,labhalo,psi,theta,istat)
**********************************************************************
***   function dsgafluxsph gives the flux of gamma-rays from DM decays
***   and/or annihilations in the halo (in the limit of zero relative
***   velocity) for a spherical dark matter profile      
***   
***   inputs:
***     egev  - gamma-ray energy [in GeV]      
***     diff  - dictates whether differential source term at egev
***             (diff=1) or integrated source term above egev (diff=0)
***             is returned      
***     xi    - factor by which the DM density should be rescaled 
***             (obtained e.g. as ratio of calculated to measured relic
***              density)
***     labhalo = halo model access label, to load the halo model via 
***         dsdmdselect_halomodel
***     psi = the angular offset between the pointing direction
***         and the direction of the center of the distributionin
***         (in rad)
***     theta = aperture of the acceptance cone (in rad; it defines
***         the solid angle within which line-of-sight-integrals are 
***         performed assuming a step-function response from the 
***         instrument; more options are available in dsomlosisph, use
***         that to compute line-of-sight-integrals in your main file
***         and compute fluxes using dscrgaflux_v0ann or
***         dscrgaflux_dec)
***
***   type : commonly used
***   desc : gamma-rays from decay/annihilation from halo specified by halo label
***
***      
***   unit of return value: cm^-2 s^-1/ GeV^diff
***   author: Piero Ullio (ullio@sissa.it)
***   date: 2017
**********************************************************************
      implicit none
      real*8 egev,xi,psi,theta
      integer diff,istat
      character(*) labhalo
      real*8 par,res,dscrgaflux_dec,dsdfactor,dscrgaflux_v0ann,
     &  dsjfactor
ccc
      par=0.d0
ccc
ccc include (eventually) a decay term decay and annihilation:        
ccc
      res=dscrgaflux_dec(egev,diff,1.d0,xi,istat)
      if(res.gt.0.d0) par=par+res*dsdfactor(labhalo,psi,theta)
ccc
ccc include (eventually) an annihilation term:        
ccc
      res=dscrgaflux_v0ann(egev,diff,1.d0,xi,istat)
      if(res.gt.0.d0) par=par+res*dsjfactor(labhalo,psi,theta)
      dsgafluxsph=par
      return
      end


