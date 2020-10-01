      subroutine dsddinit
c...  initialize the dd routines
      implicit none
      include 'dsddcom.h'
      include 'dsmpconst.h'
      integer i
      do i=1,ddnsf
        ddsf(i)=''
      enddo
      do i=1,ddnme
        ddme(i)=''
      enddo
      call dsddset('sf','default')
      call dsddset('me','default')
      call dsddDMCRquenching_set('default')
c... astro parameters, so far only used in DDCR routines
      rholocal = 0.3d0 ! GeV/cm^3 DEBUG-> extract this from halo model instead 
      vlocal = 2./sqrt(pi)*220*1.0d5 ! mean local DM velocity in cm/s
      Deff = 0.997d0 ! kpc; eff distance out to which source term is integrated
                     ! Deff = 0.997 (8.02) kpc corresponds to 1 (10) kpc in real
                     ! distance
      attenuation_how = 1 ! 1: [default] use analytic expressions from 1810.10543,
                          !    assuming a constant scattering cross section 
                          ! 2: solve differential equation for soil attenuation
                          !    numerically (important for q-dependent scattering) 
      return
      end
