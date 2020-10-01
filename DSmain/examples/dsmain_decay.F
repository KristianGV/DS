      program dsmain_decay
c
c     This program is a sample DarkSUSY main program for DECAYING DM. 
c     It asks for model parameters interactively and calculates various 
c     rates. (see dsmain_wimp.F for a similar program, adapted to
c     *annihilating* DM)
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none

      character*12 refhalolab                           ! halo label
      real*8 mwrho0,dsdmdrho0mw                         ! halo routines
      real*8 psi0,theta0,dsdfactor, dfactor             ! los routines
      real*8 fluxgacdiff,fluxgac, dsgafluxsph           ! gamma-ray flux
      real*8 cospsi0,delta                              ! gamma-ray flux
      integer istat, nlines                             ! gamma-ray flux, etc
      real*8 egam,egath                                 ! gamma-ray flux
      real*8 dscrgaflux_line_dec                        ! gamma-ray flux
      real*8 eline, fluxgaline,widthline                ! monochromatic signals
      real*8 tpbess(3),phiin,pb_a,pb_b,pb_c,dspbdphidtaxi ! pbar flux
      real*8 tp,phidb, dsdbdphidtaxi                    ! db flux
      real*8 phiep, dsepdphidpaxi                       ! positron flux
      real*8 eth, phimuhalo, dscrmuflux_dec             ! neutrino telescopes
      integer i, unphys,hwarning,acceptable
      real*8 dsmWIMP
c...model parameters for generic WIMP
      real*8 genmwimp, gendecrate
      integer genpdg
      

c
c     Here we include various "global variables" defined in common blocks
c     in src/, i.e. independent of the chosen particle physics module
c
      include 'dsio.h'      ! for I/O handling
      include 'dsidtag.h'   ! unique model ID tag for each set of model parameters
      include 'dshmcom.h'   ! for rho0 of the halo model in dfactor rescaling


c read in module definitions and compare to the one determined at compile time
c the latter file is a temporary file created during the make process
#include <module_defs.F>
#include "module_compile.F"

#if MODULE_CONFIG != MODULE_generic_decayingDM
      write(*,*) 
      write(*,*)
     &  '----------------------------------------------------------'
      write(*,*)
     &  'You have called dsmain_decay with a particle module that'
      write(*,*)
     &  'does not have a decaying DM component, so no observables'
      write(*,*)
     &  'were calculated.'     
      write(*,*)
     &  '----------------------------------------------------------'
      write(*,*)
     &  '[call ''make -B dsmain_decay DS_MODULE=<MY_MODULE>'' ', 
     &            'to try with another module <MY_MODULE>]'
      write(*,*) 
      stop 
      goto 999
#endif


c
c     A call to dsinit should be the first call in any program using DarkSUSY.  
c     It initializes some global variables and calls various other modules to 
c     set them up properly.
      call dsinit
 
c     The amount of output onto standard output can be controlled by the
c     variable prtlevel. Higher values imply more output, see 
c     for more details README.PRINTLEVEL. 
c
      prtlevel=1


**************************************
***       halo routines setup      ***
**************************************
c
c     A set of default halo models have been predefined and initialized  
c     with the call to the subroutine dsdmdinit within dsinit. Print
c     the list of available models by commenting in the following:
c
c      call dsdmdprint_halomodel('printall')
c
c     The models in this database are used via the corresponding access
c     label, i.e. 'mwnfwdef', 'mwburdef' & 'mweindef' in this case, which
c     has to be provided when calling indirect rate routines (see below). For     
c     the purpose of dstest, we will only use one refernce halo and hence label
c      
      refhalolab='mwnfwdef'
c
c     For a Milky Way profile (label name contains "MW"), the local halo 
c     density is obtained as follows (explicitly needed by DD routines)
c      
      mwrho0=dsdmdrho0mw(refhalolab)
      rho0=mwrho0     ! local halo density passed to DD via dshmcom.h
c     
c     For the fluxes of gamma-rays and neutrinos with the chosen halo
c     profile, we calculate, once and for all, the line of sight
c     integration factor D (for decaying DM!) in the direction of observation, 
c     which we define as the direction which forms an angle psi0 with respect
c     to the direction of the galactic centre (e.g. cospsi0 = 1 is the
c     galactic center direction, cospsi0 = -1 the antigalactic centre).
c
      write(*,*) 'Calculating D-factor...'

      cospsi0=1.d0
      delta=1.d-3
      psi0=dacos(cospsi0)
      theta0=dacos(1.d0-delta/(8.d0*datan(1.d0)))
      dfactor=dsdfactor(refhalolab,psi0,theta0) ! kpc sr GeV cm^-3
      write(*,*) 'D_factor = ',dfactor,
     &  ' kpc sr GeV cm^-3'



*****************************
***      model setup      ***
*****************************
c
c    First read in model parameters
c
      write(*,*) ' '
      write(*,*)
     &  'Please specify model parameters for a generic decaying ',
     &           'DM particle:'
      write(*,*)
     &  '--------------------------------------------------- '
      write(*,*) ' '
      write(*,*) 'Enter DM mass (GeV): '
      read(5,*) genmwimp
      write(*,*) 'Enter decay rate (1/s): '
      read(5,*) gendecrate
      write(*,*)
     &  'Enter PDG code for dominant annihilation channel'
      write(*,*) '(5 for bbar, 24 for W^+W^-, etc...): '
      read(5,*) genpdg
      write(*,*)
     &  '--------------------------------------------------- '
   
   
      
c... Now transfer these parameters to the DarkSUSY common blocks    
c... NB: This simple input dialogue assumes only one dominant decay channel, 
c... but you could specify any number, including invisible channels, see header
c... of dsgivemodel_decayingDM
      call dsgivemodel_decayingDM(genmwimp,gendecrate,1,1.0d0,
     &  genpdg,-genpdg)

c...For larger runs, it is also convenient to define a unique 12-character
c...idtag for each model. This one will be printed by DarkSUSY if there is
c...any problem with a specific model. For now, just set it to something.

           idtag='INTMOD000001'

c     Now we are ready to calculate the particle spectrum
c     and the 3-particle vertices. We use the routine dsmodelsetup. 
c     The flag unphys<0 means that the model is theoretically inconsistent, 
c     while the 2nd flag is more model-specifc (here, it is set to 1 if 
c     kinematically inaccessible channels are included; those are treated
c     as invisible).
c
      call dsmodelsetup(unphys,hwarning)
      acceptable=0
      if (unphys.ne.0.or.hwarning.ne.0) then
         write(*,*) unphys,hwarning
         acceptable=-1
         goto 150
      endif

      write(*,*) ' '
      write(*,*) '***** MODEL: ',idtag,' *****'
 150  write(*,*) 'acceptable =',acceptable, ' (0=OK, -1=not OK)'

c     We now set the print level to 0 to get severe error messages only.
c
        prtlevel=0
c
c     We now have all the masses and couplings calculated. Let's, e.g.,  
c     print out the DM mass (other masses, or couplings would typically be
c     particle module-dependent):
c
      write(*,*) 'Decaying dark matter mass = ', dsmWIMP()
      
      
      
***********************************************
***           Gamma-ray fluxes              ***
***********************************************
c
c     Next we compute cosmic-ray rates from neutralino annihilations 
c     in the galactic halo. We start with the gamma-ray flux, which
c     factorizes into part that depends on the particle physics 
c     and the line-of-sight integral D computed before setting up the model. 
c     There are two main ways of computing this flux, either by
c     specifying the correct label to our halo model -- or by providing
c     the D factor 'by hand' (e.g. as taken from the literature)
c
c     1) gamma-ray flux with continuum energy spectrum at a given energy
c     egam (GeV). 
c
      write(*,*)
      write (*,*) 'Calculating gamma ray fluxes...'
      egam=20.d0
      fluxgacdiff=dsgafluxsph(egam,1,1.d0,refhalolab,psi0,theta0,istat) !ph cm^-2 s^-1 GeV^-1
c     providing the previosuly calculated D factor gives the same result,
c     but require a different function call:
c      
c      fluxgacdiff=dscrgaflux_v0ann(egam,1,dfactor,1.d0,istat) !ph cm^-2 s^-1 GeV^-1

      write(*,*) '  fluxgacdiff = ',fluxgacdiff,' ph/(cm^2 s GeV)'

c     2) gamma-ray flux with continuum energy spectrum integrated above
c        some given threshold egath (GeV). E.g. :
c
      egath=1.d0
      fluxgac=dsgafluxsph(egath,0,1.d0,refhalolab,psi0,theta0,istat) !ph cm^-2 s^-1
c     again, the result is identical to:
c      
c      fluxgac=dscrgaflux_v0ann(egath,0,dfactor,1.d0,istat)  !ph cm^-2 s^-1

      write(*,*) '      fluxgac = ',fluxgac,' ph/(cm^2 s)'
      write(*,*)

c
c     3) monochromatic gamma-ray flux induced by 1-loop annihilation
c        processes into a 2-body final state containing a photon. First, we
c        need to determine the (model-dependent!) number of such lines:
      nlines=0 ! calling the routine below with nlines=0 overwrites this variable 
               ! with the total number of lines on return -- and gives us 
               ! flux, energy and width of the *first* line in the list
      fluxgaline=
     &  dscrgaflux_line_dec(nlines,eline,widthline,dfactor,1.d0,istat) 
      write(*,*) 'Total number of photon lines in module ', 
     &            moduletag(1:index(moduletag,' ')-1),':',nlines
     
c   now we can compute flux, energy and width also for the other lines: 
      do i=1,nlines
        write(*,*) 'photon flux from line No.',i,' = ',fluxgaline,' ph/(cm^2 s)'
        write(*,*) '        [at E = ',eline,' +/- ',widthline,' GeV]'
        if (i.lt.nlines) fluxgaline=
     &  dscrgaflux_line_dec(i+1,eline,widthline,dfactor,1.d0,istat) 
      end do
     
     
 
 
*********************************************************
***   Charged CR rates, interstellar, solarmodulated  ***
*********************************************************
c
c     Now we come to the antiproton routines. We get the differential
c     spectrum in units of GeV^-1 cm^-2 s^-1 sr^-1 by calling
c     dspbdphidtaxi with the antiproton kinetic energy as the first
c     argument. The second argument determines if the flux is solar
c     modulated (a la Perko). The third argument determines how the diffusion
c     equations are solved. If the argument is
c          1 = they are solved for the requested energy only
c          2 = they are tabulated for a large range of energies for the first
c              call and the table is used for subsequent calls.
c          3 = as 2, but the table is also written to disk in the current
c              directory on the first call.
c          4 = the table is read from disk on first call, and is used for
c              that and subsequent calls.
c     The tabulation can take several minutes (or hours), 
c     so if you are only interested in a few models, don't tabulate,
c     but if you want to calculate many models, do tabulate 
c     (i.e. using any of the options 2-4). Some standard tables are 
c     included in the distribution and are available in files
c     of the form data/pbtd-*.dat. You can change the default propagation
c     model with a call to dspbset (see that routine for details).

      write(*,*)
      write (*,*) 'Calculating antiproton fluxes...' ! for the following energies
      tpbess(1)=0.35d0
      tpbess(2)=1.76d0
      tpbess(3)=3.00d0

c     Here are two example calls, the first one uses the routines directly,
c     whereas the second uses tables. 
c     NB: if you use tables and change properties of the diffusion model 
c         that require a retabulation it is your responsibility to retabulate 
c         (easiest done by deleting the appropriate tabulated file 
c         data/pbtd-*.dat) !!!

      phiin=0.32d0! solar modulation parameter in GV, assuming that solar 
                  ! modulation can be treated with the force-field method
                  ! set it to below 1.d-3 and you get the interstellar

c compute the fluxes for the default NFW:
c       pb_a=dspbdphidtaxi(tpbess(1),phiin,1,refhalolab) ! no tables, slower
      pb_a=dspbdphidtaxi(tpbess(1),phiin,4,refhalolab) ! use tables, faster after slower startup
      pb_b=dspbdphidtaxi(tpbess(2),phiin,4,refhalolab) 
      pb_c=dspbdphidtaxi(tpbess(3),phiin,4,refhalolab) 
      write(*,*) 'solar modulated pbar flux at 0.35 ',
     &              'GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ', pb_a
      write(*,*) 'solar modulated pbar flux at 1.76 ',
     &              'GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ', pb_b
      write(*,*) 'solar modulated pbar flux at 3.00 ',
     &              'GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ', pb_c
      
c
c      We can also calculate the flux of antideuterons. The call is very
c      similar to that for pbar. The dbar's use the same diffusion model
c      as the pbar and can thus be changed with a call to dspbset.
c
      write (*,*) 'Calculating antideuteron fluxes...'
      tp=1.0d0

c      phidb=dsdbdphidtaxi(tp,phiin,1,refhalolab) ! no tables, slower
       phidb=dsdbdphidtaxi(tp,phiin,4,refhalolab) ! use tables, faster after slower startup

      write(*,*)
     &  '  solar modulated dbar flux at 1.00 GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ',
     &  phidb

c
c     We now continue with the rates of positrons from neutralino
c     halo annihilation. The first argument is the positron energy (GeV)
c     and the second the solar modulation parameter (GV), the third
c     argument determines if the Green's function needed for the
c     computation should be read from disk or created on the spot.
c     Default is to read from disk, and if the file does not exist,
c     it will recreate it. Some standard files are available in 
c     data/eptab-*.dat and data/epgretab-*.dat.
c
      write (*,*) 'Calculating positron fluxes at 1 GeV...'

      phiep=dsepdphidpaxi(1.0d0,0.d0,4,refhalolab) ! use of tables necessary
      write (*,*) '  phiep=',phiep,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'



***********************************
***   Neutrino telescope rates  ***
***********************************
c
c     For the decaying DM module, there is no DM scattering rate with
c     nuclei implemented, and hence DM is not captured in the sun or earth.
c     But we can calculate the neutrino-induced muon rates
c     in neutrino telescopes that would occur from neutralino
c     annihilations in our galactic halo. The neutrino-induced muon
c     flux above a threshold Eth is given by
c     

      write(*,*) 
     &  'Calculating neutrino-induced muon fluxes from the halo...'
      eth=1.0d0
      phimuhalo=dscrmuflux_dec(eth,0,dfactor,1.d0,istat) ! km^-2 yr^-1

      write(*,*) '  Muon flux from halo = ',phimuhalo,
     &  ' km^-2 yr^-1'


      write (*,*)
      write (*,*) '-------------------------------------------------------'
      write (*,*) 'The DarkSUSY example program has finished successfully.'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-------------------------------------------------------'
      write (*,*)
      write (*,*) '[call ''make -B dsmain_wimp DS_MODULE=<MY_MODULE>'' ', 
     &            'if you want to try with a module <MY_MODULE> for annihilating',
     &            ' DM instead.]'
      write (*,*)
      goto 999
      stop
 999  end


