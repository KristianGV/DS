      program dsmain_wimp
c
c     This program is a sample DarkSUSY main program. Depending on the chosen
c     particle module, it asks for model parameters interactively and 
c     calculates the relic density and various rates for each given model. 
c     This program can be used as is for quick calculations, or be copied and
c     then(!) modified to your liking for more advanced uses.
c
c-----This line is 72 columns long--------------------------------------
c
      implicit none

      real*8 oh2,xf,dsrdomega                            ! relic density
      real*8 tkd,dskdtkd,dskdmcut,mcut, dsmWIMP, mdm     ! kinetic decoupling
      real*8 sigsip,sigsin,sigsdp,sigsdn                 ! direct detection
      real*8 sigij(27,27)                                ! direct detection
      integer i,nnuc,a(10),z(10),stoich(10)              ! direct detection
      real*8 v,e,si(10),sd(10),siff(10),sdff(10),t,rsi,rsd ! direct detection
      real*8 fluxgacdiff,fluxgac                         ! gamma-ray flux
      real*8 jpsi,cospsi0,delta                          ! gamma-ray flux
      real*8 nsigvgacont,nsigvgacdiff                    ! gamma-rays
      integer istat, nlines                              ! gamma-ray flux, etc
      real*8 egam,egath, dsgafluxsph                     ! gamma-ray flux
      real*8 dscrgaflux_line_v0ann                       ! gamma-ray flux
      real*8 fluxgaline, egev, nsigvline,widthline, pdg2 ! monochromatic signals
      real*8 eline, dscrsource, dscrsource_line          ! DM source function for CR
      real*8 tpbess(3),pb_a,pb_b,pb_c                    ! pbar flux
      real*8 phiin,dspbdphidtaxi                         ! pbar flux
      real*8 tp,phidb, dsdbdphidtaxi                     ! db flux
      real*8 phiep, dsepdphidpaxi                        ! positron flux
      real*8 eth,thmax,rateea,ratesu                     ! neutrino telescopes
      integer rtype,ptype, kind                          ! neutrino telescopes
      real*8 phimuhalo, dscrmuflux_v0ann,rho0in          ! neutrino telescopes
      integer ierr,iwar,nfc, istatus,info
      character*12 refhalolab                            ! halo label
      real*8 mwrho0,dsdmdrho0mw                          ! halo routines
      real*8 psi0,theta0,dsjfactor                       ! los routines


c
c     Here we include various "global variables" defined in common blocks
c     in src/, i.e. independent of the chosen particle physics module
c
      include 'dsio.h'      ! for I/O handling
      include 'dsidtag.h'   ! unique model ID tag for each set of model parameters
      include 'dshmcom.h'   ! for rho0 of the halo model in jpsi rescaling


c read in module definitions and compare to the one determined at compile time
c the latter file is a temporary file created during the make process
#include <module_defs.F>
#include "module_compile.F"

#if MODULE_CONFIG == MODULE_generic_decayingDM
      write (*,*)
      write (*,*)
     &  '----------------------------------------------------------'
      write (*,*)
     &  'This program does not support decaying dark matter models!'
      write (*,*)
     &  '----------------------------------------------------------'
      write (*,*)
      write (*,*)
     &  '[call ''make -B dsmain_wimp DS_MODULE=<MY_MODULE>'' ', 
     &            'to try with another module <MY_MODULE>]'
      write (*,*)
      stop
      goto 999
#endif


c
c     This call should be the first call in any program using DarkSUSY. It 
c     initializes some global variables and calls various other modules to 
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
c     integration factor J (not D for WIMPs!) in the direction of observation, 
c     which we define as the direction which forms an angle psi0 with respect
c     to the direction of the galactic centre (e.g. cospsi0 = 1 is the
c     galactic center direction, cospsi0 = -1 the antigalactic centre).
c
      write(*,*) 'Calculating J-factor...'
      cospsi0=1.d0
      delta=1.d-3
      psi0=dacos(cospsi0)
      theta0=dacos(1.d0-delta/(8.d0*datan(1.d0)))
      
      jpsi=dsjfactor(refhalolab,psi0,theta0) ! kpc sr GeV^2 cm^-6
      write(*,*) 'J_factor galactic center =',jpsi,
     &  ' kpc sr GeV^2 cm^-6'


*****************************
***      model setup      ***
*****************************
c
c    Now we input the particle model parameters. The form of this subroutine 
c    depends on the chosen particle module (see bottom of this file), but
c    typically involves a call to a routine of the form dsgivemodel residing
c    in src_models/
      call dspmenterparameters

c    For larger runs, it is also convenient to define a unique 12-character
c    idtag for each model. This one will be printed by DarkSUSY if there is
c    any problem with a specific model. For now, just set it to something.
c
           idtag='INTMOD000001'

c
c     After providing model parameters with a routine like dsgivemodel_xxx,
c     the next step is always a call to the routine dsmodelsetup. This
c     actually sets up the model, by calculating the mass spectrum, 
c     relevant 3-particle vertices etc. dsmodelsetup return two flags
c     to indicate if everything is OK (both flags equal 0) or not. The 
c     flag unphys<0 means that the model is theoretically inconsistent, 
c     while the 2nd flag is more model-specific. Depending on the chosen
c     particle module, one may want to do further tests, all of which we
c     have included in a routine dspmcompute (see bottom of this file).
c     It resturns an exclusion flag "istatus" and an integer "info" 
c     explaining the reason

      call dspmcompute(istatus,info)

      write(*,*) ' '
      write(*,*) '***** MODEL: ',idtag,' *****'
      write(*,*) 'istatus =',istatus,
     &     ' (=0 OK, <0 error, >0 warning, =1 experimentally excluded)'
      if (istatus.ne.0) then
         write (*,*) 'info = ',info
      endif
      if (istatus.lt.0) goto 150


c     We now set the print level to 0 to get severe error messages only.
        prtlevel=0

c
c     We now have all the masses and couplings calculated. Let's, e.g.,  
c     print out the DM mass (other masses, or couplings would typically be
c     particle module-dependent):
c
      write(*,*) 'WIMP mass = ', dsmWIMP()


***********************************************
***   relic density + kinetic decoupling    ***
***********************************************
c
c     The relic density is calculated by the function dsrdomega.
c     The first argument determines if coannihilations should be included.
c     If the argument is 0, no coannihilations are included, if it is 1,
c     all relevant coannihilations are included (what this means is model-
c     specific: this flag is passed to the interface function dsrdparticles, 
c     and you need to look the header of that file for the respective particle 
c     module). 
c     The second argument determines the accuracy. Set it to 1 for a faster 
c     calculation (e.g. only including coannihilating particles with mass 
c     difference below 30%), or to 0 for a more accurate calculation.
c     Even faster options (fast=2, 3 for progressively faster (and less
c     accurate calculations, or fast=9 for a superfast order of magnitude
c     calculation). In practice, the fast=1 option is more than adequate
c     (better than 1% accuracy).
c     There are also other options for fast, fast=99 is a very slow option
c     that uses no tabulation and should not normally be used. fast=20
c     is a new option that uses a system with tabulation on the fly and
c     a Breit-Wigner fit to speed things up. It can be used as an alternative
c     to fast=1 and might possibly be the default in future releases.
c     The function returns the relic density, omega h^2,
c     the freeze-out temperature, xf (x=m/T), one error flag, ierr and
c     one warning flag, iwar, both of which are 0 if everything is OK.
c     nfc is the number points in momentum where the cross section had
c     to be calculated. Note that the omega calculation can be very
c     time consuming, needing up to several minutes for tricky models
c     (with many resonances, thresholds, coannihilations etc).
c     
c     The next line selects the effective degrees of freedom in the early
c     Universe. The default is Drees et al [1503.03513]. You can choose something
c     else by commenting in the line below (check the file for other options)
c
c      call dsrdset('dof','1') ! choose GG (150 MeV) dof instead of Drees default

      write(*,*)
      write(*,*) 'Calculating relic density without coannihilations,',
     &     ' please be patient...'
      oh2=dsrdomega(0,1,xf,ierr,iwar,nfc)
      write(*,*) '  Oh2 = ',oh2,ierr,iwar
      
c    the following is an example of something we can only do for modules
c    that actually include coannihilations. The precompiler directive therefore
c    makes sure that the following lines are only computed for the mssm
c
#if MODULE_CONFIG == MODULE_mssm
      write(*,*) 'Calculating omega h^2 with coannihilations,',
     &     ' please be patient...'
      oh2=dsrdomega(1,1,xf,ierr,iwar,nfc)
      write(*,*) '  with coannihilations Oh2 = ',oh2,ierr,iwar
#endif

      if (xf.ne.0.d0) then
        write(*,*) '  Chemical decoupling (freeze-out) occured at'
        write(*,*) '  T_f = ',dsmWIMP()/xf,' GeV.'
      endif


c     Now let's calculate the kinetic decoupling, and the smallest halos
c     we can have for this model.
      tkd=dskdtkd(1) ! 1=full calculation
      if (tkd.gt.1.d10) then
        if (prtlevel.gt.1) then
          write(*,*) ' '
          write(*,*) ' Note: The kinetic decoupling temperature is formally so'
          write(*,*) ' high in this model that local thermal equilibrium would'
          write(*,*) ' not have been established during the freeze-out process.'
          write(*,*) ' Likely, this is just the result of a simplified implementation'
          write(*,*) ' of the scattering process (as in the case of the '
          write(*,*) ' generic_wimp model), not taking into account scattering'
          write(*,*) ' with heavy quarks.' 
          write(*,*) ' '      
        endif  
      else
        write(*,*) ' '
        write(*,*) 'Kinetic decoupling temperature, Tkd = ',tkd, ' MeV'
        mdm=dsmwimp()
        mcut=dskdmcut(mdm,tkd,1) ! 1 = true cut-off
        write(*,*) ' The resulting cutoff in the power spectrum'//
     &     ' corresponds to a mass of',' M_cut/M_sun = ', mcut
        write(*,*) ' '
      endif


***********************************************
***           direct detection              ***
***********************************************
c
c     We are now ready to calculate rates.  Let's start with scattering
c     cross sections for direct detection experiments by calling
c     dsddsigma which returns the spin-independent and the
c     spin-dependent scattering cross sections off a target. The cross
c     sections are returned in units of cm^2.
c
c     The next lines sets the form factors to be the best available, 
c     which is also the default setting. For the spin independent form 
c     factors, this implies in order Fourier-Bessel, Sum-of-gaussians, Fermi, 
c     Lewin-Smith. For the spin dependent form factors, it means in order 
c     interacting shell model, odd group model, single particle shell model.
c     Use call dsddhelp to print out the other possibilities.
c
      call dsddset('sf_m','best')
      call dsddset('sf_sigma','best')

c     For the MSSM, there are also model-specific options implemented. Again,
c     we use a precompiler directive to make sure the following is only
c     executed in this case.
c
#if MODULE_CONFIG == MODULE_mssm
      call dsddset_mssm('dn_nopole')
#endif

      write (*,*) 'Calculating DM-nucleon scattering cross sections...'
      v=0.d0
      e=0.d0
      call dsddsigma(v,e,1,0,sigij,ierr)
      sigsin=sigij(1,1) ! SI cross section
      sigsdn=sigij(4,4) ! SD cross section
      call dsddsigma(v,e,1,1,sigij,ierr)
      sigsip=sigij(1,1) ! SI cross section
      sigsdp=sigij(4,4) ! SD cross section
c...Note: if we are only interested in the spin-independent and spin-
c...dependent scattering cross sections, we could call the wrapper
c...routine dsddsigmanucleon(v,e,sigsip,sigsin,sigsdp,sigsdn,ierr)
c...instead of calling dsddsigma above. dsddsigma is the more general
c...interface function the particle physics module should provide though.      
      write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
      write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
      write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
      write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36
      write(*,*) ' proton: sigsi (pb) = ',sigsip*1.0d36,
     &     ' sigsd (pb) = ',sigsdp*1.0d36

c
c     Now take a more complex example with scattering off a compound
c
      nnuc = 4
      a(1) = 23                 ! Na-23
      z(1) = 11                 ! Na-23
      a(2) = 127                ! I-127
      z(2) = 53                 ! I-127
      a(3) = 73                 ! Ge-73
      z(3) = 32                 ! Ge-73
      a(4) = 1                  ! p
      z(4) = 1                  ! p
c      a(5) = 1                  ! n
c      z(5) = 0                  ! n
c      a(6) = 7                  ! Li-7
c      z(6) = 3                  ! Li-7
      v = 200.d0                ! WIMP-nucleus velocity in km/s
      e = 10.0                  ! recoil energy in keV
      do i=1,nnuc
         call dsddsigma(0.d0,0.d0,a(i),z(i),sigij,ierr)
         si(i)=sigij(1,1)
         sd(i)=sigij(4,4)
         call dsddsigma(v,e,a(i),z(i),sigij,ierr)
         siff=sigij(1,1)
         sdff=sigij(4,4)
         write(*,*) '  A=',a(i),' Z=',z(i),
     &        ' sigsi (pb) = ',si(i)*1.0d36,
     &        ' sigsd (pb) = ',sd(i)*1.0d36
          write(*,*) 
     &         ' sigsi*ff(pb)=',siff(i)*1.0d36,
     &         ' sigsd*ff(pb)=',sdff(i)*1.0d36
      enddo
      stoich(1)=1
      stoich(2)=1
      t = 651.3d0
c... NB: rate calculations assume that the neutralino makes up 100%
c... of the local DM density. If this is not the case, the rate must
c... be (linearly) rescaled by the user    
      call dsdddrde(t,e,2,a,z,stoich,rsi,rsd,1)
      write(*,*) ' NaI :  dRdE [counts/kg-day-keV] (SI)=',rsi,
     &     ' (SD)=',rsd
 
 
      
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
c     the jfactor 'by hand' (e.g. as taken from the literature)
c
c     1) gamma-ray flux with continuum energy spectrum at a given energy
c     egam (GeV). 
c
      write(*,*)
      write (*,*) 'Calculating gamma ray fluxes...'
      egam=20.d0
      fluxgacdiff=dsgafluxsph(egam,1,1.d0,refhalolab,psi0,theta0,istat) !ph cm^-2 s^-1 GeV^-1
c     providing the previosuly calculated J factor gives the same result,
c     but require a different function call:
c      
c      fluxgacdiff=dscrgaflux_v0ann(egam,1,jpsi,1.d0,istat) !ph cm^-2 s^-1 GeV^-1

      write(*,*) '  fluxgacdiff = ',fluxgacdiff,' ph/(cm^2 s GeV)'

c     2) gamma-ray flux with continuum energy spectrum integrated above
c        some given threshold egath (GeV). E.g. :
c
      egath=1.d0
      fluxgac=dsgafluxsph(egath,0,1.d0,refhalolab,psi0,theta0,istat) !ph cm^-2 s^-1
c     again, the result is identical to:
c      
c      fluxgac=dscrgaflux_v0ann(egath,0,jpsi,1.d0,istat)  !ph cm^-2 s^-1

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
     &  dscrgaflux_line_v0ann(nlines,eline,widthline,jpsi,1.d0,istat) 
      write(*,*) 'Total number of photon lines in module ', 
     &            moduletag(1:index(moduletag,' ')-1),':',nlines
     
c   now we can compute flux, energy and width also for the other lines: 
      do i=1,nlines
        write(*,*) 'photon flux from line No.',i,' = ',fluxgaline,' ph/(cm^2 s)'
        write(*,*) '        [at E = ',eline,' +/- ',widthline,' GeV]'
        if (i.lt.nlines) fluxgaline=
     &  dscrgaflux_line_v0ann(i+1,eline,widthline,jpsi,1.d0,istat) 
      end do

     
c     If we want to do the gymnastics ourselves converting cross
c     sections to fluxes, we can also call dsnsigvgacont,
c     dsnsigvgacdiff and dsnsigvgaline to get the number of photons 
c     (2 for gamma gamma, 1 for other X gamma and whatever the number is for
c     continuous gammas) times the cross section and divided by a symmetry
c     factor (=2 if DM is self-conjugate, =1 if it is not).. 
c     The number returned is the dimensionless number
c        N_gammas * (sigma v) / (10^-29 cm^3 s^-1)
c     with (sigma v) in unit of cm^3 s^-1. For dsnsigvgacdiff the units
c     are of course GeV^-1 instead.

      egev=1.0d0
      nsigvgacont=dscrsource(egev,0,22,2,0d0,istat)*2d0*dsmWIMP()**2   
      nsigvgacdiff=dscrsource(egev,1,22,2,0d0,istat)*2d0*dsmWIMP()**2    
      write(*,*)
      write(*,*) ' nsigvgacont = ',nsigvgacont
      write(*,*) ' nsigvgacdiff = ',nsigvgacdiff,' GeV^-1'
      do i=1,nlines
        nsigvline=dscrsource_line(22,i,2,0d0,egev,widthline,pdg2,istat) 
     &          *dsmWIMP()**2
        write(*,*) '    nsigvline ',i,' = ',nsigvline
      enddo 



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
c     To calculate the rates in neutrino telescopes we call dssenu_rates. We
c     can either calculate the flux of neutrinos, the conversion rate
c     (per volume element) of neutrinos to muons or the muon flux. This
c     is determined by the argument rtype. We can also choose the
c     energy threshold and the maximal half-aperture angle from the
c     center of the Earth/Sun we are interested in. The rates for both
c     the earth and the sun are returned in units of km^-2 yr^-1 for
c     rtype=1,3 and in units of km^-3 yr^-1 for rtype=2. If some
c     warnings were issued, the flag istat is non-zero.
c
c     The default calculation method is to use the full expressions by
c     Gould and numerically integrate them over the velocity distribution
c     as specified by the halo profile, e.g. a gaussian for an
c     isothermal sphere. This numerical integration is rather slow and
c     there is thus an option to use tables (read from disk, recreated
c     if absent) instead. This is the default. Changing to numerical 
c     integration directly is handled by a call to dssenu_set (see dssenu_set.f
c     for details). Also other (approximate) formulae or the Damour
c     Krauss population are available and can be chosen by a call to dssenu_set.
c     Note: For the earth, which captures from a population of WIMPs bound
c     in the solar system, a new estimate (Lundberg and Edsjo, 
c     astro-ph/0401113) is used as a default. Also note that there is a
c     new option to use a cutoff for WIMPs that would reach out to Jupiter
c     after the first scatter. If you want to use this cut-off, uncomment
c     the line below. This will reduce the rates for heavy (>~ 1 TeV)
c     WIMPs.
c     
c      call dssenu_set('tabcut')  ! to take away WIMPs that reach Jupiter

      write(*,*) 'Calculating rates in neutrino telescopes'
      eth=1.0d0      ! energy threshold (of neutrino/muon), GeV
      thmax=30.0d0   ! the maximum half-aperture angle, degrees
      kind=1         ! 1=integrated
                     ! 2=differential
                     ! 3=mixed (diff in E, integrated in theta)
      rtype=3        ! 1=neutrino flux
                     ! 2=neutrino to muon conversion rate
                     ! 3=muon flux
      ptype=3        ! 1=particles only
                     ! 2=anti-particles only
                     ! 3=sum of particle and anti-particle rates

      rho0in=dsdmdrho0mw(refhalolab)    ! GeV cm^-3, local halo density
                                        ! (in the next step, we assume that 
                                        !  *all* of DM is composed of particles
                                        !  as implemented in the active particle 
                                        !  module)
      call dssenu_rates(eth,thmax,1,rtype,ptype,rho0in,rateea,ratesu,
     &  istat)

      write(*,*) '  Flux from the Earth = ',rateea, ' km^-2 yr^-1'
      write(*,*) '  Flux from the Sun =   ',ratesu, ' km^-2 yr^-1'

c     If you want differential rates, you instead call dsntdiffrates. The
c     units of the rates are then km^-2 yr^-1 GeV^-1 degree^-1 for
c     rtype=1 or rtype=3 and km^-3 yr^-1 GeV^-1 degree^-1 for
c     rtype=2. Uncomment the following lines if you want differential
c     rates.

c      energy=10.0d0  ! energy (of neutrino/muon), GeV
c      theta=30.0d0   ! angle from center of Earth/Sun, degrees
c      rtype=3        ! 1=neutrino flux
c                     ! 2=neutrino to muon conversion rate
c                     ! 3=muon flux
c      ptype=3        ! 1=particles only
c                     ! 2=anti-particles only
c                     ! 3=sum of particle and anti-particle rates
c      call dsntdiffrates(energy,theta,rtype,ptype,rateea,ratesu,istat)
c
c... Or like in the following example where we calculate the differential
c...spectrum in energy for an angle of 2 degrees.
c
c      do i=1,100
c         enu=(dble(i)-0.5d0)/100.d0*mass(kn(1))
c         call dsntdiffrates(enu,2.0d0,3,3,rateea,ratesu,istat)
c         write(67,*) enu,rateea,ratesu
c      enddo
         
c
c     Muon rates from the halo
c     
c     It is also possible to calculate the neutrino-induced muon rates
c     in neutrino telescopes that would occur from neutralino
c     annihilations in the galactic halo. The neutrino-induced muon
c     flux above a threshold Eth is given by
c     
      write(*,*) 
     &  'Calculating neutrino-induced muon fluxes from the halo...'
      eth=1.0d0

      phimuhalo=dscrmuflux_v0ann(eth,0,jpsi,1.d0,istat) ! km^-2 yr^-1
   
      write(*,*) '  Muon flux from halo = ',phimuhalo,
     &  ' km^-2 yr^-1'

 150  continue

      write (*,*)
      write (*,*) '-------------------------------------------------------'
      write (*,*) 'The DarkSUSY example program has finished successfully.'
      write (*,*) 'Particle module that was used: ', moduletag
      write (*,*) '-------------------------------------------------------'
      write (*,*)
      write (*,*) '[simply call ''make -B dsmain_wimp DS_MODULE=<MY_MODULE>'' ' 
      write (*,*) ' if you want to try with a different module <MY_MODULE>.'
      write(*,*)  ' To see which ones are available, call ',
     &            '''../scr/list_modules''.]'
      write (*,*)
      goto 999
      stop
 999  end

      
      




      
      
c============================================================================
c============================================================================
c============================================================================






ccc Here are the subroutines that depend on the particle physics model

ccc ======== MSSM =========

#if MODULE_CONFIG == MODULE_mssm

c __________ dspmenterparameters
      
      subroutine dspmenterparameters
      implicit none
      integer modchoice
      real*8 amu,am2,ama,atanbe,amsq,atm,abm,am0,amhf,aa0,asgnmu
      character slhafile*128

      write(*,*) ' '
      write(*,*) 'What kind of SUSY model do you want to look at?'
      write(*,*) '   1 = MSSM-7'
      write(*,*) '   2 = cMSSM'
      write(*,*) '   3 = as read from an SLHA2 file'
      read(5,*) modchoice
      
      if (modchoice.eq.1) then ! MSSM-7
         write(*,*) 'Enter mu (GeV): '
         read(5,*) amu
         write(*,*) 'Enter M2 (GeV): '
         read(5,*) am2
         write(*,*) 'Enter MA (mass of CP-odd Higgs) (GeV): '
         read(5,*) ama
         write(*,*) 'Enter tan(beta): '
         read(5,*) atanbe
         write(*,*) 'Enter m0 (common sfermion mass scale) (GeV): '
         read(5,*) amsq
         write(*,*) 'Enter At/m0 [-3:3]: '
         read(5,*) atm
         write(*,*) 'Enter Ab/m0 [-3:3]: '
         read(5,*) abm

c...Now transfer these variables to the DarkSUSY common blocks, properly
c...setting up all relevant low-energy parameters         
         call dsgive_model(amu,am2,ama,atanbe,amsq,atm,abm)

c...mSUGRA
      elseif (modchoice.eq.2) then

         write(*,*) 'Enter m0 (scalar mass parameter @ GUT) (GeV): '
         read(5,*) am0
         write(*,*) 'Enter m_1/2 (gaugino mass param @ GUT) (GeV): '
         read(5,*) amhf
         write(*,*) 'Enter A0 (trilinear coupling @ GUT) (GeV): '
         read(5,*) aa0
         write(*,*) 'Enter sign of mu (+1 or -1): '
         read(5,*) asgnmu
         write(*,*) 'Enter tan(beta): '
         read(5,*) atanbe

c...Now transfer these parameters to the DarkSUSY common blocks, preparing
c...for RGE running
         call dsgive_model_isasugra(am0,amhf,aa0,asgnmu,atanbe)      

      elseif (modchoice.eq.3) then
         write(*,*) 'Enter the filename for your SLHA2 file: '
         read(5,'(A)') slhafile
         call dsgive_model_SLHA(slhafile,0) ! 0=no warnings, 1=print warnings

      else
         write(*,*) 'Not a valid choice of model: ',modchoice
         stop
      endif

      return
      end

c __________ dspmcompute

      subroutine dspmcompute(istatus,info)
c_______________________________________________________________________
c  sets up and computes the particle physics model
c  and checks if it is acceptable
c    istatus = 0   acceptable
c    istatus > 0   acceptable with care (warning)
c    istatus < 0   nonacceptable (calculation should not continue)
c    info          used to pass information on kind of problem
c  in particular, istatus=1 means experimentally excluded
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c  Torsten Bringmann: unified with SLHA, earlier dssusy(_isasugra)  
c    (08/05/14)       and dsmodelsetup_isasugra
c=======================================================================
      implicit none
      integer istatus,info
      integer unphys,hwarning,excl,opt
      real*8 gm2amp,dsgm2muon   ! g-2 amplitude
      character slhaout*128
      include 'dsmssm.h'
      include 'dsidtag.h'   ! unique model ID tag for each set of model parameters

      istatus=0
      info=0

      call dsmodelsetup(unphys,hwarning)
      if (unphys.ne.0) then
         istatus=-1
         info=unphys
         return
      endif
      if (hwarning.ne.0) then
         istatus=-2
         info=hwarning
         return
      endif
    
c... Check experimental constraints
      excl=0
      call dsacbnd(excl)
      if (excl.ne.0) then
         istatus=1
         info=excl
      endif

c      write(*,*) 'Testing sparticle widths...'
c      write(*,*) 'mass, width(kh1) = ',mass(kh1), width(kh1)
c      write(*,*) 'mass, width(kh2) = ',mass(kh2), width(kh2)
c      write(*,*) 'mass, width(kh3) = ',mass(kh3), width(kh3)
c      write(*,*) 'mass, width(ksqu(1)) = ',mass(ksqu(1)), width(ksqu(1))
c      write(*,*) 'mass, width(ksqu(6)) = ',mass(ksqu(6)), width(ksqu(6))
c      write(*,*) 'mass, width(ksl(2)) = ',mass(ksl(2)), width(ksl(2))
c      write(*,*) 'mass, width(ksl(5)) = ',mass(ksl(5)), width(ksl(5))
c      write(*,*) 'mass, width(kn(3)) = ',mass(kn(3)), width(kn(3))
c      write(*,*) 'partial -> Zh2, hdwidth(12,n) = ', hdwidth(9,1),hdwidth(9,2),
c     &                                              hdwidth(9,3), hdwidth(9,4)      


c
c     Now we calculate the MSSM contribution to the g-2 amplitude.
c
      gm2amp=dsgm2muon()
      write(*,*)
      write(*,*) 'g-2 amplitude [a_mu = (g-2)/2 = ]: ',gm2amp

c
c   We also write an SLHA file
c
      write(*,*) 
      write(*,*) 
     &  'Do you want to write out an SLHA2 file for your model?'
      write(*,*) '  0 = no, but do write out particle spectrum and 3-vertices'
      write(*,*) '  1 = yes, with full 6x6 sfermion mixing'
      write(*,*) '  2 = yes, with minimal flavour violation'
      write(*,*) ' (else) = do not write out anything' 
      read(*,*) opt

      if (opt.eq.1.or.opt.eq.2) then
         write(*,*) 'Give SLHA2 file name:'
         read(5,'(A)') slhaout
         call dsSLHAwrite(slhaout,opt) 
      endif
    
      if (opt.eq.0) then
         write(*,*) 'Output file name:'
         read(5,'(A)') slhaout 
         open (unit=30,file=slhaout)
         write(30,*) ' '
         write(30,*) '***** MODEL: ',idtag,' *****'
         call dswspectrum(30)
         call dswwidth(30)
         call dswvertx(30)
         close (30)
      endif

      end

#endif


ccc ======== generic wimp =========

#if MODULE_CONFIG == MODULE_generic_wimp

c __________ dspmenterparameters
      
      subroutine dspmenterparameters
      implicit none
      integer genselfconj,genpdg
      real*8 genmwimp,gensvann,genSI

      write(*,*) ' '
      write(*,*) 'Please specify model parameters for a generic WIMP:'
      write(*,*) '--------------------------------------------------- '
      write(*,*) ' '
      write(*,*) 'Enter WIMP mass (GeV): '
      read(5,*) genmwimp
      write(*,*) 'Enter if WIMP=antiWIMP [0 no, 1 yes]: '
      read(5,*) genselfconj
      write(*,*) 'Enter annihilation cross section (cm^3/s): '
      read(5,*) gensvann
      write(*,*) 'Enter PDG code for dominant annihilation channel'
      write(*,*) '(5 for bbar, 24 for W^+W^-, etc...): '
      read(5,*) genpdg
      write(*,*) 'Enter WIMP-nucleon scattering cross section (pb)'
      read(5,*) genSI
      write(*,*) '--------------------------------------------------- '
      
c... Now transfer these parameters to the DarkSUSY common blocks      
      call dsgivemodel_generic_wimp(genmwimp,genselfconj,
     &     gensvann,genpdg,genSI)

      return
      end

c __________ dspmcompute

      subroutine dspmcompute(istatus,info)
      implicit none
      integer istatus,info
      call dsmodelsetup(istatus,info)      
      end

#endif

ccc ======== empty model =========

#if MODULE_CONFIG == MODULE_empty

c __________ dspmenterparameters
      subroutine dspmenterparameters
      implicit none
      return
      end
      
c __________ dspmcompute
      subroutine dspmcompute(istatus,info)
      implicit none
      integer istatus,info      
      call dsmodelsetup(istatus,info)
      end

#endif



ccc ======== Silveira Zee (Scalar Singlet, etc) =========

#if MODULE_CONFIG == MODULE_silveira_zee

c __________ dspmenterparameters

      subroutine dspmenterparameters
      implicit none
      real*8 myms,mylambda

      write(*,*) ' '
      write(*,*) 'Please specify model parameters for the Silveira-Zee model:'
      write(*,*) '--------------------------------------------------- '
      write(*,*) ' '
      write(*,*) 'Enter singlet Higgs mass (GeV): '
      read(5,*) myms
      write(*,*) 'Enter quartic coupling lambda: '
      read(5,*) mylambda
      write(*,*) '--------------------------------------------------- '
      
c... Now transfer these parameters to the DarkSUSY common blocks      
      call dsgivemodel_silveira_zee(mylambda,myms)

      return
      end

c __________ dspmcompute
      
      subroutine dspmcompute(istatus,info)
      implicit none
      integer istatus,info
      integer ierr,iwarn

      istatus=0
      info=0

      call dsmodelsetup(ierr,iwarn)
      if (ierr.eq.0.and.iwarn.eq.0) return

      if(ierr.ne.0) istatus=-1
      if(iwarn.ne.0) then
        istatus=1
        info=iwarn
      endif
      
      end

#endif
