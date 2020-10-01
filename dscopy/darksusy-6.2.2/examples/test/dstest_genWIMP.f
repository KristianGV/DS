      program dstest_genWIMP
c
c     This program tests the generic WIMP module, and many general DarkSUSY routines
c
c     In addition to the DarkSUSY programs and data files, this test
c     program needs the input file dstest_genWIMP.mod provided with the
c     DarkSUSY distribution. The test program compares the calculated quantities  
c     with pre-calculated values. If the reported "Total number of errors in dstest"
c     after dstest has run successfully is not equal to zero, something
c     went wrong and you should check the output more carefully.


c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none
      real*8 oh2a,oh2b,xf,dsrdomega                      ! relic density
      real*8 tkd,dskdtkd,dskdmcut,mcut                   ! kinetic decoupling
      real*8 sigsip,sigsin,sigsdp,sigsdn                 ! direct detection
      real*8 sigij(27,27)                                ! direct detection
      integer i,nnuc,a(10),z(10),stoich(10)              ! direct detection
      real*8 v,e,si(10),sd(10),siff(10),sdff(10),t,rsi,rsd ! direct detection
      real*8 fluxgacdiff,fluxgac,fluxgaline              ! gamma-ray flux
      real*8 widthline1, widthline2, eline               ! gamma-ray flux
      real*8 jpsi,cospsi0,delta                          ! gamma-ray flux
      real*8 nsigvgaga,nsigvgaz,nsigvgacont,nsigvgacdiff ! gamma-rays
      integer istat, pdg2, nlines                        ! gamma-ray flux, etc
      real*8 egam,egath,dsgafluxsph                      ! gamma-ray flux
      real*8 egev,widthline, nsigvline, dscrgaflux_line_v0ann! monochromatic signals
      real*8 tpbess(3),pb_a,pb_b,pb_c                    ! pbar flux
      real*8 phiin,dspbdphidtaxi                         ! pbar flux
      real*8 tp,phidb, dsdbdphidtaxi                     ! db flux
      real*8 phiep, dsepdphidpaxi                        ! positron flux
      real*8 dscrsource, dscrsource_line                 ! DM source function for CR
      real*8 eth,thmax,rateea,ratesu, phimuhalo          ! neutrino telescopes
      integer kind,rtype,ptype                           ! neutrino telescopes
      real*8 rho0in, dscrmuflux_v0ann                    ! neutrino telescopes
      real*8 psi0,theta0,dsjfactor,dfactor,dsdfactor     ! los routines
      real*8 gm2amp,dsgm2muon                            ! g-2 amplitude
      character*12 refhalolab                            ! halo label
      real*8 mwrho0,dsdmdrho0mw                          ! halo routines
      real*8 dsmwimp, mdm                                ! DM mass
      integer unphys,warning,hwarning,acceptable,iend,ierr,iwar,nfc 
      real*8 dsabsq,tstart, tfinish 
      character message*80, modfile*200
      logical first
      data first/.true./
      
c Set up real testing of DS, by comparing to pre-computed values
      integer testlevel, testunit
      data testlevel/1/           ! 1 -- [default] standard output to screen,
                                  !      warnings if results are different from savfile 
                                  ! 2 -- minimal output, ~only warning
                                  ! 9 -- rewrite savfile 
                                  !      [NB: ALWAYS CHANGE SAVSHORT below to
                                  !           something else than default for this option!!!
                                  !           (unless intentionally changing the default)]
      character*200 savfile, SAVSHORT
      data SAVSHORT/'DStest_genWIMP_sav.dat'/       ! default : 'DStest_genWIMP_sav.dat'
      data testunit/98/
 
c
c     Here we include various "global variables" defined in common blocks
c
      include 'dsio.h'      ! for I/O handling
      include 'dsidtag.h'   ! unique model ID tag for each set of model parameters
      include 'dshmcom.h'   ! for rho0 of the halo model in jpsi rescaling

      call CPU_TIME(tstart)

c
c     This call should be the first call in any program using DarkSUSY. It 
c     initializes some global variables and calls various other modules to 
c     set them up properly.
c
      call dsinit
      
      if (moduletag.ne.'generic_wimp') then
         write(*,*)
     &     'Sorry, dstest_genWIMP is only designed for the generic_wimp module and'
        write(*,*) 'you compiled it with module ',moduletag
        write(*,*)
        write(*,*) 'Set DS_MODULE=generic_wimp in /examples/test/makefile' //
     &              ' and try again...'
        stop
      endif
      
      write (*,*) 'Now starting dstest_genWIMP... '
      savfile=savshort
      call dsdatafile(savfile,savshort)  
      
      if (testlevel.eq.9) then
        write(*,*) 'Now WRITING numerical values to DStest savfile: ',
     &              savfile
        open(unit=testunit, file=savfile, status="unknown", err=3000)      
      elseif (testlevel.eq.1.or.testlevel.eq.2) then
        write(*,*) 
     &    'Now COMPARING to numerical values from DStest savfile: ',
     &    savfile
        open(unit=testunit, file=savfile, err=3000)      
      else
        write(*,*) 'ERROR: unsupported value of testlevel=',testlevel 
      endif

c
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
c     integration factor J (or D) in the direction of observation, which we
c     define as the direction which forms an angle psi0 with respect to
c     the direction of the galactic centre (e.g. cospsi0 = 1 is the
c     galactic center direction, cospsi0 = -1 the antigalactic centre).
c
      write(*,*) 'Calculating D-factor and J-factor...'
      cospsi0=1.d0
      delta=1.d-3
      psi0=dacos(cospsi0)
      theta0=dacos(1.d0-delta/(8.d0*datan(1.d0)))
      dfactor=dsdfactor(refhalolab,psi0,theta0) ! kpc sr GeV cm^-3
      if (testlevel.ne.2) write(*,*) 'D_factor = ',dfactor,
     &  ' kpc sr GeV cm^-3'
      jpsi=dsjfactor(refhalolab,psi0,theta0) ! kpc sr GeV^2 cm^-6
      if (testlevel.ne.2) write(*,*) 'J_factor galactic center =',jpsi,
     &  ' kpc sr GeV^2 cm^-6'
      
      call testvalue(testunit,testlevel,'j_GC',jpsi)


*******************************************************
***       read in + set up genric WIMP models       ***
*******************************************************
      modfile = 'models/dstest_genWIMP.mod'

      open (unit=11,file=modfile,err=3100)
      read (11,*, end=3100, err=3100) ! this skips the header
      iend=0

c... start loop to read in models      
 1000 call read_model_genWIMP(11,0,iend,ierr)
 
      if (ierr.ne.0) then
         write (*,*) 'Error while reading model file ',modfile
         write (*,*) 'The DarkSUSY test cannot continue'
         stop
      else if (iend.eq.1) then
         goto 2000 ! exit program
      else
         write (*,*)
         call dswrite(0,0,'Model parameters read from genric_wimp model file.')
      endif

      write(*,*) '***************************************'
      write(*,*) '***** MODEL: ',idtag,' *****'
      write(*,*) '***************************************'

c
c     After providing model parameters with a routine like dsgivemodel_xxx,
c     the next step is always a call to the routine dsmodelsetup. This
c     actually sets up the model, by calculating the mass spectrum, 
c     relevant 3-particle vertices etc. dsmodelsetup return two flags
c     to indicate if everything is OK (both flags equal 0) or not. The 
c     flag unphys<0 means that the model is theoretically inconsistent. 
c     The flag hwarning indicates the breakdown of approximations 
c     used in radiative corrections in the Higgs sector.
c
      call dsmodelsetup(unphys,hwarning)
c
c     Test summary:
      if (unphys.eq.0.and.hwarning.eq.0) then
         acceptable=0
      else 
         acceptable=-1
      endif
      write(*,*)
      write (message,*) 'acceptable =',acceptable,
     &  ' (0=OK, -1=not OK)'
      call dswrite(0,0,message)

      mdm = dsmwimp()

      write(*,*) '  Dark matter mass = ',mdm
      call testvalue(testunit,testlevel,'mdm',mdm)


***********************************************
***   relic density + kinetic decoupling    ***
***********************************************

      write(*,*)
      if (testlevel.eq.2) write(*,*) 'Calculating relic density...'
      
      if (testlevel.ne.2)  write(*,*) 'Calculating omega h^2 without coannihilations,',
     &     ' please be patient...'
      oh2a=dsrdomega(0,1,xf,ierr,iwar,nfc)
      if (testlevel.ne.2) write(*,*) '  without coannihilations Oh2 = ',
     &  oh2a,ierr,iwar

c... The generic WIMP model has no coannihilating particles
c      if (testlevel.ne.2) write(*,*) 'Calculating omega h^2 with coannihilations,',
c     &     ' please be patient...'
c      oh2b=dsrdomega(1,1,xf,ierr,iwar,nfc)
c      if (testlevel.ne.2) write(*,*) ' with coannihilations Oh2 = ',oh2b,ierr,iwar
      
      if (testlevel.ne.2) write(*,*) '  Chemical decoupling (freeze-out) occured at'
      if (testlevel.ne.2) write(*,*) '  T_f = ',mdm/xf,' GeV.'

      call testvalue(testunit,testlevel,'oh2',oh2a)
c      call testvalue(testunit,testlevel,'oh2 with coannihilation',oh2b)
      call testvalue(testunit,testlevel,'x_f',xf)

c     Now let's calculate the kinetic decoupling, and the smallest halos
c     we can have for this model.
      tkd=dskdtkd(1) ! 1=full calculation
      mcut=dskdmcut(mdm,tkd,1) ! 1 = true cut-off
      call testvalue(testunit,testlevel,'Tkd',tkd)
      call testvalue(testunit,testlevel,'Mcut',mcut)

      if (testlevel.ne.2) then
        if (tkd.gt.1.d3) then
          write(*,*) 'Kinetic decoupling happens much before the QCD phase '//
     &               'transition and cannot be determined accurately.'
        else
          write(*,*) 'Kinetic decoupling temperature, Tkd = ',tkd, ' MeV'
          write(*,*) 'The resulting cutoff in the power spectrum '//
     &               'corresponds to a mass of',' M_cut/M_sun = ', mcut
        endif  
      endif    
      if (testlevel.ne.2) write(*,*) ' '


***********************************************
***           direct detection              ***
***********************************************
c
c
      call dsddset('sf_m','best')
      call dsddset('sf_sigma','best')

      write (*,*) 'Calculating DM-nucleon scattering cross sections...'
      v=0.d0
      e=0.d0
      call dsddsigma(v,e,1,0,sigij,ierr)
      sigsin=sigij(1,1) ! SI cross section
      sigsdn=sigij(4,4) ! SD cross setion
      call dsddsigma(v,e,1,1,sigij,ierr)
      sigsip=sigij(1,1) ! SI cross section
      sigsdp=sigij(4,4) ! SD cross section

      if (testlevel.ne.2) then
c        write (*,*) ' ierr=',ierr
        write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
        write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
c Spin-dependent scattering not activated; would be switched on with call to 
c dsgivemodel_generic_wimp_opt directly after dsgivemodel_generic_wimp!
c        write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
c        write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36
      endif
      call testvalue(testunit,testlevel,'sigsip',sigsip*1.0d36)
      call testvalue(testunit,testlevel,'sigsin',sigsin*1.0d36)
c      call testvalue(testunit,testlevel,'sigsdp',sigsdp*1.0d36)
c      call testvalue(testunit,testlevel,'sigsdn',sigsdn*1.0d36)

c
c     Now take a more complex example with scattering off many
c     nuclides and/or a compound
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
        siff(i)=sigij(1,1)
        sdff(i)=sigij(4,4)
        if (testlevel.ne.2) then
          write(*,*) '  A=',a(i),' Z=',z(i),
     &         ' sigsi (pb) = ',si(i)*1.0d36 !,
c     &         ' sigsd (pb) = ',sd(i)*1.0d36
          write(*,*) 
     &         ' sigsi*ff(pb)=',siff(i)*1.0d36 !,
c     &         ' sigsd*ff(pb)=',sdff(i)*1.0d36
        endif
        call testvalue(testunit,testlevel,'sigsi(i)',si(i)*1.0d36)
c        call testvalue(testunit,testlevel,'sigsd(i)',sd(i)*1.0d36)
        call testvalue(testunit,testlevel,'sigsi*ff(i)',siff(i)*1.0d36)
c        call testvalue(testunit,testlevel,'sigsd*ff(i)',sdff(i)*1.0d36)
      enddo
      stoich(1)=1
      stoich(2)=1
      t = 651.3d0
c... NB: rate calculations assume that the Scalar Singlet makes up 100%
c... of the local DM density. If this is not the case, the rate must
c... be (linearly) rescaled by the user    
      call dsdddrde(t,e,2,a,z,stoich,rsi,rsd,1)
      if (testlevel.ne.2) write(*,*) 
     &  ' NaI :  dRdE [counts/kg-day-keV] (SI)=',rsi !,' (SD)=',rsd
      call testvalue(testunit,testlevel,'NaI (SI)',rsi)
c      call testvalue(testunit,testlevel,'naI (SD)',rsd)
      
      
      
***********************************************
***           Gamma-ray fluxes              ***
***********************************************
c
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

c     2) gamma-ray flux with continuum energy spectrum integrated above
c        some given threshold egath (GeV). E.g. :
c
      egath=1.d0
      fluxgac=dsgafluxsph(egath,0,1.d0,refhalolab,psi0,theta0,istat) !ph cm^-2 s^-1
c     again, the result is identical to:
c      
c      fluxgac=dscrgaflux_v0ann(egath,0,jpsi,1.d0,istat)  !ph cm^-2 s^-1
      
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
        if (testlevel.ne.2) then
          write(*,*) 'photon flux from line No.',i,' = ',fluxgaline,' ph/(cm^2 s)'
          write(*,*) '        [at E = ',eline,' +/- ',widthline,' GeV]'
        endif  
        call testvalue(testunit,testlevel,'fluxgaline',fluxgaline)
        call testvalue(testunit,testlevel,'eline',eline)
        if (i.lt.nlines) fluxgaline=
     &  dscrgaflux_line_v0ann(i+1,eline,widthline,jpsi,1.d0,istat) 
      end do

c
c     Lasrly, if we want to do the gymnastics ourselves converting cross
c     sections to fluxes, we can also call the source functions,
c     dscrsource and dscrsource_line, to get the number of photons 
c     (2 for gamma gamma, 1 for Z gamma and whatever the number is for
c     continuous gammas) times the cross section and divided by a symmetry
c     factor (=2 for self-conjugate DM). 
c     The numbers calculated below are given by
c        N_gammas * (sigma v) 
c     in units of cm^3 s^-1. For the differential result the unit
c     is of course multiplied by GeV^-1.
      egam=1.0d0
      nsigvgacont=dscrsource(egam,0,22,2,0d0,istat)*mdm**2   
      nsigvgacdiff=dscrsource(egam,1,22,2,0d0,istat)*mdm**2   
c     Below, we pretend not to know in which order the lines are returned
c     for our particle model, so we check explicitly which one is gg, and which gZ       
      do i=1,nlines
        nsigvline=dscrsource_line(22,i,2,0d0,egev,widthline,pdg2,istat) 
     &          *dsmWIMP()**2
        write(*,*) '    nsigvline ',i,' = ',nsigvline
        call testvalue(testunit,testlevel,'nsigvline',nsigvline)
      enddo 
      
      if (testlevel.ne.2) then
        write(*,*) '  fluxgacdiff = ',fluxgacdiff,' ph/(cm^2 s GeV)'
        write(*,*) '      fluxgac = ',fluxgac,' ph/(cm^2 s)'
        write(*,*) '  nsigvgacont = ',nsigvgacont
        write(*,*) ' nsigvgacdiff = ',nsigvgacdiff,' GeV^-1'
      endif
      
      call testvalue(testunit,testlevel,'fluxgacdiff',fluxgacdiff)
      call testvalue(testunit,testlevel,'fluxgac',fluxgac)
      call testvalue(testunit,testlevel,'nsigvgacont',nsigvgacont)
      call testvalue(testunit,testlevel,'nsigvgacdiff',nsigvgacdiff)


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
      if (testlevel.ne.2) then
         write(*,*) 'solar modulated pbar flux at 0.35 ',
     &              'GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ', pb_a
         write(*,*) 'solar modulated pbar flux at 1.76 ',
     &              'GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ', pb_b
         write(*,*) 'solar modulated pbar flux at 3.00 ',
     &              'GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ', pb_c
      endif

      call testvalue(testunit,testlevel,'pbar flux at 0.35 GeV',pb_a)
      call testvalue(testunit,testlevel,'pbar flux at 1.76 GeV',pb_b)
      call testvalue(testunit,testlevel,'pbar flux at 3.00 GeV',pb_c)
      
c
c      We can also calculate the flux of antideuterons. The call is very
c      similar to that for pbar. The dbar's use the same diffusion model
c      as the pbar and can thus be changed with a call to dspbset.
c
      write (*,*) 'Calculating antideuteron fluxes...'
      tp=1.0d0

c      phidb=dsdbdphidtaxi(tp,phiin,1,refhalolab) ! no tables, slower
      phidb=dsdbdphidtaxi(tp,phiin,4,refhalolab) ! use tables, faster after slower startup

      if (testlevel.ne.2)  write(*,*)
     &  '  solar modulated dbar flux at 1.00 GeV [GeV^-1 cm^-2 s^-1 sr^-1]: ',
     &  phidb

      call testvalue(testunit,testlevel,'dbar flux at 1.00 GeV',phidb)

c
c     We now continue with the rates of positrons from DM
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
      if (testlevel.ne.2) write (*,*) '  phiep=',phiep,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'

      call testvalue(testunit,testlevel,'pos flux at 1.00 GeV',phiep)


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
                                        !  *all* of the local DM is composed of 
                                        !  DM in the current particle module)
      call dssenu_rates(eth,thmax,1,rtype,ptype,rho0in,rateea,ratesu,
     &  istat)

      if (testlevel.ne.2) then
        write(*,*) '  Flux from the Earth = ',rateea, ' km^-2 yr^-1'
        write(*,*) '  Flux from the Sun =   ',ratesu, ' km^-2 yr^-1'
      endif

      call testvalue(testunit,testlevel,'nu flux from Earth',rateea)
      call testvalue(testunit,testlevel,'nu flux from Sun',ratesu)


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
c     in neutrino telescopes that would occur from DM
c     annihilations in the galactic halo. The neutrino-induced muon
c     flux above a threshold Eth is given by
c     
      write(*,*) 
     &  'Calculating neutrino-induced muon fluxes from the halo...'
      eth=1.0d0

      phimuhalo=dscrmuflux_v0ann(eth,0,jpsi,1.d0,istat) ! km^-2 yr^-1
      if (testlevel.ne.2) write(*,*) '  Muon flux from halo = ',phimuhalo,
     &  ' km^-2 yr^-1'

      call testvalue(testunit,testlevel,'muon flux from halo',phimuhalo)

      write(*,*) 
      write(*,*) 

c
c     end of loop to read in new models
c
      goto 1000

c exit program
2000  continue
      close (11)
      close (98)

**********************************************
***   Last lines of the main test program  ***
**********************************************

      call CPU_TIME(tfinish)
      
      write(*,*)
      write(*,*) 'The DarkSUSY test program ran successfully'
      write (*,*) 'Total time needed (in seconds): ',tfinish-tstart
      call testvalue(testunit,testlevel,'final call',gm2amp) ! NB: gm2amp is just a dummy
                                                             ! argument in this case
      write(*,*)
      write(*,*)
      stop
      
 3000 write(*,*)
      write(*,*) 'ERROR with file I/O for saved data file: ',savfile
      stop     
 3100 write(*,*) 
      write(*,*) 'ERROR with file I/O for input model file!'     
      end
********************************************************************
********************************************************************







********************************************************************
***   Subroutines used by the test program or otherwise useful.  ***
********************************************************************


      subroutine read_model_genWIMP(lunit,nmodel,iend,ierr)
c
c     To read in a model from a file
c
c     NOTE: We here show how you can read in a model file and define all
c     DarkSUSY model parameters. 

      implicit none
      include 'dsidtag.h'
      integer nmodel,lunit,iend,ierr
      real*8 mgenwimp,svann,SI
      integer pdgann
      integer i, intselfconj
      character*40 message
 200  format (1x,a12,3(1x,e14.8),2(1x,i8))
      ierr=0
      iend=0
c... When nmodel<0, skip -nmodel lines
      if (nmodel.lt.0) then
         do i=1,-nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
         return
      endif
c... If nmodel>0, read n-th model (assumes header line)
      if (nmodel.gt.0) then
         do i=1,nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
      endif
c... If nmodel=0, read next model
      read (lunit,200,end=1000,err=3000) 
     &     idtag,mgenwimp,svann,SI,intselfconj,pdgann
      call dsgivemodel_generic_wimp(mgenwimp,intselfconj,svann,pdgann,SI)
      return

 1000 continue
      iend=1
      write (message,*) 'End of model file (unit=',lunit,')'
c      call dswrite(1,0,message)
      return
 3000 continue
      ierr=1
      write(*,*) 'I/O ERROR, but not eof!'
      stop
      end


      subroutine testvalue(testunit,testlevel,msg,value)
        implicit none
        integer testunit,testlevel
        character*(*) msg
        real*8 value,tmp,eps,epscomp
        data eps/0.003/     !required accuracy for match with data file

        integer memory    ! to count total number of errors
        save memory
        data memory /0/


        if (msg.eq.'final call') then
          write(*,*) 'Total number of errors in dstest_genWIMP: ',memory
          write(*,*)
          write(*,*)
c The simpler stop + integer only works for FORTRAN >= 2008...   
c          stop memory
          open(unit=99, file="dstest_genWIMP.out", status="unknown")
          write(99,'(i3)') memory
          close(99)
          stop
c          return
        endif

        if (testlevel.eq.9) then
          write (testunit,*,err=30) value
        else
           read (testunit,*,err=30,end=30) tmp
           if (tmp.ne.0.d0) then
              epscomp=abs((value-tmp)/tmp)
           else
              epscomp=abs(value-tmp)
           endif
           
            if (epscomp.gt.eps) then
              memory=memory+1
              write(*,*)
              write(*,*) '====================================='
              write(*,*) 'ERROR in computation of ',msg,' !!!',
     &                   '   (# ',memory,'in dstest_genWIMP)'
              write(*,*) 'calculated value:      ',value
              write(*,*) 'saved reference value: ',tmp
              write(*,*) '====================================='
              write(*,*)
            endif
        endif
        return
 30     write(*,*) 'testvalue: ERROR with file I/O for savfile!'     
        stop
      end !testvalue  
