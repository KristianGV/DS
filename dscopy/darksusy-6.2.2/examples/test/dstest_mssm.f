      program dstest_mssm
c
c     This program tests the MSSM module, and many general DarkSUSY routines
c
c     In addition to the DarkSUSY programs and data files, this test
c     program needs the input files dstest_mssm.mod, dstest_mssm_long.mod,
c     dstest_slha.mod and dstest_slha_long.mod, all provided with the
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
      real*8 fluxgacdiff,fluxgac,fluxgaline1,fluxgaline2 ! gamma-ray flux
      real*8 widthline1, widthline2, eline1, eline2      ! gamma-ray flux
      real*8 jpsi,cospsi0,delta                          ! gamma-ray flux
      real*8 nsigvgaga,nsigvgaz,nsigvgacont,nsigvgacdiff ! gamma-rays
      integer istat, pdg2                                ! gamma-ray flux, etc
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
      integer unphys,warning,hwarning,acceptable,iend,ierr,iwar,nfc 
      integer modtype, testtype
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
                                  !      [NB: ALWAYS CHANGE SAVSHORT / SAVLONG below to
                                  !           something else than default for this option!!!
                                  !           (unless intentionally changing the default)]
      character*200 savfile, SAVSHORT, SAVLONG
      data SAVSHORT/'DStest_MSSM_sav.dat'/       ! default : 'DStest_MSSM_sav.dat'
      data SAVLONG/'DStest_MSSM_long_sav.dat'/   ! default : 'DStest_MSSM_long_sav.dat'
      data testunit/98/

c
c     Here we include various "global variables" defined in common blocks
c
      include 'dsmssm.h'    ! particle masses, susy parameters
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
      
      if (moduletag.ne.'MSSM') then
         write(*,*)
     &     'Sorry, dstest_mssm is only designed for the MSSM module and'
        write(*,*) 'you compiled it with module ',moduletag
        write(*,*)
        write(*,*) 'Set DS_MODULE=MSSM in /examples/test/makefile' //
     &              ' and try again...'
        stop
      endif
      
      write (*,*) 'Now starting dstest_mssm... '
      write(*,*)  'Do you want the short [0] or long [1] option?'
      read(*,*) testtype
      if (testtype.eq.1) then
        call dsdatafile(savfile,savlong)
      else
        call dsdatafile(savfile,savshort)  
      endif
      
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


***********************************************
***       read in + set up MSSM models      ***
***********************************************
      modtype=1 ! start with MSSM7 models
      if (testtype.eq.1) then
        modfile = 'models/dstest_mssm_long.mod'
      else
        modfile = 'models/dstest_mssm.mod'
      endif

 500  if (modtype.eq.2) then ! continue with reading in SLHA files
        if (testtype.eq.1) then
          modfile = 'models/dstest_mssm_slha_long.mod'
        else
          modfile = 'models/dstest_mssm_slha.mod'
        endif
      endif

      open (unit=11,file=modfile,err=3100)
      read (11,*, end=3100, err=3100) ! this skips the header
      iend=0

c... start loop to read in models      
 1000 if (modtype.eq.1) call read_model_mssm(11,0,iend,ierr)
      if (modtype.eq.2) call read_model_slha(11,0,iend,ierr)  
 
      if (ierr.ne.0) then
         write (*,*) 'Error while reading model file ',modfile
         write (*,*) 'The DarkSUSY test cannot continue'
         stop
      else if (iend.eq.1) then
         goto 2000 ! continue with slha files / exit program
      else
         write (*,*)
         if (modtype.eq.1) call dswrite(0,0,'Model parameters read from MSSM model file.')
         if (modtype.eq.2) call dswrite(0,0,'Model parameters read from SLHA model file.')
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
c     flag unphys<0 means that the model is theoretically inconsistent, and
c     (for the MSSM module) unphys=1 means that the neutralino is
c     not the LSP. The flag hwarning indicates the breakdown of approximations 
c     used in radiative corrections in the Higgs sector.
c
      call dsmodelsetup(unphys,hwarning)
c
c     For a very rough check of whether a given SUSY model is already excluded
c     by collider or flavour observables, we can call dsacbnd. This routine
c     returns a flag 'warning' which equals 0 if there is no relevant bound
c     (implemented in DarkSUSY). A call to dsacset allows to switch to 
c     previous versions of those constraints, if needed for comparison.
c
      call dsacbnd(warning)
      call testvalue(testunit,testlevel,'excl',dble(warning))
c     Test summary:
      if (warning.eq.0.and.unphys.eq.0.and.hwarning.eq.0) then
         acceptable=0
      elseif (warning.ne.0.and.(unphys.eq.0.and.hwarning.eq.0)) then
         acceptable=1
      else 
         acceptable=-1
      endif
      write(*,*)
      write (message,*) 'acceptable =',acceptable,
     &  ' (0=OK, -1=not OK, 1=excluded)'
      call dswrite(0,0,message)
c     If the model is not acceptable,  print out why
      if (unphys.ne.0) then
         call dswunph(6,unphys)  ! print meaning of flag unphys (returned by 
                                 ! dsmodelsetup) to I/O unit 6
         goto 1000               ! next model
      endif
      if (warning.ne.0) call dswexcl(6,warning)  ! print meaning warning 
      if (hwarning.ne.0) call dswhwarn(6,hwarning) ! print meaning hwarning 


c
c     We now have all the masses and couplings calculated. Let's print
c     some out to see what they are. dsabsq is just the square of the
c     absolute value of it's complex*16 argument.
c
      if (testlevel.ne.2) then
        write(*,*) '  Neutralino mass = ',mass(kn(1))
        write(*,*) '  Gaugino fraction = ',
     &    dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
        write(*,*) '  H1 mass =  ',mass(kh1),width(kh1)
        write(*,*) '  H2 mass =  ',mass(kh2),width(kh2)
        write(*,*) '  H3 mass =  ',mass(kh3),width(kh3)
        write(*,*) '  H+- mass = ',mass(khc),width(khc)
      endif
      
      call testvalue(testunit,testlevel,'mchi',mass(kn(1)))
      call testvalue(testunit,testlevel,'Z_g',
     &               dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2)))
      call testvalue(testunit,testlevel,'mh1',mass(kh1))
      call testvalue(testunit,testlevel,'Gh1',width(kh1))
      call testvalue(testunit,testlevel,'mh2',mass(kh2))
      call testvalue(testunit,testlevel,'Gh2',width(kh2))
      call testvalue(testunit,testlevel,'mh3',mass(kh3))
      call testvalue(testunit,testlevel,'Gh3',width(kh3))
      call testvalue(testunit,testlevel,'mhc',mass(khc))
      call testvalue(testunit,testlevel,'Ghc',width(khc))

c
c     Now we calculate the MSSM contribution to the g-2 amplitude
c
      gm2amp=dsgm2muon()
      if (testlevel.ne.2) write(*,*) 'g-2 amplitude [a_mu = (g-2)/2 = ]: ',
     &   gm2amp

      call testvalue(testunit,testlevel,'g-2',gm2amp)


***********************************************
***   relic density + kinetic decoupling    ***
***********************************************
c
c     The relic density is calculated by the function dsrdomega.
c     The first argument determines if coannihilations should be included.
c     If the argument is 0, no coannihilations are included, if it is 1,
c     all relevant coannihilations are included (charginos, neutralinos
c     and sfermions if light enough). One can also choose to only
c     include a partial set of coannihilations: argument=2, only chargino
c     and neutralino coannihilaitons, argument=3, only sfermion 
c     coannihiltions.
c     The second argument determines the accuracy. Set it to 1 for a faster 
c     calculation (e.g. only including coannihilating particles with mass 
c     difference below 30%), or to 0 for a more accurate calculation.
c     In practice, the fast option is more than adequate (better than
c     5% accuracy). 
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
      if (testlevel.eq.2) write(*,*) 'Calculating relic density...'
      
      if (testlevel.ne.2)  write(*,*) 'Calculating omega h^2 without coannihilations,',
     &' please be patient...'
      oh2a=dsrdomega(0,1,xf,ierr,iwar,nfc)
      if (testlevel.ne.2) write(*,*) '  without coannihilations Oh2 = ',
     &  oh2a,ierr,iwar

      if (testlevel.ne.2) write(*,*) 'Calculating omega h^2 with coannihilations,',
     &     ' please be patient...'
         oh2b=dsrdomega(1,1,xf,ierr,iwar,nfc)

      if (testlevel.ne.2) write(*,*) ' with coannihilations Oh2 = ',oh2b,ierr,iwar
      
      if (testlevel.ne.2) write(*,*) '  Chemical decoupling (freeze-out) occured at'
      if (testlevel.ne.2) write(*,*) '  T_f = ',mass(kn(1))/xf,' GeV.'

      call testvalue(testunit,testlevel,'oh2',oh2a)
      call testvalue(testunit,testlevel,'oh2 with coannihilation',oh2b)
      call testvalue(testunit,testlevel,'x_f',xf)

c     Now let's calculate the kinetic decoupling, and the smallest halos
c     we can have for this model.
      tkd=dskdtkd(1) ! 1=full calculation
      if (testlevel.ne.2) write(*,*) 'Kinetic decoupling temperature, Tkd = ',tkd, ' MeV'
      call testvalue(testunit,testlevel,'Tkd',tkd)
      
      mcut=dskdmcut(mass(kn(1)),tkd,1) ! 1 = true cut-off
      if (testlevel.ne.2) write(*,*) ' The resulting cutoff in the power spectrum'//
     &   ' corresponds to a mass of',' M_cut/M_sun = ', mcut
      call testvalue(testunit,testlevel,'Mcut',mcut)
      if (testlevel.ne.2) write(*,*) ' '


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
c
c... For the MSSM, there are also model-specific options implemented
      call dsddset_mssm('dn_nopole')

      write (*,*) 'Calculating DM-nucleon scattering cross sections...'
      v=0.d0
      e=0.d0
      call dsddsigma(v,e,1,0,sigij,ierr)
      sigsin=sigij(1,1) ! SI cross section
      sigsdn=sigij(4,4) ! SD cross setion
      call dsddsigma(v,e,1,1,sigij,ierr)
      sigsip=sigij(1,1) ! SI cross section
      sigsdp=sigij(4,4) ! SD cross section
c...Note: if we are only interested in the spin-independent and spin-
c...dependent scattering cross sections, we could call the wrapper
c...routine dsddsigmanucleon(v,e,sigsip,sigsin,sigsdp,sigsdn,ierr)
c...instead of calling dsddsigma above. dsddsigma is the more general
c...interface function the particle physics module should provide though.      

      if (testlevel.ne.2) then
c        write (*,*) ' ierr=',ierr
        write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
        write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
        write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
        write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36
        write(*,*) ' proton: sigsi (pb) = ',sigsip*1.0d36,
     &       ' sigsd (pb) = ',sigsdp*1.0d36
      endif
      call testvalue(testunit,testlevel,'sigsip',sigsip*1.0d36)
      call testvalue(testunit,testlevel,'sigsin',sigsin*1.0d36)
      call testvalue(testunit,testlevel,'sigsdp',sigsdp*1.0d36)
      call testvalue(testunit,testlevel,'sigsdn',sigsdn*1.0d36)
      call testvalue(testunit,testlevel,'proton sigsi',sigsip*1.0d36)
      call testvalue(testunit,testlevel,'proton sigsd',sigsdp*1.0d36)

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
     &         ' sigsi (pb) = ',si(i)*1.0d36,
     &         ' sigsd (pb) = ',sd(i)*1.0d36
          write(*,*) 
     &         ' sigsi*ff(pb)=',siff(i)*1.0d36,
     &         ' sigsd*ff(pb)=',sdff(i)*1.0d36
        endif
        call testvalue(testunit,testlevel,'sigsi(i)',si(i)*1.0d36)
        call testvalue(testunit,testlevel,'sigsd(i)',sd(i)*1.0d36)
        call testvalue(testunit,testlevel,'sigsi*ff(i)',siff(i)*1.0d36)
        call testvalue(testunit,testlevel,'sigsd*ff(i)',sdff(i)*1.0d36)
      enddo
      stoich(1)=1
      stoich(2)=1
      t = 651.3d0
c... NB: rate calculations assume that the neutralino makes up 100%
c... of the local DM density. If this is not the case, the rate must
c... be (linearly) rescaled by the user    
      call dsdddrde(t,e,2,a,z,stoich,rsi,rsd,1)
      if (testlevel.ne.2) write(*,*) 
     &  ' NaI :  dRdE [counts/kg-day-keV] (SI)=',rsi,' (SD)=',rsd
      call testvalue(testunit,testlevel,'NaI (SI)',rsi)
      call testvalue(testunit,testlevel,'naI (SD)',rsd)
      
      
      
***********************************************
***           Gamma-ray fluxes              ***
***********************************************
c
c     Next we compute cosmic-ray rates from neutralino annihilations 
c     in the galactic halo. We start with the gamma-ray flux, which
c     factorizes into part that depends on the particle physics (susy)
c     and the line-of-sight integral j computed before scanning in parameter
c     space. There are two main ways of computing this flux, either by
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

c     2) gamma-ray flux with continuum energy spectrum integrated above
c        some given threshold egath (GeV). E.g. :
c
      egath=1.d0
      fluxgac=dsgafluxsph(egath,0,1.d0,refhalolab,psi0,theta0,istat) !ph cm^-2 s^-1
c     again, the result is identical to:
c      
c      fluxgac=dscrgaflux_v0ann(egath,0,jpsi,1.d0,istat)  !ph cm^-2 s^-1
      
c     3) monochromatic gamma-ray flux induced by 1-loop annihilation
c        processes into a 2-body final state containing a photon. For SUSY, there
c        are two such final states: the 2 photon final state (with 
c        energy of the photons equal to the neutralino mass mx) and
c        the final state with a photon and a Z boson (with energy of 
c        the photon equal to mx * (1 - mz**2 / (4 * mx**2)) and mz
c        the mass of the Z boson). Note that asking dscrgaflux_line_v0ann for
c        a third line (first argument=3) returns zero and istat=1 
c
      fluxgaline1=
     &  dscrgaflux_line_v0ann(1,eline1,widthline1,jpsi,1.d0,istat) !ph cm^-2 s^-1 
      fluxgaline2=
     &  dscrgaflux_line_v0ann(2,eline2,widthline2,jpsi,1.d0,istat) !ph cm^-2 s^-1 
c
c     Lasrly, if we want to do the gymnastics ourselves converting cross
c     sections to fluxes, we can also call the source functions,
c     dscrsource and dscrsource_line, to get the number of photons 
c     (2 for gamma gamma, 1 for Z gamma and whatever the number is for
c     continuous gammas) times the cross section and divided by a symmetry
c     factor (=2 in the case of neutralino DM). 
c     The numbers calculated below are given by
c        N_gammas * (sigma v) 
c     in units of cm^3 s^-1. For the differential result the unit
c     is of course multiplied by GeV^-1.
      egam=1.0d0
      nsigvgacont=dscrsource(egam,0,22,2,0d0,istat)*mass(kn(1))**2   
      nsigvgacdiff=dscrsource(egam,1,22,2,0d0,istat)*mass(kn(1))**2   
c     Below, we pretend not to know in which order the lines are returned
c     for the mssm, so we check explicitly which one is gg, and which gZ       
      nsigvline=dscrsource_line(22,1,2,0d0,egev,widthline,pdg2,istat) 
     &          *mass(kn(1))**2
      if (pdg2.eq.22) nsigvgaga=nsigvline
      if (pdg2.eq.23) nsigvgaz=nsigvline
      nsigvline=dscrsource_line(22,2,2,0d0,egev,widthline,pdg2,istat)
     &          *mass(kn(1))**2
      if (pdg2.eq.22) nsigvgaga=nsigvline
      if (pdg2.eq.23) nsigvgaz=nsigvline
      
      if (testlevel.ne.2) then
        write(*,*) '  fluxgacdiff = ',fluxgacdiff,' ph/(cm^2 s GeV)'
        write(*,*) '      fluxgac = ',fluxgac,' ph/(cm^2 s)'
        write(*,*) '  fluxgaline1 = ',fluxgaline1,' ph/(cm^2 s)'
        write(*,*) '        [at E = ',eline1,' +/- ',widthline1,' GeV]'
        write(*,*) '  fluxgaline2 = ',fluxgaline2,' ph/(cm^2 s)'   
        write(*,*) '        [at E = ',eline2,' +/- ',widthline2,' GeV]'
        write(*,*) '  nsigvgacont = ',nsigvgacont
        write(*,*) ' nsigvgacdiff = ',nsigvgacdiff,' GeV^-1'
        write(*,*) '    nsigvgaga = ',nsigvgaga
        write(*,*) '     nsigvgaz = ',nsigvgaz
      endif
      
      call testvalue(testunit,testlevel,'fluxgacdiff',fluxgacdiff)
      call testvalue(testunit,testlevel,'fluxgac',fluxgac)
      call testvalue(testunit,testlevel,'fluxgaline1',fluxgaline1)
      call testvalue(testunit,testlevel,'fluxgaline2',fluxgaline2)
      call testvalue(testunit,testlevel,'nsigvgacont',nsigvgacont)
      call testvalue(testunit,testlevel,'nsigvgacdiff',nsigvgacdiff)
      call testvalue(testunit,testlevel,'nsigvgaga',nsigvgaga)
      call testvalue(testunit,testlevel,'nsigvgaz',nsigvgaz)


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
                                        !  neutralinos)
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
c     in neutrino telescopes that would occur from neutralino
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
 
c  continue with reading in SLHA files, or quit program
 2000 if (modtype.eq.1) then
        modtype = 2
        close (11)
        goto 500
      endif
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
c      
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


      subroutine read_model_mssm(lunit,nmodel,iend,ierr)
c
c     To read in a model from a file
c
c     NOTE: We here show how you can read in a model file and define all
c     DarkSUSY model parameters. For MSSM models, there exists a routine
c     src/su/dsgive_model.f that makes the definitions for you. 
c     We here read the model parameters and call that routine.
c     Caution: all of the following supersymmetric
c     parameters must be assigned: the gaugino masses m1, m2, m3; the
c     Higgs pseudoscalar mass ma; the ratio of Higgs vacuum expectation
c     values tanbe; the square of the squark and slepton mass parameters
c     mass2q(i), mass2u(i), mass2d(i), mass2l(i), mass2e(i); the
c     trilinear soft parameters asofte(i), asoftu(i), asoftd(i). Here
c     i=1,2,3 is a generation index. In the current version of DarkSUSY,
c     all these parameters are at the weak scale.

      implicit none
      include 'dsidtag.h'
      include 'dsmssm.h'
      integer nmodel,lunit,iend,ierr
      real*8 at,ab,mqtild
      integer i
      character*40 message
 200  format (1x,a12,7(1x,e14.8))
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
     &     idtag,mu,m2,ma,tanbe,mqtild,at,ab
c... modify/set additional parameters
c      higloop=5  ! 5 = Full FeynHiggs;  6 = FeynHiggsFast
      ! W mass for unitarity of tree-level annihilation amplitudes
      call dsgive_model(mu,m2,ma,tanbe,mqtild,at,ab)
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

      subroutine read_model_slha(lunit,nmodel,iend,ierr)
c
c     As read_model_mssm, but now to read in SLHA file names
c 
      implicit none
      include 'dsidtag.h'
      integer nmodel,lunit,iend,ierr
      integer i
      character*40 message
      character*300 slhafile
      include 'dsdir.h'

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
      read (lunit,*,end=1000,err=3000) slhafile
c... make sure to pass absolute path to dsgive_model_SLHA   
      slhafile=dsdatapath(1:index(dsdatapath,'data/')-1)//
     &    'examples/test/models/slha2/'//slhafile(1:index(slhafile,' ')-1)
      call dsgive_model_SLHA(slhafile,0) ! 0=no warnings, 1=print warnings
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
          write(*,*) 'Total number of errors in dstest_mssm: ',memory
          write(*,*)
          write(*,*)
c The simpler stop + integer only works for FORTRAN >= 2008...   
c          stop memory
          open(unit=99, file="dstest_mssm.out", status="unknown")
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
     &                   '   (# ',memory,'in dstest_mssm)'
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
