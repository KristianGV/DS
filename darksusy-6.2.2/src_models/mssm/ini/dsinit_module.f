*******************************************************************************
***  This is the initialization subroutine for the "MSSM" module.           ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  Last major edit by: Torsten Bringmann (torsten.bringmann@fys.uio.no)   ***
***  Date: 11/2016                                                          ***
*******************************************************************************
      subroutine dsinit_module
      implicit none
c... the following should be included by all dsinit_module versions: 
      include 'dsmpconst.h'
      include 'dssmparam.h'
      include 'dsidtag.h'
      include 'dsio.h'
c... this is the header file for particle-specific common blocks
c... (always contains dsparticles.h) 
      include 'dsmssm.h'
c... model-specific common blocks
      include 'dsandwcom.h'
      include 'dsaccom.h'
      include 'dsascom.h'
      include 'dssem_sun.h'

      integer i

      data knu,               kl,          kqu,      kqd      /
     &     knue,knumu,knutau, ke,kmu,ktau, ku,kc,kt, kd,ks,kb /     


c... naming the module is important for internal consistency checks!
      moduletag='MSSM'
   
c... We now initialize the SM physics in the "standard" way by filling
c... in particular the first 17 entries of the common block variables
c... in dsparticles.h. Note that this is a convenient, but not required way
c... of including the SM (it could be done in a completely user-defined way).
c... Additionally, any particular setting can be overridden by simply 
c... re-assigning the corresponding values after this call
      call dsinit_sm
      if (prtlevel.gt.2) write(*,*) 'AA', s2thw, sinthw, costhw

c... set electron mass in dsmpconst.h -- this is always needed by the CR routines
      m_e=mass_e_def

c... cabibbo-kobayashi-maskawa mixing matrix
      ckmdelta=0.d0       ! NB: module MSSM cannot handle cp violation!

c
c... internal fixed-for-ever values go here
c
      data kn,              kcha         /
     &     kn1,kn2,kn3,kn4, kcha1,kcha2  /
      data ksnu                          /
     &     ksnu1,ksnu2,ksnu3,0,0,0       /
      data ksl                           /
     &     ksl1,ksl3,ksl5,ksl2,ksl4,ksl6 /
      data ksqu                          /
     &     ksu1,ksu3,ksu5,ksu2,ksu4,ksu6 /
      data ksqd                          /
     &     ksd1,ksd3,ksd5,ksd2,ksd4,ksd6 /

c... Provide an *initial* assignment of flavour states, corresponding to 
c... the old convention.
c... This is re-calculated in dsmodelsetup.
      ksnu_flav(1,1)=ksnu1
      ksnu_flav(2,1)=ksnu2
      ksnu_flav(3,1)=ksnu3
      ksl_flav(1,1)=ksl1
      ksl_flav(2,1)=ksl3
      ksl_flav(3,1)=ksl5
      ksl_flav(1,2)=ksl2
      ksl_flav(2,2)=ksl4
      ksl_flav(3,2)=ksl6
      ksqu_flav(1,1)=ksu1
      ksqu_flav(2,1)=ksu3
      ksqu_flav(3,1)=ksu5
      ksqu_flav(1,2)=ksu2
      ksqu_flav(2,2)=ksu4
      ksqu_flav(3,2)=ksu6
      ksqd_flav(1,1)=ksd1
      ksqd_flav(2,1)=ksd3
      ksqd_flav(3,1)=ksd5
      ksqd_flav(1,2)=ksd2
      ksqd_flav(2,2)=ksd4
      ksqd_flav(3,2)=ksd6


c...Particle names (SM names are set in dsinit_sm)
      pname(17)='H1'
      pname(18)='H2'
      pname(19)='A'
      pname(20)='H+'
      pname(21)='s-nu_1'
      pname(22)='s-l_1'
      pname(23)='s-l_2'
      pname(24)='s-nu_2'
      pname(25)='s-l_3'
      pname(26)='s-l_4'
      pname(27)='s-nu_3'
      pname(28)='s-l_5'
      pname(29)='s-l_6'
      pname(30)='s-qu_1'
      pname(31)='s-qu_2'
      pname(32)='s-qd_1'
      pname(33)='s-qd_2'
      pname(34)='s-qu_3'
      pname(35)='s-qu_4'
      pname(36)='s-qd_3'
      pname(37)='s-qd_4'
      pname(38)='s_qu_5'
      pname(39)='s_qu_6'
      pname(40)='s-qd_5'
      pname(41)='s-qd_6'
      pname(42)='x0_1'
      pname(43)='x0_2'
      pname(44)='x0_3'
      pname(45)='x0_4'
      pname(46)='x+_1'
      pname(47)='x+_2'
      pname(48)='gluino'
      pname(49)='goldst0'
      pname(50)='goldst+'

c...Particle spins (SM spins are set in dsinit_sm)
      spin(17)=0.d0
      spin(18)=0.d0
      spin(19)=0.d0
      spin(20)=0.d0
      spin(21)=0.d0
      spin(22)=0.d0
      spin(23)=0.d0
      spin(24)=0.d0
      spin(25)=0.d0
      spin(26)=0.d0
      spin(27)=0.d0
      spin(28)=0.d0
      spin(29)=0.d0
      spin(30)=0.d0
      spin(31)=0.d0
      spin(32)=0.d0
      spin(33)=0.d0
      spin(34)=0.d0
      spin(35)=0.d0
      spin(36)=0.d0
      spin(37)=0.d0
      spin(38)=0.d0
      spin(39)=0.d0
      spin(40)=0.d0
      spin(41)=0.d0
      spin(42)=0.5d0
      spin(43)=0.5d0
      spin(44)=0.5d0
      spin(45)=0.5d0
      spin(46)=0.5d0
      spin(47)=0.5d0
      spin(48)=0.5d0
      spin(49)=0.d0
      spin(50)=0.d0

c...Degrees of freedom (SM values set in dsinit)
      kdof(17)=1
      kdof(18)=1
      kdof(19)=1
      kdof(20)=2
      kdof(21)=2
      kdof(22)=2
      kdof(23)=2
      kdof(24)=2
      kdof(25)=2
      kdof(26)=2
      kdof(27)=2
      kdof(28)=2
      kdof(29)=2
      kdof(30)=6
      kdof(31)=6
      kdof(32)=6
      kdof(33)=6
      kdof(34)=6
      kdof(35)=6
      kdof(36)=6
      kdof(37)=6
      kdof(38)=6
      kdof(39)=6
      kdof(40)=6
      kdof(41)=6
      kdof(42)=2
      kdof(43)=2
      kdof(44)=2
      kdof(45)=2
      kdof(46)=4
      kdof(47)=4
      kdof(48)=16
      kdof(49)=1
      kdof(50)=2

c...Quantum numbers
      do i=17,50
         wiso3(i)=0.d0
      enddo
      echarg(17)=0.d0
      echarg(18)=0.d0
      echarg(19)=0.d0
      echarg(20)=1.d0
      echarg(21)=0.d0
      echarg(22)=1.d0
      echarg(23)=1.d0
      echarg(24)=0.d0
      echarg(25)=1.d0
      echarg(26)=1.d0
      echarg(27)=0.d0
      echarg(28)=1.d0
      echarg(29)=1.d0
      echarg(30)=2.d0/3.d0
      echarg(31)=2.d0/3.d0
      echarg(32)=-1.d0/3.d0
      echarg(33)=-1.d0/3.d0
      echarg(34)=2.d0/3.d0
      echarg(35)=2.d0/3.d0
      echarg(36)=-1.d0/3.d0
      echarg(37)=-1.d0/3.d0
      echarg(38)=2.d0/3.d0
      echarg(39)=2.d0/3.d0
      echarg(40)=-1.d0/3.d0
      echarg(41)=-1.d0/3.d0
      echarg(42)=0.d0
      echarg(43)=0.d0
      echarg(44)=0.d0
      echarg(45)=0.d0
      echarg(46)=1.d0
      echarg(47)=1.d0
      echarg(48)=0.d0
      echarg(49)=0.d0
      echarg(50)=1.d0
      ncolor(17)=1.d0
      ncolor(18)=1.d0
      ncolor(19)=1.d0
      ncolor(20)=1.d0
      ncolor(21)=1.d0
      ncolor(22)=1.d0
      ncolor(23)=1.d0
      ncolor(24)=1.d0
      ncolor(25)=1.d0
      ncolor(26)=1.d0
      ncolor(27)=1.d0
      ncolor(28)=1.d0
      ncolor(29)=1.d0
      ncolor(30)=3.d0
      ncolor(31)=3.d0
      ncolor(32)=3.d0
      ncolor(33)=3.d0
      ncolor(34)=3.d0
      ncolor(35)=3.d0
      ncolor(36)=3.d0
      ncolor(37)=3.d0
      ncolor(38)=3.d0
      ncolor(39)=3.d0
      ncolor(40)=3.d0
      ncolor(41)=3.d0
      ncolor(42)=1.d0
      ncolor(43)=1.d0
      ncolor(44)=1.d0
      ncolor(45)=1.d0
      ncolor(46)=1.d0
      ncolor(47)=1.d0
      ncolor(48)=8.d0
      ncolor(49)=1.d0
      ncolor(50)=1.d0

c
c... default values go here
c


c... program switches
c      higloop = 3        ! Carena-Espinosa-Quiros-Wagner
      higloop = 5         ! FeynHiggs
      higwid = 5          ! Widhts from FeynHiggs as well
c      higloop = 6        ! FeynHiggsFast
      neuloop = 1
      bsgqcd = 1
      msquarks = 0.d0    ! if 0.d0, use full mass matrices (see manual)
      msleptons = 0.d0  

c... handle neutralino-neutralino annihilation as in 1510.02473
      data NLOoption/'default'/  
c     data NLOoption/'off'/       

c... supersymmetric parameters
      mu=0.d0
      m1=0.d0
      m2=0.d0
      m3=0.d0
      ma=0.d0
      tanbe=0.d0
      do i=1,3
         mass2q(i)=0.d0
         mass2u(i)=0.d0
         mass2d(i)=0.d0
         mass2l(i)=0.d0
         mass2e(i)=0.d0
         asofte(i)=0.d0
         asoftu(i)=0.d0
         asoftd(i)=0.d0
      enddo


c... added by Piero Ullio 
c... set-up fermion family codes needed for coannihilations 
c... numbers are conventional, no physical meaning
*****   itype(iii)=ivfam(ku) for up-type (s)quark
*****   itype(iii)=ivfam(kd) for down-type (s)quark
*****   itype(iii)=ivfam(iii) for (s)leptons
      do i=0,50
        ivtype(i)=0
      enddo
      ivtype(knue)=11
      ivtype(ke)=12
      ivtype(ksnu1)=11
      ivtype(ksl1)=12
      ivtype(ksl2)=12
      ivtype(knumu)=11 ! JE CORRECTION
      ivtype(kmu)=12
      ivtype(ksnu2)=11
      ivtype(ksl3)=12
      ivtype(ksl4)=12
      ivtype(knutau)=11
      ivtype(ktau)=12
      ivtype(ksnu3)=11
      ivtype(ksl5)=12
      ivtype(ksl6)=12

      ivtype(ku)=41
      ivtype(kd)=42
      ivtype(ksu1)=41
      ivtype(ksu2)=41
      ivtype(ksd1)=42
      ivtype(ksd2)=42
      ivtype(kc)=41
      ivtype(ks)=42
      ivtype(ksu3)=41
      ivtype(ksu4)=41
      ivtype(ksd3)=42
      ivtype(ksd4)=42
      ivtype(kt)=41
      ivtype(kb)=42
      ivtype(ksu5)=41
      ivtype(ksu6)=41
      ivtype(ksd5)=42
      ivtype(ksd6)=42

c... added by Piero Ullio 
c... label to get a warning statement in case of negative rate in
c... one of the partial results in dsasdwdcossfsf or dsasdwdcossfchi
c... negative results are internally reset to zero
c... no warning is printed in case aszeroprint is set to false
      aszeroprint=.false.

c... added by Erik Lundstrom (090316)
c... initialize HiggsBounds if HiggsBounds is to be used for
c... calculation of higgs boson accelerator bounds 
      if (higwid.eq.5) then

         call initialize_HiggsBounds(3,1,'LandH')
         call initialize_HiggsSignals(3,1,'latestresults')
c...HiggsSignals options         
         call setup_output_level(0) ! 0=silent,1=screen,2=more
         call setup_pdf(2)      ! 1=box, 2=gauss, 3=box+gauss
         call setup_mcmethod_dm_theory(1) ! 1=mass vari, 2=smearing

      endif

c...initialize yield routines
      call dsaninit_mssm

c... set-up defaults for modules
      call dsddset_mssm('default')
      call dsacset_mssm('default')
      call dsanset_mssm('default')


      write(*,*) 'Initialization of particle physics module MSSM complete.'

      return
      end


