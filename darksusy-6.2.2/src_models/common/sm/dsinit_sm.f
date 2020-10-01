*******************************************************************************
***  This is the recommended way to initialize SM physics, and typically    ***
***  called at the beginning of the dsinit_module subroutine of each        ***
***  particle physics module. It automatically sets basic properties of SM  ***
***  particles, running masses etc. Note that this function is just provided***
***  for convenience -- each particle physics model can either adopt its own***
***  way of representing SM physics, or override any individual settings    ***
***  made in this routine by a corresponding statement in dsinit_module     ***
***  after calling this basic initialization routine.                       ***
***                                                                         ***
***  Author: Torsten Bringmann (torsten.bringmann@fys.uio.no)               ***
***  Date: 2016-11-16                                                       ***
*******************************************************************************
      subroutine dsinit_sm
      implicit none

c... the following should be included by all dsinit_module versions: 
      include 'dssmparam.h'
      include 'dsmpconst.h'
      include 'dsparticles.h'
      include 'dssm.h'   ! note that this must be included *after* dsparticles.h
      include 'dsio.h'
      
      integer i
      real*8 dsgf2s2thw, dsmqpole4loop, mtpole


      write(*,*) 'Initializing native DarkSUSY SM routines...'

      
c... Initialize common block variables in dsparticles.h and dssm.h
      do i=0,maxnumpartspecies
        pdgcode(i) = 0
        pname(i)   = 'empty'
        mass(i)    = 0.d0
        width(i)   = 0.d0
        kdof(i)    = 0
        
        ncolor(i)  = 0.d0
        wiso3      = 0.d0
        echarg     = 0.d0
      enddo
      kdm=0
      
c...Particle names of SM particles
      pname(0)='error'
      pname(1)='nu_e'
      pname(2)='e'
      pname(3)='nu_mu'
      pname(4)='mu'
      pname(5)='nu_tau'
      pname(6)='tau'
      pname(7)='u'
      pname(8)='d'
      pname(9)='c'
      pname(10)='s'
      pname(11)='t'
      pname(12)='b'
      pname(13)='gamma'
      pname(14)='W'
      pname(15)='Z'
      pname(16)='gluon'
      pname(17)='SM_higgs'

c...Standard model masses (as read in from dssmparam.h)
      mass(1)  = mass_nue_def
      mass(2)  = mass_e_def
      mass(3)  = mass_numu_def
      mass(4)  = mass_mu_def
      mass(5)  = mass_nutau_def
      mass(6)  = mass_tau_def
      mass(13) = 0.0d0
      mass(15) = mass_Z_def       ! NB: W mass is *calculated* further down
      mass(16) = 0.0d0
      mass(17) = mass_higgs_def

c... gauge coupling constants at the z scale
      alph3mz = alph3mz_def
      alphem  = alphem_def 
      GFermi  = GFermi_def
      
      
c... determine various combinations of weinberg angle
      mtpole       = mass_top_def   
      s2thw=dsgf2s2thw(GFermi,alphem,mass(kz),mtpole,1)
      sinthw=sqrt(s2thw)
      costhw=sqrt(1.0d0-s2thw)
      mass(kw) = mass(kz)*sqrt(1.d0-s2thw) ! calc W mass from Z mass input
                                           ! (often needed in this way for unitarity of 
                                           ! tree-level annihilation amplitudes

c...When there is no gauge-invariance or unitarity issue, it might be more useful
c... to include higher-order corrections to theta_W. So we also determine it
c... at the MS-bar value at MZ
      s2wmz=dsgf2s2thw(GFermi,alphem,mass(kz),mtpole,3)
      swmz=sqrt(s2wmz)
      cwmz=sqrt(1.0d0-s2wmz)

c... Now add quark masses
c      data roption/'1loop'/
      roption = '4loop'               ! option for how to treat running of
                                      ! quark masses and alpha_s 
      first_dsralph34loop = .true.    ! make sure to compute the 4-loop results only once
      mu2gev       = mass_up_def      ! NB: ms(2GeV) rather than pole mass
      md2gev       = mass_down_def    ! NB: ms(2GeV) rather than pole mass
      mcmc         = mass_charm_def   ! NB: mc(mc) rather than  pole mass
      ms2gev       = mass_strange_def ! NB: ms(2GeV) rather than pole mass
      mbmb         = mass_bottom_def  ! NB: mb(mb) rather than  pole mass
      mass(kt)     =  mtpole
      call dsfindmtmt 
      mass(ku)     =  mu2gev
      mass(kd)     =  md2gev
      mass(ks)     =  ms2gev
      mass(kc)     =  dsmqpole4loop(kc,mcmc)
      mass(kb)     =  dsmqpole4loop(kb,mbmb)



c... standard model widths
      width(11) = width_kt_def
      width(14) = width_kw_def
      width(15) = width_kz_def
      width(17) = width_higgs_def


c... cabibbo-kobayashi-maskawa mixing matrix
      ckms12=ckms12_def   ! sin(theta_12)
      ckms23=ckms23_def   ! sin(theta_23)
      ckms13=ckms13_def   ! sin(theta_13)
      ckmdelta=0.d0

c... calc Higgs vev v0
      v0 = 2.d0*mass(14)*sqrt(s2thw/4.d0/pi/alphem_def)
      if (prtlevel.ge.3) write (*,*) 'PG dsinit_model> old v0=',v0
      v0 = 1.d0/dsqrt(dsqrt(2.d0)*GFermi_def)
      if (prtlevel.ge.3) write (*,*) 'PG dsinit_model> new v0=',v0


c...Spins of SM particles
      spin(0)=0.d0
      spin(1)=0.5d0
      spin(2)=0.5d0
      spin(3)=0.5d0
      spin(4)=0.5d0
      spin(5)=0.5d0
      spin(6)=0.5d0
      spin(7)=0.5d0
      spin(8)=0.5d0
      spin(9)=0.5d0
      spin(10)=0.5d0
      spin(11)=0.5d0
      spin(12)=0.5d0
      spin(13)=1.d0
      spin(14)=1.d0
      spin(15)=1.d0
      spin(16)=1.d0
      spin(17)=0.d0

c...Degrees of freedom
c...internal degrees of freedom of SM particles
      kdof(1)=2
      kdof(2)=4
      kdof(3)=2
      kdof(4)=4
      kdof(5)=2
      kdof(6)=4
      kdof(7)=12
      kdof(8)=12
      kdof(9)=12
      kdof(10)=12
      kdof(11)=12
      kdof(12)=12
      kdof(13)=2
      kdof(14)=6
      kdof(15)=3
      kdof(16)=16
      kdof(17)=1


c... standard model charges
      wiso3(knue)   =+0.5d0
      wiso3(ke)     =-0.5d0
      wiso3(knumu)  =+0.5d0
      wiso3(kmu)    =-0.5d0
      wiso3(knutau) =+0.5d0
      wiso3(ktau)   =-0.5d0
      wiso3(ku)     =+0.5d0
      wiso3(kd)     =-0.5d0
      wiso3(kc)     =+0.5d0
      wiso3(ks)     =-0.5d0
      wiso3(kt)     =+0.5d0
      wiso3(kb)     =-0.5d0
      wiso3(kgamma) =0.d0
      wiso3(kw)     =0.d0
      wiso3(kz)     =0.d0
      wiso3(kgluon) =0.d0
      echarg(knue)  =0.d0
      echarg(ke)    =-1.d0
      echarg(knumu) =0.d0
      echarg(kmu)   =-1.d0
      echarg(knutau)=0.d0
      echarg(ktau)  =-1.d0
      echarg(ku)    =+2.d0/3.d0
      echarg(kd)    =-1.d0/3.d0
      echarg(kc)    =+2.d0/3.d0
      echarg(ks)    =-1.d0/3.d0
      echarg(kt)    =+2.d0/3.d0
      echarg(kb)    =-1.d0/3.d0
      echarg(kgamma)=0.d0
      echarg(kw)    =1.d0
      echarg(kz)    =0.d0
      echarg(kgluon)=0.d0
      ncolor(knue)  =1.d0
      ncolor(ke)    =1.d0
      ncolor(knumu) =1.d0
      ncolor(kmu)   =1.d0
      ncolor(knutau)=1.d0
      ncolor(ktau)  =1.d0
      ncolor(ku)    =3.d0
      ncolor(kd)    =3.d0
      ncolor(kc)    =3.d0
      ncolor(ks)    =3.d0
      ncolor(kt)    =3.d0
      ncolor(kb)    =3.d0
      ncolor(kgamma)=0.d0
      ncolor(kw)    =0.d0
      ncolor(kz)    =0.d0
      ncolor(kgluon)=8.d0


c calculate SM gauge couplings from input parameters
      call dssmconst_couplings
c set CKM mixings
      call dssmconst_ckm


      return
      end


