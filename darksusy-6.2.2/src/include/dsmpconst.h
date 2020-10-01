*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                        dsmpconst.h                               ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c mathematical and physical constants, p. gondolo 2011-11-03
c additions/changes t. bringmann 2013, 2014, 2015, 2016, 2018, 2019
c additons P. Scott 2014-10-20, Cosmology      


c... math
      real*8 pi,zeta(8)
      parameter(
c     &     pi=3.1415926535897932384626433832795029D0, ! pi
     &     pi=4.0d0*datan(1.d0)) ! pi
      data zeta/1D99,1.64493D0,1.20206D0,1.20206D0,1.03693D0,1.0173D0,
     &               1.00835D0,1.00408D0/    ! Riemann Zeta for integer arguments


c... physics
      real*8 alpha_em,m_p,m_n,m_d,m_he,n_avogadro,c_light,mpl,m_e,kpc,year,
     &     gev2cm3s,fermiGeV,gev2cm2,gev3cm2g,atomicmassunit,m_p_amu,m_n_amu,
     &     GNewton,G_Fermi
      
      parameter(
     &     alpha_em= 7.29735257d-3,  ! fine-structure constant
     &     mpl=1.2209D19,  ! Planck mass [GeV]           
     &     m_p_amu=1.00727646688d0, ! proton mass [amu]
     &     m_n_amu=1.0086649156d0, ! neutron mass [amu]
     &     atomicmassunit=0.931494028d0, ! atomic mass unit [GeV/c^2]
     &     m_p=m_p_amu*atomicmassunit, ! proton mass [GeV/c^2]
     &     m_n=m_n_amu*atomicmassunit,   ! neutron mass [GeV/c^2]
c     &     m_p=0.938271998d0, ! proton mass [GeV/c^2]
c     &     m_n=0.9396d0,   ! neutron mass [GeV/c^2]
     &     m_d=1.875612762d0, ! deuteron mass [GeV/c^2]
     &     m_he=3.7273794d0, ! helium nucleus mass [GeV/c^2]
c... this is now set in dsinit_module, 
c... typically to the default value defined in src_models/include/dssmparam.h
c     &     m_e=0.000510999907d0, ! electron mass [GeV/c^2]
     &     n_avogadro=6.022d23,  ! Avogradros number [number / mol ]
     &     c_light=299792.458d0, ! speed of light [km/s]
     &     kpc=3.08567802d0,     ! 1 kpc = 3.08567802d21 cm           
     &     gev2cm3s=0.38937966d-27*c_light*1.d5, ! conversion [GeV^2 cm^3/s]
     &     fermiGeV=1.d0/0.1973269602d0, ! conversion [GeV fm]
     &     gev2cm2=(197.327053d-16)**2, ! conversion [GeV^2 cm^2]
     &     gev3cm2g= 2.1842595d-4, ! conversion [GeV^3 cm^2 g^-1]
     &     year=3.1536d7,        ! conversion [year/s]
     &     GNewton=6.6741d-8,    ! cm^3 g^-1 sec^-2 
     &     G_Fermi=1.1663787d-5  ! Fermi coupling constant [(hbar c)^3 GeV^-2]
     & )   

c...We let the module have the freedom to set m_e instead of using
c...a global parameter value      
      common /mpconst/m_e
      save /mpconst/
      
      

c...Cosmology
      double precision secperyr, cmperpc, GeVperSolarMass, SolarMass, mperkpc 
      double precision omegamh2, zeq, Mhzeq, t_0, H_0, rho_cdm, rho_conh2, kappa
      double precision t_eq, dmfrac, omegacdmh2, omegabh2, keq, M_MW
      real*8 r_earth
      parameter(r_earth=6378.14d3) ! Earth radius in m
      parameter (secperyr = 365.25d0*24.d0*60.d0*60.d0) !seconds in a yr (average)
      parameter (cmperpc = 3.0857d18)            !cm per parsec
      parameter (mperkpc = 10.d0*cmperpc)        !meters per kpc
      parameter (GeVperSolarMass = 1.1157472757094956d57)! GeV per solar mass
      parameter (SolarMass = 1.98892d33)         !g per solar mass
      parameter (M_MW = 0.94d12)                 !Milky Way Mass (in solar masses) 

      !From best fit Planck+WMAP (Planck 2018):
      parameter (t_0 = secperyr * 13.805d9)     !Current age of the universe
      parameter (H_0 = 67.4d0)                   !Current Hubble constant 
      parameter (omegamh2 = 0.1426d0)            !Present matter content of the universe * h^2
      parameter (omegabh2 = 0.02229d0)           !Present baryon content of the universe * h^2
      parameter (omegacdmh2 = 0.1196d0)          !Present CDM content of the universe * h^2

      parameter (zeq = 2.32d4*omegamh2 - 1.d0)   !redshift of matter-radiation equality (Kolb & Turner)
      parameter (Mhzeq = 6.5d15/omegamh2**2)     !Horizon mass at zeq, in solar masses (Josan, Green & Malik 2009, PRD 79:103520)
      parameter (t_eq = secperyr * 59.073d3)     !Age at matter-radiation equality (Wright 2006, PASP 118:1711)
      parameter (keq = 0.07185*omegamh2)         !Mode entering horizon at equality (derived via Friedman Eq, in Mpc^-1)
      parameter (dmfrac = omegacdmh2/(omegacdmh2+omegabh2)) !Percentage of matter that is in DM
      parameter (rho_conh2 =3.d4/(8.d0*pi*GNewton*mperkpc**2))!Critical density today / h^2 (g cm^-3)
      parameter (rho_cdm = rho_conh2*omegacdmh2) !Cosmological CDM density today (g cm^-3)
      parameter (kappa = 16.d0*pi*GNewton*rho_cdm*(2.32d4*omegamh2)**3*t_eq**2/(3.d0*dmfrac)) !Dimensionless compound quantity
      
***                                                                 ***
*********************** end of dsmpconst.h ****************************
