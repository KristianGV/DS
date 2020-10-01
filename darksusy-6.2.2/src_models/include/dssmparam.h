*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                          dssmparam.h                             ***
***  this piece of code contains SM INPUT PARAMETERS and is needed   ***
***  as a separate file that must ONLY be included in dsinit_module, ***
***  i.e. the module-specific initialization routine                 ***
c----------------------------------------------------------------------c
c  author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 04/2014
c  modified: Paolo Gondolo 20160519 added Higgs

      real*8 alphem_def, GFermi_def, alph3mz_def
      real*8 mass_Z_def, mass_W_def, mass_nue_def, mass_e_def, mass_numu_def,
     &       mass_mu_def, mass_nutau_def, mass_tau_def, mass_up_def, 
     &       mass_down_def, mass_charm_def, mass_strange_def, mass_top_def,
     &       mass_bottom_def, mass_higgs_def
      real*8 width_kz_def, width_kw_def, width_kt_def, width_higgs_def
      real*8 ckms12_def, ckms23_def, ckms13_def, ckmdelta_def

c
c    standard model default parameters
c
      parameter(
c... gauge coupling constants at the z scale
     & alphem_def  = 1d0/127.918d0, ! 0.0078186083d0
     & GFermi_def  = 1.16637d-5, ! GeV^-2
     & alph3mz_def = 0.1172d0, ! PDG 2002 value for alpha strong at the z scale
c
c... standard model masses (pole mass unless stated otherwise)
     & mass_Z_def       = 91.1876d0,
     & mass_W_def       = 80.33d0, ! don't use, will give unitarity problem
     & mass_nue_def     =  0.0d0,
     & mass_e_def       =  0.000510999907d0,
     & mass_numu_def    =  0.0d0,
     & mass_mu_def      =  0.105658389d0,
     & mass_nutau_def   =  0.0d0,
     & mass_tau_def     =  1.777d0,
     & mass_up_def      =  0.003d0,! NB: ms(2GeV) rather than pole mass
     & mass_down_def    =  0.006d0, ! NB: ms(2GeV) rather than pole mass
     & mass_charm_def   =  1.26d0,  ! NB: mc(mc) rather than  pole mass
     & mass_strange_def =  0.103d0, ! NB: ms(2GeV) rather than pole mass
     & mass_top_def     =  172.7d0,    
     & mass_bottom_def  =  4.2d0,   ! NB: mb(mb) rather than  pole mass
     & mass_higgs_def   =  125.09d0,
c
c... standard model widths
     & width_kz_def    = 2.490d0,
     & width_kw_def    = 2.07d0,
     & width_kt_def    = 2.0d0,
     & width_higgs_def = 4.21d-3,
c
c... cabibbo-kobayashi-maskawa mixing matrix
     & ckms12_def     = 0.221d0,     ! sin(theta_12)
     & ckms23_def     = 0.040d0,     ! sin(theta_23)
     & ckms13_def     = 0.0035d0,    ! sin(theta_13)
     & ckmdelta_def   = 0.d0         ! delta no cp violation
     &)


***                                                                 ***
************************ end of dsmssm.h ******************************
