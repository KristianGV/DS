      subroutine dsddhelp
c...
c...  type : commonly used
c...  desc : help with options for scattering cross sections
c...author: paolo gondolo 2000-07-07
c...modified by Gintaras Duda 2007-06-27 for new FF options
c...modified by Paolo Gondolo 2008-02-18
c...modified by Torsten Bringmann 2014-12-10: removed model-specific part
c...modified by Paolo Gondolo 2016-02-06: completely removed mssm part
c...modified by Paolo Gondolo 2016-11-20: split from dsddset
      implicit none

      write (*,*) 'dsddhelp: use "call dsddset(key,value)" where'
      write (*,*) ' key=''sf'' set all structure functions to specified value'
      write (*,*) ' key=''sf_m'' set F_M structure function'
      write (*,*) ' value=''default'' - default structure function (''best'')'
      write (*,*) '      =''q=0'' - structure function at q=0'
      write (*,*) '      =''best'' - first available in the following order'
      write (*,*) '      =''FB'' - Fourier-Bessel form factor'
      write (*,*) '      =''SOG'' - sum over gaussians form factor'
      write (*,*) '      =''Fermi'' - Helm form factor with parameters from a Fermi function'
      write (*,*) '      =''L-S'' - Helm form factor with best-fit Lewin-Smith parameters'
      write (*,*) '      =''haxton'' (not available) - Fitzpatrick et al. harmonic oscillator form factor'
      write (*,*) '      =''ds4.1'' (obsolete) - Helm form factor as in DarkSUSY 4.1'
      write (*,*) '      =''gauss'' (obsolete) - Gaussian form factor'
      write (*,*) '      =''gould'' (obsolete) - Gould''s exponential form factor'
      write (*,*) ' key=''sf_sigma'' set F_Sigma and F_Sigma'' structure functions'
      write (*,*) ' value=''default'' - default structure function (''best'')'
      write (*,*) '      =''q=0'' - structure function at q=0'
      write (*,*) '      =''best'' - first available in the following order'
      write (*,*) '      =''ISM'' - interacting shell model spin structure functions'
      write (*,*) '      =''haxton'' (not available) - Fitzpatrick et al. harmonic oscillator form factor'
      write (*,*) '      =''OddG'' - odd group model spin structure functions'
      write (*,*) '      =''SPSM'' - single particle shell model spin structure functions'
      write (*,*) '      =''SIMS'' - simplified approximate shell model for odd-odd nuclei'
      write (*,*) '      =''ISMR'' (obsolete) - Ge-73 Ressel+ 1993'
      write (*,*) ' key=''sf_delta'' set F_Delta structure function'
      write (*,*) ' value=''default'' - default structure function (''best'')'
      write (*,*) '      =''q=0'' - structure function at q=0'
      write (*,*) '      =''best'' - first available in the following order'
      write (*,*) '      =''haxton'' (not available) - Fitzpatrick et al. harmonic oscillator form factor'
      write (*,*) ' key=''sf_phi'' set F_Phi'' and F_Phi'''' structure functions'
      write (*,*) ' value=''default'' - default structure function (''best'')'
      write (*,*) '      =''q=0'' - structure function at q=0'
      write (*,*) '      =''best'' - first available in the following order'
      write (*,*) '      =''haxton'' (not available) - Fitzpatrick et al. harmonic oscillator form factor'
      write (*,*) ' key=''sf_deltasigma'' set F_DeltaSigma structure function'
      write (*,*) ' value=''default'' - default structure function (''best'')'
      write (*,*) '      =''q=0'' - structure function at q=0'
      write (*,*) '      =''best'' - first available in the following order'
      write (*,*) '      =''haxton'' (not available) - Fitzpatrick et al. harmonic oscillator form factor'
      write (*,*) ' key=''sf_phim'' set F_PhiM structure function'
      write (*,*) ' value=''default'' - default structure function (''best'')'
      write (*,*) '      =''q=0'' - structure function at q=0'
      write (*,*) '      =''best'' - first available in the following order'
      write (*,*) '      =''haxton'' (not available) - Fitzpatrick et al. harmonic oscillator form factor'
      write (*,*) ' key=''me'' set all nucleon matrix elements to specified value'
      write (*,*) ' key=''me_s'' set scalar matrix elements'
      write (*,*) ' value=''default'' - default matrix elements (''gls91'')'
      write (*,*) '      =''gls91'' - Gaisser, Leutwyler & Sainio 1991'
      write (*,*) ' key=''me_a'' set axial matrix elements'
      write (*,*) ' value=''default'' - default matrix elements (''smc'')'
      write (*,*) '      =''smc'' - SMC values'
      write (*,*) '      =''emc'' - EMC values'
      write (*,*) ' key=''me_vm'' set vector-magnetic matrix elements'
      write (*,*) ' value=(not available)'
      write (*,*) ' key=''me_am'' set axial-magnetic matrix elements'
      write (*,*) ' value=(not available)'
      write (*,*) ' key=''me_p'' set pseudoscalar matrix elements'
      write (*,*) ' value=(not available)'
      write (*,*) ' key=''me_t'' set antisymmetric tensor matrix elements'
      write (*,*) ' value=(not available)'
      write (*,*) ' key=''me_t2e'' set twist-2 C-even matrix elements'
      write (*,*) ' value=(not available)'
      write (*,*) ' key=''me_t2o'' set twist-2 C-odd matrix elements'
      write (*,*) ' value=(not available)'

      return
      end
