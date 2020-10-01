**********************************************************************
*** This file is automatically generated from the file 
*** src_models/generic_wimp/ini/dsgivemodel_generic_wimp.f
*** with the script scr/make_replaceable.pl on Aug 24, 2020.
*** The file is copied as is, but to be of any use you should of
*** course modify this file to your liking. The way the default linking
*** is set up, this file will be linked to before the corresponding
*** DarkSUSY file, meaning that this is the file that will be used,
*** i.e. it will replace the default DarkSUSY one.
*** A few lines of code are added to the executable section below
*** to remind you that you need to change this file. Delete those lines
*** (and these comments) and modify the file to your liking and you
*** are ready to go!
**********************************************************************

*******************************************************************************
***  subroutine dsgivemodel_generic_wimp reads in the parameters to         ***
***  describe a generic WIMP and transfers them to common blocks            ***
***                                                                         ***
***  input:                                                                 ***
***                                                                         ***
***    mgenwimp - WIMP mass (in GeV)                                        ***
***    genselfconj - true for self-conjugated, false for not                ***
***    svann    - constant annihilation cross section (in cm^3/s)           ***
***    pdgann   - PDG code for particle in dominant annihilation channel    ***
***               (e.g. 5 for bbar, 24 for W^+W^-)                          ***
***    SI       - cross section for spin-indep. WIMP-nucleon scattering (pb)***
***                                                                         ***
***   Note that further model details can be added / changed by a           ***
***   subsequent call to dsgivemodel_generic_wimp_opt!                      ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2014-06-20                                                         ***
*******************************************************************************
      subroutine dsgivemodel_generic_wimp(mgenwimp,genselfconj,
     &     svann,pdgann,SI)
      implicit none
      include 'dsgeneric_wimp.h'

      real*8 mgenwimp,svann,SI
      integer pdgann
      logical genselfconj
c-----------------------------------------------------------------------



      mass(kdm)=mgenwimp
      if (genselfconj) then
        selfconj=1
      else
        selfconj=2
      endif    
      kdof(kdm)=1               ! for scalar DM
      spin(kdm)=0.0d0           ! for scalar DM
      
c... write annihilatiom rate to common block
      sva=0.0d0
      svb=svann ! In principle, one can also set the p-wave contribution

c... write annihilation channel to common block
      svch=pdgann
      if (svch.lt.1.or.svch.gt.24.
     &    or.(svch.gt.6.and.svch.lt.11.).or.(svch.gt.16.and.svch.lt.21)) then
        write(*,*) 'ERROR in dsgivemodel_generic_wimp: ',pdgann, 
     &             'is not a valid particle code!'
        stop
      endif

c... write scattering cross section to common block, convert pb->cm^2
c... note that further options of specifying scattering cross sections
c... are provided by a subsequent call to dsgivemodel_generic_wimp_opt
      sigsip = SI*1.d-36 !spin-independent coupling to protons
      sigsin = SI*1.d-36 !spin-independent coupling to neutrons
      sigsdp = 0.0d0     !spin-dependent coupling to protons
      sigsdn = 0.0d0     !spin-dependent coupling to neutrons

      return
      end

