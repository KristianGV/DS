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
      sva=svann
      svb=0.0d0 ! In principle, one can also set the p-wave contribution

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


