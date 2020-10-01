      subroutine dsmodelsetup(unphys,warning)
c_______________________________________________________________________
c  set up global variables for the supersymmetric model routines.
c  uses sconst, sfesct, chasct, higsct,
c    neusct, vertx, hwidths.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c  Torsten Bringmann: unified with SLHA, earlier dssusy(_isasugra)  
c    (08/05/14)       and dsmodelsetup_isasugra
c    (24/05/19)       moved ivfam here to gurantee consistent flavour ordering
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'

      integer unphys,warning, valid, i

c-----------------------------------------------------------------------

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dsmodelsetup')

      call dsnewidnumber ! This should always be the first call in dsmodelsetup,
                         ! and makes sure the model has a new unique ID number.
                         ! Several routines in src_models/ and in src/ 
                         ! *require* this to be set anew for every new model.
      
      unphys=0
      warning=0
      mass(0)=1.d10

 
c... jump over (almost) everything if we read in an SLHA file with low-energy parameters
      if (modeltype.eq.0) goto 500 ! low energy SLHA file read

      if (modeltype.eq.3) then   ! cMSSM model
        call dsrge_isasugra(unphys,valid) !run RGEs and set masses and mixing
        if (valid.gt.0.or.unphys.ne.0) return 
      endif

      if (tanbe.eq.0.0d0) then
         write (*,*)
     &  'dsmodelsetup: susy parameters are not properly initialized'
         write (*,'(a,a)')
     &  '        pls set tanbe,ma,mu,m2,m1,m3,mass2u,mass2q,mass2d,',
     &        'mass2l,mass2e,asoftu,asoftd,asofte'
         stop
      endif

      if (modeltype.ne.3) then   ! pMSSM
        call dssmconst_couplings  ! first set global constants
        call dssmconst_ckm
        call dssuconst_yukawa
        call dssuconst_higgs
        unphys=0
        call dsspectrum(unphys,warning) ! particle spectrum and mixing
        if (unphys.ne.0.or.warning.ne.0) return      
      endif ! pMSSM
      
                        
c------------- reset running things

      call dssuconst_yukawa_running

c-------------------------------------------------- some useful vertices

      call dsvertx

c--------------------------------------------------- and particle widths

c...Add Higgs widths and QCD correction to them and to vertices

      call dshigwid

c...Add sparticle widths

      call dsspwid

 500  continue
 
c... Prepare for rate calculations
      call dsanyieldset

      kdm=kn(1)  ! In principle, one may allow for other LKPs...      

c... set-up fermion family codes needed for coannihilations 
c... (numbers are conventional, no physical meaning)
c... -> this depends on the sfermion flavour and should eventually be made obsolete
c... (once all combinations of in and out states are taken into account)
      do i=0,50
        ivfam(i)=0
      enddo
      ivfam(knue)=11
      ivfam(ke)=12
      ivfam(ksnu_flav(1,1))=11
      ivfam(ksl_flav(1,1))=12
      ivfam(ksl_flav(1,2))=12
      ivfam(knumu)=21
      ivfam(kmu)=22
      ivfam(ksnu_flav(2,1))=21
      ivfam(ksl_flav(2,1))=22
      ivfam(ksl_flav(2,2))=22
      ivfam(knutau)=31
      ivfam(ktau)=32
      ivfam(ksnu_flav(3,1))=31
      ivfam(ksl_flav(3,1))=32
      ivfam(ksl_flav(3,2))=32
      ivfam(ku)=41
      ivfam(kd)=42
      ivfam(ksqu_flav(1,1))=41
      ivfam(ksqu_flav(1,2))=41
      ivfam(ksqd_flav(1,1))=42
      ivfam(ksqd_flav(1,2))=42
      ivfam(kc)=51
      ivfam(ks)=52
      ivfam(ksqu_flav(2,1))=51
      ivfam(ksqu_flav(2,2))=51
      ivfam(ksqd_flav(2,1))=52
      ivfam(ksqd_flav(2,2))=52
      ivfam(kt)=61
      ivfam(kb)=62
      ivfam(ksqu_flav(3,1))=61
      ivfam(ksqu_flav(3,2))=61
      ivfam(ksqd_flav(3,1))=62
      ivfam(ksqd_flav(3,2))=62
 
      return
      end


