c This is a modified version of the HBandHSwithFH.F file supplied with
c HiggsSignals 1.4.0.
c The routine runs FeynHiggs to evaluate the partonic input of
c HiggsBounds/HiggsSignals for MSSM.	

c This file assumes that some FeyhHiggs subroutines have been
c called previously in dsfeynhiggs. They are marked below with 
c 'done in dsfeynhiggs'.

* Updated for FeynHiggs 2.11.3 by JE 2016-04-11
* Updated for FeynHiggs 2.13.0 by JE 2017-12-19	
* Modified for DarkSUSY by JE 2016-06-28	


      subroutine dshiggsbounds(HBresult,pvalue)
      implicit none

      include 'dsmpconst.h' ! after removing Pi below
      include 'dsio.h'

      integer error

c used by FHHiggsCorr
      double precision MHiggs(4)
      double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHCouplings:
#include "FHCouplings.h"
      double complex couplings(ncouplings), couplingsms(ncouplingsms)
      double precision gammas(ngammas), gammasms(ngammasms)
      integer fast

c used by FHHiggsProd:
        double precision sqrts, prodxs(nprodxs)

c used by FHGetPara:
      integer nmfv
      double precision MSf(2,5,3)
      double precision MASf(6,5), MCha(2), MNeu(4)
      double complex USf(2,2,5,3)
      double complex UASf(6,6,5)
      double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
      double complex DeltaMB
      double precision MGl
      double precision MHtree(4), SAtree,AlfasMT

c used by FHRetrieveSMPara:	
        double precision invAlfa0,invAlfaMZ, AlfasMZ, GF
        double precision ME, MU, MD, MM, MC, MS, ML, MB
        double precision MW, MZ, GW, GZ
        double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

c used by initialize_HiggsBounds
        integer nHiggsneut,nHiggsplus
        parameter (nHiggsneut = 3)
        parameter (nHiggsplus = 1)

c used by HiggsBounds_neutral_input_part
        double precision Mh(3),GammaTotal_hj(3)
        integer CP_value(3)
        double precision  CS_lep_hjZ_ratio(3),          
     &   CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),
     &   CS_lep_hjhi_ratio_nHbynH(3,3),               
     &   CS_gg_hj_ratio(3),CS_bb_hj_ratio(3),  
     &   CS_bg_hjb_ratio(3),                       
     &   CS_ud_hjWp_ratio(3),CS_cs_hjWp_ratio(3),
     &   CS_ud_hjWm_ratio(3),CS_cs_hjWm_ratio(3), 
     &   CS_gg_hjZ_ratio(3),
     &   CS_dd_hjZ_ratio(3),CS_uu_hjZ_ratio(3),
     &   CS_ss_hjZ_ratio(3),CS_cc_hjZ_ratio(3), 
     &   CS_bb_hjZ_ratio(3),                        
     &   CS_tev_vbf_ratio(3),CS_tev_tthj_ratio(3),
     &   CS_lhc7_vbf_ratio(3),CS_lhc7_tthj_ratio(3),
     &   CS_lhc8_vbf_ratio(3),CS_lhc8_tthj_ratio(3),
     &   BR_hjss(3),BR_hjcc(3),                         
     &   BR_hjbb(3),BR_hjmumu(3),BR_hjtautau(3),                     
     &   BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),                     
     &   BR_hjgaga(3),BR_hjgg(3),
     &   BR_hjinvisible(3),BR_hjhihi_nHbynH(3,3)

c used by HiggsBounds_charged_input
        double precision Mhplus(1),GammaTotal_Hpj(1),
     &   CS_lep_HpjHmj_ratio(1),                    
     &   BR_tWpb,BR_tHpjb(1),                     
     &   BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1) 

c used in HiggsBounds
        double precision dmhneut_hb(nHiggsneut)
        double precision dmhch_hb(nHiggsplus)
c used in HiggsSignals (can be different)
        double precision dmhneut_hs(nHiggsneut)
        
c return values of run_HiggsBounds
        integer HBresult,chan,ncombined
        double precision obsratio

c return values of run_HiggsSignals
        double precision Chisq_mu, Chisq_mh, Chisq, Pvalue
        integer nobs
        
c run options in HiggsSignals
        integer runmode
             
c misc:
        integer i,j,as,t
        double precision norm
        double precision
     &   g2hjbb(3),g2hjWW(3),g2hjZZ(3),g2hjhiZ_nHbynH(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass 

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

********************************************************************
*** NOTE. Intialization of HiggsBounds and HiggsSignals is carried
*** out in src_models/mssm/ini/dsinit_module.f
********************************************************************


*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c calls to FH subroutines:

c the subroutines setFlags,setPara,setSLHA
c are also contained in this file

c done in dsfeynhiggs
c	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
c done in dsfeynhiggs
c	call setPara
cc	call setSLHA("tmp.slha")

c...JE update for FeynHiggs 2.11.3: Msf and Usf added to arguments
c...JE update for FeynHiggs 2.13.0: AlfasMT added to arguments	
      call FHGetPara(error, nmfv, MSf,USf,MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &    MHtree, SAtree,AlfasMT)
      if( error .ne. 0 ) stop

c done in dsfeynhiggs, we do it again to get the results
        call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
        if( error .ne. 0 ) stop

c done in dsfeynhiggs
c NOTE: we are setting uzint=uzext
c        mfeff=1
c        uzint=2
c        uzext=2
c        call FHSelectUZ(error, uzint, uzext, mfeff)
c	    if( error .ne. 0 ) stop

c done in dsfeynhiggs, we do it again to get the results here
c        fast=1
        fast=0
        call FHCouplings(error,
     &    couplings, couplingsms, gammas, gammasms, fast)
        if( error .ne. 0 ) stop

c We would like FH to calculate LHC cross sections
        sqrts=8.0D0
        call FHHiggsProd(error, sqrts, prodxs)
        if( error .ne. 0 ) stop 

        call FHRetrieveSMPara(error,
     &    invAlfa0, invAlfaMZ,AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MB,
     &    MW, MZ, GW, GZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      
!--------------- Turn FeynHiggs results into HiggsBounds/Signals input ------------------!

c Set variables needed by HiggsBounds (using results from FeynHiggs).
c See HiggsBounds documentation for definition of variables used
c as arguments to HiggsBounds_neutral_input_part and run_HiggsBounds
c and FeynHiggs documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine HiggsBounds_neutral_input_part
c rather than the subroutine HiggsBounds_neutral_input_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in FeynHiggs

        do i=1,3
         Mh(i)=MHiggs(i) 
         GammaTotal_hj(i) = GammaTot(i)

         BR_hjss(i)       = BR(H0FF(i,4,2,2))
         BR_hjcc(i)       = BR(H0FF(i,3,2,2))
         BR_hjbb(i)       = BR(H0FF(i,4,3,3))
         BR_hjmumu(i)     = BR(H0FF(i,2,2,2))
         BR_hjtautau(i)   = BR(H0FF(i,2,3,3)) 

         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))
         BR_hjZga(i)    = BR(H0VV(i,2))
         BR_hjZZ(i)     = BR(H0VV(i,3))
         BR_hjgg(i)     = BR(H0VV(i,5))

         if(GammaSM(H0FF(i,4,3,3)).le.0.0D0)then
          g2hjbb(i)=0.0D0
         else
          g2hjbb(i)=Gamma(H0FF(i,4,3,3)) 
     &      /GammaSM(H0FF(i,4,3,3))
         endif

c Note that this is currently equivalent to
c         g2hjbb(i)= bbh(i)/bbhSM(i)
c         g2hjbb(i)= btagbh(i)/btagbhSM(i)
c as long as MH>80 GeV

         CS_bg_hjb_ratio(i) = g2hjbb(i)
         CS_bb_hj_ratio(i)  = g2hjbb(i)
   
         g2hjbb_s(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))+ 
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0
         g2hjbb_p(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))- 
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0

         g2hjtautau_s(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))+ 
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0
         g2hjtautau_p(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))- 
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0

         if(       g2hjbb_p(i).lt.1.0D-10)then
            CP_value(i) = 1
         elseif(   g2hjbb_s(i).lt.1.0D-10)then
            CP_value(i) = -1
         else
            CP_value(i) = 0
         endif

         CS_lep_bbhj_ratio(i)     = g2hjbb_s(i)+g2hjbb_p(i)
         CS_lep_tautauhj_ratio(i) = g2hjtautau_s(i)+g2hjtautau_p(i)

         g2hjWW(i)= dble(  Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjWW(i)= WhTev(i)/WhTevSM(i
c	  g2hjWW(i)= qqhTev(i)/qqhTevSM(i)
c as long as MH>80 GeV and uzint=uzext

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjZZ(i)= ZhTev(i)/ZhTevSM(i)
c as long as MH>80 GeV and uzint=uzext
c It is also equivalent to g2hjWW(i)

         CS_lep_hjZ_ratio(i)        = g2hjZZ(i)

         CS_gg_hjZ_ratio(i)     = 0.0D0
         CS_dd_hjZ_ratio(i)     = g2hjZZ(i)
         CS_uu_hjZ_ratio(i)     = g2hjZZ(i)
         CS_ss_hjZ_ratio(i)     = g2hjZZ(i)
         CS_cc_hjZ_ratio(i)     = g2hjZZ(i)
         CS_bb_hjZ_ratio(i)     = g2hjZZ(i)

         CS_ud_hjWp_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWp_ratio(i)    = g2hjZZ(i)
         CS_ud_hjWm_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWm_ratio(i)    = g2hjZZ(i)

         CS_tev_vbf_ratio(i)     = g2hjZZ(i)
         CS_lhc7_vbf_ratio(i)    = g2hjZZ(i)
         CS_lhc8_vbf_ratio(i)    = g2hjZZ(i)


         if(tthSM(i).gt.0.0D0)then
           CS_tev_tthj_ratio(i)    = tth(i)/tthSM(i)  
         else
           CS_tev_tthj_ratio(i)    = 0.0D0
         endif

         CS_lhc7_tthj_ratio(i) = CS_tev_tthj_ratio(i)
         CS_lhc8_tthj_ratio(i) = CS_tev_tthj_ratio(i)

c tevatron gluon fusion XS is not calculated in FH is MH<90 geV
         if(Mh(i).gt.90.0001D0)then      
           if(gghSM(i).gt.0.0D0)then
            CS_gg_hj_ratio(i) = ggh(i)/gghSM(i)
           else
            CS_gg_hj_ratio(i) = 0.0D0
           endif 
         else
           if(GammaSM(H0VV(i,5)).le.0.0D0)then
            CS_gg_hj_ratio(i)=0.0D0
           else
            CS_gg_hj_ratio(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
          endif
         endif

        enddo
  
        norm=GF*sqrt(2.0D0)*MZ**2.0D0

        do j=1,3
          do i=1,3   
            g2hjhiZ_nHbynH(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm
            
            CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)

            BR_hjhihi_nHbynH(j,i)=BR(H0HH(j,i,i))
          enddo
        enddo

c higgs->neutralino1 neutralino1 contributes the invisible Higgs decay width
c when neutralino1 or sneutrino is the LSP

        do i=1,3  
         sneutrino_lspcandidate_number=0
         invisible_lsp=.True.

c first determine whether lightest sneutrino is lighter than the lightest neutralino
c
c sneutrino_lspcandidate_number=0 indicates that lightest neutralino is 
c lighter than all the sneutrinos
         lspcandidate_mass=MNeu(1)
         do as=1,3
            if( MASf(as,1) .lt. lspcandidate_mass )then
               lspcandidate_mass=MASf(as,1)
               sneutrino_lspcandidate_number=as
            endif
         enddo

         if(     MCha(1) .lt. lspcandidate_mass )then 
           invisible_lsp=.False.
         elseif( MGl     .lt. lspcandidate_mass )then
           invisible_lsp=.False.
         else
          do as=1,6
           do t=2,4
            if( MASf(as,t) .lt. lspcandidate_mass )then
               invisible_lsp=.False.
            endif
           enddo
          enddo      
         endif

         if(invisible_lsp)then
          if(sneutrino_lspcandidate_number.eq.0)then
            BR_hjinvisible(i) = BR(H0NeuNeu(i,1,1))
          else
            BR_hjinvisible(i) = BR(H0SfSf(i,1,1,1,as))
          endif
         else
          BR_hjinvisible(i) = 0.0D0
         endif
        enddo

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c	    Charged Higgs input

        Mhplus(1)              = MHiggs(4) 
        GammaTotal_Hpj(1)      = GammaTot(4)  
        CS_lep_HpjHmj_ratio(1) = 1.0D0
        BR_tWpb                = BR( tBF(1) ) 
        BR_tHpjb(1)            = BR( tBF(2) )
        BR_Hpjcs(1)            = BR( HpFF(2,2,2) ) 
        BR_Hpjcb(1)            = BR( HpFF(2,2,3) )
        BR_Hpjtaunu(1)         = BR( HpFF(1,3,3) )

	if (prtlevel.ge.2) then
	   print*,'*** Now checking model with',
     &         ' HiggsBounds/HiggsSignals'
		print*,'The Higgs boson masses are (h, H, A, H+):'
		print*,Mh,Mhplus
        endif
	     
!--------------------------------- Run HiggsBounds --------------------------------------!

c calls to HiggsBounds_neutral_input_part,HiggsBounds_charged_input,set_mass_uncertainties
c which give input to HiggsBounds

c        print*,'calling HiggsBounds_neutral_input_part'  
        call HiggsBounds_neutral_input_part(Mh,GammaTotal_hj,CP_value, 
     &          CS_lep_hjZ_ratio,                            
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,     
     &          CS_lep_hjhi_ratio_nHbynH,                    
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,       
     &          CS_bg_hjb_ratio,                         
     &          CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,    
     &          CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,    
     &          CS_gg_hjZ_ratio,     
     &          CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,     
     &          CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,     
     &          CS_bb_hjZ_ratio,                         
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,    
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,    
     &          CS_lhc8_vbf_ratio,CS_lhc8_tthj_ratio,    
     &          BR_hjss,BR_hjcc,                             
     &          BR_hjbb,BR_hjmumu,BR_hjtautau,               
     &          BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg,  
     &          BR_hjinvisible,BR_hjhihi_nHbynH              )
     
c        print*,'calling HiggsBounds_charged_input'
        call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, 
     &          CS_lep_HpjHmj_ratio,                        
     &          BR_tWpb,BR_tHpjb,                           
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu)

c Set the theory mass uncertainties. HiggsBounds will vary the masses accordingly.
		dmhneut_hb = (/2.D0, 0.D0, 0.D0/)
		dmhch_hb   = (/0.D0/)
        call HiggsBounds_set_mass_uncertainties(dmhneut_hb,dmhch_hb)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds
        call run_HiggsBounds( HBresult,chan,obsratio,ncombined)

      if (prtlevel.ge.2) then
        print*,' '      
        print*,'*************    HiggsBounds Results  **************'
        print*,' '    
        print*,'Is this parameter point excluded by LEP, Tevatron'
        print*,'or LHC data?'         
        print*, HBresult, ',  where'
        print*,'               0 = yes, it is excluded'
        print*,'               1 = no, it has not been excluded'
        print*,'              -1 = invalid parameter set'    
        print*,' '    
        print*,'The process with the highest statistical sensitivity'
        print*,'is'
        print*, chan,'(see Key.dat)'
        print*,'This process has a theoretical rate vs. limit of'
        print*, obsratio
        print*,' '    
        print*,'The number of Higgs which have contributed to the'
        print*,'theoretical rate of this process was'
        print*, ncombined
        print*,' '    
        print*,'See HiggsBounds documentation for more information.'
        print*,'****************************************************'
        print*,' '
       endif

!--------------------------------- Run HiggsSignals -------------------------------------!

c Set the HiggsSignals runmode (1: peak-centered, 2: mass-centered, 3: both) Chi^2 method(s)
      runmode = 1

c...Number of model parameters
      call setup_nparam(0) ! consider our input model as fixed
	
c Set the theory mass uncertainties. This will affect the Chi^2 contribution
c from the mass measurements as well as the Higgs-to-peaks assignment in
c HiggsSignals .
      dmhneut_hs = (/2.D0, 0.D0, 0.D0/)
      call HiggsSignals_neutral_input_MassUncertainty(dmhneut_hs)
c Set the MSSM rate uncertainties:
!---dCS and dBR hold the model's rate uncertainties. Can be changed by user
!---with subroutine setup_rate_uncertainties. First, set to SM values.
C		dCS = (/ 0.147D0, 0.028D0, 0.037D0, 0.051D0, 0.12D0 /)  
C		dBR = (/ 0.054D0, 0.048D0, 0.048D0, 0.061D0, 0.028D0 /)
c Now, we set relative uncertainties for the processes gg->H and bb->H
c and interpolate from the effective couplings/CS-ratios (of the lightest Higgs)
c the uncertainty of single Higgs production, dCS(1).
C		dggh = 0.147D0
C		dbbh = 0.200D0		
C		call get_singleH_uncertainty(dCS(1),dggh,dbbh,CS_gg_hj_ratio(1),
C     &                                 g2hjbb_s(1)+g2hjbb_p(1),Mh(1))
C        print*,"Interpolated singleH uncertainty (in %): ",dCS(1)
C		call setup_rate_uncertainties(dCS, dBR)
!---- Run HiggsSignals																 ----!
        call run_HiggsSignals(runmode,Chisq_mu,Chisq_mh,
     &                        Chisq,nobs,Pvalue)
        if (prtlevel.ge.2) then
        print*,' '      
        print*,'*************  HiggsSignals Results  **************'
        print*,' '
        print*,'The chi^2 contr. from signal strength observables is'
        print*,Chisq_mu
        print*,'The chi^2 contr. from Higgs mass observables is'
        print*,Chisq_mh
        print*,'The total chi^2 is'
        print*,Chisq
        print*,'The number of observables is'
        print*,nobs
        print*,' '
        print*,'In the case that the number of observables is equal'
        print*,'to the number of degrees of freedom, this corresponds'
        print*,'to a probability of'
        print*,Pvalue
        print*,' '    
        print*,'See HiggsSignals documentation for more information.'
        print*,'****************************************************'
        print*,' '    
        endif
!---------------------- Finish HiggsBounds and HiggsSignals -----------------------------!

c deallocates arrays used by HiggsBounds:
c do we need to call this at the end of DarkSUSY?
c        call finish_HiggsBounds
c        call finish_HiggsSignals

      end
