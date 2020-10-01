**********************************************************************
*** subroutine dsfromSLHA transfers SLHA data in the array slhadata
*** to DarkSUSY. The SLHA things are read with with the help of
*** SLHALIB by Thomas Hahn.
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2008-07-01
*** FIXME : Better error handling
***         More options on how to treat redundant information, e.g.
***         should Yukawas from SLHA file be used instead of letting
***         DarkSUSY run them (as is currently done).
*** Modified: 2016-11-11 Joakim Edsjo, fixed lsp setting    
***           2018-12-07 Torsten Bringmann, added (s)particle widths  
**********************************************************************
      subroutine dsfromSLHA
      implicit none
c...Note: These include files need to be included with the precompiler
c...directive for include and not the fortran one. Hence, this file
c...has to end on .F instead of .f as well.


#include "SLHA.h"

#include "PDG.h"

      include 'dsmssm.h'
      include 'dsidtag.h'
#include "dsslha.h"

      integer plev
      parameter(plev=0) ! 0=only errors, 1=errors and warnings
      integer unphys,hwarning,i,opt
      integer g,n,c,n2,c2,s1,s2,g1,g2, dshiggs(4), PDGhiggs(4)
      real*8 dsgf2s2thw,tmp,dsabsq,dszarg,
     &  dsmqpole4loop, SLHAGetDecay
      complex*16 aux
c...mSUGRA parameters
      real*8 m0,mhf,a0,sgnmu,tgbeta

c...Logistics
      integer ncall
      data ncall/0/
      save ncall

      ncall=ncall+1
      if (ncall.eq.10) then
         write(*,*) 'DS Info: dsfromslha has been called 10 times,',
     &        ' non-fatal errors will no longer be printed.'
      endif
c...Zero things here
      call dsmssmzero

c...Set higwid parameter
c...-1 calculates widths, but does not update couplings, use for MSSM models    
      higwid=-1 ! like 1, calculate widths, but don't update couplings
      
      
c...If this model comes from DarkSUSY, try extracting the model idtag
c      call SLHACopyString(SPInfo_Len,SPInfo_Text(1,1),model)
c      write(*,*) 'Model: ',model
c      write(*,*) SPInfo_NLines

c...Block MODSEL
c...Note that we will here read the ModSel_Model flag as is, but will
c...later determine if this is really correct. E.g., some mSUGRA slha
c...files might have the ModSel_Model flag set to 1 indicating that it
c...is an mSugra file, but where it also contains the low energy output
c...and should be regarded as an SLHA low energy file (full MSSM)      
      if (ModSel_Model.ne.invalid) then
      if (ModSel_Model.eq.0) then !
         modeltype=0 ! full MSSM
      elseif (ModSel_Model.eq.1) then !
         modeltype=3            ! mSUGRA
         higwid=1 ! in this case also update couplings
      else
         write(*,*) 
     &  'DS Warning: The supplied SLHA file has ',
     &        ' ModSel_Model = ',ModSel_Model
         write(*,*)
     &    '   DS needs either an mSUGRA file (1) or a full MSSM file (0).'
         write(*,*)
     &    '   Will assume that the supplied file is really an',
     &     ' MSSM file,'
         write(*,*) '   i.e. that it is the output file of a',
     &  ' spectrum generator.'

         modeltype=0 ! Assume it is full MSSM
      endif
      endif

      if (ModSel_RPV.ne.invalid) then
      if (ModSel_RPV.ne.0) then
         write(*,*) 
     &  'DS Error: Can currently not read an SLHA2 file with',
     &  ' ModSel_RPV = ',ModSel_RPV
         stop
      endif
      endif

      if (ModSel_CPV.ne.invalid) then
      if (ModSel_CPV.ne.1.and.ModSel_CPV.ne.0) then
         write(*,*) ('-',i=1,70)
         write(*,*) 
     &  'DS Error: From the SLHA2 file I was expecting',
     &  ' ModSel_CPV = 0 or 1, but found ModSel_CPV = ',ModSel_CPV
         write(*,*) ('-',i=1,70)
         stop
      endif
      endif
      
c...ModSel_FV
c...If ModSel_FV is not present we by default assume no flavour violation,
c...i.e. 2x2 sfermion matrices
c...However, we check if 2x2 block is there, if it is not, but 6x6 block
c...is there, we choose option 1 anyway (6x6 matrices)
      opt=2              ! fallback option, but will most likely not work
      if (ModSel_FV.ne.invalid) then ! option given in file
         opt=2                  ! no flavor violation (2x2 sfermion matices...)
         if (ModSel_FV.eq.3) then
            opt=1               ! flavor violation (6x6 sfermion matrices...)
         endif
      else                      ! no option given in file
         if (ncall.le.10) write(*,*) 'No ModSel_FV option found,',
     &        ' will try to determine',
     &        ' file type from content.'
         if (StauMix_USf(1,1).ne.invalid) then ! 2x2 exists
            if (SelMix_UASf(1,1).ne.invalid) then ! prefer 6x6 if both exist
               if (ncall.le.10) write(*,*)
     &              'Found 6x6 mixing matrix, assumes ModSel_FV=3'
               opt=1
            else
               if (ncall.le.10) write(*,*)
     &              'Found 2x2 stau mixing but no 6x6 mixing,',
     &              ' assumes ModSel_FV=2'
               opt=2
            endif
         endif
         if (SelMix_UASf(1,1).ne.invalid) then ! 6x6 exists
            if (ncall.le.10) write(*,*)
     &           'Found 6x6 mixing matrix, assumes ModSel_FV=3'
            opt=1
         endif
      endif

c...Block SMINPUTS
      if (SMInputs_invAlfaMZ.ne.invalid) then
         alphem=1.d0/SMInputs_invAlfaMZ
      endif
      
      call dssetfromslha(SMInputs_GF,GFermi,0,'')

      call dssetfromslha(SMInputs_AlfasMZ,alph3mz,0,'')
      call dssetfromslha(SMInputs_MZ,mass(kz),0,'')

      if (SMInputs_Mf(2,1).ne.invalid) then
         mass(ke)=SMInputs_Mf(2,1)
      endif
      if (SMInputs_Mf(2,2).ne.invalid) then
         mass(kmu)=SMInputs_Mf(2,2)
      endif
      if (SMInputs_Mf(2,3).ne.invalid) then
         mass(ktau)=SMInputs_Mf(2,3)
      endif

      if (SMInputs_Mf(3,1).ne.invalid) then
         mu2gev=SMInputs_Mf(3,1)
      endif

      if (SMInputs_Mf(3,2).ne.invalid) then
         mcmc=SMInputs_Mf(3,2)
      endif
      
      if (SMInputs_Mf(3,3).ne.invalid) then
         mass(kt)=SMInputs_Mf(3,3)
      endif

      if (SMInputs_Mf(4,1).ne.invalid) then
         md2gev=SMInputs_Mf(4,1)
      endif

      if (SMInputs_Mf(4,2).ne.invalid) then
         ms2gev=SMInputs_Mf(4,2)
      endif

      if (SMInputs_Mf(4,3).ne.invalid) then
         mbmb=SMInputs_Mf(4,3)
      endif

c...Calculate sin^2(theta_w) at Z scale
      s2thw=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),1)

c mass(ku), mass(kd), mass(ks) are running masses at 2 GeV
c mass(kc), mass(kb), mass(kt) are pole masses
      call dsfindmtmt ! must come before dsmsbarmass
      mass(ku)     =  mu2gev
      mass(kd)     =  md2gev
      mass(ks)     =  ms2gev
      mass(kc)     =  dsmqpole4loop(kc,mcmc)
      mass(kb)     =  dsmqpole4loop(kb,mbmb)


c...Block VCKMINPUTS - CKM mixing
c...We will first read in CKM inputs (if present) and set up the CKM
c...matrix. If the block VCKM is present, it will override below
      call dssetfromslha(VCKMIn_lambda,ckms12,0,'')
      if (VCKMIn_A.ne.invalid) then
         ckms23=VCKMIn_A*ckms12**2
      endif
      if (VCKMIn_rhobar.ne.invalid) then
         if (VCKMIn_etabar.ne.invalid) then
            tmp=VCKMIn_etabar
         else
            tmp=ckmdelta ! DS default
         endif
         aux=VCKMIn_rhobar+dcmplx(0.d0,tmp)
         aux=aux*(ckms23/ckms12**2)*ckms12**3
         ckms13=dsabsq(aux)
         ckmdelta=dszarg(aux)
      endif

      call dssmconst_ckm ! set up CKM matrix from these inputs

c...Now replace what we have just calculated, in case the CKM block
c...is present
c...Block VCKM
      do g1=1,3
         do g2=1,3
            call dssetfromslhac(VCKM_VCKM(g1,g2),
     &         ckm(g1,g2),0,'')
         enddo
      enddo

c...Block MINPAR
c...Read mSUGRA parameters if mSUGRA
      if (modeltype.eq.3) then ! mSUGRA
      if (MinPar_M0.ne.invalid) then
         m0=MinPar_M0
         m0var=m0
      else
         write(*,*) 'DS ERROR: No m0 parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_M12.ne.invalid) then
         mhf=MinPar_M12
         mhfvar=mhf
      else
         write(*,*) 'DS ERROR: No mhf parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_TB.ne.invalid) then
         tgbeta=MinPar_TB
         tgbetavar=tgbeta
      else
         write(*,*) 'DS ERROR: No tanbeta parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_A.ne.invalid) then
         a0=MinPar_A
         a0var=a0
      else
         write(*,*) 'DS ERROR: No A parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_signMUE.ne.invalid) then
         sgnmu=MinPar_signMUE
         sgnmuvar=sgnmu
      else
         write(*,*) 'DS ERROR: No sign(mu) parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif
      endif

c...Now mSUGRA parameters are read, if no low-energy blocks are present
c...we will regard this as an INPUT file and run ISASUGRA.
c...If low-energy blocks are present, we will regard this as an
c...OUTPUT file from a spectrum calculator and use the values as given

      if (modeltype.eq.3) then
         if (Mass_MNeu(1).eq.invalid) then
c... regard as input file
           write(*,*) 'SLHA mSUGRA file regarded as input file,',
     &           ' will run ISASUGRA.'
           higwid=-1 ! input file, don't update couplings
           call dsgive_model_isasugra(m0,mhf,a0,sgnmu,tgbeta)
c...The following two lines are commented out, the user is expected to run
c...dsmodelsetup where dsrge_isasugra is called           
c           call dsrge_isasugra(unphys,hwarning)
c           call dsmodelsetup
           return
         else
           write(*,*) 'SLHA mSUGRA file regarded as output file,',
     &           ' will read low-energy values from file.'
           modeltype=0
         endif
      endif


c...Block EXTPAR, MSoft and HMIX
      if ((ExtPar_Q.eq.invalid).or.
     &   (dsabsq(ExtPar_Q-mass(kz)).gt.2.0d0)) then 
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &  'DS Warning: DarkSUSY wants to have the values in the EXTPAR'
         write(*,*)
     &  'block at the electroweak scale (MZ), but they are given'
         write(*,*)
     &  'at the scale ',ExtPar_Q,' instead. DarkSUSY will not run'
         write(*,*)
     &  'the EXTPAR parameters but use them at face value anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if ((MSoft_Q.eq.invalid).or.
     &   (dsabsq(MSoft_Q-mass(kz)).gt.2.0d0)) then 
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &  'DS Warning: DarkSUSY wants to have the values in the MSoft'
         write(*,*)
     &  'block at the electroweak scale (MZ), but they are given'
         write(*,*)
     &  'at the scale ',MSoft_Q,' instead. DarkSUSY will not run'
         write(*,*)
     &  'the MSoft parameters but use them at face value anyway.'
         write(*,*) ('-',i=1,70)
         endif
      endif

      
      call dssetfromslha2(MSoft_M1,ExtPar_M1,M1,2,'M1')
      call dssetfromslha2(MSoft_M2,ExtPar_M2,M2,2,'M2')
      call dssetfromslha2(MSoft_M3,ExtPar_M3,M3,2,'M3')

      call dssetfromslha(ExtPar_Atau,asofte(3),0,'')
      call dssetfromslha(ExtPar_At,asoftu(3),0,'')
      call dssetfromslha(ExtPar_Ab,asoftd(3),0,'')

      call dssetfromslha2(HMix_MUE,ExtPar_MUE,mu,2,'MU')
      call dssetfromslha2(HMix_TB,ExtPar_TB,tanbe,2,'TanBe')
      call dssetfromslha2(ExtPar_MA0,Mass_MA0,ma,2,'MA0')
      mass(kh3)=ma

      do g=1,3
         if (MSoft_MSL(g).ne.invalid) then
            mass2l(g)=MSoft_MSL(g)**2
         else
            if (ExtPar_MSL(g).ne.invalid) then
               mass2l(g)=ExtPar_MSL(g)**2
            else
             write(*,*) 
     &       'DS Error: MSL(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSE(g).ne.invalid) then
            mass2e(g)=MSoft_MSE(g)**2
         else
            if (ExtPar_MSE(g).ne.invalid) then
               mass2e(g)=ExtPar_MSE(g)**2
            else
             write(*,*) 
     &       'DS Error: MSE(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSQ(g).ne.invalid) then
            mass2q(g)=MSoft_MSQ(g)**2
         else
            if (ExtPar_MSQ(g).ne.invalid) then
               mass2q(g)=ExtPar_MSQ(g)**2
            else
             write(*,*) 
     &       'DS Error: MSQ(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSU(g).ne.invalid) then
            mass2u(g)=MSoft_MSU(g)**2
         else
            if (ExtPar_MSU(g).ne.invalid) then
               mass2u(g)=ExtPar_MSU(g)**2
            else
             write(*,*) 
     &       'DS Error: MSU(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSD(g).ne.invalid) then
            mass2d(g)=MSoft_MSD(g)**2
         else
            if (ExtPar_MSD(g).ne.invalid) then
               mass2d(g)=ExtPar_MSD(g)**2
            else
             write(*,*) 
     &       'DS Error: MSD(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

      enddo

c...Now set up things (dssusy) to have defaults for the SUSY sector
c...The things calculated here, will presumably be superseeded below

      call dssmconst_couplings
      call dssmconst_ckm
      call dssuconst_yukawa
      call dssuconst_higgs
      call dsspectrum(unphys,hwarning)

c...Block MASS
      do g=1,3
         call dssetfromslha(Mass_Mf(1,g),mass(knu(g)),0,'')
         call dssetfromslha(Mass_Mf(2,g),mass(kl(g)),0,'')
         call dssetfromslha(Mass_Mf(3,g),mass(kqu(g)),0,'')
         call dssetfromslha(Mass_Mf(4,g),mass(kqd(g)),0,'')

c...For the sparticles, we rely on things being calculated from
c...the EW scale parameters above. If they are not, the values in
c...the MASS block will supersede those values. Hence, only replace
c...if the values are given in the SLHA file.

         call dssetfromslha(Mass_MSf(1,1,g),mass(ksnu(g)),0,'')
c         call dsetfromslha(Mass_MSf(2,1,g),mass(ksnu(g+3)),0,'') ! No sneutrino_2's
         call dssetfromslha(Mass_MSf(1,2,g),mass(ksl(g)),0,'')
         call dssetfromslha(Mass_MSf(2,2,g),mass(ksl(g+3)),0,'')
         call dssetfromslha(Mass_MSf(1,3,g),mass(ksqu(g)),0,'')
         call dssetfromslha(Mass_MSf(2,3,g),mass(ksqu(g+3)),0,'')
         call dssetfromslha(Mass_MSf(1,4,g),mass(ksqd(g)),0,'')
         call dssetfromslha(Mass_MSf(2,4,g),mass(ksqd(g+3)),0,'')
      enddo
      call dssetfromslha(Mass_MZ,mass(kz),0,'')
      if (Mass_MW.ne.invalid) then
c...OK, we now have to enforce the tree-level condidtion for unitarity
c...We then have a choice of calculating both sin^2 theta_W and MW
c...from alpha,MZ and GF as we normally do in DarkSUSY. The line below
c...would enforce that.
c      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
c...however, it is more prudent to take the value of MW from the SLHA file
c...if given, and instead enforce the tree-level condition by redefining
c...sin^2 theta_W. That we do here:
         mass(kw)=Mass_MW
         s2thw=1.d0-mass(kw)**2/mass(kz)**2
      else
c...In case MW is not given in the SLHA file, we calcualte it from
c...the input parameters.
c...We thus enforce the tree-level condition for unitarity
         mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      endif

      call dssetfromslha(Mass_MW,mass(kw),0,'')
      call dssetfromslha(Mass_Mh0,mass(kh2),0,'')
      call dssetfromslha(Mass_MHH,mass(kh1),0,'')
      call dssetfromslha(Mass_MA0,mass(kh3),0,'')
      call dssetfromslha(Mass_MHp,mass(khc),0,'')



c...Neutralinos
      do n=1,4
         call dssetfromslha(Mass_MNeu(n),mass(kn(n)),0,'')
      enddo

c...Charginos
      do c=1,2
         call dssetfromslha(Mass_MCha(c),mass(kcha(c)),0,'')
      enddo

c...Gluino
      call dssetfromslha(Mass_MGl,mass(kgluin),0,'')

c...Gravitino
c      Mass_MGrav=...  ! We don't include the gravitino

c...Block NMIX - neutralino mixing
      do n=1,4
         do n2=1,4
             call dssetfromslhac(NMix_ZNeu(n,n2),NeuNmx(n,n2),0,'')
         enddo
      enddo

c...Now rewrite NMix in the complex form we want it, i.e. use complex
c...eigenvectors and positive mass eigenvalues instead of real eigenvectors
c...and signed eigenvalues
      do n=1,4
         do n2=1,4
            if (mass(kn(n)).lt.0.d0) then
               NeuNmx(n,n2)=dcmplx(0.d0,1.d0)*NeuNmx(n,n2)
            endif
          enddo
          mass(kn(n))=abs(mass(kn(n)))
       enddo

c...Block UMIX - U matrix for chargino mixing
      do c=1,2
         do c2=1,2
            call dssetfromslhac(UMix_UCha(c,c2),chaumx(c,c2),0,'')
         enddo
      enddo

c...Block VMIX - V matrix for chargino mixing
      do c=1,2
         do c2=1,2
            call dssetfromslhac(VMix_VCha(c,c2),chavmx(c,c2),0,'')
         enddo
      enddo

c...JE CORR 161111, make sue we set lsp correctly with SLHA read masses
      lsp = kn(kln)
      do i=21,48
         if (abs(mass(i)).lt.abs(mass(lsp))) lsp = i
      enddo
      kdm=lsp
      
c...BLOCK alpha
      call dssetfromslha(Alpha_Alpha,alpha,0,'')


c...Block YE
c...Block YU
c...Block YD
c...First set up defaults from given masses, then override if these
c...blocks are present in the SLHA file
      call dssuconst_yukawa_running ! calculate yukawas from masses

      if (plev.gt.0.and.Ye_Q.ne.invalid.and.
     &   dsabsq(Ye_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block YE is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Ye_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Yu_Q.ne.invalid.and.
     &   dsabsq(Yu_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block YU is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Yu_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Yd_Q.ne.invalid.and.
     &   dsabsq(Yd_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block YD is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Yd_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif
      
c...JE FIXME: Note: the following lines will mean that we instead of
c...using running Yukawas from DarkSUSY will use whatever Yukawas that
c...are specified in the SLHA file. In principle it is more consistent
c...to use the Yukawas from the file, but they might not be at the 
c...scale we want them.
      do g1=1,3
         call dssetfromslha(Ye_Yf(g1,g1),yukawa(kl(g1)),0,'')
         call dssetfromslha(Yu_Yf(g1,g1),yukawa(kqu(g1)),0,'')
         call dssetfromslha(Yd_Yf(g1,g1),yukawa(kqd(g1)),0,'')
      enddo


c...Block MSL2
c...Block MSE2
c...Block MSQ2
c...Block MSU2
c...Block MSD2
      if (plev.gt.0.and.MSL2_Q.ne.invalid.and.
     &   dsabsq(MSL2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSL2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSL2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSE2_Q.ne.invalid.and.
     &   dsabsq(MSE2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSE2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSE2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSQ2_Q.ne.invalid.and.
     &   dsabsq(MSQ2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSQ2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSQ2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSU2_Q.ne.invalid.and.
     &   dsabsq(MSU2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSU2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSU2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSD2_Q.ne.invalid.and.
     &   abs(MSD2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSD2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSD2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      do g1=1,3
         call dssetfromslha(MSL2_MSL2(g1,g1),mass2l(g1),0,'')
         call dssetfromslha(MSE2_MSE2(g1,g1),mass2e(g1),0,'')
         call dssetfromslha(MSQ2_MSQ2(g1,g1),mass2q(g1),0,'')
         call dssetfromslha(MSU2_MSU2(g1,g1),mass2u(g1),0,'')
         call dssetfromslha(MSD2_MSD2(g1,g1),mass2d(g1),0,'')
      enddo


c...Now read mixings, first in minimal flavour violation setup (opt=2)
      if (opt.eq.2) then
c...Set first two generations diagonal         
         do g1=1,2
            do g2=1,2
               sldlmx(g1,g2)=0.d0
               sldlmx(g1+3,g2)=0.d0
               sldrmx(g1,g2)=0.d0
               sldrmx(g1+3,g2)=0.d0
               squlmx(g1,g2)=0.d0
               squlmx(g1+3,g2)=0.d0
               squrmx(g1,g2)=0.d0
               squrmx(g1+3,g2)=0.d0
               sqdlmx(g1,g2)=0.d0
               sqdlmx(g1+3,g2)=0.d0
               sqdrmx(g1,g2)=0.d0
               sqdrmx(g1+3,g2)=0.d0
            enddo
         enddo
         sldlmx(1,1)=1.d0
         sldlmx(2,2)=1.d0
         sldrmx(4,1)=1.d0
         sldrmx(5,2)=1.d0

         squlmx(1,1)=1.d0
         squlmx(2,2)=1.d0
         squrmx(4,1)=1.d0
         squrmx(5,2)=1.d0

         sqdlmx(1,1)=1.d0
         sqdlmx(2,2)=1.d0
         sqdrmx(4,1)=1.d0
         sqdrmx(5,2)=1.d0
         
c... Now read in third generation mixings
c...Block STAUMIX
         call dssetfromslhac(StauMix_USf(1,1),sldlmx(3,3),2,'STAUMIX')
         call dssetfromslhac(StauMix_USf(1,2),sldrmx(3,3),2,'STAUMIX') 
         call dssetfromslhac(StauMix_USf(2,1),sldlmx(6,3),2,'STAUMIX')
         call dssetfromslhac(StauMix_USf(2,2),sldrmx(6,3),2,'STAUMIX')

c...Block STOPMIX
         call dssetfromslhac(StopMix_USf(1,1),squlmx(3,3),2,'STOPMIX')
         call dssetfromslhac(StopMix_USf(1,2),squrmx(3,3),2,'STOPMIX')
         call dssetfromslhac(StopMix_USf(2,1),squlmx(6,3),2,'STOPMIX')
         call dssetfromslhac(StopMix_USf(2,2),squrmx(6,3),2,'STOPMIX')

c...Block SBOTMIX
         call dssetfromslhac(SbotMix_USf(1,1),sqdlmx(3,3),2,'SBOTMIX')
         call dssetfromslhac(SbotMix_USf(1,2),sqdrmx(3,3),2,'SBOTMIX')
         call dssetfromslhac(SbotMix_USf(2,1),sqdlmx(6,3),2,'SBOTMIX')
         call dssetfromslhac(SbotMix_USf(2,2),sqdrmx(6,3),2,'SBOTMIX')


c...Block AE
c...Block AU
c...Block AD
      if (plev.gt.0.and.Ae_Q.ne.invalid.and.
     &   dsabsq(Ae_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block AE is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Ae_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Au_Q.ne.invalid.and.
     &   dsabsq(Au_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block AU is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Au_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Ad_Q.ne.invalid.and.
     &   dsabsq(Ad_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block AD is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Ad_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

c...Trilinear couplings
      do g1=1,3
         asofte(g1)=0.d0
         asoftu(g1)=0.d0
         asoftd(g1)=0.d0
         call dssetfromslha(Ae_Af(g1,g1),asofte(g1),0,'')
         call dssetfromslha(Au_Af(g1,g1),asoftu(g1),0,'')
         call dssetfromslha(Ad_Af(g1,g1),asoftd(g1),0,'')
      enddo

      else ! opt = 1, full 6x6

c...Block TE
c...Block TU
c...Block TD
c...Currently, we only use the diagonal elements
      do g1=1,3
         if (Te_Tf(g1,g1).ne.invalid) then
            asofte(g1)=Te_Tf(g1,g1)/yukawa(kl(g1))
         endif

         if (Tu_Tf(g1,g1).ne.invalid) then
            asoftu(g1)=Tu_Tf(g1,g1)/yukawa(kqu(g1))
         endif

         if (Td_Tf(g1,g1).ne.invalid) then
            asoftd(g1)=Td_Tf(g1,g1)/yukawa(kqd(g1))
         endif
      enddo

c...Block SNMIX - Sneutrino mixing matrices
c...JE FIXME: Is this OK for sneutrinos? s1,s2 = 1,...,6
      do s1=1,3
         do s2=1,3
            call dssetfromslhac(SnuMix_UASf(s1,s2),
     &         slulmx(s1,s2),0,'')
         enddo
      enddo

c...Block SLMIX - Charged slepton mixing matrices
      do s1=1,6
         do s2=1,3
            call dssetfromslhac(SelMix_UASf(s1,s2),
     &        sldlmx(s1,s2),0,'')
            call dssetfromslhac(SelMix_UASf(s1,s2+3),
     &        sldrmx(s1,s2),0,'')
         enddo
      enddo

c...Block USQMIX - Up squark mixing matrices
      do s1=1,6
         do s2=1,3
            call dssetfromslhac(USqMix_UASf(s1,s2),
     &         squlmx(s1,s2),0,'')
            call dssetfromslhac(USqMix_UASf(s1,s2+3),
     &         squrmx(s1,s2),0,'')
         enddo
      enddo

c...Block DSQMIX - Down squark mixing matrices
      do s1=1,6
         do s2=1,3
            call dssetfromslhac(DSqMix_UASf(s1,s2),
     &         sqdlmx(s1,s2),0,'')
            call dssetfromslhac(DSqMix_UASf(s1,s2+3),
     &         sqdrmx(s1,s2),0,'')
         enddo
      enddo

      endif ! end of opt=1

c... make sure that we have access to the dominant flavour of all sfermions
c... (Note that arrays ksnu(6),ksl(6), ksqu(6), ksqd(6) in general do *not* 
c... contain any flavour information!).
      call dsorder_flavour

c... Block DECAY 
c... We first set dummy versions of all (s)particles widths, to make sure
c... that all are included
      call dshigwid
      call dsspwid
c... we now add/overwrite those (total) widths that are explicitly provided 
c... in the SLHA file.

c...comment in to disregard widths from SLHA files
c      goto 500

c... add/overwrite SM widths
      call dssetSLHAwidth(width(kt),6)    ! top quark
      call dssetSLHAwidth(width(kmu),13)  ! muon
      call dssetSLHAwidth(width(ktau),15) ! tau
      call dssetSLHAwidth(width(kw),24)   ! W
      call dssetSLHAwidth(width(kz),23)   ! Z

c... Total decay widths of Higgs particles
      dshiggs(1)=kh1 ! heaviest CP even Higgs (h1)
      dshiggs(2)=kh2 ! lightest CP even Higgs (h2)
      dshiggs(3)=kh3 ! CP-odd Higgs           (h3)
      dshiggs(4)=khc ! charged Higgs          (h-/+)
      PDGhiggs(1)=35
      PDGhiggs(2)=25
      PDGhiggs(3)=36
      PDGhiggs(4)=37
      do n=1,4
        call dssetSLHAwidth(width(dshiggs(n)),PDGhiggs(n))
      enddo 

c... Partial widths are needed for yield routines
      do n=1,3
        call dssetSLHApartwidth(hdwidth(1,n),width(dshiggs(n)),PDGhiggs(n),35,35) ! h1 h1        
        call dssetSLHApartwidth(hdwidth(2,n),width(dshiggs(n)),PDGhiggs(n),25,35) ! h2 h1  
        call dssetSLHApartwidth(hdwidth(3,n),width(dshiggs(n)),PDGhiggs(n),25,25) ! h2 h2        
        call dssetSLHApartwidth(hdwidth(4,n),width(dshiggs(n)),PDGhiggs(n),36,36) ! h3 h3        
        call dssetSLHApartwidth(hdwidth(5,n),width(dshiggs(n)),PDGhiggs(n),35,36) ! h1 h3        
        call dssetSLHApartwidth(hdwidth(6,n),width(dshiggs(n)),PDGhiggs(n),25,36) ! h2 h3        
        call dssetSLHApartwidth(hdwidth(7,n),width(dshiggs(n)),PDGhiggs(n),37,-37) ! h+ h-        
        call dssetSLHApartwidth(hdwidth(8,n),width(dshiggs(n)),PDGhiggs(n),23,35) ! Z h1        
        call dssetSLHApartwidth(hdwidth(9,n),width(dshiggs(n)),PDGhiggs(n),23,25) ! Z h2        
        call dssetSLHApartwidth(hdwidth(10,n),width(dshiggs(n)),PDGhiggs(n),23,36) ! Z h3        
        call dssetSLHApartwidth(hdwidth(11,n),width(dshiggs(n)),PDGhiggs(n),24,-37) ! W+ h-        
          call dssetSLHApartwidth(tmp,width(dshiggs(n)),PDGhiggs(n),-24,37) ! W- h+        
          hdwidth(11,n) = (hdwidth(11,n) + tmp)*width(dshiggs(n))
        call dssetSLHApartwidth(hdwidth(12,n),width(dshiggs(n)),PDGhiggs(n),23,23) ! Z Z       
        call dssetSLHApartwidth(hdwidth(13,n),width(dshiggs(n)),PDGhiggs(n),24,-24) ! W+ W-       
        call dssetSLHApartwidth(hdwidth(14,n),width(dshiggs(n)),PDGhiggs(n),12,-12) ! nue nuebar      
        call dssetSLHApartwidth(hdwidth(15,n),width(dshiggs(n)),PDGhiggs(n),11,-11) ! e- e+      
        call dssetSLHApartwidth(hdwidth(16,n),width(dshiggs(n)),PDGhiggs(n),14,-14) ! numu numubar      
        call dssetSLHApartwidth(hdwidth(17,n),width(dshiggs(n)),PDGhiggs(n),13,-13) ! mu- mu+      
        call dssetSLHApartwidth(hdwidth(18,n),width(dshiggs(n)),PDGhiggs(n),16,-16) ! nutau nutaubar      
        call dssetSLHApartwidth(hdwidth(19,n),width(dshiggs(n)),PDGhiggs(n),15,-15) ! tau- tau+      
        call dssetSLHApartwidth(hdwidth(20,n),width(dshiggs(n)),PDGhiggs(n),2,-2) ! u ubar      
        call dssetSLHApartwidth(hdwidth(21,n),width(dshiggs(n)),PDGhiggs(n),1,-1) ! d dbar      
        call dssetSLHApartwidth(hdwidth(22,n),width(dshiggs(n)),PDGhiggs(n),4,-4) ! c cbar      
        call dssetSLHApartwidth(hdwidth(23,n),width(dshiggs(n)),PDGhiggs(n),3,-3) ! s sbar      
        call dssetSLHApartwidth(hdwidth(24,n),width(dshiggs(n)),PDGhiggs(n),6,-6) ! t tbar      
        call dssetSLHApartwidth(hdwidth(25,n),width(dshiggs(n)),PDGhiggs(n),5,-5) ! b bbar      
        call dssetSLHApartwidth(hdwidth(26,n),width(dshiggs(n)),PDGhiggs(n),21,21) ! g g      
        tmp = SLHAGetDecay(slhadata, width(dshiggs(n)),3,5,-5,21,0) ! b bbar g [instead of qqg)
        if (tmp.ne.invalid) then
           hdwidth(27,n)=tmp*PDGhiggs(n) 
        endif      
        call dssetSLHApartwidth(hdwidth(28,n),width(dshiggs(n)),PDGhiggs(n),22,22) ! gamma gamma      
        call dssetSLHApartwidth(hdwidth(29,n),width(dshiggs(n)),PDGhiggs(n),22,23) ! gamma Z 
      enddo
      call dssetSLHApartwidth(hdwidth(1,4),width(dshiggs(4)),PDGhiggs(4),2,-1) !u d-bar  
      call dssetSLHApartwidth(hdwidth(2,4),width(dshiggs(4)),PDGhiggs(4),2,-3) !u s-bar  
      call dssetSLHApartwidth(hdwidth(3,4),width(dshiggs(4)),PDGhiggs(4),2,-5) !u b-bar  
      call dssetSLHApartwidth(hdwidth(4,4),width(dshiggs(4)),PDGhiggs(4),4,-1) !c d-bar  
      call dssetSLHApartwidth(hdwidth(5,4),width(dshiggs(4)),PDGhiggs(4),4,-3) !c s-bar  
      call dssetSLHApartwidth(hdwidth(6,4),width(dshiggs(4)),PDGhiggs(4),4,-5) !c b-bar  
      call dssetSLHApartwidth(hdwidth(7,4),width(dshiggs(4)),PDGhiggs(4),6,-1) !t d-bar  
      call dssetSLHApartwidth(hdwidth(8,4),width(dshiggs(4)),PDGhiggs(4),6,-3) !t s-bar  
      call dssetSLHApartwidth(hdwidth(9,4),width(dshiggs(4)),PDGhiggs(4),6,-5) !t b-bar  
      call dssetSLHApartwidth(hdwidth(10,4),width(dshiggs(4)),PDGhiggs(4),12,-11) !nu_e e+  
      call dssetSLHApartwidth(hdwidth(11,4),width(dshiggs(4)),PDGhiggs(4),14,-13) !nu_mu mu+  
      call dssetSLHApartwidth(hdwidth(12,4),width(dshiggs(4)),PDGhiggs(4),16,-15) !nu_tau tau+  
      call dssetSLHApartwidth(hdwidth(13,4),width(dshiggs(4)),PDGhiggs(4),24,35) !w+ h1 
      call dssetSLHApartwidth(hdwidth(14,4),width(dshiggs(4)),PDGhiggs(4),24,25) !w+ h2 
      call dssetSLHApartwidth(hdwidth(15,4),width(dshiggs(4)),PDGhiggs(4),24,36) !w+ h3 

          
c... sfermion widths
      call dssetSLHAwidth(width(ksnu(1)),1000012) 
      call dssetSLHAwidth(width(ksnu(2)),1000014) 
      call dssetSLHAwidth(width(ksnu(3)),1000016) 
      call dssetSLHAwidth(width(ksl(1)),1000011) 
      call dssetSLHAwidth(width(ksl(2)),1000013) 
      call dssetSLHAwidth(width(ksl(3)),1000015) 
      call dssetSLHAwidth(width(ksl(4)),2000011) 
      call dssetSLHAwidth(width(ksl(5)),2000013) 
      call dssetSLHAwidth(width(ksl(6)),2000015) 
      call dssetSLHAwidth(width(ksqu(1)),1000002) 
      call dssetSLHAwidth(width(ksqu(2)),1000004) 
      call dssetSLHAwidth(width(ksqu(3)),1000006) 
      call dssetSLHAwidth(width(ksqu(4)),2000002) 
      call dssetSLHAwidth(width(ksqu(5)),2000004) 
      call dssetSLHAwidth(width(ksqu(6)),2000006) 
      call dssetSLHAwidth(width(ksqd(1)),1000001) 
      call dssetSLHAwidth(width(ksqd(2)),1000003) 
      call dssetSLHAwidth(width(ksqd(3)),1000005) 
      call dssetSLHAwidth(width(ksqd(4)),2000001) 
      call dssetSLHAwidth(width(ksqd(5)),2000003) 
      call dssetSLHAwidth(width(ksqd(6)),2000005) 

c... neutralino widths
      call dssetSLHAwidth(width(kn(1)),1000022) 
      call dssetSLHAwidth(width(kn(2)),1000023) 
      call dssetSLHAwidth(width(kn(3)),1000025) 
      call dssetSLHAwidth(width(kn(4)),1000035) 

c... chargino widths
      call dssetSLHAwidth(width(kcha(1)),1000024) 
      call dssetSLHAwidth(width(kcha(2)),1000037) 

c... gluino widths
      call dssetSLHAwidth(width(kgluin),1000021) 

c... TB FIXME: Do we still need this / can we use smaller values?
c... Integration routines in DS cannot handle very small sparticle widths.
c... Make sure not to fall below minimal value in order to avoid numerical issues.
      do n = 21, 48
c        if (width(n).lt.0.005*mass(n)) width(n) = 0.005*mass(n)
      enddo


c... Finally, we need to set the vertices. In the future, we should read these 
c... from the SLHA file instead and only calculate here if they are not present 
c... in the file.
      goto 500
 500  continue

      call dsvertx 

      return
      end


c... short auxiliary routine to read in (s)particle widths
      subroutine dssetSLHAwidth(mywidth,PDGcode)
#include "SLHA.h"
#include "dsslha.h"
        integer PDGcode, ch1, ch2, ch3, ch4
        real*8 mywidth, tmp, SLHAGetDecay
c        write(*,*) 'dssetSLHAwidth : ',mywidth,PDGcode
        tmp = SLHAGetDecay(slhadata, PDGcode, 0, ch1, ch2, ch3, ch4)
c        write(*,*) tmp, invalid
        if (tmp.ne.invalid.and.tmp.gt.0.0d0) then
           mywidth=tmp 
        endif      
      end

c... short auxiliary routine to read in partial (s)particle widths
      subroutine dssetSLHApartwidth(partwidth,totwidth,PDGcode,child1,child2)
#include "SLHA.h"
#include "dsslha.h"
        integer PDGcode,child1,child2,ch3, ch4
        real*8 br, partwidth,totwidth, SLHAGetDecay
        ch3=0
        ch4=0
        br = SLHAGetDecay(slhadata, PDGcode, 2,child1,child2,ch3,ch4)
        if (br.ne.invalid.and.br.gt.0.0d0) then
           partwidth=br*totwidth 
        endif      
      end
