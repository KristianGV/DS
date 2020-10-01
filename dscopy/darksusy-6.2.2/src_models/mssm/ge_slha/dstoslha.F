**********************************************************************
*** subroutine dstoSLHA transfers DarkSUSY data to the SUSY Les
*** Houches Accord format (into the array slhadata)
*** with the help of SLHALIB by Thomas Hahn.
*** Input: opt = 1 - an SLHA2 file with full 6x6 squark and slepton
***                  mixings and 3x3 trilinear couplings is written.
***                  For minimal flavour violation
***                  this is a bit of an overkill, but more consistent
***                  with how things are stored internally in DarkSUSY
***              2 - a minimal flavour violation SLHA2 file is written
***                  with only third generation mixings and trilinear
***                  couplings. This option will in principle be an
***                  approximation to the internal setup in DarkSUSY,
***                  but most likely good enough for most scenarios.
***              3 - an mSUGRA input file (i.e. only with the SM and
***                  mSUGRA input parameters). This option could be
***                  of interest when the RGEs should be solved with
***                  an external program.
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: 2008-06-24
*** mod:  2019-04-26 torsten.bringmann@fys.uio.no (added block DECAY)
**********************************************************************
      subroutine dstoSLHA(opt)
      implicit none
c...Note: These include files need to be included with the precompiler
c...directive for include and not the fortran one. Hence, this file
c...has to end on .F instead of .f as well.

#include "SLHA.h"

#include "PDG.h"

      include 'dsmssm.h'
      include 'dsidtag.h'
      include 'dsver.h'
#include "dsslha.h"

c      complex*16 aux
      integer opt
      integer g,n,c,n2,c2,s1,s2,g1,g2
      integer decayID, dshiggs(4),PDGhiggs(4)
      real*8 tmp, br

c...Clear slhadata
      call SLHAClear(slhadata)

c...Block SPINFO
      if (opt.ne.3) then
         call SLHAPutInfo(BlockSPInfo(1),1,"DarkSUSY")
         call SLHAPutInfo(BlockSPInfo(1),2,dsversion)
         call SLHAPutInfo(BlockSPInfo(1),3,"Model: "//idtag)
      endif


C...Block MODSEL
c      ModSel_Model=...
      ModSel_Model=0 ! 0=MSSM
      if (modeltype.eq.3) then
         ModSel_Model=1 ! 1=mSUGRA
      endif
      if (opt.ne.3) then
         ModSel_RPV = 0         ! no R-partity violation
         if (abs(ckmdelta).lt.1.d-10) then
            ModSel_CPV = 0      ! 0=No CP violation
         else
            ModSel_CPV = 1  ! 1=Only CKM-CP violation
         endif
      endif
c...For the flavour violation flag, we set the flag to 3, even in the MFV
c...scenario. The reason for this is that we want to give the full
c...6x6 slepton and squark mixings below
      if (opt.eq.2) then
         ModSel_FV = 0  ! 0=no flavour violation, 6x6 SQMIX etc should be absent
      elseif (opt.eq.1) then
         ModSel_FV = 3  ! 3=flavour violation in quarks and leptons
      elseif (opt.eq.3) then
         continue ! do nothing in this case
      else
         write(*,*) 'DS Error in dstoslha: invalid option: ',opt
         ModSel_FV = 3  ! 3=flavour violation in quarks and leptons
         opt=1
      endif


c...Block SMINPUTS
      SMInputs_invAlfaMZ=1.d0/alphem
      SMInputs_GF=GFermi
      SMInputs_AlfasMZ = alph3mz
      SMInputs_MZ = mass(kz)

c...JE FIXME: running masses at 2 GeV
      SMInputs_Mf(3,3) = mass(kt)
      SMInputs_Mf(4,3) = mbmb
      SMInputs_Mf(2,3) = mass(ktau)
      if (opt.ne.3) then ! add more stuff
        SMInputs_Mf(2,1) = mass(ke)
        SMInputs_Mf(2,2) = mass(kmu)
        SMInputs_Mf(3,1) = mass(ku)
        SMInputs_Mf(3,2) = mcmc
        SMInputs_Mf(4,1) = mass(kd)
        SMInputs_Mf(4,2) = mass(ks)
      endif


c...Give MINPAR block if we started with mSUGRA
      if (modeltype.eq.3) then
         MinPar_M0=m0var
         MinPar_M12=mhfvar
         MinPar_A=a0var
         MinPar_signMUE=sgnmuvar
         MinPar_TB=tgbetavar
      endif

      if (opt.eq.3) return ! no extra blocks needed for mSUGRA input SLHA file

c...Block EXTPAR
c      ExtPar_Q = -1 ! this gives scale sqrt(m_stop1*m_stop2), OK?
      ExtPar_Q = mass(kz) ! Input scale is mZ.
      ExtPar_M1 = M1
      ExtPar_M2 = M2
      ExtPar_M3 = M3
      ExtPar_Atau=asofte(3)
      ExtPar_At=asoftu(3)
      ExtPar_Ab=asoftd(3)
      ExtPar_MUE = mu
c      ExtPar_MA02 = ...  ! this is tree-level mass, not pole mass. 
      ExtPar_TB = tanbe
      ExtPar_MA0=mass(kh3) ! this is pole mass
      do g=1,3
          ExtPar_MSL(g)=sqrt(mass2l(g))
          ExtPar_MSE(g)=sqrt(mass2e(g))
          ExtPar_MSQ(g)=sqrt(mass2q(g))
          ExtPar_MSU(g)=sqrt(mass2u(g))
          ExtPar_MSD(g)=sqrt(mass2d(g))
      enddo

c...Block MSoft
      MSoft_Q=mass(kz) ! Input scale is mz
      MSoft_M1 = M1
      MSoft_M2 = M2
      MSoft_M3 = M3
      do g=1,3
         MSoft_MSL(g)=sqrt(mass2l(g))
         MSoft_MSE(g)=sqrt(mass2e(g))
         MSoft_MSQ(g)=sqrt(mass2q(g))
         MSoft_MSU(g)=sqrt(mass2u(g))
         MSoft_MSD(g)=sqrt(mass2d(g))
      enddo

c...Block MASS
      do g=1,3
         Mass_Mf(1,g) = mass(knu(g))
         Mass_Mf(2,g) = mass(kl(g))
         Mass_Mf(3,g) = mass(kqu(g))
         Mass_Mf(4,g) = mass(kqd(g))
         Mass_MSf(1,1,g) = mass(ksnu(g))
c         Mass_MSf(2,1,g) = mass(ksnu(g+3)) ! No sneutrino_2's
         Mass_MSf(1,2,g) = mass(ksl(g))   
         Mass_MSf(2,2,g) = mass(ksl(g+3))
         Mass_MSf(1,3,g) = mass(ksqu(g))
         Mass_MSf(2,3,g) = mass(ksqu(g+3))
         Mass_MSf(1,4,g) = mass(ksqd(g))
         Mass_MSf(2,4,g) = mass(ksqd(g+3))
      enddo
c...Note, in DS we might have the first two generations of suarks
c...and sleptons mass-ordered instead of left-right-ordered. If that
c...is the case and we have opt=2 (i.e. if we do not write out the full
c...mixing matrices, change order of sfermion masses)
      if (opt.eq.2) then
         do g=1,2
            if (abs(sldlmx(g,g)).lt.0.5d0) then ! charged sleptons
               tmp=Mass_MSf(1,2,g)
               Mass_MSf(1,2,g)=Mass_MSf(2,2,g)
               Mass_MSf(2,2,g)=tmp
            endif
c JE debug
c            write(*,*) squlmx(4,1)
            if (abs(squlmx(g,g)).lt.0.5d0) then ! up squarks
               tmp=Mass_MSf(1,3,g)
               Mass_MSf(1,3,g)=Mass_MSf(2,3,g)
               Mass_MSf(2,3,g)=tmp
            endif
            if (abs(sqdlmx(g,g)).lt.0.5d0) then ! down squarks
               tmp=Mass_MSf(1,4,g)
               Mass_MSf(1,4,g)=Mass_MSf(2,4,g)
               Mass_MSf(2,4,g)=tmp
            endif
         enddo
      endif
               
      Mass_MZ = mass(kz)
      Mass_MW = mass(kw)
      Mass_Mh0 = mass(kh2)
      Mass_MHH = mass(kh1)
      Mass_MA0 = mass(kh3)
      Mass_MHp = mass(khc)

c...Neutralinos
      do n=1,4
         if (abs(dimag(NeuNmx(n,1))).gt.0.d0) then ! complex
            Mass_MNeu(n)=-mass(kn(n))
         else ! real
            Mass_MNeu(n)=mass(kn(n))
         endif
      enddo

c...Charginos
      do c=1,2
         Mass_MCha(c)=mass(kcha(c))
      enddo

c...Gluino
      Mass_MGl=mass(kgluin)

c...Gravitino
c      Mass_MGrav=...  ! We don't include the gravitino

c...Block NMIX - neutralino mixing
      do n=1,4
         if (abs(dimag(NeuNmx(n,1))).gt.0.d0) then ! complex, pick complex part
            do n2=1,4
               NMix_ZNeu(n,n2)=dcmplx(0.d0,-1.d0)*NeuNmx(n,n2)
            enddo
         else ! pick real part
            do n2=1,4
               NMix_ZNeu(n,n2)=NeuNmx(n,n2)
            enddo
         endif
      enddo

c...Block UMIX - U matrix for chargino mixing
      do c=1,2
         do c2=1,2
            UMix_UCha(c,c2)=chaumx(c,c2)
         enddo
      enddo

c...Block VMIX - V matrix for chargino mixing
      do c=1,2
         do c2=1,2
            VMix_VCha(c,c2)=chavmx(c,c2)
         enddo
      enddo

c...BLOCK alpha
      Alpha_Alpha = alpha

      if (opt.eq.2) then ! simplified MFV output

c...Block STAUMIX
         StauMix_USf(1,1)=sldlmx(3,3) 
         StauMix_USf(1,2)=sldrmx(3,3) 
         StauMix_USf(2,1)=sldlmx(6,3) 
         StauMix_USf(2,2)=sldrmx(6,3)

c...Block STOPMIX
         StopMix_USf(1,1)=squlmx(3,3) 
         StopMix_USf(1,2)=squrmx(3,3) 
         StopMix_USf(2,1)=squlmx(6,3) 
         StopMix_USf(2,2)=squrmx(6,3)

c...Block SBOTMIX
         SbotMix_USf(1,1)=sqdlmx(3,3) 
         SbotMix_USf(1,2)=sqdrmx(3,3) 
         SbotMix_USf(2,1)=sqdlmx(6,3) 
         SbotMix_USf(2,2)=sqdrmx(6,3)

c...NOTE: For the blocks AE, AU and AD, these things are 
c...in principle defined in the block ExtPar. We here give
c...them anyway for consistency. For the full thing
c...we should use the blocks TE, TU and TD below instead.
c...Block AE
c...Block AU
c...Block AD
         Ae_Q = mass(kz)
         Au_Q = mass(kz)
         Ad_Q = mass(kz)

cc      do g1=1,3
cc         do g2=1,3
cc            Ae_Af(g1,g2)=dcmplx(0.0d0,0.0d0)
cc            Au_Af(g1,g2)=dcmplx(0.0d0,0.0d0)
cc            Ad_Af(g1,g2)=dcmplx(0.0d0,0.0d0)
cc         enddo
cc         Ae_Af(g1,g1)=asofte(g1)
cc         Au_Af(g1,g1)=asoftu(g1)
cc         Ad_Af(g1,g1)=asoftd(g1)
cc      enddo
cc...Note: only (3,3) component should be given, hence do not use above      
         Ae_Af(3,3)=asofte(3)
         Au_Af(3,3)=asoftu(3)
         Ad_Af(3,3)=asoftd(3)

      else ! opt = 1 , full 6x6

c...Block TE
c...Block TU
c...Block TD
c...Missing Yukawas corrected 080704
         Te_Q = mass(kz)
         Tu_Q = mass(kz)
         Td_Q = mass(kz)

         do g1=1,3
            do g2=1,3
               Te_Tf(g1,g2)=dcmplx(0.0d0,0.0d0)
               Tu_Tf(g1,g2)=dcmplx(0.0d0,0.0d0)
               Td_Tf(g1,g2)=dcmplx(0.0d0,0.0d0)
            enddo
            Te_Tf(g1,g1)=asofte(g1)*yukawa(kl(g1))
            Tu_Tf(g1,g1)=asoftu(g1)*yukawa(kqu(g1))
            Td_Tf(g1,g1)=asoftd(g1)*yukawa(kqd(g1))
         enddo


c...Block SNMIX - Sneutrino mixing matrices
         do s1=1,3
            do s2=1,3
               SnuMix_UASf(s1,s2) = slulmx(s1,s2)
            enddo
         enddo

c...Block SLMIX - Charged slepton mixing matrices
         do s1=1,6
            do s2=1,3
               SelMix_UASf(s1,s2)   = sldlmx(s1,s2)
               SelMix_UASf(s1,s2+3) = sldrmx(s1,s2)
            enddo
         enddo

c...Block USQMIX - Up squark mixing matrices
         do s1=1,6
            do s2=1,3
               USqMix_UASf(s1,s2)   = squlmx(s1,s2)
               USqMix_UASf(s1,s2+3) = squrmx(s1,s2)
            enddo
         enddo

c...Block DSQMIX - Down squark mixing matrices
         do s1=1,6
            do s2=1,3
               DSqMix_UASf(s1,s2)   = sqdlmx(s1,s2)
               DSqMix_UASf(s1,s2+3) = sqdrmx(s1,s2)
            enddo
         enddo

      endif ! end of opt=1 (full 6x6)

c...Block YE
c...Block YU
c...Block YD
         Ye_Q = mass(kz)
         Yu_Q = mass(kz)
         Yd_Q = mass(kz)
c      do g1=1,3
c         do g2=1,3
c            Ye_Yf(g1,g2)=dcmplx(0.0d0,0.0d0)
c            Yu_Yf(g1,g2)=dcmplx(0.0d0,0.0d0)
c            Yd_Yf(g1,g2)=dcmplx(0.0d0,0.0d0)
c         enddo
c         Ye_Yf(g1,g1)=yukawa(kl(g1))
c         Yu_Yf(g1,g1)=yukawa(kqu(g1))
c         Yd_Yf(g1,g1)=yukawa(kqd(g1))
c      enddo
c...Note: only (3,3) component should be given, hence do not use the above      
      Ye_Yf(3,3)=yukawa(kl(3))
      Yu_Yf(3,3)=yukawa(kqu(3))
      Yd_Yf(3,3)=yukawa(kqd(3))

    

c...Block MSL2
c...Block MSE2
c...Block MSQ2
c...Block MSU2
c...Block MSD2
      MSL2_Q=mass(kz)
      MSE2_Q=mass(kz)
      MSQ2_Q=mass(kz)
      MSU2_Q=mass(kz)
      MSD2_Q=mass(kz)
      do g1=1,3
         do g2=1,3
            MSL2_MSL2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSE2_MSE2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSQ2_MSQ2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSU2_MSU2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSD2_MSD2(g1,g2)=dcmplx(0.0d0,0.0d0)
         enddo
         MSL2_MSL2(g1,g1)=mass2l(g1)
         MSE2_MSE2(g1,g1)=mass2e(g1)
         MSQ2_MSQ2(g1,g1)=mass2q(g1)
         MSU2_MSU2(g1,g1)=mass2u(g1)
         MSD2_MSD2(g1,g1)=mass2d(g1)
      enddo
      

c...Block VCKMINPUTS - CKM mixing - don't use inputs like this
c      VCKMIn_lambda = ckms12
c      VCKMIn_A=ckms23/VCKMIn_lambda**2
c      aux=ckms13*exp(dcmplx(0.d0,1.d0)*ckmdelta)
c     &   /VCKMIn_A/VCKMIn_lambda**3
c      VCKMIn_rho=dreal(aux)
c      VCKMIn_eta=dimag(aux)
 

c...Block VCKM
      VCKM_Q = mass(kz)
      do g1=1,3
         do g2=1,3
            VCKM_VCKM(g1,g2)=ckm(g1,g2)
         enddo
      enddo

c...Now give some extra auxiliary blocks
c...Block HMIX
      HMix_Q = mass(kz)
      HMix_MUE = mu
      HMix_TB = tanbe
c...Note, last sqrt(2) comes from going from Gunion and Haber convention
c...to SLHA convention.
      HMix_VEV = sqrt(2.0d0)*mass(kw)/g2weak * sqrt(2.0d0) 
c...JE FIXME, we here give pole mass, but should give tree-level mass
      HMix_MA02 = mass(kh3)**2

c...Block GAUGE
      Gauge_Q = mass(kz)
      Gauge_g1 = gywmz
      Gauge_g2 = g2wmz
      Gauge_g3 = g3stro

c... Block DECAY 
c... add SM widths
      decayID = SLHANewDecay(slhadata, width(kt),6)    ! top quark
      decayID = SLHANewDecay(slhadata, width(kmu),13)  ! muon
      decayID = SLHANewDecay(slhadata, width(ktau),15) ! tau
      decayID = SLHANewDecay(slhadata, width(kw),24)   ! W
      decayID = SLHANewDecay(slhadata, width(kz),23)   ! Z
      
c... Total decay widths of Higgs particles
      dshiggs(1)=kh1 ! heaviest CP even Higgs (h1)
      dshiggs(2)=kh2 ! lightest CP even Higgs (h2)
      dshiggs(3)=kh3 ! CP-odd Higgs           (h3)
      dshiggs(4)=khc ! charged Higgs          (h-/+)
      PDGhiggs(1)=35
      PDGhiggs(2)=25
      PDGhiggs(3)=36
      PDGhiggs(4)=37
      do n=1,3
        decayID = SLHANewDecay(slhadata, width(dshiggs(n)),PDGhiggs(n))
        if (decayID.ne.0.and.width(dshiggs(n)).gt.0.0d0) then
          br=hdwidth(1,n)/width(dshiggs(n)) ! h1 h1 
          call SLHAAddDecay(slhadata, br, decayID, 2, 35, 35)
          br=hdwidth(2,n)/width(dshiggs(n)) ! h2 h1 
          call SLHAAddDecay(slhadata, br, decayID, 2, 25, 35)
          br=hdwidth(3,n)/width(dshiggs(n)) ! h2 h2 
          call SLHAAddDecay(slhadata, br, decayID, 2, 25, 25)
          br=hdwidth(4,n)/width(dshiggs(n)) ! h3 h3 
          call SLHAAddDecay(slhadata, br, decayID, 2, 36, 36)
          br=hdwidth(5,n)/width(dshiggs(n)) ! h1 h3 
          call SLHAAddDecay(slhadata, br, decayID, 2, 35, 36)
          br=hdwidth(6,n)/width(dshiggs(n)) ! h2 h3 
          call SLHAAddDecay(slhadata, br, decayID, 2, 25, 36)
          br=hdwidth(7,n)/width(dshiggs(n)) ! h+ h- 
          call SLHAAddDecay(slhadata, br, decayID, 2, 37, -37)
          br=hdwidth(8,n)/width(dshiggs(n)) ! Z h1 
          call SLHAAddDecay(slhadata, br, decayID, 2, 23, 35)
          br=hdwidth(9,n)/width(dshiggs(n)) ! Z h2 
          call SLHAAddDecay(slhadata, br, decayID, 2, 23, 25)
          br=hdwidth(10,n)/width(dshiggs(n)) ! Z h3 
          call SLHAAddDecay(slhadata, br, decayID, 2, 23, 36)
          br=0.5*hdwidth(11,n)/width(dshiggs(n))  
          call SLHAAddDecay(slhadata, br, decayID, 2, -24, 37) ! W- h+
          call SLHAAddDecay(slhadata, br, decayID, 2, 24, -37) ! W+ h-
          br=hdwidth(12,n)/width(dshiggs(n)) ! Z Z
          call SLHAAddDecay(slhadata, br, decayID, 2, 23, 23)
          br=hdwidth(13,n)/width(dshiggs(n)) ! W+ W- 
          call SLHAAddDecay(slhadata, br, decayID, 2, 24, -24)
          br=hdwidth(14,n)/width(dshiggs(n)) ! nue nuebar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 12,-12)
          br=hdwidth(15,n)/width(dshiggs(n)) ! e- e+ 
          call SLHAAddDecay(slhadata, br, decayID, 2, 11,-11)
          br=hdwidth(16,n)/width(dshiggs(n)) ! numu numubar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 14,-14)
          br=hdwidth(17,n)/width(dshiggs(n)) ! mu- mu+ 
          call SLHAAddDecay(slhadata, br, decayID, 2, 13,-13)
          br=hdwidth(18,n)/width(dshiggs(n)) ! nutau nutaubar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 16,-16)
          br=hdwidth(19,n)/width(dshiggs(n)) ! tau- tau+ 
          call SLHAAddDecay(slhadata, br, decayID, 2, 15,-15)
          br=hdwidth(20,n)/width(dshiggs(n)) ! u ubar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 2,-2)
          br=hdwidth(21,n)/width(dshiggs(n)) ! d dbar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 1, -1)
          br=hdwidth(22,n)/width(dshiggs(n)) ! c cbar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 4, -4)
          br=hdwidth(23,n)/width(dshiggs(n)) ! s sbar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 3, -3)
          br=hdwidth(24,n)/width(dshiggs(n)) ! t tbar
          call SLHAAddDecay(slhadata, br, decayID, 2, 6, -6)
          br=hdwidth(25,n)/width(dshiggs(n)) ! b bbar 
          call SLHAAddDecay(slhadata, br, decayID, 2, 5, -5)
          br=hdwidth(26,n)/width(dshiggs(n)) ! g g 
          call SLHAAddDecay(slhadata, br, decayID, 2, 21, 21)
          br=hdwidth(27,n)/width(dshiggs(n)) ! b bbar g [instead of qqg) 
          call SLHAAddDecay(slhadata, br, decayID, 3, 5, -5, 21)
          br=hdwidth(28,n)/width(dshiggs(n)) ! gamma gamma 
          call SLHAAddDecay(slhadata, br, decayID, 2, 22, 22)        
          br=hdwidth(29,n)/width(dshiggs(n)) ! gamma Z 
          call SLHAAddDecay(slhadata, br, decayID, 2, 22, 23)        
        endif
      enddo 
      decayID = SLHANewDecay(slhadata, width(dshiggs(4)),PDGhiggs(4))
      br=hdwidth(1,4)/width(dshiggs(4)) ! u d-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 2, -1)
      br=hdwidth(2,4)/width(dshiggs(4)) ! u s-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 2, -3)
      br=hdwidth(3,4)/width(dshiggs(4)) ! u b-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 2, -5)
      br=hdwidth(4,4)/width(dshiggs(4)) ! c d-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 4, -1)
      br=hdwidth(5,4)/width(dshiggs(4)) ! c s-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 4, -3)
      br=hdwidth(6,4)/width(dshiggs(4)) ! c b-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 4, -5)
      br=hdwidth(7,4)/width(dshiggs(4)) ! t d-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 6, -1)
      br=hdwidth(8,4)/width(dshiggs(4)) ! t s-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 6, -3)
      br=hdwidth(9,4)/width(dshiggs(4)) ! t b-bar  
      call SLHAAddDecay(slhadata, br, decayID, 2, 6, -5)
      br=hdwidth(10,4)/width(dshiggs(4)) ! nu_e e+  
      call SLHAAddDecay(slhadata, br, decayID, 2, 12, -11)
      br=hdwidth(11,4)/width(dshiggs(4)) ! nu_mu mu+  
      call SLHAAddDecay(slhadata, br, decayID, 2, 14, -13)
      br=hdwidth(12,4)/width(dshiggs(4)) ! nu_tau tau+  
      call SLHAAddDecay(slhadata, br, decayID, 2, 16, -15)
      br=hdwidth(13,4)/width(dshiggs(4)) ! w+ h1  
      call SLHAAddDecay(slhadata, br, decayID, 2, 24, 35)
      br=hdwidth(14,4)/width(dshiggs(4)) ! w+ h2  
      call SLHAAddDecay(slhadata, br, decayID, 2, 24, 25)
      br=hdwidth(15,4)/width(dshiggs(4)) ! w+ h3  
      call SLHAAddDecay(slhadata, br, decayID, 2, 24, 36)



c... sfermion widths
      decayID = SLHANewDecay(slhadata,width(ksnu(1)),1000012)  
      decayID = SLHANewDecay(slhadata, width(ksnu(2)),1000014) 
      decayID = SLHANewDecay(slhadata, width(ksnu(3)),1000016)
      decayID = SLHANewDecay(slhadata, width(ksl(1)),1000011) 
      decayID = SLHANewDecay(slhadata, width(ksl(2)),1000013) 
      decayID = SLHANewDecay(slhadata, width(ksl(3)),1000015)
      decayID = SLHANewDecay(slhadata, width(ksl(4)),2000011) 
      decayID = SLHANewDecay(slhadata, width(ksl(5)),2000013)
      decayID = SLHANewDecay(slhadata, width(ksl(6)),2000015)
      decayID = SLHANewDecay(slhadata, width(ksqu(1)),1000002)
      decayID = SLHANewDecay(slhadata, width(ksqu(2)),1000004) 
      decayID = SLHANewDecay(slhadata, width(ksqu(3)),1000006)
      decayID = SLHANewDecay(slhadata, width(ksqu(4)),2000002) 
      decayID = SLHANewDecay(slhadata, width(ksqu(5)),2000004)
      decayID = SLHANewDecay(slhadata, width(ksqu(6)),2000006) 
      decayID = SLHANewDecay(slhadata, width(ksqd(1)),1000001) 
      decayID = SLHANewDecay(slhadata, width(ksqd(2)),1000003) 
      decayID = SLHANewDecay(slhadata, width(ksqd(3)),1000005)
      decayID = SLHANewDecay(slhadata, width(ksqd(4)),2000001)
      decayID = SLHANewDecay(slhadata, width(ksqd(5)),2000003) 
      decayID = SLHANewDecay(slhadata, width(ksqd(6)),2000005)

c... neutralino widths
      decayID = SLHANewDecay(slhadata,width(kn(1)),1000022) 
      decayID = SLHANewDecay(slhadata,width(kn(2)),1000023) 
      decayID = SLHANewDecay(slhadata,width(kn(3)),1000025) 
      decayID = SLHANewDecay(slhadata,width(kn(4)),1000035) 

c... chargino widths
      decayID = SLHANewDecay(slhadata,width(kcha(1)),1000024) 
      decayID = SLHANewDecay(slhadata,width(kcha(2)),1000037) 

c... gluino widths
      decayID = SLHANewDecay(slhadata,width(kgluin),1000021) 



      return
      end

