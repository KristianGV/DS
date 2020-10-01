      subroutine dsacbnd13(warning)
c_______________________________________________________________________
c  check if accelerator data warningude the present point.
c  This routine uses simplified models and approximate constraints.      
c  output:
c    warning - code of the reason for warningusion (integer); 0 if allowed
c    if not allowed, the reasons are coded as follows
c      bit set   dec.   oct.   reason
c      -------   ----   ----   ------
c            0      1      1   chargino mass
c            1      2      2   gluino mass
c            2      4      4   squark mass
c            3      8     10   slepton mass
c            4     16     20   invisible z width
c            5     32     40   higgs mass in excluded region
c            6     64    100   neutralino mass
c            7    128    200   b -> s gamma
c            8    256    400   rho parameter
c            9    512   1000   (g-2)_mu  
c           10   1024   2000   B_s -> mu+ mu-
c           11   2048   4000   squark-gluino
c           12   4096  10000   Higgs mass does not fit observed Higgs        
c  author: paolo gondolo 1994-1999
c  history:
c     940407 first version paolo gondolo
c     950323 update paolo gondolo
c     971200 partial update joakim edsjo
c     980428 update joakim edsjo
c     990719 update paolo gondolo
c     000310 update piero ullio
c     000424 added delrho joakim edsjo
c     000904 update according to pdg2000 lars bergstrom      
c     010214 mh2 limits corrected, joakim edsjo
c     020927 higgs limits update according to pdg2002 mia schelke
c     021001 susy part. mass limits update to pdg2002 je/ms
c     031204 standard model higgs like mh2 limit for msugra models
c     070529 calling new bsg routine+new experim bsg value, ms
c     090105 corrected do-loop for mh3 mass bounds
c     090316 higgs limits by HiggsBounds, erik lundstrom
c     110905 new higgs limits with HiggsBounds, paolo gondolo
c     130412 Added B_s -> mu+ mu- (calculated with SuperIso), Joakim Edsjo
c     160209 Added LHC 8 TeV mgluino msquark bounds, lars bergstrom
c     171212 Updated b -> s\gamma, B_s -> mu+ mu-, (g-2), torsten bringmann
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsaccom.h'
      include 'dsmpconst.h'
      integer i,j,warning,error
      logical warninguded
      real*8 temp,mi,mj,ei,ej,p2,mz,mz2,dsabsq
      real*8 dsgm2muon,gm2mu,bsg_MSSM,
     &    bsmm_mssm,bsmm_untag_mssm,delta0_mssm,
     &    brbtaunu_mssm,brbdtaunu_mssm,amu_mssm,rmu23_mssm
c      real*8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6 ! JE TMP
c      integer itmp ! JE TMP
      complex*16 gzij
      real*8 diffneucha,absmacha,arrmh3min(6),arrh3mintanbe(6),
     & arrmh3(18),arrh3lowtanbe(18),arrmh3cdf(10),arrh3hightanbe(10),
     & arrmh2(14),arrh2lowtanbe(14),
     & arrmh2cdf(10),arrh2hightanbe(10),ylim
       real*8 msq,mgl,mchi,boundx(0:10),boundy(0:10)
       real*8 test
       integer npoints,ii
       integer iiikkk
c...Declaration of HiggsBounds variables (el 090316)
       integer HBresult
       real*8 HSpvalue
c...(PG: no longer needed with HiggsBounds 3.4.0)
c      integer anH,achan,ancombined,anprodxs
c      parameter (anH=3,anprodxs=49)
c      character*5 awhichexpt
c      parameter (awhichexpt="onlyL")
c      double precision aMh(anH)
c      double precision aCS_lep_hjZ_ratio(anH)
c      double precision aCS_lep_hjhi_ratio(anH,anH)
c      double precision aCS_tev_pp_hj_ratio(anH)
c      double precision aCS_tev_pp_hjb_ratio(anH)
c      double precision aCS_tev_pp_hjW_ratio(anH)
c      double precision aCS_tev_pp_hjZ_ratio(anH)
c      double precision aCS_tev_pp_vbf_ratio(anH)
c      double precision aBR_hjbb(anH),aBR_hjtautau(anH)
c      double precision aBR_hjWW(anH),aBR_hjgaga(anH)
c      double precision aBR_hjhihi(anH,anH)
c      double precision aobsratio
c      double precision SMBR_Hbb,SMBR_Htautau
c      double precision SMBR_HWW,SMBR_Hgamgam
c...


      data (arrmh3min(i),i=1,6)/3.8800d+02,3.9300d+02,4.0800d+02,
     & 4.2000d+02,4.6900d+02,5.0000d+02/

      data (arrh3mintanbe(i),i=1,6)/0.400d0,0.420d0,0.440d0,
     & 0.450d0,0.479d0,0.480d0/

      data (arrmh3(i),i=1,18)/0.9190d+02,0.9189d+02,0.9188d+02,
     & 0.9300d+02,
     & 0.9800d+02,1.0300d+02,1.0500d+02,1.1000d+02,1.1500d+02,
     & 1.2500d+02,1.3900d+02,1.5400d+02,1.6900d+02,1.9500d+02,
     & 2.3600d+02,2.8600d+02,3.5900d+02,5.0000d+02/

      data (arrh3lowtanbe(i),i=1,18)/47.500d0,30.000d0,15.240d0,
     & 10.000d0,
     & 6.660d0,7.040d0,8.000d0,8.570d0,8.000d0,   
     & 6.120d0,4.820d0,4.070d0,3.540d0,3.210d0,
     & 2.910d0,2.670d0,2.490d0,2.400d0/

      data (arrmh3cdf(i),i=1,10)/0.9190d+02,0.9800d+02,1.0100d+02,
     & 1.1000d+02,1.2000d+02,1.3000d+02,1.4000d+02,1.5100d+02,
     & 2.0100d+02,2.4700d+02/

      data (arrh3hightanbe(i),i=1,10)/47.500d0,53.000d0,54.500d0,
     & 59.000d0,68.500d0,76.500d0,87.000d0,83.500d0,
     & 96.000d0,100.000d0/



      data (arrmh2(i),i=1,14)/
     & 0.9000d+02,0.8999d+02,0.9270d+02,0.9230d+02,
     & 0.9450d+02,0.9680d+02,1.0000d+02,1.0320d+02,
     & 1.0820d+02,1.1090d+02,1.1180d+02,1.1319d+02,
     & 1.1320d+02,1.1321d+02/

      data (arrh2lowtanbe(i),i=1,14)/
     & 45.500d0,30.001d0,30.000d0,9.000d0,
     &  7.000d0,6.660d0,7.450d0,8.330d0,
     &  8.690d0,8.000d0,7.000d0,5.000d0,
     &  2.400d0,0.500d0/

      data (arrmh2cdf(i),i=1,10)/0.9000d+02,0.9720d+02,1.0060d+02,
     & 1.1500d+02,1.1830d+02,1.2000d+02,1.2220d+02,1.2380d+02,
     & 1.2500d+02,1.2560d+02/

      data (arrh2hightanbe(i),i=1,10)/45.500d0,53.000d0,55.500d0,
     & 61.500d0,65.000d0,70.000d0,83.500d0,90.000d0,
     & 91.000d0,100.000d0/

c...Zero warningusion flag      
      warning=0

c----------------------------------------------bounds from b -> s gamma
c...Heavy Flavour Averaging Group 
c...http://www.slac.stanford.edu/xorg/hfag/ or arXiv:0704:3575v1
c...Barberio et al.:
c...The experimental world average
c...BR(b->s gamma) = (3.55 +- 0.24 +0.09-0.1 +- 0.03) e-4
c...               = (3.55 +- 0.26) e-4
c...For the theoretical error:
c...SM calculation, Misiak et al (hep-ph/0609232 and hep-ph/0609241)
c...+ M. Misiak and M. Steinhauser, private communication:
c...SM error: +- 0.23 e-4
c...Let's assume theoretical SUSY error is also +- 0.23 e-4
c...The total error (experiment+theory) 0.42 e-4. 
c...The two sigma confidence level is then (2.71 -> 4.39 )e-4 

c...Obsolete b->s gamma calculation      
c      call dsbsg2007full(bsg,1)  ! bsg w/ SUSY contribution

c...Call SuperIso to get flavour observables      
      call dsflavour(bsg_MSSM,delta0_MSSM,bsmm_MSSM,
     &  bsmm_untag_MSSM,
     &     brbtaunu_mssm,brbdtaunu_mssm,amu_mssm,rmu23_mssm,error)
c      write(*,*) 'bsmm_mssm, bsg_MSSM = ',bsmm_mssm,bsg_MSSM 
      bsg=bsg_MSSM
      bsgsi=bsg_MSSM ! SuperIso b->s gamma calculation JE FIX: remove
      delta0=delta0_MSSM
      bsmumu=bsmm_mssm
      bsmumu_untag=bsmm_untag_mssm
      brbtaunu=brbtaunu_mssm
      brbdtaunu=brbdtaunu_mssm
      amusi=amu_mssm
      rmu23=rmu23_mssm
c      call FHFlavour(itmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6) ! JE TMP
c      write(55,*) bsg,tmp1,tmp2 ! JE TMP
c      if (bsg.lt.2.71d-4.or.bsg.gt.4.39d-4) then
c         warning=ibset(warning,7)
c      endif

c... TB 2017-12-19: update to 2s limits from 1705.07933, based on BarBar & Belle
      if (bsg.lt.2.99d-4.or.bsg.gt.3.65d-4) then
         warning=ibset(warning,7)
      endif
 
 
      
c--------------------------------------------bounds from B_s -> mu+ mu-
c...We get B_s -> mu+ mu- from SuperIso. One could also get it from
c...FeynHiggs, but we don't use that option. The calling sequence is
c...kept for reference though.      
c
c...Obsolete FH call.
c...We need to set up FH for the calculation unless it is already
c...done in the spectrum calculation
c      call FHFlavour(error,
c     &    bsg_MSSM, bsg_SM,
c     &    deltaM_sMSSM, deltaM_sSM,
c     &    bsmm_MSSM, bsmm_SM)
c      if (error.ne.0) then 
c        write(*,*) 'ERROR from FeynHiggs (called from dsacbnd11)'
c      endif


c...SuperIso already called earlier      
c...Limit comes from combined limit given in 
c...LHCb-CONF-2012-017
c      if (bsmumu.gt.4.2d-9) then 
c         warning=ibset(warning,10)
c      endif

c... TB update: LHCb measurement, 2s combined theory and exp. error.
      if (bsmumu.gt.4.3d-9.or.bsmumu.lt.1.7d-9) then 
         warning=ibset(warning,10)
      endif



c------------------------------------------lower bound on chargino mass
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
c...the following limits still hold 
c  ... lbe update 000904 to pdg2000    
      absmacha=min(mass(kcha(1)),mass(kcha(2)))
      diffneucha=absmacha-mass(kn(kln))
c ...  lbe, pdg2000, l3 limit:
      if (absmacha.lt.67.7d0) warning=ibset(warning,0) 
c ... the following is from pdg2000, abreu et al (delphi):      
      if (absmacha.lt.88.4d0.and.diffneucha.ge.3.d0
     & .and.tanbe.ge.1.d0.and.mass(ksnu(1)).ge.absmacha)
     &  warning=ibset(warning,0)
c        pdg 94 p 1801 iv hidaka91 cdf 
c      if (abs(mass(kcha(1))).lt.99.0d0) warning=ibset(warning,0)

c--------------------------------------------lower bound on gluino mass
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
      if (mass(kgluin).lt.195.d0) warning=ibset(warning,1)
c  ... lbe update 000904 to pdg2000; use only limit with cascade decays:
c      if (mass(kgluin).lt.173.d0) warning=ibset(warning,1)
c      warninguded=.false.
c      do i=1,6
c        warninguded = warninguded .or. (mass(ksqu(i)).lt.mass(kgluin)) .or.
c     &        (mass(ksqd(i)).lt.mass(kgluin))
c      enddo
c      if (warninguded .and. mass(kgluin).lt.260.0d0) warning=ibset(warning,1)

c------------------------------------------lower bound on squark masses
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)

      warninguded=.false.
      do i=1,2 ! DO limit pdg2002 
        warninguded = warninguded .or.
     &    mass(ksqu(i)).lt.250.d0 .or.
     &    mass(ksqd(i)).lt.250.d0        !corrected u->d, ms 021001
      enddo      
      do i=1,2 ! DO limit pdg2002 
        warninguded = warninguded .or.
     &    mass(ksqu(i+3)).lt.250.d0 .or.
     &    mass(ksqd(i+3)).lt.250.d0      !corrected u->d, ms 021001
      enddo                                                 
      if (warninguded) warning=ibset(warning,2)
c...je/ms update 02-10-01 news:stop/sbottom limits included
      if (mass(ksqd_flav(3,1)).lt.91.0d0
     &  .and.(mass(ksqd_flav(3,1))-mass(kn(kln))).gt.8.0d0)
     &  warning=ibset(warning,2)
      if (mass(ksqd_flav(3,2)).lt.91.0d0
     &  .and.(mass(ksqd_flav(3,2))-mass(kn(kln))).gt.8.0d0)
     &  warning=ibset(warning,2)
      if (mass(ksqu_flav(3,1)).lt.86.4d0
     &  .and.(mass(ksqu_flav(3,1))-mass(kn(kln))).gt.5.0d0)
     &  warning=ibset(warning,2)
      if (mass(ksqu_flav(3,2)).lt.86.4d0
     &  .and.(mass(ksqu_flav(3,2))-mass(kn(kln))).gt.5.0d0)
     &  warning=ibset(warning,2)
      warninguded=.false.
      do i=1,6 ! lep, pdg2000
        warninguded = warninguded .or. (mass(ksl(i)).lt.41.d0)
      enddo
c...021001update:ksl numbers corrected,limits+massdiff.updated,ALEP added 
c ... opal limits, pdg2002:        
      if (mass(ksl(1)).lt.87.1d0
     &     .and.(mass(ksl(1))-mass(kn(kln))).gt.5.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(4)).lt.87.1d0
     &     .and.(mass(ksl(4))-mass(kn(kln))).gt.5.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(2)).lt.82.3d0
     &     .and.(mass(ksl(2))-mass(kn(kln))).gt.3.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(5)).lt.82.3d0
     &     .and.(mass(ksl(5))-mass(kn(kln))).gt.3.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(3)).lt.73.0d0
     &     .and.(mass(ksl(3))-mass(kn(kln))).gt.10.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(6)).lt.73.0d0
     &     .and.(mass(ksl(6))-mass(kn(kln))).gt.10.0d0)
     &     warning=ibset(warning,3)
c ... alep limits, pdg2002:        
      if (mass(ksl(1)).lt.95.0d0
     &   .and.(mass(ksl(1))-mass(kn(kln))).gt.15.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(4)).lt.95.0d0
     &     .and.(mass(ksl(4))-mass(kn(kln))).gt.15.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(2)).lt.88.0d0
     &     .and.(mass(ksl(2))-mass(kn(kln))).gt.15.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(5)).lt.88.0d0
     &     .and.(mass(ksl(5))-mass(kn(kln))).gt.15.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(3)).lt.76.0d0
     &     .and.(mass(ksl(3))-mass(kn(kln))).gt.15.0d0)
     &     warning=ibset(warning,3)
      if (mass(ksl(6)).lt.76.0d0
     &     .and.(mass(ksl(6))-mass(kn(kln))).gt.15.0d0)
     &     warning=ibset(warning,3)
c---------------------------------------lower bound on sneutrino masses
c...je/ms update 021001 to pdg2002(see ref above);the following still holds
c ... lbe update 000904 to pdg2000 "unpublished lep limits":
      if (mass(ksnu(1)).lt.43.7d0) warning=ibset(warning,3)
      if (mass(ksnu(2)).lt.43.7d0) warning=ibset(warning,3)
      if (mass(ksnu(3)).lt.43.7d0) warning=ibset(warning,3)

c-------------------------------------lower bounds on neutralino masses
c     assumes neutralinos are ordered by increasing mass
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
      if (abs(mass(kn(1))).lt.37.0d0) warning=ibset(warning,6)
      if (abs(mass(kn(2))).lt.62.4d0.and.tanbe.ge.1.0d0
     &    .and.tanbe.le.40.0d0)
     &    warning=ibset(warning,6)
      if (abs(mass(kn(3))).lt.99.9d0.and.tanbe.ge.1.0d0
     &    .and.tanbe.le.40.0d0)
     &    warning=ibset(warning,6)
      if (abs(mass(kn(4))).lt.116.0d0.and.tanbe.ge.1.0d0
     &    .and.tanbe.le.40.0d0)
     &    warning=ibset(warning,6)

c------------------------------------------------bounds on higgs masses
c... el update 09-03-16
c... for default settings (higwid=5, in dsinit) higgs limits
c....are extracted from HiggsBounds.
c... if higwid=1, darksusy internal limits (i.e. the same 
c... as in dsacbnd7) are instead used.

c...Using limits from various literature (same as in dsacbnd7)
c...for higwid=1 (or mSUGRA models)
      if (abs(higwid).eq.1.or.modeltype.eq.3) then

c h3:
c ... lbe update, pdg2000 delphi limit:
c      if(mass(kh3).lt.84.1d0) warning=ibset(warning,5)

c...ms update 02-09-27
c...using limits stated for the maximal mixing scenario
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings
c...see sec. 4 and fig 7 in this
c...(citation:Phys.Rev.D66(2002)010001)
c...for more detailed figures see
c...low tanbe:LEP preliminary limits LHWG Note 2001-04 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
c...we use fig 4_3
c...high tanbe:CDF limits, Phys.Rev.Lett. 86(2001)4472
c...(hep-ex/0010052) we use fig 3 (p.18)

      if(tanbe.lt.arrh3mintanbe(6)) then
       if(mass(kh3).lt.arrmh3min(1)) then
         warning=ibset(warning,5)
         goto 110
       elseif(mass(kh3).gt.arrmh3min(6)) then
         goto 110
       else
         ylim=0.0d0
         do iiikkk=1,6
           if(mass(kh3).lt.arrmh3min(iiikkk+1).and
     &          .mass(kh3).ge.arrmh3min(iiikkk)) then
             ylim=arrh3mintanbe(iiikkk)+
     &          (arrh3mintanbe(iiikkk+1)-arrh3mintanbe(iiikkk))
     &          *(mass(kh3)-arrmh3min(iiikkk))
     &          /(arrmh3min(iiikkk+1)-arrmh3min(iiikkk))
             goto 111
           endif
         enddo
 111     if(tanbe.gt.ylim) warning=ibset(warning,5)
       endif
      endif
 110  continue
         
      if(tanbe.ge.arrh3mintanbe(6).and.tanbe.le.arrh3lowtanbe(18)) then
         warning=ibset(warning,5)
      endif
           
      if(tanbe.lt.arrh3lowtanbe(1).and.tanbe.gt.arrh3lowtanbe(18)) then
       if(mass(kh3).lt.arrmh3(1)) then
         warning=ibset(warning,5)
         goto 220
       elseif(mass(kh3).gt.arrmh3(18)) then
         goto 220
       else
         ylim=0.0d0
         do iiikkk=1,17
           if(mass(kh3).lt.arrmh3(iiikkk+1).and
     &          .mass(kh3).ge.arrmh3(iiikkk)) then
             ylim=arrh3lowtanbe(iiikkk)+
     &          (arrh3lowtanbe(iiikkk+1)-arrh3lowtanbe(iiikkk))
     &          *(mass(kh3)-arrmh3(iiikkk))
     &          /(arrmh3(iiikkk+1)-arrmh3(iiikkk))
             goto 222
           endif
         enddo
 222     if(tanbe.lt.ylim) warning=ibset(warning,5)
       endif
      endif
 220  continue
         
       
      if(tanbe.ge.arrh3hightanbe(1)
     &   .and.tanbe.lt.arrh3hightanbe(10)) then
       if(mass(kh3).lt.arrmh3cdf(1)) then
         warning=ibset(warning,5)
         goto 330
       elseif(mass(kh3).gt.arrmh3cdf(10)) then
         goto 330
       else
         ylim=0.0d0  
         do iiikkk=1,10
           if(mass(kh3).lt.arrmh3cdf(iiikkk+1).and
     &          .mass(kh3).ge.arrmh3cdf(iiikkk)) then
             ylim=arrh3hightanbe(iiikkk)+
     &          (arrh3hightanbe(iiikkk+1)-arrh3hightanbe(iiikkk))
     &          *(mass(kh3)-arrmh3cdf(iiikkk))
     &          /(arrmh3cdf(iiikkk+1)-arrmh3cdf(iiikkk))
             goto 333
           endif
         enddo
 333     if(tanbe.gt.ylim) warning=ibset(warning,5)
       endif
      endif
 330  continue
         
       
c h2:
c...pu update 00-03-10 limit implemented from aleph 2000-006
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
c...Corrected, JE 2001-02-14
             
c...ms update 02-09-27
c...using limits stated for the maximal mixing scenario
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings
c...see sec 4 and fig 7 in this (use this fig for tanbe~30)
c...(citation:Phys.Rev.D66(2002)010001)
c...for more detailed figures see
c...for low tanbe:LEP preliminary limits LHWG Note 2001-04 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
c...we use fig 4_2
c...for high tanbe:CDF limits, Phys.Rev.Lett. 86(2001)4472
c...(hep-ex/0010052) we use fig 2 (p.17)



      if(tanbe.lt.arrh2lowtanbe(1)) then
       if(mass(kh2).lt.arrmh2(1)) then
         warning=ibset(warning,5)
         goto 440
       elseif(mass(kh2).gt.arrmh2(14)) then
         goto 440
       else
         ylim=0.0d0
         do iiikkk=1,14
           if(mass(kh2).lt.arrmh2(iiikkk+1).and
     &          .mass(kh2).ge.arrmh2(iiikkk)) then
             ylim=arrh2lowtanbe(iiikkk)+
     &          (arrh2lowtanbe(iiikkk+1)-arrh2lowtanbe(iiikkk))
     &          *(mass(kh2)-arrmh2(iiikkk))
     &          /(arrmh2(iiikkk+1)-arrmh2(iiikkk))
             goto 444
           endif
         enddo
 444     if(tanbe.lt.ylim) warning=ibset(warning,5)
       endif
      endif
 440  continue
         
       
      if(tanbe.ge.arrh2hightanbe(1)
     &   .and.tanbe.lt.arrh2hightanbe(10)) then
       if(mass(kh2).lt.arrmh2cdf(1)) then
         warning=ibset(warning,5)
         goto 550
       elseif(mass(kh2).gt.arrmh2cdf(10)) then
         goto 550
       else
         ylim=0.0d0  
         do iiikkk=1,10
           if(mass(kh2).lt.arrmh2cdf(iiikkk+1).and
     &          .mass(kh2).ge.arrmh2cdf(iiikkk)) then
             ylim=arrh2hightanbe(iiikkk)+
     &          (arrh2hightanbe(iiikkk+1)-arrh2hightanbe(iiikkk))
     &          *(mass(kh2)-arrmh2cdf(iiikkk))
     &          /(arrmh2cdf(iiikkk+1)-arrmh2cdf(iiikkk))
             goto 555
           endif
         enddo
 555     if(tanbe.gt.ylim) warning=ibset(warning,5)
       endif
      endif
 550  continue

c change this when we set a global variable to say whether we are in the
c Mssm or in mSUGRA
c limit on a standard model like h2, standard model limit from
c CERN-EP-2003-011 at: http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
c

      if(dabs(dabs(sgnmuvar)-1.d0).lt.1.d-10) then
        if(mass(kh2).lt.114.4d0) warning=ibset(warning,5)
      endif

      
c hc:
c...ms update 02-09-27
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings
c...see sec 5 (citation:Phys.Rev.D66(2002)010001)
c...for more details see
c...LEP preliminary limits LHWG Note 2001-05 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
      if (mass(khc).lt.78.6d0) warning=ibset(warning,5)

c...pu update 00-03-10. limit implemented from aleph 2000-011
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
c      if (mass(khc).lt.77.7d0) warning=ibset(warning,5)
c ... lbe update, pdg2000, d0 limit:
c      if (mass(khc).lt.82.8d0) warning=ibset(warning,5)


      elseif (higwid.eq.5) then

c...using HiggsBounds
         call dshiggsbounds(HBresult,HSpvalue)
         if (HBresult.ne.1) warning=ibset(warning,5)
         if (HSpvalue.lt.1.d-2) warning=ibset(warning,12) ! 99%

c...limits on charged higgs taken from literature (same as for higwid=1)
c...ms update 02-09-27
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings
c...see sec 5 (citation:Phys.Rev.D66(2002)010001)
c...for more details see
c...LEP preliminary limits LHWG Note 2001-05 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
         if (mass(khc).lt.78.6d0) warning=ibset(warning,5)

c...pu update 00-03-10. limit implemented from aleph 2000-011
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
c      if (mass(khc).lt.77.7d0) warning=ibset(warning,5)
c ... lbe update, pdg2000, d0 limit:
c      if (mass(khc).lt.82.8d0) warning=ibset(warning,5)


c...Invalid higwid option
      else
         write(*,*) 'DS ERROR in dsacbnd12:'
         write(*,*) 'Invalid higwid=',higwid
         write(*,*) 'Stopping...'
         stop
      endif

c-----------------------------------------bound from gamma_z(invisible)
      mz = mass(kz)
      mz2 = mz*mz
      ! 3 neutrinos, corrected 2002-10-01, factor of 3 added
      ! (Thanks to Dan Hooper)
      gzinv = 3.0d0 * (g2wmz/2/cwmz)**2 * mz/(24.d0*pi)
      ! 3 sneutrinos
      do i=1,3
         temp = 1.d0-4.d0*mass(ksnu(i))**2/mz2
         if (temp.gt.0.d0) then
            temp = temp**(1.5d0) *
     &           (g2wmz/2/cwmz)**2 * mz/(48.d0*pi)
            gzinv = gzinv + temp
         endif
      enddo
      ! neutralinos
c...pg corrected june 30, 2000
      do i=1,4
         mi = abs(mass(kn(i)))
         do j=i,4
            mj = abs(mass(kn(j)))
            p2 = (mz2-(mi+mj)**2)*(mz2-(mi-mj)**2)/(4*mz2)
            ei = (mz2-mj**2+mi**2)/(2*mz)
            ej = (mz2-mi**2+mj**2)/(2*mz)
            if (p2.gt.0.d0.and.ei.gt.0.0d0.and.ej.gt.0.0d0) then
               gzij = g2wmz/(2.d0*cwmz)*
     &              (conjg(neunmx(i,3))*neunmx(j,3)-
     &               conjg(neunmx(i,4))*neunmx(j,4))
               temp=sqrt(p2)/(2.d0*pi*mz2) * (
     &              dsabsq(gzij) * (ei*ej+p2/3.0d0)-
     &              dreal(gzij**2)*mi*mj)
               if (i.eq.j) temp=0.5d0*temp
               gzinv = gzinv + temp
            endif
         enddo
      enddo
c        pdg 94 p 1358 (498.2+4.2)
c        pdg 2002 (499 +- 1.5) => 2 sigma limit, 499 +- 3
      if (gzinv.gt.0.502d0) warning=ibset(warning,4)

c--------------------------------------------limit on the rho parameter
c...je addition 2000-04-24
c...uses delta rho as calculated by feynhiggs
      if (delrho.gt.1.3d-3) then
        warning=ibset(warning,8)
      endif


c-------------------- bounds from magnetic moment of the muon, (g-2)_mu
c...We will here use a rather large SM theoretical bound, both using
c...the hadronic and tau estimates
c...e+ e-: Teubner et al, arXiv:1001.5401,
c...a_mu^exp-a_mu^SM=(31.3 +- 7.9) x 10^-10
c...e+ e- (w/Babar), Davier et al, arXiv:1001.2243
c...a_mu^exp-a_mu^SM=(25.5 +- 8.0) x 10^-10
c...tau, Davier et al, arXiv:0906.5443
c...a_mu^exp-a_mu^SM=(15.7 +- 8.2) x 10^-10
c...To be conservative, the 2sigma overlapping limit from these are
c...a_mu^exp-a_mu^SM = (-0.7 --  47.1) x 10^-10
c...To allow for a theoretical uncertainty, let us take
c...a_mu^exp-a_mu^SM = (-2.7 --  49.1) x 10^-10
c...as our limit

      gm2mu=dsgm2muon() ! this is a_mu = (g-2)/2
c      if (gm2mu.lt.-2.7d-10.or.gm2mu.gt.49.1d-10) then
c         warning=ibset(warning,9)
c      endif

c... TB 2017-12-19 update to following roughly the treatment in PrecisionBit
c... (11659180.2 ± 4.9)*1d-10 from e+- data
c... (11659189.4 ± 5.4)*1d-10 from tau+tau-data
c... (11659208.9 ± 6.3)*1d-10 observed
c... ------------------ [2s errors in quad for largest diff]
c... 
      if (gm2mu.lt.-2.9d-10.or.gm2mu.gt.44.7d-10) then
         warning=ibset(warning,9)
      endif


c----------------------------------------
c... Bounds from 8 TeV run at CERN, combined squark-gluino for neutralino
c... masses smaller than 395 GeV, doi:10.1007/JHEP09(2014)176, lbe 160212
c... simple, piecewise linear approx of LHC 8 TeV data, fig. 9
c----------------------------------------
c       write(*,*) 'START OF NEW ACCELERATOR BBOUNDS'
       mchi=mass(kn(1))
       mgl=mass(kgluin)
       msq=mass(ksqd(1)) ! assume all have same mass
       boundx(1)=1150.   ! GeV 
       boundy(1)=2800.
       boundx(2)=1180.
       boundy(2)=2430.
       boundx(3)=1250.
       boundy(3)=2300.
       boundx(4)=1350.
       boundy(4)=2310.
       boundx(5)=1450.
       boundy(5)=2250.
       boundx(6)=1670.
       boundy(6)=1870.
       boundx(7)=1690.
       boundy(7)=1800.
       boundx(8)=1760.
       boundy(8)=1600.
       boundx(9)=2400.
       boundy(9)=1400.
       npoints=9
       test=0.
       if ((mgl.le.1150..or.msq.le.1400.).and.mchi.le.395.) then
            warning=ibset(warning,11)
c        write(*,*) 'warningUDED! '
c        stop         
        goto 20  
       endif
       do 10 ii=1,npoints-1
          if (mgl.lt.boundx(ii+1).and.mgl.gt.boundx(ii)) then 
            test=boundy(ii)+(mgl-boundx(ii))/(boundx(ii+1)-boundx(ii))*
     &         (boundy(ii+1)-boundy(ii))   
          endif
          if (msq.lt.test) then 
            warning=ibset(warning,11)
            goto 20
          endif
 10    continue
 20    continue


c   reason for warningusion
c... should be handled by main program!
      continue
c      if (prtlevel.ge.1) call dswexcl(6,warning)
      return

      end

