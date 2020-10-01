************************************************************************
*** auxiliary function that provides integrand for convolution of 3-body 
*** final state particle distribution and MC fragmentation functions
***
*** Author: Torsten Bringmann, 2012-12-03 (2015-11_19: added Higgs final states)
***         2016-02-05 updated to DS6 conventions
************************************************************************

      real*8 function dsIB2convint(zint) 
      implicit none
      include 'dsib2com.h'
      include 'dsmssm.h'
      include 'dsio.h'

c------------------------ functions ------------------------------------
 
           real*8 dsIB2dsde_aux,dsanyield_ch,dsanyield_sim
      
c------------------------ variables ------------------------------------
 
      real*8 zint
      real*8 mpp, epp, mppVS, mppVS2, mppVV, dsanyieldcorr
      integer istat,ch_ha, chpdg, pdgyield, diff
      
c-----------------------------------------------------------------------

      dsIB2convint=0.d0
      ch_ha=0

      if (c_xstable.gt.zint) return


      pdgyield=0
      diff=c_yieldk/100   ! diff=0: integrated; diff=1: differential
      if (mod(c_yieldk,100).eq.51) pdgyield=-11           ! positrons
      if (mod(c_yieldk,100).eq.52) pdgyield=22            ! cont. gamma rays
      if (mod(c_yieldk,100).eq.54) pdgyield=-2212         ! antiprotons
      if (mod(c_yieldk,100).eq.56) pdgyield=12            ! nu_e and nu_e-bar
      if (mod(c_yieldk,100).eq.71.or.mod(c_yieldk,100).eq.53) pdgyield=14      ! nu_mu and nu_mu-bar
      if (mod(c_yieldk,100).eq.57) pdgyield=16            ! nu_tau and nu_tau-bar
      if (mod(c_yieldk,100).eq.72) pdgyield=130072        ! muons from nu at creation
      if (mod(c_yieldk,100).eq.73) pdgyield=130073        ! muons from nu, as seen by a detector in ice
                                                        ! (i.e. integrating 130072 over the mean muon path)
      if (pdgyield.eq.0.or.diff.lt.0.or.diff.gt.1) then
        if (prtlevel.ge.2) then
          write(*,*) 'WARNING in dsib2yieldone (dsIB2convint): called with unsopported ',
     &               'parameter c_yieldk = ',c_yieldk
          write(*,*) 'Returning with zero yield...'
        endif  
        return
      endif

      chpdg = 0
c... determine channel for dsanyield routines
      if (c_pfinal.eq.'B') then
        if (c_Btype.eq.1) then
          ch_ha = 12   ! ZZ
          chpdg = 23
        elseif (c_Btype.eq.2) then
          ch_ha = 13   ! WW
          chpdg = 24
        endif  
        if (c_Btype.eq.3) ch_ha = 3    ! hh
        if (c_Btype.eq.4) ch_ha = 1    ! HH
        if (c_Btype.eq.5) ch_ha = 4    ! AA
        if (c_Btype.eq.6) ch_ha = 7    ! H+H-
      elseif (c_pfinal.eq.'f') then
        ch_ha=13+c_ftype
        chpdg=c_ftype - 7 + 2*mod(c_fbartype,2)
      elseif (c_pfinal.eq.'fbar') then
        ch_ha=13+c_fbartype      
        chpdg=-c_fbartype + 7 - 2*mod(c_fbartype,2)
      endif

      if (ib2dnhow.eq.'full') then
        mpp=zint*mx
        epp=c_xstable*mx

        if (ch_ha.ge.8.and.ch_ha.le.11) then
           write(*,*) 'ERROR in dsIB2convint: ch_ha = ',ch_ha
        elseif (ch_ha.gt.11) then
          dsIB2convint = 0.5d0*
     &                   dsIB2dsde_aux(zint)*
     &                   dsanyield_sim(mpp,epp,chpdg,0,pdgyield,diff,istat)
        else

c... temporary FIX: dsanyield_sim returns zero for 2 identical Higgs final states.
c... Approximate yield by Z(W)+Higgs, and subtract yield from Z (W) = ZZ (WW)/2. 
c... FIXME: this can be replaced by calls to dshayields... !  
          if (ch_ha.eq.1) then
            mppVV = sqrt(mpp**2+mass(kz)**2-mass(kh1)**2)
            mppVS = (mppVV + mpp)/2.d0 
            dsanyieldcorr=dsanyield_ch(mppVS,epp,23,35,c_yieldk,istat)       ! Zh channel
     &                    -0.5*dsanyield_sim(mppVV,epp,23,0,pdgyield,diff,istat)  ! ZZ channel
          elseif (ch_ha.eq.3) then
            mppVV = sqrt(mpp**2+mass(kz)**2-mass(kh2)**2)
            mppVS = (mppVV + mpp)/2.d0 
            dsanyieldcorr=dsanyield_ch(mppVS,epp,23,25,c_yieldk,istat)       ! ZH channel
     &                    -0.5*dsanyield_sim(mppVV,epp,23,0,pdgyield,diff,istat) ! ZZ channel
          elseif (ch_ha.eq.4) then
            mppVV = sqrt(mpp**2+mass(kz)**2-mass(kh3)**2)
            mppVS = (sqrt(mpp**2+mass(kh1)**2-mass(kh3)**2) + mpp)/2.d0 
            mppVS2 = (2.0d0*mppVS + mppVV - mpp)/2.d0        
            dsanyieldcorr=dsanyield_ch(mppVS,epp,35,36,c_yieldk,istat)        ! Ah channel
     &                    - dsanyield_ch(mppVS2,epp,23,35,c_yieldk,istat)     ! Zh channel
     &                    + 0.5*dsanyield_sim(mppVV,epp,23,0,pdgyield,diff,istat) ! ZZ channel
          elseif (ch_ha.eq.7) then
            mppVV = sqrt(mpp**2+mass(kw)**2-mass(khc)**2)
            mppVS = (mppVV + mpp)/2.d0 
            dsanyieldcorr=dsanyield_ch(mppVS,epp,24,-37,c_yieldk,istat)       ! WH channel (no factor 0.5!)
     &                    -0.5*dsanyield_sim(mppVV,epp,24,0,pdgyield,diff,istat) ! WW channel
          else
            write(*,*) 'ERROR in dsib2convint: wrong channel ch_ha = ',ch_ha
            stop
          endif
          dsIB2convint = dsIB2dsde_aux(zint)*dsanyieldcorr
        endif

      elseif (ib2dnhow.eq.'FSR') then
        return  ! FIXME not yet implemented...
c        dsIB2convint = 0.5d0*
c     &        dsib2dsde_auxfsr(zint,c_pfinal,c_Btype,c_ftype)*
c     &          dshayield(zint*mx,c_xstable*mx,ch_ha,c_yieldk,istat)

      else
        write(*,*) 'ERROR in dsIB2convint: undefined value of ib2dnhow = ',
     &      ib2dnhow
      
      endif

      return 
      end



