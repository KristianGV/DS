*****************************************************************************
*** function dsIB3yieldone calculates the IB3 yield from one given 
*** annihilation channel, i.e. the *difference* to the corresponding qq
*** yield when taking into account qqg final states. 
***
*** Currently included are:
***   yieldk =  54 - antiproton yield above threshold emuthr
***            154 - differential antiproton yield at emuthr
***             52 - photon yield above threshold emuthr
***            152 - differential photon yield at emuthr
***
*** The annihilation channels are:
***  qch  =  7 - u u-bar
***          8 - d d-bar
***          9 - c c-bar
***         10 - s s-bar
***         11 - t t-bar
***         12 - b b-bar
***
*** See arXiv: 1510.02473 for more details on the scheme implemented here.
***
*** the units are (annihilation into IBch)**-1
*** for the differential yields, the units are the same times gev**-1.
***
*** istat will set upon return in case of errors
***   bit  decimal  reason
***     0        1  dsIBf_intdxdy failed
***     1        2  dsIBf_intdy failed
*** Author: torsten.bringmann@fys.uio.no
*** Date: 2015-06-01
*****************************************************************************

      real*8 function dsIB3yieldone(emuthr,qch,yieldk,istat)
      implicit none
      include 'dsmssm.h'


c------------------------ variables ------------------------------------

      real*8 emuthr,mDM,tmpresult, y
      integer qch,istat,yieldk, mix, pdg, yieldpdg, diff

c------------------------ functions ------------------------------------

      real*8 dsanyield_sim, dsIB3yieldtab, dsib3svqqgratio, dsib3svqqratio

c-----------------------------------------------------------------------

      istat=-10 ! channel not implemented
      dsIB3yieldone=0d0
      
      if (yieldk.ne.54.and.yieldk.ne.154.and.yieldk.ne.52.and.yieldk.ne.152)
     &  return      

      if (qch.lt.7.or.qch.gt.12) then
        write(*,*) 'ERROR in dsIB3yieldone: unknown channel IBch = ', qch
        istat=-20
        return
      endif

      istat=0
      tmpresult=0d0
      mDM=mass(kn(1))
  
c...only if kinematically allowed go on to compute yields
      if (mdm.gt.mass(qch).and.emuthr.le.(0.9999*mdm*(1.-mass(qch)**2/mDM**2))) then
        
        call dsIB3yieldfit(qch,yieldk,mix,y)
        if (y.gt.1.5.or.y.lt.-0.5) then ! something went wrong, return zero
          return
        endif

c... interpolate 3-body spectra between the two extreme cases        
        tmpresult = y*dsIB3yieldtab(emuthr,mDM,qch,2-mix,yieldk) ! VIB-type spectrum
        
        if (y.lt.0.999) ! add heavy-squark-type spectrum
     &     tmpresult= tmpresult+(1-y)*dsIB3yieldtab(emuthr,mDM,qch,4-mix,yieldk)
     
        tmpresult = tmpresult*dsib3svqqgratio(kn(1),qch)

        if (NLOoption.eq.'default') ! correct for the fact that the qq cross section
                                    ! already contains QCD corrections
     &    tmpresult = tmpresult/dsib3svqqratio(2*mDM,kn(1),qch)
        
c... add change in normalization of 2-body spectrum (if not already done in dssigmav0)
        if (NLOoption.ne.'default') then
          pdg  = qch-7+2*mod(qch,2) 
          diff = yieldk/100
          if (mod(yieldk,100).eq.54) yieldpdg = -2212 ! pbar
          if (mod(yieldk,100).eq.52) yieldpdg = 22    ! gamma        
          tmpresult = tmpresult +  (dsib3svqqratio(2*mDM,kn(1),qch)-1.)*
     &                      dsanyield_sim(mDM,emuthr,pdg,0,yieldpdg,diff,istat)
        endif
 
      endif

      dsIB3yieldone=tmpresult

      return

      end


