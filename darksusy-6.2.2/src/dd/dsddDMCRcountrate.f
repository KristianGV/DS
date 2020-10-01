*******************************************************************************
*** Function dsddDMCRcountrate compares the total event rate from nuclear   ***
*** recoils caused by high-energy DM particles (from their interaction with ***
*** cosmic rays, see Bringmann & Pospelov (2018)) to limits reported by     ***
*** various experiments.                                                    ***
***                                                                         ***
***  Output: *ratio* of expected count rates to reported sensitivity        ***
***                                                                         ***
*** In a sense, this is a collection of examples for how to call            ***
*** dsddDMCRdgammadt in order to get the differential countrate for a       ***
*** number of possible configurations. See each of the implementd           ***
*** options further down for further details.                               ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-02-02                                                         ***
*** mod  2019-03-15 merged with previous 'dsddDMCRcountrate_full', added    ***
***                 dsddDMCRsigtarget with options set here                 ***
*******************************************************************************
      real*8 function dsddDMCRcountrate(option)
      implicit none
      include 'dsio.h'
      include 'dsmpconst.h'
      include 'dsddcom.h'
      character*(*) option

      real*8 Tmin,Tmax, ntargets, mtarget,depth 
      real*8 rhoDeff, ratelimit, sigtlimit, res, kappaSHM, corr, mdm
      real*8 sip, sdp, sin, sdn, sigsi
      integer ierr, targetoption
      common /sigtargetcom/ targetoption
      save /sigtargetcom/

      real*8 mNcom,depthcom,rhoDeffcom
      common /DMDRcountaux/ mNcom,depthcom,rhoDeffcom 
      save /DMDRcountaux/ 

c... functions
      real*8 dsddDMCRcountrate_aux, dsmwimp
      external dsddDMCRcountrate_aux
      
*******************************************************************************
c... default options for the prediction of the CR-induced DM flux
c... that are independent pf the experimental settings (see further down for those)
*******************************************************************************

c... the nucleon cross sections here only enter for the calculation of the mean 
c... free path. These settings do NOT affect the DMCR flux or the scattering rate
c... in the experiment. For those, see dsddDMCRsigCR and dsddDMCRsigtarget!
      call dsddsigmanucleon(0.0d0,0.0d0,sip,sin,sdp,sdn,ierr)
      sigsi = (sip+sin)/2.  ! take *average* to compute simplified mean free path      
                                           
c... astro parameters -> now set in dsddinit instead!
c      rholocal = 0.3d0 ! GeV/cm^3  
c      vlocal = 2./sqrt(pi)*220*1.0d5 ! mean local DM velocity in cm/s
c      Deff = 1.0d0 ! kpc; distance out to which source term is integrated
      rhoDeff = rholocal*Deff

*******************************************************************************
c... settings specific to the respective experimental configurations    
*******************************************************************************
      if (option.eq.'Xenon1t') then
         targetoption=10 ! for dsddDMCRsigtarget
         call dsddDMCRquenching_set('no_quenching') 
         Tmin = 4.9d-6        ! this is the Xenon1t analysis window for signals
         Tmax = 40.8d-6       ! (in GeV)
         mtarget = 122.05     ! mXe [GeV]
         mdm = dsmwimp()      ! DM mass [GeV]
         kappaSHM = 0.23      ! fraction of SHM that is sampled, see (16) in 1810.10543
         depth =  1400.0d2    ! detector located 1400m underground
         sigtlimit = 8.3d-46  ! Xenon 1t limit *per nucleon* for 1 TeV DM particles 
                              ! [1805.12562]
         corr = (mdm+mtarget)**2/(mdm+m_p)**2 ! correct for difference in coherence  
                                              ! factor at 1 TeV and mdm.  
                                              ! For a more accurate estimate, initialize
                                              ! model with 1 Tev, and divide result
                                              ! from dsddg2sigma with that for the
                                              ! actual model (in each case averaging
                                              ! over all isotopes)

         ratelimit = kappaSHM*rholocal/1.0d3 ! because we chose to take sigtlimit at 1 TeV
     &               *corr*sigtlimit*vlocal  ! this should now (roughly) be the  
                                             ! maximal rate [in 1/s] per Xe nucleus 
                                             ! that is compatible with the published limit
                                              
      elseif (option.eq.'Borexino'.or.option.eq.'Borexino_SD') then
         targetoption=20 ! for dsddDMCRsigtarget, dsddTDMattenuation
         if (option.eq.'Borexino_SD') targetoption=21
         call dsddDMCRquenching_set('borexino') 
         Tmin = 12.5d-3       ! lower energy for seeing no signal in Borexino
         Tmax = 50*Tmin       ! (in GeV)
         mtarget = m_p        ! at these energies, targets are free protons [GeV]
         depth =  1400.0d2    ! detector located 1400m underground
         ratelimit = 2.44             ! 90% C.L. limit for
     &               / (1.282*year)   ! data taking period 
     &               / 3.2d31         ! per proton in detector volume          

                                          
      elseif (option.eq.'MiniBoone') then
         targetoption=30 ! for dsddDMCRsigtarget
         call dsddDMCRquenching_set('no_quenching') 
         Tmin = 35.0d-3       ! this is the lower energy for the counting rate in MiniBoone
         Tmax = 40*Tmin       ! (in GeV)
         mtarget = m_p        ! at these energies, targets are free protons [GeV]
         depth =  2*3.0d2     ! detector located 3m underground
                              ! times 2 for effective thickness of atmosphere
         ratelimit = 800.0          ! counting limit for 
     &               / 893.0        ! data taking period (in seconds)
     &               / 5.83d31      ! per proton in detector volume          

      else
         write(*,*) 'FATAL: dsddDMCRcountrate has been called'
         write(*,*) 'with an unsupported option: ', option
         write(*,*) 'I do not know what to do... '
         stop
      endif


c... tramsfer input values to common blocks for integration over recoil energy
      mNcom = mtarget
      depthcom = depth
      rhoDeffcom = rhoDeff
      Tmin=log(Tmin)  ! to better sample fast falling integrand
      Tmax=log(Tmax)      
      call dgadap(Tmin,Tmax,dsddDMCRcountrate_aux,1.0d-3,res)

c... Finally, we can compare the rate to the respective limit
c... (undoing the rescaling in the integrand)
      dsddDMCRcountrate = res*1.0d-30 / ratelimit 

c      write(*,*) 'In dsddDMCRcountrate: sigsi, zlfree, res, ratelimit = ',sigsi,zlfree, res, ratelimit

      return
      end



*******************************************************************************
*** auxiliary routine just for integration
*******************************************************************************
      real*8 function dsddDMCRcountrate_aux(Tquench)
      implicit none
      real*8 Tquench, Tr, dTrdTq, dsddDMCRdgammadt
      real*8 mNcom,depthcom,rhoDeffcom
      common /DMDRcountaux/ mNcom,depthcom,rhoDeffcom 
      real*8 dsddDMCRsigtarget
      external dsddDMCRsigtarget

      call dsddDMCRquenching(exp(Tquench), Tr, dTrdTq)

c... rescale with a large number to improve convergence
      dsddDMCRcountrate_aux = 1.0d30*dTrdTq
     &     *exp(Tquench) !from transformation to log integrand
     &     *dsddDMCRdgammadt(Tr,mNcom,dsddDMCRsigtarget,depthcom,rhoDeffcom)


      return
      end
