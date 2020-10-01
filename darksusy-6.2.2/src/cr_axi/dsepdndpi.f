************************************************************************
*** initial (i.e. prior propagation) electron/positron spectrum 
*** according to the the default in ds, i.e. dsanyield. 
***
***  type : REPLACEABLE                                                   
***                                                                       
***  input: pp = momentum (GeV)
***         power = integer selecting annihilations (=2) or decays (=1)
***  output: GeV^-1  / [(cm^3 s) / (GeV / cm^3)^crdmpow]
***
*** author: Piero Ullio
*** modified: Torsten Bringmann, 12/06/2015 
***           (made replaceable & linked to dscrsource)
************************************************************************
      real*8 function dsepdndpi(pp,power)
      implicit none
      include 'dsanyieldcom.h'  ! for ansmooth
      include 'dsmpconst.h'  
      real*8 pp,ee,ratio,dedp,eemin
      integer power
      real*8 dscrsource
      integer istat,ansmoothbkp
ccc  
      ratio=m_e/pp
      if(ratio.gt.1.d-8) then
        ee=dsqrt(pp**2+m_e**2) ! energy in GeV
        dedp=pp/ee
      else
        ee=pp*(1.d0+ratio**2/2.d0)
        dedp=1.d0-ratio**2/2.d0
      endif
      eemin=0.05d0 ! minimum energy in GeV at which the function
                   ! dsanyield is trusted, below dn/dE is assumed 
                   ! to be constant, is the chosen value ok?
      ee=max(ee,eemin)
      ansmoothbkp=ansmooth
      ansmooth=2  ! 2 to smooth well
c... Note extra factor 1d30 required by integration routines! taken out
c and reintroduced in dsepdndpaxi_int      
c... FIXME: consider taking out the dependence on pspower in dscrsource
      dsepdndpi=dedp*dscrsource(ee,1,-11,power,0.0d0,istat)
      ansmooth=ansmoothbkp
      return
      end

      
