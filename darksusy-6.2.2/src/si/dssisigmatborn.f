*******************************************************************************
*** Function dssisigmatborn provides the momentum-transfer cross section    ***
*** for a Yukawa potential in the Born regime.                              ***
***                                                                         ***
***  Input:                                                                 ***
***    mmed - mediator mass (in GeV)                                        ***
***    beta - ratio of potential to kinetic energy,                         ***
***           beta = 2*alpha*mmed/mdm/vrel**2                               ***
***    R    - ratio of interaction range to de Broglie wavelength           ***
***           R = vrel*mDM/mmed                                             ***
***                                                                         ***
***  Output: momentum-transfer cross section in units of GeV**-2.           ***
***                                                                         ***
***  This implements Eq. (5) in 0911.0422, and provides the correct         ***
***  description of \sigma_T for \alpha mdm/mmed << 1.                      ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-05-17                                                         ***
*******************************************************************************
      real*8 function dssisigmatborn(mmed,beta,R)
      implicit none
      include 'dsmpconst.h'

      real*8 mmed,beta,R, tmp

      dssisigmatborn = 0.0d0

      if (mmed.le.0.d0.or.beta.le.0d0.or.R.le.0.0d0) then
        write(*,*) 'FATAL error: dssisigmatborn called with' 
        write(*,*) 'negative mass or coupling strength!'
        write(*,*) 'mmed, beta, R = ', mmed,beta,R
        stop        
      endif
      
      tmp = 2*pi*beta**2/mmed**2
      
      if (R.gt.1.0d-3) then
        dssisigmatborn =  tmp*(dlog(1.+R**2) - R**2/(1.+R**2))           
      else ! evaluation of log is not accurate enough -> take limit of small R    
        dssisigmatborn = tmp*(R**4/2.-R**6*2./3.)
      endif
            
      return
      end

