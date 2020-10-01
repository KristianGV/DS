*******************************************************************************
*** Function dssisigmatclassical provides the momentum-transfer cross       ***
*** section for a Yukawa potential in the classical regime.                 ***
***                                                                         ***
***  Input:                                                                 ***
***    mmed - mediator mass (in GeV)                                        ***
***    beta - ratio of potential to kinetic energy,                         ***
***           beta = 2*alpha*mmed/mdm/vrel**2                               ***
***    rep  - specifies whether the potential is attractive (rep=0)         ***
***           or repulsive (rep=1)                                          ***
***                                                                         ***
***  Output: momentum-transfer cross section in units of GeV**-2.           ***
***                                                                         ***
***  This implements Eqs (45,46) in 1512.05344, based on numerical results  ***   
***  by Khrapak et al. (2003,2004), and provides the correct description of ***
***  \sigma_T for vrel mdm/mmed >> 1.                                       ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-05-17                                                         ***
*******************************************************************************
      real*8 function dssisigmatclassical(mmed,beta,rep)
      implicit none
      include 'dsmpconst.h'

      real*8 mmed,beta, tmp
      integer rep

      dssisigmatclassical = 0.0d0
      tmp=0.0d0

      if (mmed.le.0.d0.or.beta.le.0d0) then
        write(*,*) 'FATAL error: dssisigmatclassical called with' 
        write(*,*) 'negative mass or coupling strength!'
        write(*,*) 'mmed, beta = ',mmed,beta 
        stop        
      endif

      if (rep.eq.0) then ! attractive potential

        if (beta.lt.0.01) then
          tmp = 2*beta**2*log(1.+1./beta**2)
        elseif (beta.lt.100.0) then
          tmp = 7*(beta**1.8 + 280*(beta/1.0d1)**10.3)
     &          /(1. + 1.4*beta + 0.006*beta**4 + 160.*(beta/1.0d1)**10)
        else ! beta > 100
          tmp = 0.81*(1. + log(beta) - 1./(2.*log(beta)))**2
        endif

      elseif (rep.eq.1) then ! repulsive potential

        if (beta.lt.0.01) then
          tmp = 2*beta**2*log(1.+1./beta**2)
        elseif (beta.lt.1.0d4) then
          tmp = 8*beta**1.8 / (1. + 5*beta**0.9 + 0.85*beta**1.6)
        else ! beta > 10000
          tmp = (log(2.*beta) - log(log(2.*beta)))**2
        endif

      else 
        write(*,*) 'FATAL error: dssisigmatclassical called with' 
        write(*,*) 'unrecognized parameter rep = ',rep 
        stop
      endif
            
      dssisigmatclassical = tmp*pi/mmed**2      
            
      return
      end

