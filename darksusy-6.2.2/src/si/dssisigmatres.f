*******************************************************************************
*** Function dssisigmatres provides the momentum-transfer cross section     ***
*** for a Hulthén potential in the resonant regime.                         ***
***                                                                         ***
***  Input:                                                                 ***
***    mdm   - dark matter mass (in GeV)                                    ***
***    mmed  - mediator mass (in GeV)                                       ***
***    vrel  - relative velocity of DM particles (in units of c)            ***
***    alpha - DM-mediator coupling (alpha=g**2/4/pi)                       ***
***            beta = 2*alpha*mmed/mdm/vrel**2                              ***
***    rep   - specifies whether the potential is attractive (rep=0)        ***
***            or repulsive (rep=1)                                         ***
***                                                                         ***
***  Output: momentum-transfer cross section in units of GeV**-2.           ***
***                                                                         ***
***  This implements Eq. (A4,A5) in 1302.3898, and provides the correct     ***
***  description of \sigma_T for vrel mdm/mmed < 1, i.e. outside the.       ***
***  classical regime.                                                      ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-05-17                                                         ***
*******************************************************************************
      real*8 function dssisigmatres(mdm,mmed,vrel,alpha,rep)
      implicit none
      include 'dsmpconst.h'

      real*8 mdm,mmed,vrel,alpha
      integer rep

      complex*16 wgamma, wlgama ! complex gamma function, from CERN library
    
      real*8 kappares, a, b, delta
      complex*16 tmp, lambdaplus, lambdaminus


      dssisigmatres = 0.0d0
c      kappares = 1.6  ! see 1302.3898 for a discussion
      kappares = 1.5505205 ! =sqrt(2*zeta(3))

      if (mmed.le.0.d0.or.mdm.le.0.d0.or.vrel.le.0.d0.or.alpha.le.0.d0) then
        write(*,*) 'FATAL error: dssisigmatres called with' 
        write(*,*) 'negative mass, relative velocity or coupling strength!'
        write(*,*) 'mdm,mmed,vrel,alpha = ', mdm,mmed,vrel,alpha
        stop        
      endif

      a = mdm*vrel/2.d0/kappares/mmed
      b = mdm*alpha/kappares/mmed
            
      if (rep.eq.0) then ! attractive potential

        if (b.gt.a**2) then
          lambdaplus  = dcmplx(1.0d0 + sqrt(b-a**2), a)
          lambdaminus = dcmplx(1.0d0 - sqrt(b-a**2), a)
        else 
          lambdaplus  = dcmplx(1.0d0, a + sqrt(a**2-b))
          lambdaminus = dcmplx(1.0d0, a - sqrt(a**2-b))
        endif

      elseif (rep.eq.1) then ! repulsive potential

          lambdaplus  = dcmplx(1.0d0, a + sqrt(b+a**2))
          lambdaminus = dcmplx(1.0d0, a - sqrt(b+a**2))

      else 
        write(*,*) 'FATAL error: dssisigmatres called with' 
        write(*,*) 'unrecognized parameter rep = ',rep 
        stop
      endif

      if (a.gt.10.0.or.sqrt(b).gt.10.0) then
        tmp=dcmplx(0.0d0,1.0d0)*exp(wlgama(dcmplx(0.0d0,2.0d0*a))
     &       - wlgama(lambdaplus) - wlgama(lambdaminus))
      else ! use gamma itself
        tmp = dcmplx(0.0d0,1.0d0)*wgamma(dcmplx(0.0d0,2*a))
     &        /wgamma(lambdaplus)/wgamma(lambdaminus)
      endif

      
      delta= datan(dimag(tmp)/dble(tmp))
                  
      dssisigmatres = 16.*pi/mdm**2/vrel**2 * sin(delta)**2      
            
c      write(*,*) delta,dssisigmatres
            
      return
      end

