!Integrand for calculating alpha^2 for power-law cosmological power spectra (Eq B9 in arXiv:1110.2484)
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2011

      double precision function alphaIntegrand(x,n,alpha,kratio)

      implicit none

      double precision, intent(IN) :: x, n, alpha, kratio
      double precision :: dsucmh_transFunc, sinc, logx    

      logx = log(x)
      alphaIntegrand = dsucmh_transFunc(x)
      alphaIntegrand = alphaIntegrand*(sinc(x)-cos(x)) 
      alphaIntegrand = x**(n-2.d0+alpha*(logx+log(kratio))) * 
     &                  kratio**(alpha*logx) * alphaIntegrand*alphaIntegrand
      
      end function alphaIntegrand

