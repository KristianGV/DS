c_______________________________________________________________________
c  Function dsfidecint.
c
c  
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfidecint(lnT,M,eta,stat)
      implicit none

      include 'dsmpconst.h'

      real*8 tmp,lnT,M,eta,H,s,sqrtgstar,heff, geff,HPrime,zero,
     &k1
      integer stat
      tmp=exp(lnT)
      zero=0
      call dsrdset('dof','default')
      call dskdgeff(tmp,geff)
      call dsrddof(tmp,sqrtgstar,heff)
      call dsfik1(M/tmp,zero,zero,eta,zero,zero,stat,k1)
      s=heff*2*pi*pi/45*tmp*tmp*tmp
      H=sqrt(4*pi**3*geff/45)*tmp*tmp/mpl
      HPrime=H*heff/sqrt(geff)/sqrtgstar 
      dsfidecint= tmp*k1/HPrime/s
      return
      end

      real*8 function dsfidecint_simp(lnT)
      implicit none

      real*8 lnT,M,eta,dsfidecint,x1,x2,x3,eta1,eta2,eta3
      integer stat
      common /dsfidecvar/ M,eta
      common /dsfik1var/ x1,x2,x3,eta1,eta2,eta3,stat

      dsfidecint_simp=dsfidecint(lnT,M,eta,stat)
      return
      end