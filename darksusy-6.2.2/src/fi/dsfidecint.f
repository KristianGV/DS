c_______________________________________________________________________
c  Function dsfidecint gives the integrand use to calculate freeze-in
c     decay abundance.
c
c  
c    
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-17
c=======================================================================

      real*8 function dsfidecint(lnT,M,eta)
      implicit none

      include 'dsmpconst.h'
      include 'dsficom.h'

      real*8 tmp,lnT,M,eta,H,s,sqrtgstar,heff, geff,HPrime,zero,
     &k1
      tmp=exp(lnT); 
      zero=0


      call dsrdset('dof','default')
      call dskdgeff(tmp,geff)
      call dsrddof(tmp,sqrtgstar,heff)
      call dsfik1(M/tmp,zero,zero,eta,zero,zero,k1)

      s=heff*2*pi*pi/45*tmp*tmp*tmp
      H=sqrt(4*pi**3*geff/45)*tmp*tmp/mpl
      HPrime=H*heff/sqrt(geff)/sqrtgstar 

      dsfidecint= k1/HPrime/s*tmp
      return
      end

c     Made integrant into function of one variable to use in inegration routine
      real*8 function dsfidecint_simp(tmp)
      implicit none

      include 'dsficom.h'

      real*8 tmp,M,eta,dsfidecint,x1,x2,x3,eta1,eta2,eta3
      if(stat.eq.0) then
            eta_dec=0
      end if
      dsfidecint_simp=dsfidecint(tmp,M_dec,eta_dec)
      return
      end