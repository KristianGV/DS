c_______________________________________________________________________
c  Function dsfi2to2rhs gives the right hand side of the Boltzmann equation
c   s*dY/dt=RHS. Done by integrating dsfi2to2int over s, and multiplying by
c   temperature dependent factors in the end.
c   
c   input:
c     tmp - temperature of heath bath
c
c
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-09-07
c=======================================================================
      real*8 function dsfi2to2rhs(lnT)
      implicit none
      include 'dsmpconst.h'
      include 'dsficom.h'
      include 'dsidtag.h'



      real*8 lnT,tmp,a,b,eps,sum,geff,sqrtgstar,heff,s_ent,H,HPrime,
     &dsf_int2,abserr,resabs,resasc
      real*8, external :: dsfi2to2int_simp

      tmp=exp(lnT)
      call dsrdset('dof','default')
      call dskdgeff(tmp,geff)
      call dsrddof(tmp,sqrtgstar,heff)

      s_ent=heff*2*pi*pi/45*tmp*tmp*tmp
      H=sqrt(4*pi**3*geff/45)*tmp*tmp/mpl
      HPrime=H*heff/sqrt(geff)/sqrtgstar 

      T_22=tmp;x1_22=m1_22/tmp;x2_22=m1_22/tmp
c     Integration limits
      a=1.d-10;eps=1.d-5
      b=1/((2*mdm)**2)
      ! b=1

      ! a=4*mdm**2; b=1.d10

c     Integrating over s using adaptive Gaussian integration routine    
      call dgadap(a,b,dsfi2to2int_simp,eps,sum)
      ! call dqk21(dsfi2to2int_simp,a,b,sum,eps,resabs,resasc)
      dsfi2to2rhs=sum/HPrime/s_ent*tmp
      return
      end
      
