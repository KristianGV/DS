      real*8 function dsrdbreit_wigner(mres,gres,alpha,mwimp,pwimp)
c_______________________________________________________________________
c  
c  input:
c    - mres - mass of resonance (GeV)
c    - gres - width of resonance (GeV)
c    - alpha, exponent of correction factor, pwimp**alpha      
c    - mwimp - mass of WIMP (GeV)
c    - pwimp - momentum of WIMP (GeV)
c  output:
c    the Breit-Wigner distribution, unnormalized for a resonance of mass
c    mres and width gres, calculated at centre of mass
c    energy 2*sqrt(mwimp**2+pwimp**2)      
c  author: joakim edsjo (edsjo@fysik.su.se) 2018-10-25
c=======================================================================
      implicit none
      real*8 mres,gres,alpha,mwimp,pwimp
      real*8 E
      real*8 k
      real*8 corr
      parameter(k=1.d0) ! normalization factor

      E=sqrt(mwimp**2+pwimp**2)*2.d0
      corr=pwimp**alpha  ! correction 
      dsrdbreit_wigner=k/((E**2-mres**2)**2+mres**2*gres**2)*corr
c      dsrdbreit_wigner=k/((E-mres)**2+gres**2/4.d0)*corr      
            
      return
      end








