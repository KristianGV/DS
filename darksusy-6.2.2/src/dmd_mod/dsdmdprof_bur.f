************************************************************************
***  dark matter density profiles written in the form:
***
***    rho(r) = rhos * dsdmprof(r/rs)
***
***  namely dsdmdprof is the function g(x), e.g., in Eq. (10) of 
***  astro-ph/0207125.
***      
***  NOTE: version valid for a Burkert profile only      
***      
***  input:
***    x [1]
************************************************************************
      real*8 function dsdmdprof_bur(x)
      implicit none
      real*8 x
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprof_bur=1.d0/(1.d0+x)/(1.d0+x**2)
      end



************************************************************************
***  first derivative with respect to x of the function dsdmdprof(x) as 
***  implemented in the definition of the density profile:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***      
***  NOTE: version valid for a Burkert profile only     
***      
***  input:
***    x = r/rs [1]
************************************************************************
      real*8 function dsdmdprofdx_bur(x)
      implicit none
      real*8 x
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprofdx_bur=-((2.d0*x)/((1.d0+x)*(1.d0+x**2)**2))
     &  -1.d0/((1.d0+x)**2*(1.d0+x**2))
      end


************************************************************************
***  second derivative with respect to x of the function dsdmdprof(x) as 
***  implemented in the definition of the density profile:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***      
***  NOTE: version valid for a Burkert profile only      
***      
***  input:
***    x = r/rs [1]
************************************************************************
      real*8 function dsdmd2profdx2_bur(x)
      implicit none
      real*8 x
      if(x.lt.1.d-16) x=1.d-16
      dsdmd2profdx2_bur=(4.d0*x**2*(3.d0+4.d0*x+3.d0*x**2))
     &  /((1.d0+x)**3*(1.d0+x**2)**3)
      end
      
      
************************************************************************
***  volume integral I1 of the dark matter density profile:
***
***    dsI1(x) = \int_0^x dxp xp**2 dsdmprof(xp)
***
***  namely, e.g., Eq. (22) of astro-ph/0207125 with n=1.
***      
***  NOTE: version valid for a Burkert profile only     
***      
***  input:
***    x [1]
************************************************************************
      real*8 function dsI1_bur(x)
      implicit none
      real*8 x,I1
      if(x.lt.1.d-5) then
        I1=x**3/3.d0-x**4/4.d0
      else  
        I1=0.5d0*(dlog(1.d0+x)+0.5d0*dlog(1.d0+x**2)-datan(x))
      endif
      dsI1_bur=I1
      end
      
