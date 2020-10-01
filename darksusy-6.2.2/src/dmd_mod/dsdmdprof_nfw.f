************************************************************************
***  dark matter density profiles written in the form:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***
***  namely dsdmdprof is the function g(x), e.g., in Eq. (10) of 
***  astro-ph/0207125.
***      
***  NOTE: version valid for a NFW profile only      
***      
***  input:
***    x [1]
************************************************************************
      real*8 function dsdmdprof_nfw(x)
      implicit none
      real*8 x
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprof_nfw=1.d0/x/(1.d0+x)**2
      end


************************************************************************
***  first derivative with respect to x of the function dsdmdprof(x) as 
***  implemented in the definition of the density profile:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***      
***  NOTE: version valid for a NFW profile only      
***      
***  input:
***    x = r/rs [1]
************************************************************************
      real*8 function dsdmdprofdx_nfw(x)
      implicit none
      real*8 x
      if(x.lt.1.d-16) x=1.d-16
      dsdmdprofdx_nfw=-1.d0/x**2/(1.d0+x)**2-2.d0/x/(1.d0+x)**3
      end


************************************************************************
***  second derivative with respect to x of the function dsdmdprof(x) as 
***  implemented in the definition of the density profile:
***
***    rho(r) = rhos * dsdmdprof(r/rs)
***      
***  NOTE: version valid for a NFW profile only      
***      
***  input:
***    x = r/rs [1]
************************************************************************
      real*8 function dsdmd2profdx2_nfw(x)
      implicit none
      real*8 x
      if(x.lt.1.d-16) x=1.d-16
      dsdmd2profdx2_nfw=6.d0/(x*(1.d0+x)**4) + 4.d0/(x**2*(1.d0+x)**3)+ 
     &        2.d0/(x**3*(1.d0+x)**2)
      end
      


************************************************************************
***  volume integral I1 of the dark matter density profile:
***
***    dsI1(x) = \int_0^x dxp xp**2 dsdmprof(xp)
***
***  namely, e.g., Eq. (22) of astro-ph/0207125 with n=1.
***      
***  NOTE: version valid for a NFW profile only      
***      
***  input:
***    x [1]
************************************************************************
      real*8 function dsI1_nfw(x)
      implicit none
      real*8 x,I1
      if(x.lt.1.d-5) then
        I1=x**2/2.d0-x**3*2.d0/3.d0+x**4*3.d0/4.d0
      else  
        I1=dlog(1.d0+x)+1.d0/(1.d0+x)-1.d0
      endif
      dsI1_nfw=I1
      end
      
