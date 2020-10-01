*******************************************************************************
*** Function dssisigtmav provides the velocity-averaged momentum-transfer   ***
*** cross section per DM mass, assuming a Maxwellian velocity distribution  ***
*** of the DM particles.                                                    ***
***                                                                         ***
***  type : commonly used                                                   ***
***  desc : velocity-averaged momentum-transfer cross section               ***
***                                                                         ***
***  Input:                                                                 ***
***    v0 - most probable velocity of individual DM particles               ***
***         (i.e. the mean speed is given by 2/sqrt(pi)*v0)                 ***
***         units: km/s                                                     ***
***                                                                         ***
***  Output:                                                                ***
***  BeginTex
***    \begin{displaymath} 
***    \langle \sigma_T \rangle = 
***       \frac{4 \pi}{(\sqrt{2 \pi} v_0)^3} 
***        \int_0^\infty dv\,v^2 \exp{-v^2/(2 v_0^2)} \sigma_T (v) \,,
***    \end{displaymath}
***    where $v$ is the {\it relative} (i.e. not individual) DM velocity.
***  EndTex
***                                                                         ***
***  This function requires the particle module to provide the interface    ***
***  function dssisigtm.                                                    ***
***                                                                         ***
***  Units of output: cm**2 / g.                                            ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-05-17                                                         ***
*******************************************************************************
      real*8 function dssisigtmav(v0)
      implicit none
      include 'dsmpconst.h'

      real*8 v0
      real*8 intres,a,b, dsauxxtMBint
      external dsauxxtMBint
      real*8 v0MB
      common /v0MB/ v0MB

      dssisigtmav = 0.0d0
      
      a = 0.0d0
      b = 7.0d0 ! "infinity" for the MB distribution 
                ! (warning: too large values result in wrong results from dgadap!!)
      v0MB = v0
      call dgadap(a,b,dsauxxtMBint,1.0d-4,intres)
      
      dssisigtmav = 4./sqrt(pi) * intres
            
      return
      end


*******************************************************************************
*******************************************************************************
c... auxiliary function to integrate \sigma_T. xv = v/(sqrt(2)*v0)
      real*8 function dsauxxtMBint(xv)
      implicit none
      real*8 xv,vrel,dssisigtm    
      external dssisigtm
      real*8 v0MB
      common /v0MB/ v0MB

      vrel = 1.41421356*v0MB*xv
      dsauxxtMBint = xv**2*exp(-xv**2)*dssisigtm(vrel)

      return 
      end
