*******************************************************************************
*** Function dsanwx provides the  WIMP self-annihilation invariant rate.    ***
***                                                                         ***
***  type : INTERFACE                                                        ***
***                                                                         ***
***  Input:                                                                 ***
***    p - initial cm momentum (real) for DM annihilations                  ***
***  Output:                                                                ***
***  BeginTex
***    \begin{displaymath}
***    W_{\rm{eff}} = \sum_{ij}\frac{p_{ij}}{p_{11}}
***    \frac{g_ig_j}{g_1^2} W_{ij} = 
***    \sum_{ij} \sqrt{\frac{[s-(m_{i}-m_{j})^2][s-(m_{i}+m_{j})^2]}
***    {s(s-4m_1^2)}} \frac{g_ig_j}{g_1^2} W_{ij}.
***    \end{displaymath}
***    where the $p$'s are the momenta, the $g$'s are the internal
***    degrees of freedom, the $m$'s are the masses and $W_{ij}$ is
***    the invariant annihilation rate for the included subprocess.
***  EndTex
***  passed to dsrdens by dsrdomega.                                        ***
***                                                                         ***
***  author: Paolo Gondolo 2016                                             ***
*******************************************************************************
      function dsanwx(p)
      implicit none
      include 'dssilveira_zee.h'

      real*8 p,dsanwx, dsmwimp
      real*8 mx,p1p2,sv,dssigmavpartial
      integer ich

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dsanwx')

      dsanwx=0.0d0
c      write (*,*) 'PG dsanwx> entered'

c... however, if v denotes the relative velocity of one particle in the frame of the other,
c... then the invariant rate W is
c...  W = 4 p1.p2 sigma v
c...    = 4 (p^2+sqrt(p^2+m1^2)*sqrt(p^2+m2^2)) sigma v
c... and for equal mass m1=m2==mx
c...  W = 4 (mx^2+2*p^2) sigma v

      mx=dsmwimp()
      p1p2=mx**2+2.d0*p**2
      sv=0.d0
      do ich=1,18
        sv=sv+dssigmavpartial(ich,p)  ! note that v.ne.vCMS here!
      enddo
      dsanwx=4.d0*p1p2*sv/gev2cm3s
c      write (*,*) 'PG dsanwx> exited'
      return
      end

      
