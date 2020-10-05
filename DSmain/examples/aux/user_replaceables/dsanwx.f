*******************************************************************************
*** Function dsanwx provides the  WIMP self-annihilation invariant rate.    ***
***                                                                         ***
***  type : interface                                                       ***
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
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2015-06-11                                                         ***
*** mod  2020-04-12 TB changed to Møller velocity as reference              ***
*******************************************************************************
      real*8 function dsanwx(p)
      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsmpconst.h'

      real*8 p
      real*8 s, vmoeller, dsmwimp, vcm2, vratio, dsmass

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('generic_wimp','dsanwx')


      dsanwx=0.0d0
      
      s = 4.0d0|*(dsmwimp()**2+p**2)

      vmoeller = 2.0d0*p*sqrt(s)/(s-2.0d0*dsmwimp()**2)
      
c      if (4*dsmass(svch)**2.gt.s) return ! CMS energy must be large enough to
                                         ! produce final states!
      
c... the expression below assumes that sva gives the cross section
c... times the *Møller* velocity
       dsanwx = 2.0d0*(s-2.0d0*dsmwimp()**2)*(sva + svb*vmoeller**2)/gev2cm3s

c... useful relations to express this in a different way:
c       vmoeller = 2.0d0*p*sqrt(s)/(s-2.0d0*dsmwimp()**2)
c       vcm2     = 1.0d0-4*dsmwimp()**2/s                  ! vcm    = p/(2 sqrt(s))
c       vratio   = 1.0d0/(1.0d0+vcm2)                      ! vratio = vmoeller / (2vcms)

      return
      end
