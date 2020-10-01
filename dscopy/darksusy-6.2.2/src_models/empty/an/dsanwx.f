*******************************************************************************
*** Function dsanwx provides the  WIMP self-annihilation invariant rate.    ***
***                                                                         ***
***  type : interface                                                       ***
***                                                                         ***
***  desc : Self-annihilation invariant rate                                ***
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
*** author: Torsten.Bringmann@fys.uio.no                                    ***
*** date 2015-06-11                                                         ***
*******************************************************************************
      real*8 function dsanwx(p)
      implicit none
      include 'dsempty.h'

      real*8 p

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('empty','dsanwx')

      dsanwx=0.0d0 ! This produces NAN for the relic density: OK

      return

      end
