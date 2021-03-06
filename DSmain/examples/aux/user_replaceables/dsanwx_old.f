**********************************************************************
*** This file is automatically generated from the file 
*** src_models/generic_wimp/an/dsanwx.f
*** with the script scr/make_replaceable.pl on Aug 27, 2020.
*** The file is copied as is, but to be of any use you should of
*** course modify this file to your liking. The way the default linking
*** is set up, this file will be linked to before the corresponding
*** DarkSUSY file, meaning that this is the file that will be used,
*** i.e. it will replace the default DarkSUSY one.
*** A few lines of code are added to the executable section below
*** to remind you that you need to change this file. Delete those lines
*** (and these comments) and modify the file to your liking and you
*** are ready to go!
**********************************************************************

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
*******************************************************************************
      real*8 function dsanwx(p)
      implicit none
      include 'dsgeneric_wimp.h'
      include 'dsmpconst.h'

      real*8 p
      real*8 s, vmoeller, dsmwimp, pdv, dsmass,vcm

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions


      call dscheckmodule('generic_wimp','dsanwx')


      dsanwx=0.0d0
      
      s = 4.0d0*(dsmwimp()**2+p**2)
      
      if (4*dsmass(svch)**2.gt.s) return ! CMS energy must be large enough to
                                        ! produce final states!
      
      vmoeller = 4.*p/sqrt(s) ! Moeller velocity in the CMS


      pdv=sqrt(p**2+dsmwimp()**2)/2.0d0 ! p/vmoeller
      
c... the expression below assumes that sva gives the cross section
c... times twice the CMS velocity of one of the initial DM particles      
      dsanwx= 4.0d0*pdv*sqrt(s)*(sva + svb*vmoeller**2)/gev2cm3s
      
      return
      end
