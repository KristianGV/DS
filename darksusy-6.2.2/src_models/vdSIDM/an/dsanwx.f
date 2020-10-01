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
      include 'dsvdSIDM.h'
      include 'dsmpconst.h'

      real*8 p
      real*8 s, sigmavmoeller, sv0, vrelcms, SE, alpha, dsansommerfeld
      integer l ! s-wave or p-wave
      

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('vdSIDM','dsanwx')


      dsanwx=0.0d0
      sv0 = 0.0d0 ! this is the contribution in the v->0 limits (s- or p-wave)
      l=-1
      
      s = 4.0d0*(mass(kdm)**2+p**2)
      
      if (4.*mass(kmed)**2.gt.s) return ! CMS energy must be large enough to
                                        ! produce final states!

      vrelcms = 4.0d0*p/sqrt(s) ! 2*v_cms 

      sigmavmoeller=0.0d0

c... DM fermion, vector mediator:
      if (DMtype.eq.1.and.Mediatortype.eq.2) then
        
         l = 0 ! s-wave
      
c... t- and u-channel from mathematica
         sv0 = gDM**4*(mass(kdm)*(mass(kdm)**2 - mass(kmed)**2)**1.5)/
     -         (4.*Pi*(-2*mass(kdm)**3 + mass(kdm)*mass(kmed)**2)**2)
c... This corrects leading order result (28) in 1302.3898      
c         sigmavmoeller = 1./16./pi*gDM**4/mass(kdm)**2
c     &                   *sqrt(1.0d0-mass(kmed)**2/mass(kdm)**2) 
         if (p/mass(kdm).lt.5.d-5) then ! sv = a+ bv^2
            sigmavmoeller = sv0 + gDM**4*  
     -         (p**2*mass(kdm)*(76*mass(kdm)**6*mass(kmed)**2 - 
     -         142*mass(kdm)**4*mass(kmed)**4 + 
     -         89*mass(kdm)**2*mass(kmed)**6 - 23*mass(kmed)**8))/
     -         (24.*Pi*Sqrt(mass(kdm)**2 - mass(kmed)**2)*
     -        (-2*mass(kdm)**3 + mass(kdm)*mass(kmed)**2)**4)     
         else
            sigmavmoeller = gDM**4*  
     -        (-2*Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2)) - 
     -        (Log((s - 2*mass(kmed)**2 + 
     -        Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2)))**2/
     -        (-s + 2*mass(kmed)**2 + 
     -        Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2)))**2)*
     -        (s**2 - 8*mass(kdm)**4 + 4*mass(kmed)**4 + 
     -         mass(kdm)**2*(4*s - 8*mass(kmed)**2)))/
     -         (-s + 2*mass(kmed)**2) - 
     -        (2*Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2))*
     -         (2*mass(kdm)**2 + mass(kmed)**2)**2)/
     -       (mass(kmed)**4 + mass(kdm)**2*(s - 4*mass(kmed)**2)))/
     -       (16.*Pi*Sqrt(s*(s - 4*mass(kdm)**2))*(s - 2*mass(kdm)**2))
         endif
c... now add s-channel (annihilation to DR)
         sv0 = sv0 + gDM**2*gDR**2* 
     -          (mass(kdm)**2/  (Pi*(16*mass(kdm)**4 - 8*mass(kdm)**2*
     -          mass(kmed)**2 + mass(kmed)**4 + mass(kmed)**2*width(kmed)**2)))
         if (p/mass(kdm).lt.5.d-5) then ! sv = a+ bv^2
            sigmavmoeller = sigmavmoeller + gDM**2*gDR**2* 
     -          (-p**2*(112*mass(kdm)**4 - 32*mass(kdm)**2*mass(kmed)**2 + 
     -          mass(kmed)**4 + mass(kmed)**2*width(kmed)**2))/
     -          (3.*Pi*(16*mass(kdm)**4 - 8*mass(kdm)**2*mass(kmed)**2 + 
     -          mass(kmed)**4 + mass(kmed)**2*width(kmed)**2)**2) 
         else
            sigmavmoeller = sigmavmoeller + gDM**2*gDR**2* 
     -          (s*(s + 2*mass(kdm)**2))/(12.*Pi*(s - 2*mass(kdm)**2)*
     -          ((s - mass(kmed)**2)**2 + mass(kmed)**2*width(kmed)**2))            
         endif


c... DM fermion, scalar mediator
      elseif (DMtype.eq.1.and.Mediatortype.eq.3) then
         l = 1 ! p-wave

c... DM fermion, vector mediator: t- and u-channel from mathematica
         sv0 = gDM**4*(p**2*(-2*mass(kmed)**6 + 10*mass(kmed)**4*mass(kdm)**2 - 
     -         17*mass(kmed)**2*mass(kdm)**4 + 9*mass(kdm)**6))/
     -         (6.*Pi*mass(kdm)*(mass(kmed)**2 - 2*mass(kdm)**2)**4*
     -         Sqrt(-mass(kmed)**2 + mass(kdm)**2))
c... This corrects leading order result (28) in 1302.3898      
c          sv0 = 3.*gDM**4/64.0d0/pi*vrelcms**2/mass(kdm)**2
c     &                   *sqrt(1.0d0-mass(kmed)**2/mass(kdm)**2)  
         if (p/mass(kdm).lt.5.d-5) then ! sv = b v^2 + c v^4
            sigmavmoeller = sv0 + gDM**4*  
     -          (p**4*(-1248*mass(kdm)**10 + 3396*mass(kdm)**8*mass(kmed)**2 - 
     -          3536*mass(kdm)**6*mass(kmed)**4 + 
     -          1817*mass(kdm)**4*mass(kmed)**6 - 
     -          460*mass(kdm)**2*mass(kmed)**8 + 46*mass(kmed)**10))/
     -          (60.*Pi*mass(kdm)**3*Sqrt(mass(kdm)**2 - mass(kmed)**2)*
     -          (-2*mass(kdm)**2 + mass(kmed)**2)**6)
         else
            sigmavmoeller = gDM**4*  
     -          (-Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2)) + 
     -          (-4*mass(kdm)**2 + mass(kmed)**2)**2/(s - 2*mass(kmed)**2 + 
     -          Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2))) + 
     -          (-4*mass(kdm)**2 + mass(kmed)**2)**2/(-s + 2*mass(kmed)**2 + 
     -          Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2))) + 
     -          (Log((s - 2*mass(kmed)**2 + 
     -          Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2)))**2/
     -          (-s + 2*mass(kmed)**2 + 
     -          Sqrt((s - 4*mass(kdm)**2)*(s - 4*mass(kmed)**2)))**2)*
     -          (s**2/4. - 8*mass(kdm)**4 - s*mass(kmed)**2 + 
     -          (3*mass(kmed)**4)/2. + 4*mass(kdm)**2*(s - mass(kmed)**2)
     -          ))/(s - 2*mass(kmed)**2))/
     -          (16.*Pi*Sqrt(s*(s - 4*mass(kdm)**2))*(s - 2*mass(kdm)**2))
         endif
c... now add s-channel (annihilation to DR)
         sv0 = sv0 + gDM**2*gDR**2* 
     -          p**2/(32*Pi*mass(kdm)**4 - 16*Pi*mass(kdm)**2*mass(kmed)**2 + 
     -          2*Pi*mass(kmed)**4 + 2*Pi*mass(kmed)**2*width(kmed)**2)
         if (p/mass(kdm).lt.5.d-5) then ! sv = bv^2 + cv^4
            sigmavmoeller = sigmavmoeller + gDM**2*gDR**2* 
     -          (-p**4*(48*mass(kdm)**4 - 16*mass(kdm)**2*mass(kmed)**2 + 
     -          mass(kmed)**4 + mass(kmed)**2*width(kmed)**2))/
     -          (2.*Pi*mass(kdm)**2*(16*mass(kdm)**4 - 
     -          8*mass(kdm)**2*mass(kmed)**2 + mass(kmed)**4 + 
     -          mass(kmed)**2*width(kmed)**2)**2)
         else
            sigmavmoeller = sigmavmoeller + gDM**2*gDR**2* 
     -          (s*(s - 4*mass(kdm)**2))/(16.*Pi*(s - 2*mass(kdm)**2)*
     -          ((s - mass(kmed)**2)**2 + mass(kmed)**2*width(kmed)**2))            
         endif

         

            
      else
        write(*,*) 'ERROR in dsanwx: called with unsupported '
        write(*,*) 'DM/DR/mediator spin combination: ',  
     &              DMtype, Mediatortype, DRtype   
        write(*,*) 'program stopping...'
        stop
      endif
      
c... calulcate Sommerfeld factor      
      alpha = gDM**2/4.0d0/pi
      SE = dsansommerfeld(mass(kdm),mass(kmed),alpha,vrelcms,l)


c... final result (apply Sommerfeld factor only to leading term)     
      dsanwx= 2.*(s-2*mass(kdm)**2)*
     &   (sigmavmoeller + sv0*(SE - 1.0d0))

c TB debug.  
c      write(*,*) s,sigmavmoeller, mass(kdm),mass(kmed), alpha,SE

            
      return
      end
