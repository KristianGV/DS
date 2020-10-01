**********************************************************************
*** function dssigmav returns the partial annihilation cross section
*** sigma v for WIMP-WIMP annihilation into channel ichannel at
*** center-of-mass momentum p.
*** Here v is defined to be the relative velocity of one particle in
*** the frame of the other.
***
*** Notice that in 1306.4710, their sigma v_rel is related to the invariant rate W via
***  W = s sigma v_rel, from comparison with the invariant rate in their formula for <sigma v>
*** From the general relation
***  W = 4 p1.p2 sigma v = 2 (s-m1^2-m2^2) sigma v = 4 (p^2+sqrt(p^2+m1^2)*sqrt(p^2+m2^2)) sigma v
*** one obtains
***  sigma v = s/4/(p^2+sqrt(p^2+m1^2)*sqrt(p^2+m2^2)) sigma v_rel
*** and for m1=m2==mx
***  sigma v = W/4/(mx^2+2*p^2)
***          = (mx^2+p^2)/(mx^2+2*p^2) sigma v_rel
***                                                                         
***  type : INTERFACE                                                        
***                                                                         
*** Units of returned cross section: cm^3 s^-1
***
*** List of channels:
***   ichannel = 1  :  nue + anti-nue (outside 90GeV-300GeV only)
***   ichannel = 2  :  e+ + e- (outside 90GeV-300GeV only)
***   ichannel = 3  :  numu + anti-numu (outside 90GeV-300GeV only)
***   ichannel = 4  :  mu+ + mu-
***   ichannel = 5  :  nutau + anti-nutau (outside 90GeV-300GeV only)
***   ichannel = 6  :  tau+ + tau-
***   ichannel = 7  :  u + ubar (outside 90GeV-300GeV only)
***   ichannel = 8  :  d + dbar (outside 90GeV-300GeV only)
***   ichannel = 9  :  c + cbar
***   ichannel = 10  :  s + sbar
***   ichannel = 11  :  t + tbar
***   ichannel = 12  :  b + bbar
***   ichannel = 13  :  gamma + gamma (90GeV-300GeV only)
***   ichannel = 14  :  W+ + W-
***   ichannel = 15  :  Z + Z
***   ichannel = 16  :  g + g (90GeV-300GeV only)
***   ichannel = 17  :  H + H
***   ichannel = 18  :  Z + gamma (90GeV-300GeV only)
***
*** author: Paolo Gondolo
*** date: 2016
**********************************************************************

c... FIXME: eventually, it would be nice to take PDG codes rather than uchannel as input here...

      function dssigmavpartial(ichannel,p)
      implicit none
      include 'dssilveira_zee.h'

      real*8 dssigmavpartial,p
      integer ichannel
      real*8 kin,mx,s,sqrts,sv,dh,dsgammah,dssmgammahpartial,dsmwimp
      real*8 mh,vs,vh,ar,ai,aux,x,xx,y,auxl
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dssigmavpartial')

c      write (*,*) 'PG dssigmavpartial> entered'
c      write (*,*) 'PG dssigmavpartial> ichannel=',ichannel
c      write (*,*) 'PG dssigmavpartial> sqrts=',sqrts

      dssigmavpartial=0.0d0
      if (p.lt.0.d0) return

      sv=0.d0
      mx=dsmwimp()
      mh=mass(khsm)
      s=4.d0*(mx**2+p**2)
      kin=(mx**2+p**2)/(mx**2+2.d0*p**2)
      sqrts=sqrt(s)
c      dh=1.d0/((s-mh**2)**2+mh**2*width(khsm))
      dh=1.d0/((s-mh**2)**2+mh**2*width(khsm)**2) !TB corr
! SM channel except H+H
      if (ichannel.ge.1.and.ichannel.le.18.and.ichannel.ne.17) then
        sv=2.d0/sqrts*(lambda*v0)**2*dh*dssmgammahpartial(ichannel,sqrts)
! H + H
      else if (ichannel.eq.17) then
        vh=(p**2+mx**2-mh**2)/(mx**2+p**2)
        if (vh.gt.0.d0) then
          vh=sqrt(vh)
          vs=p/sqrt(mx**2+p**2)
          ar=1.d0+3.d0*mh**2*(s-mh**2)*dh
          ai=3.d0*mh**2*sqrts*width(khsm)*dh
          x=2.d0*vh*vs/(1.d0+vh**2) ! 0 <= x < 1
          y=lambda*v0**2/(s-2.d0*mh**2)
          if (x.gt.1.d-5) then
            auxl=log((1.d0+x)/(1.d0-x))/(2.d0*x)
          else
            xx=x*x
            auxl=1.d0+xx*(1.d0/3.d0+xx/5.d0)
          endif
          aux=ar**2+ai**2+8.d0*y**2/(1.d0-x**2)-8.d0*y*(ar-y)*auxl
c          sv=lambda**2*vh/8.d0/pi/s*aux
          sv=lambda**2*vh/16.d0/pi/s*aux !TB corr
        endif
! non-existent ichannel
      else
        write (*,*) 'Warning in dssigmav : non-existent ichannel = ',ichannel
        write (*,*) '   sigma v will be set to 0'
      endif

      dssigmavpartial = gev2cm3s*kin*sv
      
c      write (*,*) 'PG dssigmavpartial> dssigmavpartial=',dssigmavpartial
c      write (*,*) 'PG dssigmavpartial> exited'
      
      return
      end
