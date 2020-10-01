*******************************************************************************
*** function dsansommerfeld returns the Sommerfeld correction factor S      ***
*** for a p-wave and s-wave process that results from the multiple          ***
*** exchange of a single vector or scalar mediator type. The total          ***
*** crosss ection is thus given as sv = S * (sv)_0, where (sv)_0 is the     ***
*** tree level result.                                                      ***
*** This implements Eq. (34) in 1302.3898, based on S. Cassel,              ***
*** J. Phys. G 37, 105009 (2010) [arXiv:0903.5307 [hep-ph]]                 ***
***                                                                         ***
*** input:                                                                  ***
***                                                                         ***
***    mDM   - DM mass (in GeV)                                             ***
***    mmed  - vector mediator mass (in GeV)                                ***
***    alpha - (DM-mediator coupling)^2 /(4pi)                              ***
***    vrel  - relatove velocity between DM particles /in CMS)              ***
***    l     - determines partial wave [0:s-wave, 1:p-wave]                 ***
***                                                                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date: 2018-02-21                                                        ***
*******************************************************************************

      real*8 function dsansommerfeld(mDM,mmed,alpha,vrel,l)
      implicit none
      include 'dsmpconst.h'
      include 'dsio.h'

      real*8 mDM,mmed,alpha,vrel
      integer l
      logical dsisnan
      
      real*8 a,c, S
     
      dsansommerfeld=1.0d0
      
      if (vrel.le.0.or.vrel.ge.1.d0.or.alpha.gt.1d2) return
      
      a = vrel / (2.*alpha)
      c = 6.*alpha*mDM/pi**2/mmed

c... original, analytic form
c      S = (Pi*Sinh(2.*a*c*Pi))/
c     -    (a*(-Cos(2.*Sqrt(c - a**2*c**2)*Pi) + Cosh(2.*a*c*Pi)))         
c... rewrite this for numerical stability

      if (2.*a*c*Pi.lt.4.0d1) then
         S = Pi*Tanh(2.*a*c*Pi)/a
      else
         S = Pi/a
      endif      
            
      if (c.ge.(a*c)**2) then
        if (2.*a*c*Pi.lt.4.0d1) ! catches too large arguments
     &    S = S / (1.0d0 - Cos(2.*Sqrt(c - (a*c)**2)*Pi) / Cosh(2.*a*c*Pi)) 
      else
         if (2.*a*c*Pi.lt.3.0d1.and.(a*c)**2.lt.1.0d2*c) then
            S = S / (1.0d0 - Cosh(2.*Sqrt((a*c)**2-c)*Pi) / Cosh(2.*a*c*Pi)) 
            if (dsisnan(S)) write(*,*) 'B'
         else ! here we are in the clasical regime -- with many extremely
              ! closely spaced resonances. This is numerically hard to 
              ! handle, so we revert to the limit of a massless mediator. 
            if (a.gt.1.0d3) then
              S = 1.0d0 + 0.5*pi/a
            else   
              S = pi/a / (1.0d0 - exp(-pi/a) )
            endif  
         endif          
      endif 

c... TB debug
c            S = pi/a / (1.0d0 - exp(-pi/a) )

c... TB debug
c    MISSING: unitarity bound !


      if (dsisnan(S)) then
c TB debug. 
        write(*,*) ' dsansommerfeld returned NaN !'
        write(*,*) mDM,mmed,alpha,vrel
        write(*,*) a,c,Sinh(2.*a*c*Pi),Cosh(2.*a*c*Pi)
        write(*,*) (c - a**2*c**2), Cos(2.*Sqrt(c - a**2*c**2)*Pi), 
     &             Cosh(2.*Sqrt((a*c)**2-c)*Pi) / Cosh(2.*a*c*Pi), 1.0d0 - exp(-pi/a)
        write(*,*) '(hit ENTER to set Sommerfeld factor to 1 and continue)'
        read(*,*)
        S=1.0d0
      endif

      if (l.eq.0) then ! s-wave
      
         dsansommerfeld = S
      
      elseif (l.eq.1) then ! p-wave
      
         dsansommerfeld= S*(((-1.0d0 + c)**2 + 4.*a**2*c**2))/(1.0d0 + 4.*a**2*c**2)
      
      else 
        if (prtlevel.gt.0) then
          write(*,*)
     &    'WARNING in dsansommerfeld: called with unsupported option l =',l        
          write(*,*) 'Setting Sommerfeld factor to 1.0'
        endif
      endif
      
     
      return
      end

      
