      real*8 function dsgammahpartial(ichannel,sqrts)
***
*** Returns the partial Higgs decay width in GeV for channel ichannel
***
*** List of channels:
***   ichannel = 1  :  nue + anti-nue
***   ichannel = 2  :  e+ + e-
***   ichannel = 3  :  numu + anti-numu
***   ichannel = 4  :  mu+ + mu-
***   ichannel = 5  :  nutau + anti-nutau
***   ichannel = 6  :  tau+ + tau-
***   ichannel = 7  :  u + ubar
***   ichannel = 8  :  d + dbar
***   ichannel = 9  :  c + cbar
***   ichannel = 10  :  s + sbar
***   ichannel = 11  :  t + tbar
***   ichannel = 12  :  b + bbar
***   ichannel = 13  :  gamma + gamma
***   ichannel = 14  :  W+ + W-
***   ichannel = 15  :  Z + Z
***   ichannel = 16  :  g + g
***   ichannel = 17  :  H + H
***   ichannel = 18  :  gamma + Z
***   ichannel = 19  :  S + S
***
*** Author: Paolo Gondolo 2016
***
      implicit none
      real*8 sqrts
      integer ichannel
      include 'dssilveira_zee.h'

      real*8 ans,dssmgammahpartial,aux, dsmwimp
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dsgammahpartial')

c      write (*,*) 'PG dssmgammahpartial> entered'
c      write (*,*) 'PG dssmgammahpartial> ichannel=',ichannel
c      write (*,*) 'PG dssmgammahpartial> sqrts=',sqrts

      ans=0.0d0
! SM channel
      if (ichannel.ge.1.and.ichannel.le.18) then
        ans=dssmgammahpartial(ichannel,sqrts)
      else if (ichannel.eq.19) then
! H -> S + S
        aux=1.d0-4.d0*(dsmwimp()/mass(khsm))**2
        if (aux.gt.0.d0) then
          ans = (lambda*v0)**2/(32.d0*pi*mass(khsm))*sqrt(aux)        
        endif
      else
! non-existent ichannel
        write (*,*) 'Warning in dsgammahpartial : non-existent ichannel = ',ichannel
        write (*,*) '   dsgammahpartial will be set to 0'
        ans=0.d0
      endif

      dsgammahpartial = ans
c      write (*,*) 'PG dsgammahpartial> dsgammahpartial=',dsgammahpartial
c      write (*,*) 'PG dsgammahpartial> exited'
      
      return
      end

