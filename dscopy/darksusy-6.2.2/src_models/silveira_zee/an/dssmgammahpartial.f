      function dssmgammahpartial(ichannel,sqrts)
***
*** Returns the standard model partial Higgs decay width in GeV for channel ichannel
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
***
*** Tables are calculated with HDecay v. 6.51 
*** Analytic expressions from Cline et al 1306.4710, Ilisie 2011
***
*** Author: Paolo Gondolo 2016
*** Modified by Joakim Edsjo, Torsten Bringmann 2016 to allow using 
*** Hdecay tables (dssmgammah_hdecay_tab) instead of Dittmaier tables
***   (dssmgammahpartialtab). The former is set as default
***
      implicit none
      real*8 dssmgammahpartial,sqrts
      integer ichannel
      include 'dssilveira_zee.h'

      real*8 ans,x
      real*8 dssmgammah_options
      complex*16 cli2,dd
      integer i
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dssmgammahpartial')

c      write (*,*) 'PG dssmgammahpartial> entered'
c      write (*,*) 'PG dssmgammahpartial> ichannel=',ichannel
c      write (*,*) 'PG dssmgammahpartial> sqrts=',sqrts
      
      ans=0.0d0
! ll
      if (ichannel.ge.1.and.ichannel.le.6) then
        if ((ichannel.eq.4.or.ichannel.eq.6).and.sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
          ans=dssmgammah_options(ichannel,sqrts)
        else
          i=ichannel
          x=(mass(i)/sqrts)**2
          if (1.d0-4.d0*x.gt.0.d0) then
            ans=sqrts*mass(i)**2/8.d0/pi/v0**2*sqrt(1.d0-4.d0*x)**3
          endif
c          write (*,*) 'PG dssmgammahpartial> inside channel ',ichannel,' : x=',x,' ans=',ans
        endif
! qq
      else if (ichannel.ge.7.and.ichannel.le.12) then
        if (ichannel.ge.9.and.sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
          ans=dssmgammah_options(ichannel,sqrts)
        else
          i=ichannel
          x=(mass(i)/sqrts)**2
          if (1.d0-4.d0*x.gt.0.d0) then
            ans=sqrts*mass(i)**2/8.d0/pi/v0**2*sqrt(1.d0-4.d0*x)**3*3.d0
            if (1.5d0*log(x)+2.25d0.gt.0.d0) then
              ans=ans*(1.d0+(1.5d0*log(x)+2.25d0)*4.d0*alph3mz_def/3.d0/pi)
            endif
          endif
c          write (*,*) 'PG dssmgammahpartial> inside channel ',ichannel,' : x=',x,' ans=',ans
        endif
! gamma+gamma and Z+gamma
      else if (ichannel.eq.13.or.ichannel.eq.18) then
        if (sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
          ans=dssmgammah_options(ichannel,sqrts)
        else
c          write (*,*) 'Warning in dssmgammahpartial: H->gamma+gamma,Z+gamma unavailable at sqrts=',sqrts
c          write (*,*) '   dssmgammahpartial will be set to 0'
          ans=0.d0
        endif
! WW
      else if (ichannel.eq.14) then
        if (sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
          ans=dssmgammah_options(ichannel,sqrts)
        else
          x=(mass(kw)/sqrts)**2
          if (1.d0-4.d0*x.gt.0.d0) then
            ans=sqrts**3/16.d0/pi/v0**2*sqrt(1.d0-4.d0*x)*(1.d0+x*(-4.d0+x*12.d0))
          endif
        endif
! ZZ
      else if (ichannel.eq.15) then
        if (sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
          ans=dssmgammah_options(ichannel,sqrts)
        else
          x=(mass(kz)/sqrts)**2
          if (1.d0-4.d0*x.gt.0.d0) then
            ans=sqrts**3/32.d0/pi/v0**2*sqrt(1.d0-4.d0*x)*(1.d0+x*(-4.d0+x*12.d0))
          endif
        endif
! gg
      else if (ichannel.eq.16) then ! JE FIX. Do we want to use Ilisie or Hdecay
        !if (sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
        !  ans=dssmgammah_options(ichannel,sqrts)
        !else
          ! Ilisie
          ans=0.d0
          do i=7,12
            x=(mass(i)/sqrts)**2
            dd=-2.d0+(4.d0*x-1.d0)*(cli2(-2.d0/(zsqrt(dcmplx(1.d0-4.d0*x))-1.d0))+cli2(2.d0/(zsqrt(dcmplx(1.d0-4.d0*x))+1.d0)))
            ans=ans+sqrts**3/8.d0/pi/v0**2*(alph3mz_def/pi)**2*x**2*(dreal(dd)**2+dimag(dd)**2)
          enddo
        !endif
! hh
      else if (ichannel.eq.17) then
          write (*,*) 'Warning in dssmgammahpartial : H -> H + H unavailable.'
          write (*,*) '   dssmgammahpartial will be set to 0'
          ans=0.d0
c! Zgamma
c      else if (ichannel.eq.18) then
c        if (sqrts.ge.1.d0.and.sqrts.le.1020.d0) then
c          ans=dssmgammah_options(ichannel,sqrts)
c        else
c          write (*,*) 'Warning in dssmgammahpartial : H -> Z + gamma unavailable outside 1GeV<sqrts<1020GecV.'
cc          write (*,*) '   dssmgammahpartial will be set to 0'
c          ans=0.d0
c        endif
! non-existent ichannel
      else
        write (*,*) 'Warning in dssmgammahpartial : non-existent ichannel = ',ichannel
        write (*,*) '   dssmgammahpartial will be set to 0'
        ans=0.d0
      endif

      dssmgammahpartial = ans
c      write (*,*) 'PG dssmgammahpartial> dssmgammahpartial=',dssmgammahpartial
c      write (*,*) 'PG dssmgammahpartial> exited'
      
      return
      end

      
c... auxiliary function to return Higgs widths based on global settings, TB 21/12/2016      
      real*8 function dssmgammah_options(ichannel,sqrts)
        implicit none
        real*8 dssmgammahpartial,sqrts
        integer ichannel
        include 'dssilveira_zee.h'
      
        real*8 dssmgammah_hdecay_tab, dssmgammahpartialtab
        
        dssmgammah_options = 0.d0
        
        if (gammahhow.eq.1) then
          dssmgammah_options = dssmgammah_hdecay_tab(ichannel,sqrts)
        elseif (gammahhow.eq.2) then
          dssmgammah_options = dssmgammahpartialtab(ichannel,sqrts)
        else
          write(*,*) 'ERROR in dssmgammahpartial: undefined global option gammahhow = ',gammahhow
          stop 'Program stopping...'
        endif
      

      end ! dssmgammah_options




c$$$      function dsmqrunmsbar(i,mu)
c$$$c... Running quark mass in the MSbar scheme at 2 loops from Gorishny et al,
c$$$c... Kataev et al, Surguladze, and Chetyrkin as quoted in Djouadi 9712334
c$$$      implicit none
c$$$      real*8 dsmqrunmsbar,mu
c$$$      integer i
c$$$      data 
c$$$      if (i.lt.7.or.i.gt.12) then
c$$$        write (*,*) 'ERROR in dsmqrunmsbar : called with wrong quark index i=',i
c$$$      else if (i.eq.7) then
c$$$        mq=2.d0
c$$$        mqmq=mass_up_def
c$$$      else if (i.eq.8) then
c$$$        mq=2.d0
c$$$        mqmq=mass_down_def
c$$$      else if (i.eq.9) then
c$$$        mq=mcmc
c$$$        mqmq=mass_charm_def
c$$$      else if (i.eq.10) then
c$$$        mq=2.d0
c$$$        mqmq=mass_strange_def
c$$$      else if (i.eq.11) then
c$$$        mq=mtmt
c$$$        mqmq=mass_top_def
c$$$      else if (i.eq.12) then
c$$$        mq=mbmb
c$$$        mqmq=mass_bottom_def
c$$$      endif
c$$$      xmu=dsa
c$$$      dsmqrunmsbar=mqmq*dsmqrunc(xmu)/dsmqrunc(xmq)
c$$$      return
c$$$      end
      
