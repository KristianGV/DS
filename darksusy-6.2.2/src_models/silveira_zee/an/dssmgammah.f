      function dssmgammah(sqrts)
***
*** Returns the total Standard-Model Higgs decay width in GeV
***
*** Tabulated width from Dittmaier et al 1101.0539
*** Analytic expressions from Cline et al 1306.4710, Ilisie 2011
***
*** Author: Paolo Gondolo 2016
***
      implicit none
      real*8 dssmgammah,sqrts
      include 'dssilveira_zee.h'

      real*8 ans,dssmgammahtab,dssmgammahpartial
      integer ich
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dssmgammah')

c      write (*,*) 'PG dssmgammah> entered with sqrts=',sqrts

      ans=0.0d0
      if (sqrts.ge.90.d0.and.sqrts.le.300.d0) then
        ans=dssmgammahtab(sqrts)
      else
        ans=0.d0
        do ich=1,18
          ans = ans + dssmgammahpartial(ich,sqrts)
        enddo
      endif

      dssmgammah = ans
c      write (*,*) 'PG dssmgammah> dssmgammah=',dssmgammah
c      write (*,*) 'PG dssmgammah> exited'
      
      return
      end
