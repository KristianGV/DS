      function dsgammah(sqrts)
***
*** Returns the total Higgs decay width in GeV
***
*** Tabulated width from Dittmaier et al 1101.0539
*** Analytic expressions from Cline et al 1306.4710, Ilisie 2011
***
*** Author: Paolo Gondolo 2016
***
      implicit none
      real*8 dsgammah,sqrts
      real*8 dssmgammah, dsmwimp
      include 'dssilveira_zee.h'

      real*8 ans,aux
      
c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('silveira_zee','dsgammah')

c      write (*,*) 'PG dsgammah> entered with sqrts=',sqrts

      ans=dssmgammah(sqrts)

! add H -> S + S
      aux=1.d0-4.d0*(dsmwimp()/mass(khsm))**2
      if (aux.gt.0.d0) then
        ans = ans + (lambda*v0)**2/(32.d0*pi*mass(khsm))*sqrt(aux)        
      endif

      dsgammah = ans
c      write (*,*) 'PG dsgammah> dsgammah=',dsgammah
c      write (*,*) 'PG dsgammah> exited'
      
      return
      end
