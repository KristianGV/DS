      subroutine dsrddof(t,sqrtgstar,heff)
c_______________________________________________________________________
c     return effective degrees of freedom at temperature t (in GeV)
c
c         sqrtgstar - g_*^0.5 as defined in (2.24) of Gelmini&Gondolo 1991   
c         heff - entropy d.o.f: s=heff*(2*pi^2/45)*T^3 
c
c     common:
c       'dsrdcom.h' - included common blocks
c  author: paolo gondolo (paolo@physics.utah.edu) 2005
c  modified: paolo gondolo (paolo@physics.utah.edu) 2008
c  modified: torsten bringmann 2017 -- extended header+moved geff to separate routine in /kd
c=======================================================================
      implicit none
      real*8 t,sqrtgstar,heff
      include 'dsrdcom.h'
      integer k,kk

      if (t.ge.tgev(1)) then
         sqrtgstar = fg(1)
         heff = fh(1)

      elseif (t.le.tgev(nf)) then ! TB add to get correct results for small T!
         sqrtgstar = fg(nf)
         heff = fh(nf)

      else

         kk=1
 100     if (tgev(khi).gt.t) then
            khi=khi+kk
            kk=2*kk
            if (khi.le.nf) goto 100
            khi=nf
         endif
 110     if (tgev(klo).lt.t) then
            klo=klo-kk
            kk=2*kk
            if (klo.ge.1) goto 110
            klo=1
         endif
 120     if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if (tgev(k).lt.t) then
               khi=k
            else
               klo=k
            endif
            goto 120
         endif
         
         sqrtgstar = (t*(fg(khi)-fg(klo))+
     &        fg(klo)*tgev(khi)-fg(khi)*tgev(klo))/(tgev(khi)-tgev(klo))
         heff = (t*(fh(khi)-fh(klo))+
     &        fh(klo)*tgev(khi)-fh(khi)*tgev(klo))/(tgev(khi)-tgev(klo))
      endif

      end
