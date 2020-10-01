      real*8 function dsansumaa()
c_______________________________________________________________________
c  sum the amplitude matrix
c  author: joakim edsjo (edsjo@physto.se) 96-02-02
c          paolo gondolo 99-1-15 factor of 3 faster
c  called by: dwdcos
c mod TB 2016-02-08: added OMP compatibility
c=======================================================================
      implicit none
c     include 'dsandiacom.h'

      real*8 aaaa(108)
      real*8 mp1,mp2,mk,mp3,mp4,ppl,pmi,epl,emi,mw,mz,pp,kpl,
     &  kmi,efpl,efmi,kk,e1,e2,e3,e4,s,t,u,wd(2,-2:2,-2:2)
      common /diacom/ aaaa,
     &  mp1,mp2,mk,mp3,mp4,ppl,pmi,epl,emi,mw,mz,pp,kpl,kmi,
     &  Efpl,efmi,kk,e1,e2,e3,e4,s,t,u,wd
      save /diacom/

*$OMP THREADPRIVATE (/diacom/)


      real*8 sumf
      integer i
c      equivalence (aa,aaaa) ! TB: not compatible with !OMP !

      sumf=0.0d0

      do i=1,108
         sumf = sumf + aaaa(i)**2
      enddo
      
      dsansumaa = sumf

      end
