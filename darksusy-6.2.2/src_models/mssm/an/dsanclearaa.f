      subroutine dsanclearaa
c_______________________________________________________________________
c  clear the amplitude matrix
c  author: joakim edsjo (edsjo@physto.se) 95-10-25
c          paolo gondolo 99-1-15 factor of 3.7 faster
c          torsten bringmann: added OMP compatibility
c  called by: dwdcos
c=======================================================================
      implicit none
c      include 'dsandiacom.h'

      real*8 aaaa(108)
      real*8 mp1,mp2,mk,mp3,mp4,ppl,pmi,epl,emi,mw,mz,pp,kpl,
     &  kmi,efpl,efmi,kk,e1,e2,e3,e4,s,t,u,wd(2,-2:2,-2:2)
      common /diacom/ aaaa,
     &  mp1,mp2,mk,mp3,mp4,ppl,pmi,epl,emi,mw,mz,pp,kpl,kmi,
     &  Efpl,efmi,kk,e1,e2,e3,e4,s,t,u,wd
      save /diacom/

*$OMP THREADPRIVATE (/diacom/)


      integer i
c      equivalence (aa,aaaa) ! TB: not compatible with $OMP !

      do i=1,108
         aaaa(i)=0.0d0
      enddo

      end
