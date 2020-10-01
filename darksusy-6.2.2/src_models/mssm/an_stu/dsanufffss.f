****************************************************
*** subroutine dsanufffss                        ***
*** fermion + fermion -> scalar + scalar in      ***
*** u-channel fermion exchange (index k)         ***
*** 1 - arrow in, 2 - arrow out, k intermediate  ***
*** this code is computer generated by reduce    ***
*** and gentran.                                 ***
*** author: joakim edsjo, edsjo@fysik.su.se        ***
****************************************************

      subroutine dsanufffss(p,costheta,kp1,kp2,kpk,kp3,kp4)
      implicit none
      include 'dsmssm.h'
      include 'dsandiacom.h'
      integer kp1,kp2,kpk,kp3,kp4
      real*8 p,costheta
      complex*16 du,gp2,gm2,gx2,gy2

      call dsankinvar(p,costheta,kp1,kp2,kpk,kp3,kp4)
      if (s.lt.(mp3+mp4)**2) return
      du=1.0d0/dcmplx(mk**2-u,-width(kpk)*mk)
      gp2=gl(kp3,kpk,kp2)*gl(kp4,kp1,kpk)+
     &   gr(kp3,kpk,kp2)*gr(kp4,kp1,kpk)
      gm2=gl(kp3,kpk,kp2)*gl(kp4,kp1,kpk)-
     &   gr(kp3,kpk,kp2)*gr(kp4,kp1,kpk)
      gx2=gl(kp3,kpk,kp2)*gr(kp4,kp1,kpk)+
     &   gr(kp3,kpk,kp2)*gl(kp4,kp1,kpk)
      gy2=gl(kp3,kpk,kp2)*gr(kp4,kp1,kpk)-
     &   gr(kp3,kpk,kp2)*gl(kp4,kp1,kpk)


      aa(0,0,0,0)=aa(0,0,0,0)-(dsqrt(2.0d0)*(wd(1,0,0)*gy2*kk*pmi-(e2*
     . emi*gy2)+e3*emi*gy2+epl*gm2*mk+gy2*pmi*pp)*du)
      aa(1,-1,0,0)=aa(1,-1,0,0)+dsqrt(2.0d0)*(epl*gx2-(gy2*ppl))*wd(1,
     . -1,0)*du*kk
      aa(1,0,0,0)=aa(1,0,0,0)+dsqrt(2.0d0)*(wd(1,0,0)*emi*gx2*kk-(e2*
     . gx2*pmi)+e3*gx2*pmi+emi*gx2*pp+gp2*mk*ppl)*du
      aa(1,1,0,0)=aa(1,1,0,0)+dsqrt(2.0d0)*(epl*gx2+gy2*ppl)*wd(1,1,0)
     . *du*kk

      return

      end
