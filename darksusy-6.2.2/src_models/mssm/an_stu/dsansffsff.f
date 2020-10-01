****************************************************
*** subroutine dsansffsff                        ***
*** fermion + fermion -> fermion + fermion in    ***
*** s-channel scalar exchange (index k)          ***
*** 1 - arrow in, 2 - arrow out, k intermediate  ***
*** this code is computer generated by reduce    ***
*** and gentran.                                 ***
*** author: joakim edsjo, edsjo@fysik.su.se        ***
****************************************************

      subroutine dsansffsff(p,costheta,kp1,kp2,kpk,kp3,kp4)
      implicit none
      include 'dsmssm.h'
      include 'dsandiacom.h'
      integer kp1,kp2,kpk,kp3,kp4
      real*8 p,costheta
      complex*16 dh,gpp2,gmm2,gpm2,gmp2

      call dsankinvar(p,costheta,kp1,kp2,kpk,kp3,kp4)
      if (s.lt.(mp3+mp4)**2) return
      dh=1.0d0/dcmplx(mk**2-s,-width(kpk)*mk)
      gpp2=(gl(kpk,kp2,kp1)+gr(kpk,kp2,kp1))*
     &     (gl(kpk,kp3,kp4)+gr(kpk,kp3,kp4))
      gmm2=(gl(kpk,kp2,kp1)-gr(kpk,kp2,kp1))*
     &     (gl(kpk,kp3,kp4)-gr(kpk,kp3,kp4))
      gpm2=(gl(kpk,kp2,kp1)+gr(kpk,kp2,kp1))*
     &     (gl(kpk,kp3,kp4)-gr(kpk,kp3,kp4))
      gmp2=(gl(kpk,kp2,kp1)-gr(kpk,kp2,kp1))*
     &     (gl(kpk,kp3,kp4)+gr(kpk,kp3,kp4))

      aa(0,0,0,0)=aa(0,0,0,0)-(2.0d0*dh*efpl*epl*gmm2)
      aa(0,0,1,0)=aa(0,0,1,0)-(2.0d0*dh*epl*gmp2*kpl)
      aa(1,0,0,0)=aa(1,0,0,0)+2.0d0*dh*efpl*gpm2*ppl
      aa(1,0,1,0)=aa(1,0,1,0)+2.0d0*dh*gpp2*kpl*ppl

      return

      end
