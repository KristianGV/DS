********************************************************
*** subroutine dsansffvvs                            ***
*** fermion + fermion -> scalar + gauge boson in     ***
*** s-channel gauge boson exchange (index k)         ***
*** 1 - arrow in, 2 - arrow out, k intermediate      ***
*** this code is computer generated by reduce        ***
*** and gentran.                                     ***
*** author: joakim edsjo, edsjo@fysik.su.se            ***
********************************************************

      subroutine dsansffvvs(p,costheta,kp1,kp2,kpk,kp3,kp4)
      implicit none
      include 'dsmssm.h'
      include 'dsandiacom.h'
      integer kp1,kp2,kpk,kp3,kp4
      real*8 p,costheta
      complex*16 dv,gv2,ga2

      call dsankinvar(p,costheta,kp1,kp2,kpk,kp3,kp4)
      if (s.lt.(mp3+mp4)**2) return
      dv=1.0d0/dcmplx(mk**2-s,-width(kpk)*mk)
      gv2=gl(kp4,kp3,kpk)*
     &  (gl(kpk,kp2,kp1)+gr(kpk,kp2,kp1))
      ga2=gl(kp4,kp3,kpk)*
     &  (gl(kpk,kp2,kp1)-gr(kpk,kp2,kp1))

      if (kpk.eq.kgamma) then
        if (kp3.eq.kgamma) then

      aa(0,0,-1,0)=aa(0,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*ga2*mw*pmi
     . )
      aa(0,0,1,0)=aa(0,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*ga2*mw*pmi)
      aa(1,-1,-1,0)=aa(1,-1,-1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1
     . ,-1,-1)*dv*mw)
      aa(1,-1,1,0)=aa(1,-1,1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1,
     . -1,1)*dv*mw)
      aa(1,0,-1,0)=aa(1,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*emi*gv2*mw
     . )
      aa(1,0,1,0)=aa(1,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*emi*gv2*mw)
      aa(1,1,-1,0)=aa(1,1,-1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1
     . ,1,-1)*dv*mw)
      aa(1,1,1,0)=aa(1,1,1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1,1
     . ,1)*dv*mw)

        else

      aa(0,0,-1,0)=aa(0,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*ga2*mw*pmi
     . )
      aa(0,0,0,0)=aa(0,0,0,0)-(dsqrt(2.0d0)*(wd(1,0,0)*e3*pmi+emi*kk)
     . *dv*ga2*mw)/mp3
      aa(0,0,1,0)=aa(0,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*ga2*mw*pmi)
      aa(1,-1,-1,0)=aa(1,-1,-1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1
     . ,-1,-1)*dv*mw)
      aa(1,-1,0,0)=aa(1,-1,0,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1,
     . -1,0)*dv*e3*mw)/mp3
      aa(1,-1,1,0)=aa(1,-1,1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1,
     . -1,1)*dv*mw)
      aa(1,0,-1,0)=aa(1,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*emi*gv2*mw
     . )
      aa(1,0,0,0)=aa(1,0,0,0)-(dsqrt(2.0d0)*(wd(1,0,0)*e3*emi+kk*pmi)
     . *dv*gv2*mw)/mp3
      aa(1,0,1,0)=aa(1,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*emi*gv2*mw)
      aa(1,1,-1,0)=aa(1,1,-1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1
     . ,1,-1)*dv*mw)
      aa(1,1,0,0)=aa(1,1,0,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1,
     . 1,0)*dv*e3*mw)/mp3
      aa(1,1,1,0)=aa(1,1,1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1,1
     . ,1)*dv*mw)

        endif
      else
        if (kp3.eq.kgamma) then

      aa(0,0,-1,0)=aa(0,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*ga2*mw*pmi
     . )
      aa(0,0,1,0)=aa(0,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*ga2*mw*pmi)
      aa(1,-1,-1,0)=aa(1,-1,-1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1
     . ,-1,-1)*dv*mw)
      aa(1,-1,1,0)=aa(1,-1,1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1,
     . -1,1)*dv*mw)
      aa(1,0,-1,0)=aa(1,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*emi*gv2*mw
     . )
      aa(1,0,1,0)=aa(1,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*emi*gv2*mw)
      aa(1,1,-1,0)=aa(1,1,-1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1
     . ,1,-1)*dv*mw)
      aa(1,1,1,0)=aa(1,1,1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1,1
     . ,1)*dv*mw)

        else

      aa(0,0,-1,0)=aa(0,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*ga2*mw*pmi
     . )
      aa(0,0,0,0)=aa(0,0,0,0)+(dsqrt(2.0d0)*((e1+e2+mk)*(e1+e2-mk)*emi
     . *kk-(wd(1,0,0)*e3*mk**2*pmi))*dv*ga2*mw)/(mk**2*mp3)
      aa(0,0,1,0)=aa(0,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*ga2*mw*pmi)
      aa(1,-1,-1,0)=aa(1,-1,-1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1
     . ,-1,-1)*dv*mw)
      aa(1,-1,0,0)=aa(1,-1,0,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1,
     . -1,0)*dv*e3*mw)/mp3
      aa(1,-1,1,0)=aa(1,-1,1,0)-(dsqrt(2.0d0)*(epl*gv2+ga2*ppl)*wd(1,
     . -1,1)*dv*mw)
      aa(1,0,-1,0)=aa(1,0,-1,0)-(dsqrt(2.0d0)*wd(1,0,-1)*dv*emi*gv2*mw
     . )
      aa(1,0,0,0)=aa(1,0,0,0)+(dsqrt(2.0d0)*((e1+e2+mk)*(e1+e2-mk)*kk*
     . pmi-(wd(1,0,0)*e3*emi*mk**2))*dv*gv2*mw)/(mk**2*mp3)
      aa(1,0,1,0)=aa(1,0,1,0)-(dsqrt(2.0d0)*wd(1,0,1)*dv*emi*gv2*mw)
      aa(1,1,-1,0)=aa(1,1,-1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1
     . ,1,-1)*dv*mw)
      aa(1,1,0,0)=aa(1,1,0,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1,
     . 1,0)*dv*e3*mw)/mp3
      aa(1,1,1,0)=aa(1,1,1,0)-(dsqrt(2.0d0)*(epl*gv2-(ga2*ppl))*wd(1,1
     . ,1)*dv*mw)

        endif
      endif

      return

      end














