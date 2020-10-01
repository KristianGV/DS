*******************************************************************
*** subroutine dsib2ffZtuFSR add to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel FSR diagrams. 
*** The Z is emitted by final legs 1 (FSR1) or 2 (FSR2) of the t+u channel diagrams.
*** isf1 is the index if the exchanged sfermion in th t-channel
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffZtuFSR(f, leg, isf1)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf1, leg
      
      integer sf ! (2)
      real*8 Msf ! (2)
      complex*16 delta
      complex*16 CFSR1L, CFSR1R, CFSR2L, CFSR2R, CFSR3L, CFSR3R
      complex*16 amp00, amp01, amp02, amp03, 
     &     amp10, amp11, amp12, amp13, amp20, amp21, amp22, amp23

      real*8 q1 
      integer dsib2getsfermion
      
      sf  = dsib2getsfermion(f, isf1)
      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics
 
c TB: obsolete code after introducing dsib2getsfermion    
c      if(f.EQ.knue) sf(1) = ksnue
c      if(f.EQ.knumu) sf(1) = ksnumu
c      if(f.EQ.knutau) sf(1) = ksnutau
c      if(f.EQ.ke) sf(1) = kse(1)
c      if(f.EQ.kmu) sf(1) = ksmu(1)
c      if(f.EQ.ktau) sf(1) = kstau(1)
c      if(f.EQ.ku) sf(1) = ksu(1)
c      if(f.EQ.kd) sf(1) = ksd(1)
c      if(f.EQ.kc) sf(1) = ksc(1)
c      if(f.EQ.ks) sf(1) = kss(1)
c      if(f.EQ.kb) sf(1) = ksb(1)
c      if(f.EQ.kt) sf(1) = kst(1)
c      if(f.eq.knue.or.f.eq.knumu.or.f.eq.knutau) then
c         sf(2) = sf(1)
c      else
c         sf(2)=sf(1)+1
c      endif
c... masses of exchanged particle are divided by mx in order to be consistent with the kinematics
c      Msf(1) = mass(sf(1))/mx
c      Msf(2) = mass(sf(2))/mx

c... Defines FSR1 and FSR2 amplitudes through their simmetry relations via q1 exchange
c      if (leg.eq.1) then
         q1 = kJ
c      else 
      if (leg.eq.2) q1 = -kJ

c... denominators
      delta = (
     -     dcmplx(2*E1J*EvJ + MB**2 - 2*CW*kv*q1,Mf*width(f)/mx)*
     -     dcmplx(-E1J**2 + Mf**2 + MB**2/4. - CW*kv*q1 - Msf**2,
     -     Msf*width(sf)/mx*0d0))

c... couplings
      CFSR1L = Conjg(gl(sf,f,kn(1)))*gl(kz,f,f)*
     -     gl(sf,f,kn(1))
      CFSR2L = Conjg(gl(sf,f,kn(1)))*gl(kz,f,f)*
     -     gr(sf,f,kn(1)) ! same as Conjg(gr(sf,f,kn(1)))*gl(kz,f,f)*gl(sf,f,kn(1))
      CFSR3L = Conjg(gr(sf,f,kn(1)))*gl(kz,f,f)*
     -     gr(sf,f,kn(1))
      CFSR1R = Conjg(gl(sf,f,kn(1)))*gl(sf,f,kn(1))*
     -     gr(kz,f,f)
      CFSR2R = Conjg(gl(sf,f,kn(1)))*gr(kz,f,f)*
     -     gr(sf,f,kn(1)) ! same as Conjg(gr(sf,f,kn(1)))*gl(sf,f,kn(1))*gr(kz,f,f)
      CFSR3R = Conjg(gr(sf,f,kn(1)))*gr(kz,f,f)*
     -     gr(sf,f,kn(1))


c... helicity amplitudes (entries as for FSR1). (factor of 2 to get the t+u channel)
      amp00 = (
     -     ((0,-0.5)*(2*CFSR2L*(-2*E1J**2*kv+CW*(2*E1J*EvJ + MB**2)*q1)+ 
     -     2*CFSR2R*(-2*E1J**2*kv + CW*(2*E1J*EvJ + MB**2)*q1) + 
     -     Mf*(-2*CFSR1L*E1J**2*kv - 2*CFSR1R*E1J**2*kv
     -     - 2*CFSR3L*E1J**2*kv - 2*CFSR3R*E1J**2*kv 
     -     - CFSR1R*kv*MB**2 - CFSR3L*kv*MB**2 + 
     -     2*(CFSR1R + CFSR3L)*CW*EvJ**2*q1 + CFSR1L*CW*MB**2*q1 - 
     -     CFSR1R*CW*MB**2*q1 - CFSR3L*CW*MB**2*q1+CFSR3R*CW*MB**2*q1 + 
     -     2*E1J*EvJ*(-((CFSR1R + CFSR3L)*kv) + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*CW*q1))))/MB)/delta

      amp01 = (
     -     ((Mf*(-2*CFSR2L + 2*CFSR2R 
     -     + (-CFSR1L + CFSR1R - CFSR3L + CFSR3R)*Mf)*
     -     MB**2 - 2*E1J*EvJ**2*(CFSR1R*(-E1J + q1)+CFSR3L*(E1J + q1)) - 
     -     EvJ*((CFSR1R + CFSR3L)*q1*(MB**2 - 2*CW*kv*q1) + 
     -     E1J*(4*CFSR2L*Mf - 4*CFSR2R*Mf + 2*CFSR1L*Mf**2 - 
     -     2*CFSR1R*Mf**2 + 2*CFSR3L*Mf**2 - 2*CFSR3R*Mf**2 - 
     -     CFSR1R*MB**2 + CFSR3L*MB**2 
     -     + 2*(CFSR1R - CFSR3L)*CW*kv*q1)))*SW)/(2.*Sqrt(2.)*MB))/delta

      amp02 = (
     -     ((0,0.5)*(2*(-CFSR1R + CFSR3L)*CW*E1J*EvJ**2*Mf + 
     -     E1J*(2*CFSR2L - 2*CFSR2R
     -     + (CFSR1L - CFSR1R + CFSR3L - CFSR3R)*Mf)*
     -     (CW*MB**2 - 2*kv*q1) + CW*EvJ*
     -     (2*E1J**2*(2*CFSR2L - 2*CFSR2R + 
     -     (CFSR1L - CFSR1R + CFSR3L - CFSR3R)*Mf) - 
     -     (CFSR1R - CFSR3L)*Mf*(MB**2 - 2*CW*kv*q1))))/MB)/delta

      amp03 = (
     -     ((Mf*(-2*CFSR2L + 2*CFSR2R 
     -     + (-CFSR1L + CFSR1R - CFSR3L + CFSR3R)*Mf)*
     -     MB**2 + 2*E1J*EvJ**2*(CFSR3L*(-E1J + q1)+CFSR1R*(E1J + q1)) + 
     -     EvJ*((CFSR1R + CFSR3L)*q1*(MB**2 - 2*CW*kv*q1) + 
     -     E1J*(-4*CFSR2L*Mf + 4*CFSR2R*Mf - 2*CFSR1L*Mf**2 + 
     -     2*CFSR1R*Mf**2 - 2*CFSR3L*Mf**2 + 2*CFSR3R*Mf**2 + 
     -     CFSR1R*MB**2 - CFSR3L*MB**2 
     -     +2*(-CFSR1R + CFSR3L)*CW*kv*q1)))*SW)/(2.*Sqrt(2.)*MB))/delta

      amp10 = (
     -     ((2*CFSR2L*(-2*E1J - EvJ + kv)-2*CFSR2R*(2*E1J + EvJ + kv) - 
     -     (2*CFSR1L*E1J + 2*CFSR1R*E1J + 2*CFSR3L*E1J + 2*CFSR3R*E1J + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*EvJ - CFSR1L*kv + 
     -     CFSR1R*kv - CFSR3L*kv + CFSR3R*kv)*Mf)*q1*SW)
     -     /(2.*Sqrt(2.)))/delta

      amp11 = (
     -     (0,0.25)*(1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf) + 
     -     Mf*(CFSR3R*(2*E1J + EvJ + kv)*Mf - 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR1R*((2*E1J + EvJ + kv)*Mf**2 - 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) - 
     -     CFSR3L*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta

      amp12 = (
     -     ((2*E1J**2*(2*CFSR2L - 2*CFSR2R + 
     -     (CFSR1L - CFSR1R + CFSR3L - CFSR3R)*Mf) + 
     -     E1J*EvJ*(2*CFSR2L - 2*CFSR2R + 
     -     (CFSR1L - 3*CFSR1R + 3*CFSR3L - CFSR3R)*Mf) - 
     -     E1J*kv*(2*(CFSR2L + CFSR2R) + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*Mf) - 
     -     (CFSR1R - CFSR3L)*Mf*(MB**2 - 2*CW*kv*q1))*SW)
     -     /(2.*Sqrt(2.)))/delta

      amp13 = (
     -     (0,0.25)*(-1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf) + 
     -     Mf*(CFSR3R*(2*E1J + EvJ + kv)*Mf - 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR3L*((-2*E1J - EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) + 
     -     CFSR1R*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta

      amp20 = (
     -     -((2*CFSR2R*(2*E1J + EvJ - kv)+2*CFSR2L*(2*E1J + EvJ + kv) + 
     -     (2*CFSR1L*E1J + 2*CFSR1R*E1J + 2*CFSR3L*E1J + 2*CFSR3R*E1J + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*EvJ + CFSR1L*kv - 
     -     CFSR1R*kv + CFSR3L*kv - CFSR3R*kv)*Mf)*q1*SW)
     -     /(2.*Sqrt(2.)))/delta

      amp21 = (
     -     (0,-0.25)*(-1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf)+ 
     -     Mf*(CFSR3R*(-2*E1J - EvJ + kv)*Mf + 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR1R*((-2*E1J - EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) + 
     -     CFSR3L*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta
       
      amp22 = (
     -     ((E1J*(2*CFSR2R*(-2*E1J - EvJ + kv) 
     -     + 2*CFSR2L*(2*E1J + EvJ + kv) + 
     -     (CFSR3R*(-2*E1J - EvJ + kv)+CFSR1L*(2*E1J + EvJ + kv))*Mf) - 
     -     CFSR1R*Mf*(E1J*(2*E1J + 3*EvJ - kv)+MB**2 - 2*CW*kv*q1) + 
     -     CFSR3L*Mf*(E1J*(2*E1J + 3*EvJ + kv)+MB**2 - 2*CW*kv*q1))*SW)/
     -     (2.*Sqrt(2.)))/delta
       
      amp23 = (
     -     (0,-0.25)*(1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf) + 
     -     Mf*(CFSR3R*(-2*E1J - EvJ + kv)*Mf + 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR3L*((2*E1J + EvJ + kv)*Mf**2 - 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) - 
     -     CFSR1R*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta

      if (leg.eq.1) then
c... FSR1
         amp(0, 0) = amp(0, 0) + 2*amp00
         amp(0, 1) = amp(0, 1) + 2*amp01
         amp(0, 2) = amp(0, 2) + 2*amp02
         amp(0, 3) = amp(0, 3) + 2*amp03
         amp(1, 0) = amp(1, 0) + 2*amp10
         amp(1, 1) = amp(1, 1) + 2*amp11
         amp(1, 2) = amp(1, 2) + 2*amp12
         amp(1, 3) = amp(1, 3) + 2*amp13
         amp(-1, 0) = amp(-1, 0) + 2*amp20
         amp(-1, 1) = amp(-1, 1) + 2*amp21
         amp(-1, 2) = amp(-1, 2) + 2*amp22
         amp(-1, 3) = amp(-1, 3) + 2*amp23
      else if (leg.eq.2) then
c... FSR2
         amp(0, 0) = amp(0, 0) - 2*amp00
         amp(0, 1) = amp(0, 1) - 2*amp03
         amp(0, 2) = amp(0, 2) - 2*amp02
         amp(0, 3) = amp(0, 3) - 2*amp01
         amp(1, 0) = amp(1, 0) - 2*amp20
         amp(1, 1) = amp(1, 1) - 2*amp23
         amp(1, 2) = amp(1, 2) - 2*amp22
         amp(1, 3) = amp(1, 3) - 2*amp21
         amp(-1, 0) = amp(-1, 0) - 2*amp10
         amp(-1, 1) = amp(-1, 1) - 2*amp13
         amp(-1, 2) = amp(-1, 2) - 2*amp12
         amp(-1, 3) = amp(-1, 3) - 2*amp11
      endif

      return
      end

c==================== TEMPORARY ===================================
******************************************************************
***   TEMPORARY SUBROUTINE: for IB check
***   1) w/o width propagator
***   2) isf2 -> isf1
***   3) kz -> kgamma
********************************************************************
      subroutine dsib2ffZtFSRIB(f, leg, isf1)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf1, leg
      
      integer sf !(2)
      real*8 Msf ! (2)
      complex*16 delta
      complex*16 CFSR1L, CFSR1R, CFSR2L, CFSR2R, CFSR3L, CFSR3R
      complex*16 amp10, amp11, amp12, amp13, amp20, amp21, amp22, amp23

      real*8 q1
      integer dsib2getsfermion
      
      sf  = dsib2getsfermion(f, isf1)
      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics
 
c TB: obsolete code after introducing dsib2getsfermion    
c      if(f.EQ.knue) sf(1) = ksnue
c      if(f.EQ.knumu) sf(1) = ksnumu
c      if(f.EQ.knutau) sf(1) = ksnutau
c      if(f.EQ.ke) sf(1) = kse(1)
c      if(f.EQ.kmu) sf(1) = ksmu(1)
c      if(f.EQ.ktau) sf(1) = kstau(1)
c      if(f.EQ.ku) sf(1) = ksu(1)
c      if(f.EQ.kd) sf(1) = ksd(1)
c      if(f.EQ.kc) sf(1) = ksc(1)
c      if(f.EQ.ks) sf(1) = kss(1)
c      if(f.EQ.kb) sf(1) = ksb(1)
c      if(f.EQ.kt) sf(1) = kst(1)
c      if(f.eq.knue.or.f.eq.knumu.or.f.eq.knutau) then
c         sf(2) = sf(1)
c      else
c         sf(2)=sf(1)+1
c      endif
c... masses of exchanged particle are divided by mx in order to be consistent with the kinematics
c      Msf(1) = mass(sf(1))/mx
c      Msf(2) = mass(sf(2))/mx

c... Defines FSR1 and FSR2 amplitudes through their simmetry relations via q1 exchange
c      if (leg.eq.1) then
         q1 = kJ
c      else 
       if (leg.eq.2) q1 = -kJ

c... denominators
      delta = (
     -     (2*E1J*EvJ + MB**2 - 2*CW*kv*q1)*
     -     (-E1J**2 + Mf**2 + MB**2/4. - CW*kv*q1 - Msf**2))


c... couplings 
      CFSR1L = Conjg(gl(sf,f,kn(1)))*gl(kgamma,f,f)*
     -     gl(sf,f,kn(1))
      CFSR2L = Conjg(gl(sf,f,kn(1)))*gl(kgamma,f,f)*
     -     gr(sf,f,kn(1)) ! same as Conjg(gr(sf,f,kn(1)))*gl(kz,f,f)*gl(sf,f,kn(1))
      CFSR3L = Conjg(gr(sf,f,kn(1)))*gl(kgamma,f,f)*
     -     gr(sf,f,kn(1))
      CFSR1R = Conjg(gl(sf,f,kn(1)))*gl(sf,f,kn(1))*
     -     gr(kgamma,f,f)
      CFSR2R = Conjg(gl(sf,f,kn(1)))*gr(kgamma,f,f)*
     -     gr(sf,f,kn(1)) ! same as Conjg(gr(sf,f,kn(1)))*gl(sf,f,kn(1))*gr(kz,f,f)
      CFSR3R = Conjg(gr(sf,f,kn(1)))*gr(kgamma,f,f)*
     -     gr(sf,f,kn(1))

c$$$      CFSR1L = Conjg(gl(sf,f,kn(1)))*gl(kgamma,f,f)*
c$$$     -     gl(sf,f,kn(1))
c$$$      CFSR2L =0.d0
c$$$      CFSR3L = 0.d0
c$$$      CFSR1R = Conjg(gl(sf,f,kn(1)))*gl(sf,f,kn(1))*
c$$$     -     gr(kgamma,f,f)
c$$$      CFSR2R =0.d0
c$$$      CFSR3R =0.d0

c... helicity amplitudes (entries as for FSR1). (factor of 2 to get the t+u channel)
      amp10 = (
     -     ((2*CFSR2L*(-2*E1J - EvJ + kv)-2*CFSR2R*(2*E1J + EvJ + kv) - 
     -     (2*CFSR1L*E1J + 2*CFSR1R*E1J + 2*CFSR3L*E1J + 2*CFSR3R*E1J + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*EvJ - CFSR1L*kv + 
     -     CFSR1R*kv - CFSR3L*kv + CFSR3R*kv)*Mf)*q1*SW)
     -     /(2.*Sqrt(2.)))/delta

      amp11 = (
     -     (0,0.25)*(1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf) + 
     -     Mf*(CFSR3R*(2*E1J + EvJ + kv)*Mf - 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR1R*((2*E1J + EvJ + kv)*Mf**2 - 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) - 
     -     CFSR3L*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta

      amp12 = (
     -     ((2*E1J**2*(2*CFSR2L - 2*CFSR2R + 
     -     (CFSR1L - CFSR1R + CFSR3L - CFSR3R)*Mf) + 
     -     E1J*EvJ*(2*CFSR2L - 2*CFSR2R + 
     -     (CFSR1L - 3*CFSR1R + 3*CFSR3L - CFSR3R)*Mf) - 
     -     E1J*kv*(2*(CFSR2L + CFSR2R) + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*Mf) - 
     -     (CFSR1R - CFSR3L)*Mf*(MB**2 - 2*CW*kv*q1))*SW)
     -     /(2.*Sqrt(2.)))/delta

      amp13 = (
     -     (0,0.25)*(-1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf) + 
     -     Mf*(CFSR3R*(2*E1J + EvJ + kv)*Mf - 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR3L*((-2*E1J - EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) + 
     -     CFSR1R*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta

      amp20 = (
     -     -((2*CFSR2R*(2*E1J + EvJ - kv)+2*CFSR2L*(2*E1J + EvJ + kv) + 
     -     (2*CFSR1L*E1J + 2*CFSR1R*E1J + 2*CFSR3L*E1J + 2*CFSR3R*E1J + 
     -     (CFSR1L + CFSR1R + CFSR3L + CFSR3R)*EvJ + CFSR1L*kv - 
     -     CFSR1R*kv + CFSR3L*kv - CFSR3R*kv)*Mf)*q1*SW)
     -     /(2.*Sqrt(2.)))/delta

      amp21 = (
     -     (0,-0.25)*(-1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf)+ 
     -     Mf*(CFSR3R*(-2*E1J - EvJ + kv)*Mf + 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR1R*((-2*E1J - EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) + 
     -     CFSR3L*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta
       
      amp22 = (
     -     ((E1J*(2*CFSR2R*(-2*E1J - EvJ + kv) 
     -     + 2*CFSR2L*(2*E1J + EvJ + kv) + 
     -     (CFSR3R*(-2*E1J - EvJ + kv)+CFSR1L*(2*E1J + EvJ + kv))*Mf) - 
     -     CFSR1R*Mf*(E1J*(2*E1J + 3*EvJ - kv)+MB**2 - 2*CW*kv*q1) + 
     -     CFSR3L*Mf*(E1J*(2*E1J + 3*EvJ + kv)+MB**2 - 2*CW*kv*q1))*SW)/
     -     (2.*Sqrt(2.)))/delta
       
      amp23 = (
     -     (0,-0.25)*(1 + CW)*(kv*Mf*(2*(CFSR2L + CFSR2R) + CFSR1L*Mf) + 
     -     Mf*(CFSR3R*(-2*E1J - EvJ + kv)*Mf + 
     -     (2*E1J + EvJ)*(2*CFSR2L - 2*CFSR2R + CFSR1L*Mf)) + 
     -     CFSR3L*((2*E1J + EvJ + kv)*Mf**2 - 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) - 
     -     CFSR1R*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))))/delta

      if (leg.eq.1) then
c... FSR1
         amp(1, 0) = amp(1, 0) + 2*amp10
         amp(1, 1) = amp(1, 1) + 2*amp11
         amp(1, 2) = amp(1, 2) + 2*amp12
         amp(1, 3) = amp(1, 3) + 2*amp13
         amp(-1, 0) = amp(-1, 0) + 2*amp20
         amp(-1, 1) = amp(-1, 1) + 2*amp21
         amp(-1, 2) = amp(-1, 2) + 2*amp22
         amp(-1, 3) = amp(-1, 3) + 2*amp23
      else if (leg.eq.2) then
c... FSR2
         amp(1, 0) = amp(1, 0) - 2*amp20
         amp(1, 1) = amp(1, 1) - 2*amp23
         amp(1, 2) = amp(1, 2) - 2*amp22
         amp(1, 3) = amp(1, 3) - 2*amp21
         amp(-1, 0) = amp(-1, 0) - 2*amp10
         amp(-1, 1) = amp(-1, 1) - 2*amp13
         amp(-1, 2) = amp(-1, 2) - 2*amp12
         amp(-1, 3) = amp(-1, 3) - 2*amp11
      endif

      return
      end
