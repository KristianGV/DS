*******************************************************************
*** subroutine dsib2ffZtuISR add to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel ISR diagrams. 
*** Z emitted by initial legs 1 (ISR1) or 2 (ISR2) the neutralino 
*** radiates a Z by becoming a "second" neutralino - ineu -
*** isf1 is the index if the exchanged sfermion in the t-channel
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffZtuISR(f, leg, isf1, ineu)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg, isf1, ineu
      
      integer sf !(2)
      real*8 Mkn, Msf ! (4), (2)
      complex*16 delta

      real*8 q1, q3
      complex*16 CISR1, CISR2, CISR3, CISR4, CISR5, CISR6, CISR7, CISR8,
     &     C1, C2, C3, C4, C5, C6, C7, C8
      complex*16 amp00, amp01, amp02, amp03, amp10, amp11, amp12, amp13,
     &     amp20, amp21, amp22, amp23
      integer dsib2getsfermion
       
      sf  = dsib2getsfermion(f, isf1)
      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics
      Mkn = mass(kn(ineu))/mx
 
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
c
c      do j = 1,4
c         Mkn(j) = mass(kn(j))/mx
c      end do

c... Defines ISR1 and ISR2 amplitudes through their simmetry relations via q3 exchange and conjugation of couplings. Further internal simmetry in ISRx via q1 & q3 exchange.
c      if (leg.eq.1) then
         q1 = kJ
         q3 = kv
c      else 
      if (leg.eq.2) then
         q1 = kJ
         q3 = -kv
      endif

c... denominators
      delta = (
     -     dcmplx(E1J*(E1J - EvJ) + MB**2/4. - Mkn**2,
     -     Mkn*width(kn(ineu))/mx*0d0)*
     -     dcmplx(-E1J**2 + Mf**2 + MB**2/4. - 
     -     CW*q1*q3 - Msf**2,Msf*width(sf)/mx*0d0))

c... couplings
      CISR1 = Conjg(gl(sf,f,kn(1)))*gl(kz,kn(ineu),kn(1))*
     -     gl(sf,f,kn(ineu))
      CISR2 = Conjg(gl(sf,f,kn(1)))*gl(kz,kn(ineu),kn(1))*
     -     gr(sf,f,kn(ineu))
      CISR3 = Conjg(gr(sf,f,kn(1)))*gl(kz,kn(ineu),kn(1))*
     -     gl(sf,f,kn(ineu))
      CISR4 = Conjg(gr(sf,f,kn(1)))*gl(kz,kn(ineu),kn(1))*
     -     gr(sf,f,kn(ineu))
      CISR5 = Conjg(gl(kz,kn(ineu),kn(1)))*Conjg(gl(sf,f,kn(1)))*
     -     gl(sf,f,kn(ineu))
      CISR6 = Conjg(gl(kz,kn(ineu),kn(1)))*Conjg(gl(sf,f,kn(1)))*
     -     gr(sf,f,kn(ineu))
      CISR7 = Conjg(gl(kz,kn(ineu),kn(1)))*Conjg(gr(sf,f,kn(1)))*
     -     gl(sf,f,kn(ineu))
      CISR8 = Conjg(gl(kz,kn(ineu),kn(1)))*Conjg(gr(sf,f,kn(1)))*
     -     gr(sf,f,kn(ineu))

c      if (leg.eq.1) then
         C1 = CISR1
         C2 = CISR2
         C3 = CISR3
         C4 = CISR4
         C5 = CISR5
         C6 = CISR6
         C7 = CISR7
         C8 = CISR8
c      else 
      if (leg.eq.2) then
         C1 = Conjg( CISR1 )
         C2 = Conjg( CISR2 )
         C3 = Conjg( CISR3 )
         C4 = Conjg( CISR4 )
         C5 = Conjg( CISR5 )
         C6 = Conjg( CISR6 )
         C7 = Conjg( CISR7 )
         C8 = Conjg( CISR8 )    
      endif 

c... helicity amplitudes (entries as for ISR1)
      amp00 = (
     -     ((0,0.25)*(-4*C6*CW*E1J*EvJ*q1 + 2*C6*CW*MB**2*q1 - 
     -     4*C6*E1J**2*q3 + 4*C1*E1J**2*Mf*q3 - 4*C8*E1J**2*Mf*q3 + 
     -     C1*Mf*MB**2*q3 - C8*Mf*MB**2*q3 + 
     -     C3*(-2*CW*(-2*E1J*EvJ + MB**2)*q1 + 4*E1J**2*q3) + 
     -     (4*(C4 - C5)*Mf*q3 - 
     -     2*C2*(CW*(2*E1J*EvJ + MB**2)*q1 - 2*E1J**2*q3) + 
     -     2*C7*(CW*(2*E1J*EvJ + MB**2)*q1 - 2*E1J**2*q3)
     -     )*Mkn))/MB)/delta

      amp01 = (
     -     (SW*(-4*C1*E1J**3*EvJ - 4*C8*E1J**3*EvJ - 4*C6*E1J*EvJ*Mf + 
     -     C1*E1J*EvJ*MB**2 + C8*E1J*EvJ*MB**2 + 2*C6*Mf*MB**2 + 
     -     2*C3*Mf*(-2*E1J*EvJ + MB**2) + 4*C1*E1J**2*EvJ*q1 - 
     -     4*C8*E1J**2*EvJ*q1 - C1*EvJ*MB**2*q1 + C8*EvJ*MB**2*q1 - 
     -     2*(C2*Mf*(2*E1J*EvJ + MB**2) + C7*Mf*(2*E1J*EvJ + MB**2) + 
     -     2*EvJ*(C5*(E1J - q1) + C4*(E1J + q1)))*Mkn))/
     -     (4.*Sqrt(2.)*MB))/delta

      amp02 = (
     -     ((0,-0.25)*(-4*C6*CW*E1J**2*EvJ - 4*C1*CW*E1J**2*EvJ*Mf - 
     -     4*C8*CW*E1J**2*EvJ*Mf + 2*C6*CW*E1J*MB**2 + 
     -     C1*CW*EvJ*Mf*MB**2 + C8*CW*EvJ*Mf*MB**2 - 4*C6*E1J*q1*q3 + 
     -     2*C3*E1J*(CW*(-2*E1J*EvJ + MB**2) - 2*q1*q3) - 
     -     2*(2*(C4 + C5)*CW*EvJ*Mf + 
     -     C7*E1J*(2*CW*E1J*EvJ + CW*MB**2 - 2*q1*q3) + 
     -     C2*E1J*(CW*(2*E1J*EvJ + MB**2) - 2*q1*q3)
     -     )*Mkn))/MB)/delta
      
      amp03 = (
     -     (SW*(-4*C1*E1J**3*EvJ - 4*C8*E1J**3*EvJ - 4*C6*E1J*EvJ*Mf + 
     -     C1*E1J*EvJ*MB**2 + C8*E1J*EvJ*MB**2 + 2*C6*Mf*MB**2 + 
     -     2*C3*Mf*(-2*E1J*EvJ + MB**2) - 4*C1*E1J**2*EvJ*q1 + 
     -     4*C8*E1J**2*EvJ*q1 + C1*EvJ*MB**2*q1 - C8*EvJ*MB**2*q1 - 
     -     2*(C2*Mf*(2*E1J*EvJ + MB**2) + C7*Mf*(2*E1J*EvJ + MB**2) + 
     -     2*EvJ*(C4*(E1J - q1) + C5*(E1J + q1)))*Mkn))/
     -     (4.*Sqrt(2.)*MB))/delta

      amp10 = (
     -     -(q1*SW*(C6*(2*E1J - EvJ + q3) + C3*(-2*E1J + EvJ + q3) - 
     -     (C2*(-2*E1J - EvJ + q3) + C7*(2*E1J + EvJ + q3))*Mkn))/
     -     (2.*Sqrt(2.)))/delta

      amp11 = (
     -     (0,-0.125)*(1 + CW)*(2*C6*Mf*(2*E1J - EvJ + q3) - 
     -     2*C3*Mf*(-2*E1J + EvJ + q3) + 
     -     C8*(4*Mf**2*q3 + (E1J + q1)*(-MB**2 + 4*(E1J**2 + q1*q3))) + 
     -     C1*(MB**2*(-E1J + q1) - 
     -     4*(Mf**2*q3 + (-E1J + q1)*(E1J**2 + q1*q3))) - 
     -     4*C5*(-E1J + q1)*Mkn + 4*C4*(E1J + q1)*Mkn - 
     -     2*C2*Mf*(-2*E1J - EvJ + q3)*Mkn + 
     -     2*C7*Mf*(2*E1J + EvJ + q3)*Mkn))/delta

      amp12 = (
     -     (SW*(2*C3*E1J*(2*E1J - EvJ - q3)+2*C6*E1J*(2*E1J - EvJ + q3)- 
     -     C1*Mf*(MB**2 + 4*E1J*(-E1J + q3)) - 
     -     C8*Mf*(MB**2 - 4*E1J*(E1J + q3)) + 4*C4*Mf*Mkn + 
     -     4*C5*Mf*Mkn + 2*C2*E1J*(2*E1J + EvJ - q3)*Mkn + 
     -     2*C7*E1J*(2*E1J + EvJ + q3)*Mkn))/(4.*Sqrt(2.)))/delta

      amp13 = (
     -     (0,-0.125)*(-1 + CW)*(2*C6*Mf*(2*E1J - EvJ + q3) - 
     -     2*C3*Mf*(-2*E1J + EvJ + q3) + 
     -     C8*(4*Mf**2*q3 + (-E1J + q1)*(-4*E1J**2 + MB**2 + 4*q1*q3)) - 
     -     C1*(MB**2*(E1J + q1) + 4*(Mf**2*q3
     -     + (E1J + q1)*(-E1J**2 + q1*q3))) - 
     -     4*C4*(-E1J + q1)*Mkn + 4*C5*(E1J + q1)*Mkn - 
     -     2*C2*Mf*(-2*E1J - EvJ + q3)*Mkn + 
     -     2*C7*Mf*(2*E1J + EvJ + q3)*Mkn))/delta

      amp20 = (
     -     (q1*SW*(C3*(2*E1J - EvJ + q3) + C6*(-2*E1J + EvJ + q3) - 
     -     (C7*(-2*E1J - EvJ + q3) + C2*(2*E1J + EvJ + q3))*Mkn))/
     -     (2.*Sqrt(2.)))/delta

      amp21 = (
     -     (0,-0.125)*(-1 + CW)*(2*C3*Mf*(2*E1J - EvJ + q3) - 
     -     2*C6*Mf*(-2*E1J + EvJ + q3) + 
     -     C8*(-4*Mf**2*q3 - (E1J + q1)*(-4*E1J**2 + MB**2 + 4*q1*q3)) + 
     -     C1*(MB**2*(-E1J + q1) + 
     -     4*(Mf**2*q3 + (-E1J + q1)*(-E1J**2 + q1*q3))) - 
     -     4*C5*(-E1J + q1)*Mkn + 4*C4*(E1J + q1)*Mkn - 
     -     2*C7*Mf*(-2*E1J - EvJ + q3)*Mkn + 
     -     2*C2*Mf*(2*E1J + EvJ + q3)*Mkn))/delta
      
      amp22 = (
     -     (SW*(2*C3*E1J*(2*E1J - EvJ + q3)
     -     - 2*C6*E1J*(-2*E1J + EvJ + q3) - 
     -     C8*Mf*(MB**2 + 4*E1J*(-E1J + q3)) - 
     -     C1*Mf*(MB**2 - 4*E1J*(E1J + q3)) + 4*C4*Mf*Mkn + 
     -     4*C5*Mf*Mkn + 2*C7*E1J*(2*E1J + EvJ - q3)*Mkn + 
     -     2*C2*E1J*(2*E1J + EvJ + q3)*Mkn))/(4.*Sqrt(2.)))/delta
      
      amp23 = (
     -     (0,0.125)*(1 + CW)*(-2*C3*Mf*(2*E1J - EvJ + q3) + 
     -     2*C6*Mf*(-2*E1J + EvJ + q3) + 
     -     C8*(MB**2*(E1J - q1) + 4*(Mf**2*q3
     -     + (-E1J + q1)*(E1J**2 + q1*q3))) + 
     -     C1*(MB**2*(E1J + q1) - 4*(Mf**2*q3 
     -     + (E1J + q1)*(E1J**2 + q1*q3))) + 
     -     4*C4*(-E1J + q1)*Mkn - 4*C5*(E1J + q1)*Mkn + 
     -     2*C7*Mf*(-2*E1J - EvJ + q3)*Mkn - 
     -     2*C2*Mf*(2*E1J + EvJ + q3)*Mkn))/delta

      if (leg.eq.1) then
c... ISR1 (factor of 2 to get the t+u channel)
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
c... ISR2 (factor of 2 to get the t+u channel)
         amp(0, 0) = amp(0, 0) + 2*amp00
         amp(0, 1) = amp(0, 1) - 2*amp01
         amp(0, 2) = amp(0, 2) - 2*amp02
         amp(0, 3) = amp(0, 3) - 2*amp03
         amp(1, 0) = amp(1, 0) + 2*amp10
         amp(1, 1) = amp(1, 1) - 2*amp11
         amp(1, 2) = amp(1, 2) - 2*amp12
         amp(1, 3) = amp(1, 3) - 2*amp13
         amp(-1, 0) = amp(-1, 0) + 2*amp20
         amp(-1, 1) = amp(-1, 1) - 2*amp21
         amp(-1, 2) = amp(-1, 2) - 2*amp22
         amp(-1, 3) = amp(-1, 3) - 2*amp23
      endif

      return
      end
