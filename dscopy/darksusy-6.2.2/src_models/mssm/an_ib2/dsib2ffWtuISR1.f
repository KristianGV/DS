*******************************************************************
*** subroutine dsib2ffWtuISRX adds to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel ISR diagrams. 
*** W emitted by initial legs 1 (ISR1) or 2 (ISR2). isf is the index of the exchanged sfermion.
***
*** Subroutines: dsib2ffWtuISR1:  W emitted from p1 leg
***              dsib2ffWtuISR2:  W emitted from p2 leg
***
*** IMPORTANT: do *NOT* use separately dsib2ffWtuISR1 and dsib2ffWtuISR2.
***            The correct t+u channel for the 4 involved diagrams is given by ISR1 + ISR2. 
***            Separately, the result is NOT the same as ISR1_t + ISR1_u (ISR2_t + ISR2_u), 
***            being ISR1(t) \equiv ISR2(u) and ISR2(t) \equiv ISR1(u)
***
***            ISR1 and ISR2 are not invariant under CP symmetry serately, but
***            their sum is!
***
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
***         Francesca Calore, 2015-01-30
*******************************************************************
      subroutine dsib2ffWtuISR1(f, isf, icha)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf, icha

      integer ff, sff
      real*8 Msff, Mkcha
      integer dsib2getsfermion

      complex*16 delta
      complex*16 CISR1, CISR2, CISR3, CISR4, CISR5, CISR6, CISR7, CISR8
      complex*16 CCISR1, CCISR2, CCISR3, CCISR4, CCISR5,
     &     CCISR6, CCISR7, CCISR8

c      sf  = dsib2getsfermion(f, isf)
c      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics

      ff   = c_fbartype
      sff  = dsib2getsfermion(ff, isf)
      Msff = mass(sff)/mx ! dimensionless in order to be consistent with kinematics

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
c
c      if (f.EQ.kd.or.f.EQ.ks.or.f.EQ.kb) then
c         sf(2) = sf(1)+1
c         sff(1) = sf(1)-2
c         sff(2)= sff(1)+1
c      elseif(f.EQ.ku.or.f.EQ.kc.or.f.EQ.kt) then
c         sf(2) = sf(1)+1
c         sff(1) = sf(1)+2
c         sff(2)= sff(1)+1
c      elseif(f.EQ.ke.or.f.EQ.kmu.or.f.EQ.ktau) then
c         sf(2) = sf(1)+1
c         sff(1) = sf(1)-1
c         sff(2)= sff(1)
c      elseif(f.EQ.knue.or.f.EQ.knumu.or.f.EQ.knutau) then 
c         sf(2) = sf(1)
c         sff(1) = sf(1)+1
c         sff(2) = sff(1)+1
c      endif
c      
c      Msf(1) = mass(sf(1))/mx
c      Msf(2) = mass(sf(2))/mx
c      Msff(1) = mass(sff(1))/mx
c      Msff(2) = mass(sff(2))/mx 
 
c      do j = 1,2
         Mkcha = mass(kcha(icha))/mx
c      end do 

c... couplings
      if(mod(c_ftype,2).EQ.0) then
         CCISR1 = Conjg(gl(sff,ff,kn(1)))*
     -    gl(sff,f,kcha(icha))*gr(kw,kn(1),kcha(icha))
         CCISR2 = Conjg(gl(sff,ff,kn(1)))*gr(kw,kn(1),kcha(icha))*
     -    gr(sff,f,kcha(icha))
         CCISR3 =Conjg(gr(sff,ff,kn(1)))*gl(sff,f,kcha(icha))*
     -    gr(kw,kn(1),kcha(icha))
         CCISR4 =Conjg(gr(sff,ff,kn(1)))*gr(kw,kn(1),kcha(icha))*
     -    gr(sff,f,kcha(icha))
         CCISR5 = Conjg(gl(sff,ff,kn(1)))*gl(kw,kn(1),kcha(icha))*
     -    gl(sff,f,kcha(icha))
         CCISR6 = Conjg(gl(sff,ff,kn(1)))*gl(kw,kn(1),kcha(icha))*
     -    gr(sff,f,kcha(icha))
         CCISR7 = Conjg(gr(sff,ff,kn(1)))*gl(kw,kn(1),kcha(icha))*
     -    gl(sff,f,kcha(icha))
         CCISR8 =Conjg(gr(sff,ff,kn(1)))*gl(kw,kn(1),kcha(icha))*
     -    gr(sff,f,kcha(icha))

         CISR1  = - CCISR1
         CISR2  = - CCISR2 
         CISR3  = - CCISR3
         CISR4  = - CCISR4 
         CISR5  = - CCISR5 
         CISR6  = - CCISR6 
         CISR7  = - CCISR7 
         CISR8  = - CCISR8 

      else
         CISR1 = Conjg(gl(kw,kn(1),kcha(icha)))*
     -        Conjg(gl(sff,ff,kn(1)))*gl(sff,f,kcha(icha))
         CISR2 = Conjg(gl(kw,kn(1),kcha(icha)))*
     -        Conjg(gl(sff,ff,kn(1)))*gr(sff,f,kcha(icha))
         CISR3 = Conjg(gl(kw,kn(1),kcha(icha)))*
     -        Conjg(gr(sff,ff,kn(1)))*gl(sff,f,kcha(icha))
         CISR4 = Conjg(gl(kw,kn(1),kcha(icha)))*
     -        Conjg(gr(sff,ff,kn(1)))*gr(sff,f,kcha(icha))
         CISR5 = Conjg(gl(sff,ff,kn(1)))*
     -        Conjg(gr(kw,kn(1),kcha(icha)))*gl(sff,f,kcha(icha))
         CISR6 = Conjg(gl(sff,ff,kn(1)))*
     -        Conjg(gr(kw,kn(1),kcha(icha)))*gr(sff,f,kcha(icha))
         CISR7 = Conjg(gr(kw,kn(1),kcha(icha)))*
     -        Conjg(gr(sff,ff,kn(1)))*gl(sff,f,kcha(icha))
         CISR8 = Conjg(gr(kw,kn(1),kcha(icha)))*
     -        Conjg(gr(sff,ff,kn(1)))*gr(sff,f,kcha(icha))
      endif
      
c... denominators
      delta = (
     -     dcmplx((2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ + 2*kJ**2 + 
     -      MB**2 + Mf**2 + Mff**2 - 4*Mkcha**2)/4.,
     -     Mkcha*width(kcha(icha))/mx*0d0)*
     -     dcmplx((-2*E1J*E2J + 2*E1J*EvJ - 2*E2J*EvJ - 2*kJ**2 - 
     -      4*CW*kJ*kv + MB**2 + Mf**2 + Mff**2 - 4*Msff**2)/4.
     -     ,Msff*width(sff)/mx*0d0))

c... Helicity amplitudes
      amp(0, 0) = amp(0, 0) + 2*(
     - ((0,0.25)*(2*CISR3*(2*E2J*Epl*kv + Emi*kv*Mf - Emi*kv*Mff + 
     - CW*EvJ*(2*Epl*kJ + Mf*pmi - Mff*pmi) - CW*MB**2*ppl) + 
     - 2*CISR6*(2*E2J*Epl*kv + Emi*kv*Mf - Emi*kv*Mff + 
     -  CW*EvJ*(2*Epl*kJ + Mf*pmi - Mff*pmi) - CW*MB**2*ppl) + 
     -  (CISR1 + CISR8)*(kv*(2*(E1J + E2J)*Epl*(Mf + Mff) + 
     -  Emi*(-2*E1J*(E1J + E2J) + MB**2 + Mf**2 - Mff**2)) + 
     -  CW*EvJ*(-2*E1J*(E1J + E2J) + MB**2 + Mf**2 - Mff**2)*pmi - 
     -    2*CW*MB**2*(-2*Emi*kJ + (-E1J + E2J)*pmi + (Mf + Mff)*ppl)) - 
     -  2*(-2*(CISR4 + CISR5)*Emi*kv + 
     -   CW*EvJ*(-2*CISR4*pmi - 2*CISR5*pmi + 
     -   (CISR2 + CISR7)*(2*Epl*kJ + Mf*pmi - Mff*pmi)) + 
     -   (CISR2 + CISR7)*(kv*(-2*E1J*Epl +
     -   Emi*(Mf - Mff)) + CW*MB**2*ppl))*
     -   Mkcha))/MB
     -   )/delta

      amp(0, 1) = amp(0, 1) + 2*(
     - (SW*(2*CISR3*MB**2*(Emi + pmi) - 
     - 2*CISR3*EvJ*(Epl*(Mf + Mff) + (Mf - Mff)*ppl) + 
     - EvJ*(2*CISR6*(Epl*(Mf + Mff) + (-Mf + Mff)*ppl) + 
     - (-2*E1J*(E1J + E2J) + MB**2 + Mf**2 - Mff**2)*
     - (CISR1*(Epl - ppl) - CISR8*(Epl + ppl))) + 
     - 2*MB**2*(CISR6*(-Emi + pmi) + 
     -  CISR1*(Emi*(-Mf + Mff) + (Mf + Mff)*pmi + 
     -  (E1J - E2J)*(Epl - ppl)) + 
     -  CISR8*(Emi*(Mf - Mff) + (Mf + Mff)*pmi -
     -  (E1J - E2J)*(Epl + ppl))) 
     -  - 2*(MB**2*(CISR2*(Emi - pmi) - CISR7*(Emi + pmi)) + 
     -  EvJ*((CISR2 - CISR7)*Epl*(Mf + Mff) - 
     -  (CISR2 + CISR7)*(Mf - Mff)*ppl + 2*CISR5*(-Epl + ppl) + 
     -  2*CISR4*(Epl + ppl)))*Mkcha))/(4.*Sqrt(2.)*MB)
     -  )/delta
     
      amp(0, 2) = amp(0, 2) + 2*(
     - ((0,0.25)*(2*CISR1*CW*E1J*Emi*kv**2 - 2*CISR1*CW*E2J*Emi*kv**2 + 
     - 2*CISR8*CW*E1J*Emi*MB**2 - 2*CISR8*CW*E2J*Emi*MB**2 + 
     - 2*CISR1*CW*Epl*MB**2*Mf - 2*CISR8*CW*Epl*MB**2*Mf - 
     - 2*CISR1*CW*Epl*MB**2*Mff + 2*CISR8*CW*Epl*MB**2*Mff - 
     - (CISR1 - CISR8)*CW*EvJ*Emi*
     -  (-2*E1J*(E1J + E2J) + MB**2 + Mf**2 - Mff**2) + 
     - 2*CISR1*E1J**2*kv*pmi - 2*CISR8*E1J**2*kv*pmi + 
     - 2*CISR1*E1J*E2J*kv*pmi - 2*CISR8*E1J*E2J*kv*pmi + 
     - 4*CISR1*CW*kJ*kv**2*pmi + 4*CISR8*CW*kJ*MB**2*pmi - 
     - CISR1*kv*MB**2*pmi + CISR8*kv*MB**2*pmi - CISR1*kv*Mf**2*pmi + 
     - CISR8*kv*Mf**2*pmi + CISR1*kv*Mff**2*pmi - CISR8*kv*Mff**2*pmi - 
     - 2*CISR1*CW*EvJ**2*((E1J - E2J)*Emi + 2*kJ*pmi) - 
     - 2*CISR1*E1J*kv*Mf*ppl + 2*CISR8*E1J*kv*Mf*ppl - 
     - 2*CISR1*E2J*kv*Mf*ppl + 2*CISR8*E2J*kv*Mf*ppl + 
     - 2*CISR1*E1J*kv*Mff*ppl - 2*CISR8*E1J*kv*Mff*ppl + 
     - 2*CISR1*E2J*kv*Mff*ppl - 2*CISR8*E2J*kv*Mff*ppl + 
     - 2*CISR3*(kv*(Mf + Mff)*pmi + 2*E2J*kv*ppl + 
     -    CW*(-(Epl*MB**2) + EvJ*(Emi*(Mf + Mff) + 2*kJ*ppl))) - 
     - 2*CISR6*(kv*(Mf + Mff)*pmi + 2*E2J*kv*ppl + 
     -    CW*(-(Epl*MB**2) + EvJ*(Emi*(Mf + Mff) + 2*kJ*ppl))) + 
     - 2*(2*(CISR4 - CISR5)*kv*pmi + 
     -  (CISR2 - CISR7)*(CW*Epl*MB**2 + kv*(Mf + Mff)*pmi - 
     -       2*E1J*kv*ppl) + CW*EvJ*
     -  (2*CISR4*Emi - 2*CISR5*Emi + 
     -  (CISR2 - CISR7)*(Emi*(Mf + Mff) + 2*kJ*ppl)))*Mkcha))/MB
     -     )/delta

      amp(0, 3) = amp(0, 3) + 2*(
     - -(SW*(2*CISR3*(MB**2*(-Emi + pmi) + 
     - EvJ*(Epl*(Mf + Mff) + (-Mf + Mff)*ppl)) - 
     - EvJ*(2*CISR6*(Epl*(Mf + Mff) + (Mf - Mff)*ppl) + 
     - (-2*E1J*(E1J + E2J) + MB**2 + Mf**2 - Mff**2)*
     - (CISR8*(-Epl + ppl) + CISR1*(Epl + ppl))) + 
     - 2*MB**2*(CISR6*(Emi + pmi) + 
     -  CISR8*(Emi*(-Mf + Mff) + (Mf + Mff)*pmi + 
     -  (E1J - E2J)*(Epl - ppl)) + 
     -  CISR1*(Emi*(Mf - Mff) + (Mf + Mff)*pmi - (E1J - E2J)*(Epl + ppl))
     -  ) + 2*(MB**2*(CISR7*(-Emi + pmi) + CISR2*(Emi + pmi)) + 
     -  EvJ*((CISR2 - CISR7)*Epl*(Mf + Mff) + 2*CISR4*(Epl - ppl) + 
     -  (CISR2 + CISR7)*(Mf - Mff)*ppl - 2*CISR5*(Epl + ppl)))*
     -  Mkcha))/(4.*Sqrt(2.)*MB)
     -     )/delta

      amp(1, 0) = amp(1, 0) + 2*(
     - (SW*(-2*CISR1*E1J**2*pmi 
     - - 2*CISR8*E1J**2*pmi - 2*CISR1*E1J*E2J*pmi - 
     - 2*CISR8*E1J*E2J*pmi + CISR1*MB**2*pmi + CISR8*MB**2*pmi + 
     - CISR1*Mf**2*pmi + CISR8*Mf**2*pmi - CISR1*Mff**2*pmi - 
     - CISR8*Mff**2*pmi + 2*CISR1*kv*Mf*ppl - 2*CISR8*kv*Mf*ppl - 
     - 2*CISR1*kv*Mff*ppl + 2*CISR8*kv*Mff*ppl + 
     - 2*CISR6*(2*Epl*kJ + Mf*pmi - Mff*pmi - EvJ*ppl + kv*ppl) + 
     - 2*CISR3*(2*Epl*kJ + Mf*pmi - Mff*pmi - (EvJ + kv)*ppl) - 
     - 2*(CISR1 + CISR8)*EvJ*(-2*Emi*kJ - E1J*pmi + E2J*pmi + 
     -    (Mf + Mff)*ppl) - 2*(2*CISR2*Epl*kJ + 2*CISR7*Epl*kJ - 
     -    2*CISR4*pmi - 2*CISR5*pmi + CISR2*Mf*pmi + CISR7*Mf*pmi - 
     -    CISR2*Mff*pmi - CISR7*Mff*pmi + (CISR2 + CISR7)*EvJ*ppl - 
     -    CISR2*kv*ppl + CISR7*kv*ppl)*Mkcha))/(4.*Sqrt(2.))
     -     )/delta

      amp(1, 1) = amp(1, 1) + 2*(
     - (0,0.125)*(1 + CW)*(-2*CISR1*E1J**2*Epl + 2*CISR8*E1J**2*Epl - 
     - 2*CISR1*E1J*E2J*Epl + 2*CISR8*E1J*E2J*Epl - 4*CISR1*Epl*kJ*kv + 
     - 4*CISR8*Epl*kJ*kv + CISR1*Epl*MB**2 - CISR8*Epl*MB**2 + 
     - 2*CISR1*Emi*kv*Mf + 2*CISR8*Emi*kv*Mf + CISR1*Epl*Mf**2 - 
     - CISR8*Epl*Mf**2 + 2*CISR1*Emi*kv*Mff + 2*CISR8*Emi*kv*Mff - 
     - CISR1*Epl*Mff**2 + CISR8*Epl*Mff**2 - 2*CISR1*kv*Mf*pmi + 
     - 2*CISR8*kv*Mf*pmi + 2*CISR1*kv*Mff*pmi - 2*CISR8*kv*Mff*pmi + 
     - 2*CISR1*E1J**2*ppl + 2*CISR8*E1J**2*ppl + 2*CISR1*E1J*E2J*ppl + 
     - 2*CISR8*E1J*E2J*ppl + 4*CISR1*kJ*kv*ppl + 4*CISR8*kJ*kv*ppl - 
     -  CISR1*MB**2*ppl - CISR8*MB**2*ppl 
     -  - CISR1*Mf**2*ppl - CISR8*Mf**2*ppl + 
     -  CISR1*Mff**2*ppl + CISR8*Mff**2*ppl + 
     - 2*CISR6*(Epl*(Mf + Mff) + (-EvJ + kv)*(Emi - pmi) 
     -  - Mf*ppl + Mff*ppl) + 
     - 2*CISR3*(-(Epl*(Mf + Mff)) + (EvJ + kv)*(Emi + pmi) - Mf*ppl + 
     -  Mff*ppl) + 2*EvJ*(CISR1*
     -  (Emi*(-Mf + Mff) + (Mf + Mff)*pmi + (E1J - E2J)*(Epl - ppl)) + 
     -  CISR8*(Emi*(Mf - Mff) + (Mf + Mff)*pmi -
     -  (E1J - E2J)*(Epl + ppl))) + 
     -  2*(-(CISR2*EvJ*Emi) + CISR7*EvJ*Emi + CISR2*Emi*kv +CISR7*Emi*kv- 
     -   CISR2*Epl*Mf + CISR7*Epl*Mf - CISR2*Epl*Mff + CISR7*Epl*Mff + 
     -   CISR2*EvJ*pmi + CISR7*EvJ*pmi - CISR2*kv*pmi + CISR7*kv*pmi + 
     -  2*CISR5*(Epl - ppl) + CISR2*Mf*ppl + CISR7*Mf*ppl -CISR2*Mff*ppl- 
     -   CISR7*Mff*ppl - 2*CISR4*(Epl + ppl))*Mkcha)
     -     )/delta

      amp(1, 2) = amp(1, 2) + 2*(
     - (SW*(CISR1*(-2*Epl*((-EvJ + kv)*Mf + (EvJ + kv)*Mff) + 
     - Emi*(2*(E1J**2 + E1J*E2J - E1J*EvJ + E2J*EvJ) - MB**2 - Mf**2 + 
     - Mff**2) - 4*EvJ*kJ*pmi) - 
     - CISR8*(Emi*(2*(E1J**2 + E1J*E2J - E1J*EvJ + E2J*EvJ) - MB**2 - 
     - Mf**2 + Mff**2) + 2*Epl*(EvJ*(Mf - Mff) + kv*(Mf + Mff)) - 
     -    4*EvJ*kJ*pmi) - 2*CISR6*
     -     (Epl*(-EvJ + kv) + Emi*(Mf + Mff) + 2*kJ*ppl) + 
     - 2*CISR3*(-(Epl*(EvJ + kv)) + Emi*(Mf + Mff) + 2*kJ*ppl) + 
     - 4*CISR4*Emi*Mkcha - 4*CISR5*Emi*Mkcha + 
     - 2*CISR2*(Epl*(EvJ - kv) + Emi*(Mf + Mff) 
     - + 2*kJ*ppl)*Mkcha - 
     - 2*CISR7*(Epl*(EvJ + kv) + Emi*(Mf + Mff) 
     -  + 2*kJ*ppl)*Mkcha))/
     -  (4.*Sqrt(2.))
     -     )/delta

      amp(1, 3) = amp(1, 3) + 2*(
     - (0,0.125)*(-1 + CW)*(-2*CISR1*E1J**2*Epl + 2*CISR8*E1J**2*Epl - 
     - 2*CISR1*E1J*E2J*Epl + 2*CISR8*E1J*E2J*Epl + 4*CISR1*Epl*kJ*kv - 
     - 4*CISR8*Epl*kJ*kv + CISR1*Epl*MB**2 - CISR8*Epl*MB**2 + 
     - 2*CISR1*Emi*kv*Mf + 2*CISR8*Emi*kv*Mf + CISR1*Epl*Mf**2 - 
     -  CISR8*Epl*Mf**2 + 2*CISR1*Emi*kv*Mff + 2*CISR8*Emi*kv*Mff - 
     -  CISR1*Epl*Mff**2 + CISR8*Epl*Mff**2 + 2*CISR1*kv*Mf*pmi - 
     - 2*CISR8*kv*Mf*pmi - 2*CISR1*kv*Mff*pmi + 2*CISR8*kv*Mff*pmi - 
     - 2*CISR1*E1J**2*ppl - 2*CISR8*E1J**2*ppl - 2*CISR1*E1J*E2J*ppl - 
     - 2*CISR8*E1J*E2J*ppl + 4*CISR1*kJ*kv*ppl + 4*CISR8*kJ*kv*ppl + 
     -  CISR1*MB**2*ppl + CISR8*MB**2*ppl
     -  + CISR1*Mf**2*ppl + CISR8*Mf**2*ppl - 
     -  CISR1*Mff**2*ppl - CISR8*Mff**2*ppl + 
     - 2*CISR6*(Epl*(Mf + Mff) + (-EvJ + kv)*
     - (Emi + pmi) + Mf*ppl - Mff*ppl) - 
     - 2*CISR3*(Epl*(Mf + Mff) - (EvJ + kv)*
     - (Emi - pmi) - Mf*ppl + Mff*ppl) - 
     - 2*EvJ*(CISR8*(Emi*(-Mf + Mff) + (Mf + Mff)*pmi + 
     -  (E1J - E2J)*(Epl - ppl)) + 
     -  CISR1*(Emi*(Mf - Mff) + (Mf + Mff)*pmi 
     - - (E1J - E2J)*(Epl + ppl))) + 
     - 2*(CISR2*Emi*kv + CISR7*Emi*kv - CISR2*Epl*Mf + CISR7*Epl*Mf - 
     -  CISR2*Epl*Mff + CISR7*Epl*Mff + CISR2*kv*pmi - CISR7*kv*pmi - 
     -  EvJ*((CISR2 - CISR7)*Emi + (CISR2 + CISR7)*pmi) - CISR2*Mf*ppl - 
     -  CISR7*Mf*ppl + CISR2*Mff*ppl + CISR7*Mff*ppl + 
     - 2*CISR4*(-Epl + ppl) + 2*CISR5*(Epl + ppl))*Mkcha)
     -     )/delta
      
      amp(-1, 0) = amp(-1, 0) + 2*(
     - (SW*(-2*CISR1*E1J**2*pmi - 2*CISR8*E1J**2*pmi
     - - 2*CISR1*E1J*E2J*pmi - 
     - 2*CISR8*E1J*E2J*pmi + CISR1*MB**2*pmi + CISR8*MB**2*pmi + 
     - CISR1*Mf**2*pmi + CISR8*Mf**2*pmi - CISR1*Mff**2*pmi - 
     - CISR8*Mff**2*pmi - 2*CISR1*kv*Mf*ppl + 2*CISR8*kv*Mf*ppl + 
     - 2*CISR1*kv*Mff*ppl - 2*CISR8*kv*Mff*ppl + 
     - 2*CISR3*(2*Epl*kJ + Mf*pmi - Mff*pmi - EvJ*ppl + kv*ppl) + 
     - 2*CISR6*(2*Epl*kJ + Mf*pmi - Mff*pmi - (EvJ + kv)*ppl) - 
     - 2*(CISR1 + CISR8)*EvJ*(-2*Emi*kJ - E1J*pmi + E2J*pmi + 
     -  (Mf + Mff)*ppl) - 2*(2*CISR2*Epl*kJ + 2*CISR7*Epl*kJ - 
     - 2*CISR4*pmi - 2*CISR5*pmi + CISR2*Mf*pmi + CISR7*Mf*pmi - 
     -  CISR2*Mff*pmi - CISR7*Mff*pmi + (CISR2 + CISR7)*EvJ*ppl + 
     -  CISR2*kv*ppl - CISR7*kv*ppl)*Mkcha))/(4.*Sqrt(2.))
     -     )/delta
      
      amp(-1, 1) = amp(-1, 1) + 2*(
     - (0,-0.125)*(-1 + CW)*(2*CISR1*E1J**2*Epl - 2*CISR8*E1J**2*Epl + 
     - 2*CISR1*E1J*E2J*Epl - 2*CISR8*E1J*E2J*Epl - 4*CISR1*Epl*kJ*kv + 
     - 4*CISR8*Epl*kJ*kv - CISR1*Epl*MB**2 + CISR8*Epl*MB**2 + 
     - 2*CISR1*Emi*kv*Mf + 2*CISR8*Emi*kv*Mf - CISR1*Epl*Mf**2 + 
     -  CISR8*Epl*Mf**2 + 2*CISR1*Emi*kv*Mff + 2*CISR8*Emi*kv*Mff + 
     -  CISR1*Epl*Mff**2 - CISR8*Epl*Mff**2 - 2*CISR1*kv*Mf*pmi + 
     - 2*CISR8*kv*Mf*pmi + 2*CISR1*kv*Mff*pmi - 2*CISR8*kv*Mff*pmi - 
     - 2*CISR1*E1J**2*ppl - 2*CISR8*E1J**2*ppl - 2*CISR1*E1J*E2J*ppl - 
     - 2*CISR8*E1J*E2J*ppl + 4*CISR1*kJ*kv*ppl + 4*CISR8*kJ*kv*ppl + 
     -  CISR1*MB**2*ppl + CISR8*MB**2*ppl 
     -  + CISR1*Mf**2*ppl + CISR8*Mf**2*ppl - 
     -  CISR1*Mff**2*ppl - CISR8*Mff**2*ppl + 
     - 2*CISR3*(Epl*(Mf + Mff) + (-EvJ + kv)*
     -  (Emi + pmi) + Mf*ppl - Mff*ppl) - 
     - 2*CISR6*(Epl*(Mf + Mff) - (EvJ + kv)*
     -  (Emi - pmi) - Mf*ppl + Mff*ppl) + 
     - 2*EvJ*(CISR1*(Emi*(Mf - Mff) - (Mf + Mff)*pmi - 
     -  (E1J - E2J)*(Epl - ppl)) - 
     -  CISR8*(Emi*(Mf - Mff) + (Mf + Mff)*pmi - 
     -  (E1J - E2J)*(Epl + ppl))) + 
     - 2*(CISR2*Emi*kv + CISR7*Emi*kv + CISR2*Epl*Mf - CISR7*Epl*Mf + 
     -  CISR2*Epl*Mff - CISR7*Epl*Mff - CISR2*kv*pmi + CISR7*kv*pmi + 
     -  EvJ*((CISR2 - CISR7)*Emi - (CISR2 + CISR7)*pmi) - CISR2*Mf*ppl - 
     -  CISR7*Mf*ppl + CISR2*Mff*ppl + CISR7*Mff*ppl + 
     -  2*CISR5*(-Epl + ppl) + 2*CISR4*(Epl + ppl))*Mkcha)
     -     )/delta

      amp(-1, 2) = amp(-1, 2) + 2*(
     - (SW*(CISR1*(Emi*(2*(E1J**2 + E1J*E2J - 
     - E1J*EvJ + E2J*EvJ) - MB**2 - 
     - Mf**2 + Mff**2) + 2*Epl*(EvJ*(Mf - Mff) + kv*(Mf + Mff)) - 
     - 4*EvJ*kJ*pmi) + CISR8*
     - (Emi*(-2*E1J*(E1J + E2J) + 2*(E1J - E2J)*EvJ + MB**2 + Mf**2 - 
     - Mff**2) + 2*(Epl*kv*(Mf + Mff) + 
     - EvJ*(Epl*(-Mf + Mff) + 2*kJ*pmi))) + 
     - 2*CISR3*(Epl*(-EvJ + kv) + Emi*(Mf + Mff) + 2*kJ*ppl) - 
     - 2*CISR6*(-(Epl*(EvJ + kv)) + Emi*(Mf + Mff) + 2*kJ*ppl) + 
     - 4*CISR4*Emi*Mkcha - 4*CISR5*Emi*Mkcha - 
     - 2*CISR7*(Epl*(EvJ - kv) + Emi*(Mf + Mff) 
     - + 2*kJ*ppl)*Mkcha + 
     - 2*CISR2*(Epl*(EvJ + kv) + Emi*(Mf + Mff) 
     - + 2*kJ*ppl)*Mkcha))/
     -  (4.*Sqrt(2.))
     -     )/delta
      
      amp(-1, 3) = amp(-1, 3) + 2*(
     - (0,-0.125)*(1 + CW)*(2*CISR1*E1J**2*Epl - 2*CISR8*E1J**2*Epl + 
     - 2*CISR1*E1J*E2J*Epl - 2*CISR8*E1J*E2J*Epl + 4*CISR1*Epl*kJ*kv - 
     - 4*CISR8*Epl*kJ*kv - CISR1*Epl*MB**2 + CISR8*Epl*MB**2 + 
     - 2*CISR1*Emi*kv*Mf + 2*CISR8*Emi*kv*Mf - CISR1*Epl*Mf**2 + 
     -  CISR8*Epl*Mf**2 + 2*CISR1*Emi*kv*Mff + 2*CISR8*Emi*kv*Mff + 
     -  CISR1*Epl*Mff**2 - CISR8*Epl*Mff**2 + 2*CISR1*kv*Mf*pmi - 
     - 2*CISR8*kv*Mf*pmi - 2*CISR1*kv*Mff*pmi + 2*CISR8*kv*Mff*pmi + 
     - 2*CISR1*E1J**2*ppl + 2*CISR8*E1J**2*ppl + 2*CISR1*E1J*E2J*ppl + 
     - 2*CISR8*E1J*E2J*ppl + 4*CISR1*kJ*kv*ppl + 4*CISR8*kJ*kv*ppl - 
     -  CISR1*MB**2*ppl - CISR8*MB**2*ppl 
     - - CISR1*Mf**2*ppl - CISR8*Mf**2*ppl + 
     -  CISR1*Mff**2*ppl + CISR8*Mff**2*ppl + 
     - 2*CISR3*(Epl*(Mf + Mff) + (-EvJ + kv)*
     - (Emi - pmi) - Mf*ppl + Mff*ppl) + 
     - 2*CISR6*(-(Epl*(Mf + Mff)) + (EvJ + kv)*(Emi + pmi) - Mf*ppl + 
     -  Mff*ppl) + 2*EvJ*(CISR8*
     -   (Emi*(-Mf + Mff) + (Mf + Mff)*pmi + (E1J - E2J)*(Epl - ppl)) + 
     -  CISR1*(Emi*(Mf - Mff) + (Mf + Mff)*pmi - 
     - (E1J - E2J)*(Epl + ppl))) + 
     - 2*(CISR2*Emi*kv + CISR7*Emi*kv + CISR2*Epl*Mf - CISR7*Epl*Mf + 
     -  CISR2*Epl*Mff - CISR7*Epl*Mff + CISR2*kv*pmi - CISR7*kv*pmi + 
     -  EvJ*((CISR2 - CISR7)*Emi + (CISR2 + CISR7)*pmi) + 
     -  2*CISR4*(Epl - ppl) + CISR2*Mf*ppl 
     - + CISR7*Mf*ppl - CISR2*Mff*ppl - 
     -  CISR7*Mff*ppl - 2*CISR5*(Epl + ppl))*Mkcha)
     -     )/delta
      
      
      return
      end

