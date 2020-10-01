*******************************************************************
*** subroutine dsib2ffHPtuISRX adds to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel ISR diagrams. 
*** HP emitted by initial legs 1 (ISR1) or 2 (ISR2). isf is the index of the exchanged sfermion.
***
*** Subroutines: dsib2ffHPtuISR1:  HP emitted from p1 leg
***              dsib2ffHPtuISR2:  HP emitted from p2 leg
***
*** IMPORTANT: do NOT use separately dsib2ffHPtuISR1 and dsib2ffHPtuISR2.
***            The correct t+u channel for the 4 involved diagrams is given by ISR1 + ISR2. 
***            Separately, the result is NOT the same as ISR1_t + ISR1_u (ISR2_t + ISR2_u), 
***            being ISR1(t) \equiv ISR2(u) and ISR2(t) \equiv ISR1(u)
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
***         Francesca Calore, 2015-02-15 (bug solved cf with CH)
*******************************************************************
      subroutine dsib2ffHPtuISR2(f, isf, icha)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf, icha

      integer ff, sf
      real*8 Msf, Mkcha
      integer dsib2getsfermion

      complex*16 delta
      complex*16 CISR1, CISR2, CISR3, CISR4, CISR5, CISR6, CISR7, CISR8
      complex*16 CCISR1, CCISR2, CCISR3, CCISR4, CCISR5,
     &     CCISR6, CCISR7, CCISR8

      sf  = dsib2getsfermion(f, isf)
      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics

      ff   = c_fbartype

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
c
c      do j = 1,2
         Mkcha = mass(kcha(icha))/mx
c      end do   

c... couplings
      if(mod(c_ftype,2).EQ.0) then
         CCISR1 = Conjg(gl(sf,ff,kcha(icha)))*
     -    gl(sf,f,kn(1))*gr(khc,kn(1),kcha(icha))
         CCISR2 = Conjg(gl(sf,ff,kcha(icha)))*
     -    gr(khc,kn(1),kcha(icha))*gr(sf,f,kn(1))
         CCISR3 = Conjg(gr(sf,ff,kcha(icha)))*
     -    gl(sf,f,kn(1))*gr(khc,kn(1),kcha(icha))
         CCISR4 =Conjg(gr(sf,ff,kcha(icha)))*
     -    gr(khc,kn(1),kcha(icha))*gr(sf,f,kn(1))
         CCISR5 = Conjg(gl(sf,ff,kcha(icha)))*
     -    gl(khc,kn(1),kcha(icha))*gl(sf,f,kn(1))
         CCISR6 = Conjg(gl(sf,ff,kcha(icha)))*
     -    gl(khc,kn(1),kcha(icha))*gr(sf,f,kn(1))
         CCISR7 = Conjg(gr(sf,ff,kcha(icha)))*
     -    gl(khc,kn(1),kcha(icha))*gl(sf,f,kn(1))
         CCISR8 =Conjg(gr(sf,ff,kcha(icha)))*
     -    gl(khc,kn(1),kcha(icha))*gr(sf,f,kn(1))

         CISR1  = - CCISR5
         CISR2  = - CCISR7 
         CISR3  = - CCISR6
         CISR4  = - CCISR8 
         CISR5  = - CCISR1 
         CISR6  = - CCISR3 
         CISR7  = - CCISR2 
         CISR8  = - CCISR4 
      else
         CISR1 =Conjg(gl(sf,ff,kcha(icha)))*
     -    Conjg(gr(khc,kn(1),kcha(icha)))*gl(sf,f,kn(1))
         CISR2 =Conjg(gr(khc,kn(1),kcha(icha)))*
     -    Conjg(gr(sf,ff,kcha(icha)))*gl(sf,f,kn(1))
         CISR3 = Conjg(gl(sf,ff,kcha(icha)))*
     -    Conjg(gr(khc,kn(1),kcha(icha)))*gr(sf,f,kn(1))
         CISR4 =Conjg(gr(khc,kn(1),kcha(icha)))*
     -    Conjg(gr(sf,ff,kcha(icha)))*gr(sf,f,kn(1))
         CISR5 = Conjg(gl(khc,kn(1),kcha(icha)))*
     -    Conjg(gl(sf,ff,kcha(icha)))*gl(sf,f,kn(1))
         CISR6 =Conjg(gl(khc,kn(1),kcha(icha)))*
     -    Conjg(gr(sf,ff,kcha(icha)))*gl(sf,f,kn(1))
         CISR7 =Conjg(gl(khc,kn(1),kcha(icha)))*
     -    Conjg(gl(sf,ff,kcha(icha)))*gr(sf,f,kn(1))
         CISR8 =Conjg(gl(khc,kn(1),kcha(icha)))*
     -    Conjg(gr(sf,ff,kcha(icha)))*gr(sf,f,kn(1))
      endif

c... denominators
      delta = (
     -     dcmplx((2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ + 2*kJ**2 + MB**2 + 
     -     Mf**2 + Mff**2 - 4*Mkcha**2)/4.,
     -     Mkcha*width(kcha(icha))/mx*0d0)*
     -     dcmplx((-2*E1J*E2J - 2*E1J*EvJ + 2*E2J*EvJ
     -     - 2*kJ**2 + 4*CW*kJ*kv + 
     -     MB**2 + Mf**2 + Mff**2 - 4*Msf**2)/4.,
     -     Msf*width(sf)/mx*0d0))

      amp(0, 0) = amp(0, 0) + 2*((0,0.25)*
     - (2*CISR1*(-(Emi*EvJ) + Epl*(Mf + Mff) - CW*kv*pmi) +
     -  2*CISR8*(-(Emi*EvJ) + Epl*(Mf + Mff) - CW*kv*pmi) - 
     -  (CISR3 + CISR6)*(Epl*(-2*(E1J*(E1J + E2J) + (E1J - E2J)*EvJ) + 
     -  4*CW*kJ*kv + MB**2 + Mf**2 - Mff**2) + 
     -  2*(Mf - Mff)*(Emi*EvJ + CW*kv*pmi)) + 
     -  2*(2*CISR2*Epl + 2*CISR7*Epl + 
     -  (CISR4 + CISR5)*(Emi*EvJ + Epl*(Mf + Mff) 
     -  + CW*kv*pmi))*Mkcha))/delta

      amp(0, 1) = amp(0, 1) + 2*((kv*SW*(Epl*(-CISR1 + CISR8 -
     - (CISR3 - CISR6)*(Mf + Mff)) + 
     - (CISR1 + CISR8 + (CISR3 + CISR6)*(Mf - Mff))*ppl - 
     - (CISR5*(-Epl + ppl) + CISR4*(Epl + ppl))*
     - Mkcha))/(2.*Sqrt(2.)))/delta

      amp(0, 2) = amp(0, 2) + 2*((0,0.25)*(2*CISR1*(CW*Emi*kv 
     - + EvJ*pmi - Mf*ppl + Mff*ppl) - 
     - 2*CISR8*(CW*Emi*kv + EvJ*pmi - Mf*ppl + Mff*ppl) + 
     - (CISR3 - CISR6)*(2*Mf*(CW*Emi*kv + EvJ*pmi) + 
     - 2*Mff*(CW*Emi*kv + EvJ*pmi) + 
     - (-2*(E1J*(E1J + E2J) + (E1J - E2J)*EvJ) + 4*CW*kJ*kv + MB**2)*ppl + 
     - Mf**2*ppl - Mff**2*ppl) + 
     - 2*(2*CISR2*ppl - 2*CISR7*ppl + 
     - (CISR4 - CISR5)*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl))*
     - Mkcha))/delta

      amp(0, 3) = amp(0, 3) + 2*((kv*SW*(Epl*(-CISR1 + CISR8 - 
     - (CISR3 - CISR6)*(Mf + Mff)) - 
     - (CISR1 + CISR8 + (CISR3 + CISR6)*(Mf - Mff))*ppl + 
     - (CISR4*(-Epl + ppl) + CISR5*(Epl + ppl))*
     - Mkcha))/(2.*Sqrt(2.)))/delta

      return
      end
