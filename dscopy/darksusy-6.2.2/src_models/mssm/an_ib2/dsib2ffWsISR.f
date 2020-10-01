*******************************************************************
*** subroutine dsib2ffWsISRX add to the common block array amp(i,j) 
*** the helicity amplitudes of the s channel ISR diagrams. 
*** W emitted from the initial legs of the s channel diagrams.
*** Exchanged particle index icha stays for chargino.
***
*** Subroutines: dsib2ffWsISRH: Hc exchanged; ISR1 & ISR2
***              dsib2ffWsISRW: W wxhanged; ISR1 & ISR2
*** 
*** Author: Francesca Calore, 2013-04-10
***         Francesca Calore, 2015-01-28
*******************************************************************
c     Hc exchanged ISR1 & ISR2
      subroutine dsib2ffWsISRH(f, icha) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, icha

      integer ff, i
      real*8 Mkhc, Mkcha(2)
      complex*16 delta
      complex*16 CISR1, CISR2, CISR3, CISR4, CISR5, CISR6, CISR7, CISR8
      complex*16 amp00ISR1, amp02ISR1, amp00ISR2, amp02ISR2

      ff = c_fbartype

c... masses/mx exchanged particles
      Mkhc = mass(khc)/mx
      do i = 1,2
         Mkcha(i) = mass(kcha(i))/mx
      end do

c... denominators
      delta = (dcmplx(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - Mkhc**2
     -     ,Mkhc*width(khc)/mx)*
     -     dcmplx((2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ 
     -     + 2*kJ**2 + Mf**2 + Mff**2 + MB**2 - 4*Mkcha(icha)**2)/4.
     -     ,Mkcha(icha)*width(kcha(icha))/mx*0d0))

c... couplings (no simmetries)
      CISR1 = Conjg(gr(kw,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -     gl(khc,kn(1),kcha(icha))
      CISR2 = Conjg(gr(kw,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -     gr(khc,kn(1),kcha(icha))
      CISR3 = Conjg(gl(kw,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -     gl(khc,kn(1),kcha(icha))
      CISR4 = Conjg(gl(kw,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -     gr(khc,kn(1),kcha(icha))
      CISR5 = Conjg(gr(kw,kn(1),kcha(icha)))*gl(khc,kn(1),kcha(icha))*
     -     gr(khc,f,ff)
      CISR6 = Conjg(gr(kw,kn(1),kcha(icha)))*gr(khc,f,ff)*
     -     gr(khc,kn(1),kcha(icha))
      CISR7 = Conjg(gl(kw,kn(1),kcha(icha)))*gl(khc,kn(1),kcha(icha))*
     -     gr(khc,f,ff)
      CISR8 = Conjg(gl(kw,kn(1),kcha(icha)))*gr(khc,f,ff)*
     -     gr(khc,kn(1),kcha(icha))

      amp00ISR1 = (
     -     ((0,1)*(E1J + E2J)*Epl*kv*(CISR2 - CISR3 - CISR6 + CISR7 + 
     -     (-CISR1 + CISR4 + CISR5 - CISR8)*Mkcha(icha)))/MB)/delta

      amp02ISR1 = (
     -     ((0,1)*(E1J + E2J)*kv*ppl*(CISR2 - CISR3 + CISR6 - CISR7 + 
     -     (-CISR1 + CISR4 - CISR5 + CISR8)*Mkcha(icha)))/MB)/delta

      amp00ISR2 = (
     -    ((0,1)*(E1J + E2J)*Epl*kv*(CISR2 - CISR3 - CISR6 + CISR7 + 
     -    (-CISR1 + CISR4 + CISR5 - CISR8)*Mkcha(icha)))/MB)/delta

      amp02ISR2 = (
     -     ((0,1)*(E1J + E2J)*kv*ppl*(CISR2 - CISR3 + CISR6 - CISR7 + 
     -     (-CISR1 + CISR4 - CISR5 + CISR8)*Mkcha(icha)))/MB)/delta

      amp(0, 0) = amp(0, 0) + amp00ISR1 + amp00ISR2
      amp(0, 2) = amp(0, 2) + amp02ISR1 + amp02ISR2

      return
      end

c     W exchanged ISR1 & ISR2
      subroutine dsib2ffWsISRW(f, leg, icha) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg, icha
      integer ff, i

      real*8 m11, m22, q1, tmpppl
      real*8 Mw, Mkcha(2)
      complex*16 delta
      complex*16 CISR1, CISR2, CISR3
      complex*16 amp00, amp01, amp02, amp03, amp10, amp11, amp12, amp13,
     &     amp20, amp21, amp22, amp23

      ff = c_fbartype

c... Defines ISR1 and ISR2 amplitudes through their simmetry relations: only ppl exchange
c      if (leg.eq.1) then
         m11 = Mf
         m22 = Mff
         q1 = kJ
         tmpppl = ppl
c      else 
      if (leg.eq.2) then
         m11 = Mf
         m22 = Mff
         q1 = kJ
         tmpppl = ppl
      endif

c... masses/mx exchanged particles
      Mw = mass(kw)/mx
      do i = 1, 2
         Mkcha(i) = mass(kcha(i))/mx
      end do 

c... denominators
      delta = (dcmplx(2*E1J*E2J + m11**2 + m22**2 - Mw**2 + 2*q1**2
     -     ,Mw*width(kw)/mx)*
     -     dcmplx((2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ 
     -     + m11**2 + m22**2 + MB**2 + 2*q1**2 - 4*Mkcha(icha)**2)/4.
     -     ,Mkcha(icha)*width(kcha(icha))/mx*0d0))

c... couplings
      CISR1 = Conjg(gr(kw,kn(1),kcha(icha)))*gl(kw,f,ff)*
     -     gl(kw,kn(1),kcha(icha)) ! same as (checked numerically) Conjg(gl(kw,kn(1),kcha(icha)))*gl(kw,f,ff)*gr(kw,kn(1),kcha(icha))
      CISR2 = Conjg(gr(kw,kn(1),kcha(icha)))*gl(kw,f,ff)*
     -     gr(kw,kn(1),kcha(icha))
      CISR3 = Conjg(gl(kw,kn(1),kcha(icha)))*gl(kw,f,ff)*
     -     gl(kw,kn(1),kcha(icha))

c... helicity amplitudes
      amp00 = (
     -     ((0,-0.5)*(CISR2 - CISR3)*(CW*EvJ*Mw**2*pmi*
     -     (2*E1J*E2J + m11**2 + m22**2 - MB**2 + 2*q1**2) + 
     -     Emi*kv*(2*E1J*E2J**3 - 2*E1J**2*Mw**2 + 
     -     m11**2*((E1J + E2J)**2 + Mw**2) + 
     -     m22**2*((E1J + E2J)**2 + Mw**2) + E1J**2*MB**2-Mw**2*MB**2 + 
     -     2*E1J**2*q1**2 + 2*Mw**2*q1**2 + 
     -     E2J**2*(4*E1J**2 - 2*Mw**2 + MB**2 + 2*q1**2) + 
     -     2*E1J*E2J*(E1J**2 - Mw**2 + MB**2 + 2*q1**2)))
     -     )/(Mw**2*MB))/delta
      
      amp01 = (
     -     ((CISR2 - CISR3)*EvJ*(2*E1J*E2J 
     -     + m11**2 + m22**2 - MB**2 + 2*q1**2)*SW*
     -     (Epl + tmpppl))/(2.*Sqrt(2.)*MB))/delta
      
      amp02 = (
     -     ((0,-0.5)*(CISR2 - CISR3)*(m11**2*
     -     (CW*EvJ*Emi*Mw**2 + kv*((E1J + E2J)**2 + Mw**2)*pmi) + 
     -     m22**2*(CW*EvJ*Emi*Mw**2 +kv*((E1J + E2J)**2 + Mw**2)*pmi) + 
     -     CW*EvJ*Emi*Mw**2*(-MB**2 + 2*(E1J*E2J + q1**2)) + 
     -     kv*pmi*((E1J + E2J)**2*(2*E1J*E2J + MB**2 + 2*q1**2) - 
     -     Mw**2*(MB**2 + 2*(E1J**2 + E1J*E2J + E2J**2 - q1**2)))))/
     -     (Mw**2*MB))/delta

      amp03 = (
     -     ((CISR2 - CISR3)*EvJ*(2*E1J*E2J
     -     + m11**2 + m22**2 - MB**2 + 2*q1**2)*SW*
     -     (Epl - tmpppl))/(2.*Sqrt(2.)*MB))/delta
      
      amp10 = (
     -     -(pmi*(-(CISR3*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2))) + 
     -     CISR2*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     SW)/(2.*Sqrt(2.)))/delta
      
      amp11 =(
     -     (0,0.25)*(1 + CW)*(-(CISR3*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2))) + 
     -     CISR2*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     (Epl + tmpppl))/delta
      
      amp12= (
     -     -(Emi*(-(CISR3*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2))) + 
     -     CISR2*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     SW)/(2.*Sqrt(2.)))/delta
      
      amp13 = (
     -     (0,0.25)*(-1 + CW)*(-(CISR3*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2))) + 
     -     CISR2*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     (Epl - tmpppl))/delta
      
      amp20 = (
     -     (pmi*(-(CISR2*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2))) + 
     -     CISR3*(m11**2 + m22**2 - MB**2
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     SW)/(2.*Sqrt(2.)))/delta
      
      amp21 = (
     -     (0,0.25)*(-1 + CW)*(CISR2*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2)) - 
     -     CISR3*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     (Epl + tmpppl))/delta
      
      amp22 = (
     -     (Emi*(-(CISR2*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2))) + 
     -     CISR3*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     SW)/(2.*Sqrt(2.)))/delta
      
      amp23 = (
     -     (0,0.25)*(1 + CW)*(CISR2*(m11**2 + m22**2 - MB**2 + 
     -     2*(E1J*E2J - (E1J + E2J)*kv + q1**2)) - 
     -     CISR3*(m11**2 + m22**2 - MB**2 
     -     + 2*(E1J*kv + E2J*(E1J + kv) + q1**2)))*
     -     (Epl - tmpppl))/delta
 
      if (leg.eq.1) then
c... ISR1
         amp(0, 0) = amp(0, 0) + amp00
         amp(0, 1) = amp(0, 1) + amp01
         amp(0, 2) = amp(0, 2) + amp02
         amp(0, 3) = amp(0, 3) + amp03
         amp(1, 0) = amp(1, 0) + amp10
         amp(1, 1) = amp(1, 1) + amp11
         amp(1, 2) = amp(1, 2) + amp12
         amp(1, 3) = amp(1, 3) + amp13
         amp(-1, 0) = amp(-1, 0) + amp20
         amp(-1, 1) = amp(-1, 1) + amp21
         amp(-1, 2) = amp(-1, 2) + amp22
         amp(-1, 3) = amp(-1, 3) + amp23
      else if (leg.eq.2) then
c... ISR2
         amp(0, 0) = amp(0, 0) + amp00
         amp(0, 1) = amp(0, 1) + amp01
         amp(0, 2) = amp(0, 2) + amp02
         amp(0, 3) = amp(0, 3) + amp03
         amp(1, 0) = amp(1, 0) + amp10
         amp(1, 1) = amp(1, 1) + amp11
         amp(1, 2) = amp(1, 2) + amp12
         amp(1, 3) = amp(1, 3) + amp13
         amp(-1, 0) = amp(-1, 0) + amp20
         amp(-1, 1) = amp(-1, 1) + amp21
         amp(-1, 2) = amp(-1, 2) + amp22
         amp(-1, 3) = amp(-1, 3) + amp23
      endif    
      

      return
      end
