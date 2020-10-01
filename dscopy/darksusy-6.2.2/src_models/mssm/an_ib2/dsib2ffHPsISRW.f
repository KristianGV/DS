*******************************************************************
*** subroutine dsib2ffWsISRX add to the common block array amp(i,j), i = 0,
*** the helicity amplitudes of the s channel ISR diagrams. 
*** W emitted from the initial legs of the s channel diagrams.
*** Exchanged particle index icha stays for chargino.
***
*** Subroutines: dsib2ffHPsISRH: Hc exchanged; ISR1 & ISR2
***              dsib2ffHPsISRW: W wxhanged; ISR1 & ISR2
*** 
*** Author: Francesca Calore, 2014-02-22
***         Francesca Calore, 2015-02-15 (bug solved cf with CH)
*******************************************************************
c     W exchanged ISR1 & ISR2
      subroutine dsib2ffHPsISRW(f, icha) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, icha
      integer ff, i

      real*8 Mkw, Mkcha(2)
      complex*16 delta
      complex*16 CISR1, CISR2, CISR3, CISR4
      complex*16 CCISR1, CCISR2, CCISR3, CCISR4
      complex*16 amp00, amp01, amp02, amp03

      ff = c_fbartype

c... masses/mx exchanged particles
      Mkw = mass(kw)/mx
      do i = 1, 2
         Mkcha(i) = mass(kcha(i))/mx
      end do 

c... denominators
      delta = (dcmplx(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - Mkw**2,
     -     Mkw*width(kw)/mx)* dcmplx((2*E1J*E2J - 2*E1J*EvJ 
     -     - 2*E2J*EvJ + 2*kJ**2 + MB**2 + Mf**2 + 
     -     Mff**2 - 4*Mkcha(icha)**2)/4.,
     -     Mkcha(icha)*width(kcha(icha))/mx*0d0))

c... couplings
      if(mod(c_ftype,2).EQ.0) then
         CCISR1 = Conjg(gr(kw,kn(1),kcha(icha)))*
     -        gl(khc,kn(1),kcha(icha))*gl(kw,f,ff)
         CCISR2 = Conjg(gl(kw,kn(1),kcha(icha)))*
     -        gl(khc,kn(1),kcha(icha))*gl(kw,f,ff)
         CCISR3 = Conjg(gr(kw,kn(1),kcha(icha)))*gl(kw,f,ff)*
     -        gr(khc,kn(1),kcha(icha))
         CCISR4 = Conjg(gl(kw,kn(1),kcha(icha)))*gl(kw,f,ff)*
     -        gr(khc,kn(1),kcha(icha))

         CISR1 = CCISR1
         CISR2 = CCISR2
         CISR3 = CCISR3
         CISR4 = CCISR4
         
      else
         CISR1 = Conjg(gr(khc,kn(1),kcha(icha)))*
     -        gl(kw,f,ff)*gl(kw,kn(1),kcha(icha))
         CISR2 = Conjg(gr(khc,kn(1),kcha(icha)))*
     -        gl(kw,f,ff)*gr(kw,kn(1),kcha(icha))         
         CISR3 = Conjg(gl(khc,kn(1),kcha(icha)))*
     -        gl(kw,f,ff)*gl(kw,kn(1),kcha(icha))
         CISR4 =Conjg(gl(khc,kn(1),kcha(icha)))*
     -        gl(kw,f,ff)*gr(kw,kn(1),kcha(icha))
      endif
 
c... helicity amplitudes
         
      amp00 = (((0,-1)*((CISR2 - CISR3)*(-(Emi*EvJ*(E1J + E2J - Mkw)*
     - (E1J + E2J + Mkw)) + (E1J + E2J)*Emi*(2*E1J*E2J + 2*kJ**2 
     - + Mf**2 + Mff**2 - Mkw**2) + 
     - CW*kv*Mkw**2*pmi) + (CISR1 - CISR4)*
     - (-(Emi*EvJ*(E1J + E2J - Mkw)*(E1J + E2J + Mkw)) - 
     - (E1J + E2J)*Emi*(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - Mkw**2) + 
     - CW*kv*Mkw**2*pmi)*Mkcha(icha)))/Mkw**2)/delta
         
      amp01 = ((kv*(Epl + ppl)*SW*(CISR2 - CISR3 +
     -     (CISR1 - CISR4)*Mkcha(icha)))/Sqrt(2.))/delta
         
      amp02 = (((0,-1)*((CISR2 - CISR3)*(CW*Emi*kv*Mkw**2 - 
     - EvJ*(E1J + E2J - Mkw)*(E1J + E2J + Mkw)*pmi + 
     - (E1J + E2J)*(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 -Mkw**2)*pmi) + 
     - (CISR1 - CISR4)*(CW*Emi*kv*Mkw**2 - 
     - EvJ*(E1J + E2J - Mkw)*(E1J + E2J + Mkw)*pmi - 
     - (E1J + E2J)*(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - Mkw**2)*pmi)*
     - Mkcha(icha)))/Mkw**2)/delta
         
      amp03 = ((kv*(Epl - ppl)*SW*(CISR2 - CISR3 + 
     -     (CISR1 - CISR4)*Mkcha(icha)))/Sqrt(2.))/delta


      amp(0, 0) = amp(0, 0) + 2 * amp00
      amp(0, 1) = amp(0, 1) + 2 * amp01
      amp(0, 2) = amp(0, 2) + 2 * amp02
      amp(0, 3) = amp(0, 3) + 2 * amp03

      return
      end
