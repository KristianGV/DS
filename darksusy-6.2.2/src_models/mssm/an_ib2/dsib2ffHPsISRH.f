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
c     Hc exchanged ISR1 & ISR2
      subroutine dsib2ffHPsISRH(f, icha) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, icha

      integer ff, i
      real*8 Mkhc, Mkcha(2)
      complex*16 delta
      complex*16 CISR1, CISR2, CISR3, CISR4, CISR5, CISR6, CISR7, CISR8
      complex*16 CCISR1, CCISR2, CCISR3, CCISR4, 
     &     CCISR5, CCISR6, CCISR7, CCISR8
      complex*16 amp00ISR2, amp02ISR2

      ff = c_fbartype

c... masses/mx exchanged particles
      Mkhc = mass(khc)/mx
      do i = 1,2
         Mkcha(i) = mass(kcha(i))/mx
      end do

c... denominators
      delta = dcmplx(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2
     -     - Mkhc**2,Mkhc*width(khc)/mx)*
     -     dcmplx((2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ
     -     + 2*kJ**2 + Mf**2 + Mff**2 + 
     -     MB**2 - 4*Mkcha(icha)**2)/4.,
     -     Mkcha(icha)*width(kcha(icha))/mx*0d0)

c... couplings (no simmetries)
      if(mod(c_ftype,2).EQ.0) then

         CCISR1 = Conjg(gr(khc,kn(1),kcha(icha)))*
     -    gl(khc,kn(1),kcha(icha))*gr(khc,f,ff)
         CCISR2 = Conjg(gl(khc,kn(1),kcha(icha)))*
     -    gl(khc,kn(1),kcha(icha))*gr(khc,f,ff)
         CCISR3 = Conjg(gr(khc,kn(1),kcha(icha)))*gr(khc,f,ff)*
     -    gr(khc,kn(1),kcha(icha))
         CCISR4 = Conjg(gl(khc,kn(1),kcha(icha)))*gr(khc,f,ff)*
     -    gr(khc,kn(1),kcha(icha))
         CCISR5 = Conjg(gr(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -    gl(khc,kn(1),kcha(icha))
         CCISR6 = Conjg(gl(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -    gl(khc,kn(1),kcha(icha))
         CCISR7 = Conjg(gr(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -    gr(khc,kn(1),kcha(icha))
         CCISR8 = Conjg(gl(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -    gr(khc,kn(1),kcha(icha))

         CISR1 = - CCISR1
         CISR2 = - CCISR2
         CISR3 = - CCISR3
         CISR4 = - CCISR4
         CISR5 = - CCISR5
         CISR6 = - CCISR6
         CISR7 = - CCISR7
         CISR8 = - CCISR8

      else
         CISR1 = Conjg(gr(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -        gl(khc,kn(1),kcha(icha))
         CISR2 = Conjg(gr(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -        gr(khc,kn(1),kcha(icha))
         CISR3 =Conjg(gl(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -        gl(khc,kn(1),kcha(icha))
         CISR4 = Conjg(gl(khc,kn(1),kcha(icha)))*gl(khc,f,ff)*
     -        gr(khc,kn(1),kcha(icha))
         CISR5 =Conjg(gr(khc,kn(1),kcha(icha)))*
     -        gl(khc,kn(1),kcha(icha))*gr(khc,f,ff)
         CISR6 = Conjg(gr(khc,kn(1),kcha(icha)))*gr(khc,f,ff)*
     -        gr(khc,kn(1),kcha(icha))
         CISR7 =Conjg(gl(khc,kn(1),kcha(icha)))*
     -        gl(khc,kn(1),kcha(icha))*gr(khc,f,ff)
      CISR8 =Conjg(gl(khc,kn(1),kcha(icha)))*gr(khc,f,ff)*
     -        gr(khc,kn(1),kcha(icha))

      endif

c... ISR1 == ISR2
      amp00ISR2 = ((0,0.5)*Epl*((CISR2 - CISR3 - CISR6 + CISR7)*
     -     (2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - MB**2) + 
     -     4*(-CISR1 + CISR4 + CISR5 - CISR8)*Mkcha(icha)))/delta

      amp02ISR2 = ((0,0.5)*ppl*((CISR2 - CISR3 + CISR6 - CISR7)*
     -  (2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - MB**2) + 
     -    4*(-CISR1 + CISR4 - CISR5 + CISR8)*Mkcha(icha)))/delta

      amp(0, 0) = amp(0, 0) + 2 * amp00ISR2
      amp(0, 2) = amp(0, 2) + 2 * amp02ISR2

      return
      end
