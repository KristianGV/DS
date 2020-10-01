*******************************************************************
*** subroutine dsib2ffHPsFSR add to the common block array amp(i,j), i = 0,
*** the helicity amplitudes of the s channel FSR diagrams. 
*** The Hp is emitted from the final legs of the s channel diagrams. 
*** 
*** Author: Francesca Calore, 2014-02-22
***         Francesca Calore, 2015-02-15 (bug solved cf with CH)
*******************************************************************
c...  H3 mediated diagrams
      subroutine dsib2ffHPsFSRH3(f, leg) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg
      integer ff

      real*8 Mh3
      complex*16 CH31, CH32
      complex*16 deltaH3
      complex*16 ampH00, ampH01, ampH02, ampH03 

      ff = c_fbartype

c... masses/mx exchanged particles
      Mh3 = mass(kh3)/mx
c... couplings
c... H3      
c      if (leg.eq.1) then
         CH31 =gl(khc,f,ff)*gr(kh3,ff,ff)*
     -    (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))
         CH32 =gr(khc,f,ff)*gr(kh3,ff,ff)*
     -    (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))


         deltaH3 = (dcmplx(2*E1J*EvJ 
     -        - 2*CW*kJ*kv + Mf**2 - Mff**2 + MB**2,
     -        Mff*width(ff)/mx)*dcmplx(2*E1J*E2J
     -        + 2*E1J*EvJ + 2*E2J*EvJ + 
     -        2*kJ**2 + Mf**2 + Mff**2 - Mh3**2
     -        + MB**2,Mh3*width(kh3)/mx))

c      else 
       if (leg.eq.2) then

         CH31 =gl(khc,f,ff)*gr(kh3,f,f)*
     -    (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))
         CH32 =gr(khc,f,ff)*gr(kh3,f,f)*
     -    (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))

         deltaH3 = (dcmplx(2*E2J*EvJ
     -        + 2*CW*kJ*kv - Mf**2 + Mff**2 + MB**2,
     -        Mf*width(f)/mx)*dcmplx(2*E1J*E2J + 2*E1J*EvJ
     -        + 2*E2J*EvJ + 2*kJ**2 + Mf**2 + Mff**2
     -        - Mh3**2 + MB**2,Mh3*width(kh3)/mx))

      endif


c... helicity amplitudes
c... H3 
      ampH00 =((0,2)*(CH31 + CH32)*(Emi*EvJ
     -     + Epl*(Mf + Mff) + CW*kv*pmi))/deltaH3

      ampH01 =(Sqrt(2.)*kv*(CH31*(Epl - ppl) 
     -     - CH32*(Epl + ppl))*SW)/deltaH3

      ampH02 =((0,-2)*(CH31 - CH32)*(CW*Emi*kv 
     -     + EvJ*pmi + Mf*ppl - Mff*ppl))/deltaH3

      ampH03 =(Sqrt(2.)*kv*(CH32*(-Epl + ppl)
     -     + CH31*(Epl + ppl))*SW)/deltaH3

      if (leg.eq.1) then
c... FSR1
         amp(0, 0) = amp(0, 0) + ampH00 
         amp(0, 1) = amp(0, 1) + ampH01 
         amp(0, 2) = amp(0, 2) + ampH02 
         amp(0, 3) = amp(0, 3) + ampH03 
      else if (leg.eq.2) then
c... FSR2
         amp(0, 0) = amp(0, 0) + ampH00 
         amp(0, 1) = amp(0, 1) - ampH01 
         amp(0, 2) = amp(0, 2) - ampH02 
         amp(0, 3) = amp(0, 3) - ampH03 
      endif

      return
      end
