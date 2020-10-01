*******************************************************************
*** subroutine dsib2ffWsVIB add to the common block array amp(i,j) 
*** the helicity amplitudes of the s channel VIB diagrams. 
*** The W is emitted from the mediators (H3, Z), becoming Hc or W respectively.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffWsVIB(f) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f
      integer ff

      real*8 Mkhc, Mh3, Mw, Mz
      complex*16 deltaH3, deltaZ
      complex*16 ampH300, ampH302, ampZ00, ampZ02
      complex*16 CIBZ
      complex*16 amp10, amp12 ! Z channel
      real*8 epsilon

      ff = c_fbartype
      if(mod(c_ftype,2).EQ.0) then
         epsilon = -1d0 ! epsilon sign defines the sign of the s VIB Z diagram
      else
         epsilon = +1d0
      endif 

c..  masses/mx of exchanged particles
      Mkhc = mass(khc)/mx
      Mh3 = mass(kh3)/mx
      Mw = mass(kw)/mx
      Mz = mass(kz)/mx

c... denominators
c... H3
      deltaH3 = (dcmplx(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - Mkhc**2,
     -     Mkhc*width(khc)/mx)*
     -     dcmplx(2*E1J*E2J + 2*E1J*EvJ + 2*E2J*EvJ + 2*kJ**2 
     -     + Mf**2 + Mff**2 - Mh3**2 + MB**2,Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(2*E1J*E2J + 2*kJ**2 + Mf**2 + Mff**2 - Mw**2,
     -     Mw*width(kw)/mx)*
     -     dcmplx(2*E1J*E2J + 2*E1J*EvJ + 2*E2J*EvJ + 2*kJ**2 
     -     + Mf**2 + Mff**2 - Mz**2 + MB**2,Mz*width(kz)/mx))

c... couplings
c... Z
      CIBZ = gl(kz,kw,kw)*gl(kw,f,ff)*
     -     (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))

c... helicity amplitudes
c... firstly H3 + Z (amp00, amp02)
      ampH300 = 0.d0
      ampH302 = 0.d0
      ampZ00 = 0.d0
      ampZ02 = 0.d0
c... H3
      ampH300 = + (
     -     ((0,-4)*(E1J + E2J)*Epl*kv*gl(kw,khc,kh3)*
     -     (gl(khc,f,ff) - gr(khc,f,ff))*(Conjg(gr(kh3,kn(1),kn(1)))
     -     - gr(kh3,kn(1),kn(1))))/MB)/deltaH3
      
      ampH302 = + (
     -     ((0,-4)*(E1J + E2J)*kv*ppl*gl(kw,khc,kh3)*
     -     (gl(khc,f,ff) + gr(khc,f,ff))*(Conjg(gr(kh3,kn(1),kn(1)))
     -     - gr(kh3,kn(1),kn(1))))/MB)/deltaH3
c... Z 
      ampZ00 = + (
     -     ((0,-1)*CIBZ*(-2*E1J**2 - 2*E1J*E2J  
     -     - 2*E1J*EvJ - 2*E2J*EvJ + Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Emi*kv + CW*EvJ*pmi))/(Mz**2*MB))/deltaZ

      ampZ02 = + (
     -     ((0,-1)*CIBZ*(-2*E1J**2 - 2*E1J*E2J 
     -     - 2*E1J*EvJ - 2*E2J*EvJ + Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (CW*EvJ*Emi + kv*pmi))/(Mz**2*MB))/deltaZ


      amp(0, 0) = amp(0, 0) + ampH300 + (epsilon*ampZ00)
      amp(0, 2) = amp(0, 2) + ampH302 + ( epsilon*ampZ02)

      amp(0, 1) = amp(0, 1) +  (epsilon*(
     -     (CIBZ*EvJ*(-2*E1J**2 - 2*E1J*E2J
     -      - 2*E1J*EvJ - 2*E2J*EvJ + Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Epl + ppl)*SW)/(Sqrt(2.)*Mz**2*MB))/deltaZ)

      amp(0, 3) = amp(0, 3) + (epsilon*(
     -     (CIBZ*EvJ*(-2*E1J**2 - 2*E1J*E2J  
     -     - 2*E1J*EvJ - 2*E2J*EvJ + Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Epl - ppl)*SW)/(Sqrt(2.)*Mz**2*MB))/deltaZ)

      amp10 = (
     -     -((CIBZ*(-2*E1J**2 - 2*E1J*E2J
     -     - 2*E1J*EvJ - 2*E2J*EvJ + Mf**2 - Mff**2 + 
     -     Mz**2 - MB**2)*(-2*E1J**2 - 2*E1J*E2J + Mf**2 
     -     - Mff**2 + MB**2)*pmi*SW)/(Sqrt(2.)*Mz**2)))/deltaZ

      amp(1, 0) = amp(1, 0) + (epsilon*amp10)
      amp(-1, 0) = amp(-1, 0) + (epsilon*amp10)

      amp(1, 1) = amp(1, 1) + (epsilon*(
     -     ((0,0.5)*CIBZ*(1 + CW)*(-2*E1J**2 
     -     - 2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ + 
     -     Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Epl + ppl))/Mz**2)/deltaZ)

      amp12 = (
     -     -((CIBZ*Emi*(-2*E1J**2 - 2*E1J*E2J  
     -     - 2*E1J*EvJ - 2*E2J*EvJ + Mf**2 -Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*SW)
     -     /(Sqrt(2.)*Mz**2)))/deltaZ

      amp(1, 2) = amp(1, 2) + (epsilon*amp12)
      amp(-1, 2) = amp(-1, 2) + (epsilon*amp12)

      amp(1, 3) = amp(1, 3) + (epsilon*(
     -     ((0,0.5)*CIBZ*(-1 + CW)*(-2*E1J**2
     -     - 2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ + 
     -     Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Epl - ppl))/Mz**2)/deltaZ)

      amp(-1, 1) = amp(-1, 1) + (epsilon*(
     -     ((0,0.5)*CIBZ*(-1 + CW)*(-2*E1J**2 
     -     - 2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ + 
     -     Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Epl + ppl))/Mz**2)/deltaZ)

      amp(-1, 3) = amp(-1, 3) + (epsilon*(
     -     ((0,0.5)*CIBZ*(1 + CW)*(-2*E1J**2
     -     - 2*E1J*E2J - 2*E1J*EvJ - 2*E2J*EvJ + 
     -     Mf**2 - Mff**2 + Mz**2 - MB**2)*
     -     (-2*E1J**2 - 2*E1J*E2J + Mf**2 - Mff**2 + MB**2)*
     -     (Epl - ppl))/Mz**2)/deltaZ)

      return 
      end
