*******************************************************************
*** subroutine dsib2ffHPsVIB add to the common block array amp(i,j), i = 0, 
*** the helicity amplitudes of the s channel VIB diagrams. 
*** The HP is emitted from the mediators (H3, Z).
*** 
*** Author: Francesca Calore, 2014-02-22
***         Francesca Calore, 2015-02-15 (bug solved cf with CH)
*******************************************************************
      subroutine dsib2ffHPsVIB(f) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f
      integer ff

      real*8 Mkhc, Mh3, Mkw, Mkz
      complex*16 deltaH3, deltaZ
      complex*16 ampH300, ampH301, ampH302, ampH303, 
     -     ampZ00, ampZ01, ampZ02, ampZ03
      complex*16 CZ1, CZ2
      complex*16 CH00, CH01, CH02, CH03
      real*8 epsilon

      ff = c_fbartype

      if(mod(c_ftype,2).EQ.0) then
         epsilon = -1d0 ! epsilon sign defines the sign of the s VIB H3 diagram under CP transformation
      else
         epsilon = +1d0
      endif 

c..  masses/mx of exchanged particles
      Mkhc = mass(khc)/mx
      Mh3 = mass(kh3)/mx
      Mkw = mass(kw)/mx
      Mkz = mass(kz)/mx

c... denominators
c... H3
      deltaH3 = (dcmplx(2*E1J*E2J + 2*kJ**2
     -     + Mf**2 + Mff**2 - Mkw**2,Mkw*width(kw)/mx)*
     -     dcmplx(2*E1J*E2J + 2*E1J*EvJ + 2*E2J*EvJ
     -     + 2*kJ**2 + MB**2 + Mf**2 + 
     -     Mff**2 - Mh3**2,Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(2*E1J*E2J + 2*kJ**2 + Mf**2 
     -     + Mff**2 - Mkhc**2,Mkhc*width(khc)/mx)*
     -     dcmplx(2*E1J*E2J + 2*E1J*EvJ + 2*E2J*EvJ
     -     + 2*kJ**2 + MB**2 + Mf**2 + 
     -     Mff**2 - Mkz**2,Mkz*width(kz)/mx))

c... couplings
c... H3
      CH00 =  gl(kw,khc,kh3)*gl(kw,f,ff)*
     -      (Conjg(gr(kh3,kn(1),kn(1))) - 
     -      gr(kh3,kn(1),kn(1)))
      CH01 = gl(kw,khc,kh3)*gl(kw,f,ff)*
     -     (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))
      CH02 = gl(kw,khc,kh3)*
     -     gl(kw,f,ff)*(Conjg(gr(kh3,kn(1),kn(1))) - 
     -     gr(kh3,kn(1),kn(1)))
      CH03 = gl(kw,khc,kh3)*gl(kw,f,ff)*
     -     (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))

c... Z
      CZ1 =gl(khc,f,ff)*gl(kz,khc,khc)*
     - (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))
      CZ2 = gl(kz,khc,khc)*
     - (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*
     - gr(khc,f,ff)

c... helicity amplitudes
c... H3
      ampH300 = + ((0,2)*(-(Emi*(2*E1J**2*E2J
     -      + 2*E1J*E2J**2 + 2*E1J**2*EvJ + 
     -      4*E1J*E2J*EvJ + 2*E2J**2*EvJ + 2*E1J*kJ**2 + 2*E2J*kJ**2 + 
     -      (E1J + E2J)*Mf**2 + (E1J + E2J)*Mff**2 - E1J*Mkw**2 - 
     -      E2J*Mkw**2 - 2*EvJ*Mkw**2)) + 2*CW*kv*Mkw**2*pmi)*CH00
     -     )/Mkw**2/deltaH3
      
      ampH301 = + ((-2*Sqrt(2.)*kv*(Epl + ppl)*SW*
     -     CH01))/deltaH3
      
      ampH302 = (((0,2)*(2*CW*Emi*kv*Mkw**2 - 
     -     (2*E1J**2*E2J + 2*E1J*E2J**2 + 
     -    2*E1J**2*EvJ + 4*E1J*E2J*EvJ + 2*E2J**2*EvJ + 2*E1J*kJ**2 + 
     -    2*E2J*kJ**2 + (E1J + E2J)*Mf**2 + (E1J + E2J)*Mff**2 - 
     -    E1J*Mkw**2 - E2J*Mkw**2 - 2*EvJ*Mkw**2)*pmi)*CH02
     -     )/Mkw**2)/deltaH3

      ampH303 = (-2*Sqrt(2.)*kv*(Epl - ppl)*SW*
     -     CH03)/deltaH3

c... Z 
      ampZ00 = + (((0,-1)*(CZ1 - CZ2)*Epl*
     -     (2*E1J*E2J + 2*kJ**2 - MB**2 + Mf**2 + Mff**2)*
     -     (2*(E1J*E2J + (E1J + E2J)*EvJ + kJ**2) 
     -     + MB**2 + Mf**2 + Mff**2 - 
     -     Mkz**2))/Mkz**2)/deltaZ

      ampZ01 = 0.d0

      ampZ02 = + (((0,-1)*(CZ1 + CZ2)*(2*E1J*E2J
     -     + 2*kJ**2 - MB**2 + Mf**2 + Mff**2)*
     -     (2*(E1J*E2J + (E1J + E2J)*EvJ + kJ**2) 
     -     + MB**2 + Mf**2 + Mff**2 - 
     -     Mkz**2)*ppl)/Mkz**2)/deltaZ

      ampZ03 = 0d0

      amp(0, 0) = amp(0, 0) + (epsilon*ampH300) + ( ampZ00)
      amp(0, 1) = amp(0, 1) + (epsilon*ampH301) + ( ampZ01)
      amp(0, 2) = amp(0, 2) + (epsilon*ampH302) + ( ampZ02)
      amp(0, 3) = amp(0, 3) + (epsilon*ampH303) + ( ampZ03)

      return 
      end
