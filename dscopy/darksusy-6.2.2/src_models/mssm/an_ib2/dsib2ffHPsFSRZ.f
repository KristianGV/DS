*******************************************************************
*** subroutine dsib2ffHPsFSR add to the common block array amp(i,j), i = 0,
*** the helicity amplitudes of the s channel FSR diagrams. 
*** The Hp is emitted from the final legs of the s channel diagrams. 
*** 
*** Author: Francesca Calore, 2014-02-22
***         Francesca Calore, 2015-02-15 (bug solved cf with CH)
*******************************************************************
c...  Z mediated diagrams
      subroutine dsib2ffHPsFSRZ(f, leg) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg
      integer ff

      real*8 Mkz
      complex*16 CZ1, CZ2, CZ3, CZ4
      complex*16 deltaZ

      ff = c_fbartype

c... masses/mx exchanged particles
      Mkz = mass(kz)/mx

c... Z FSR1
      if (leg.eq.1) then

         CZ1 = gl(khc,f,ff)*gl(kz,ff,ff)*gl(kz,kn(1),kn(1)) ! RealPart(gl(kz,kn(1),kn(1)) == gl(kz,kn(1),kn(1))
         CZ2 = gl(khc,f,ff)*gr(kz,ff,ff)*gl(kz,kn(1),kn(1))
         CZ3 = gl(kz,ff,ff)*gr(khc,f,ff)*gl(kz,kn(1),kn(1))
         CZ4 = gr(khc,f,ff)*gr(kz,ff,ff)*gl(kz,kn(1),kn(1))

         deltaZ = dcmplx(2*E1J*EvJ - 2*CW*kJ*kv
     -    + MB**2 + Mf**2 - Mff**2,Mff*width(ff)/mx)*
     -    dcmplx(2*E1J*E2J + 2*E1J*EvJ + 2*E2J*EvJ
     -    + 2*kJ**2 + MB**2 + Mf**2 + 
     -    Mff**2 - Mkz**2,Mkz*width(kz)/mx)

         amp(0, 0) = amp(0, 0) + (((0,-2)*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    (-((CZ2 - CZ3)*Mff*(Emi*(E1J + E2J + EvJ) + CW*kv*pmi)) + 
     -    CZ1*(Epl*((E1J + E2J)*EvJ + MB**2)
     -    + Emi*(E1J + E2J + EvJ)*Mf + 
     -    CW*kv*(Mf*pmi - (E1J + E2J)*ppl)) - 
     -    CZ4*(Epl*((E1J + E2J)*EvJ + MB**2)
     -    + Emi*(E1J + E2J + EvJ)*Mf + 
     -    CW*kv*(Mf*pmi - (E1J + E2J)*ppl))))/Mkz**2)/deltaZ

         amp(0, 1) = amp(0, 1) + ((Sqrt(2.)*kv*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    (CZ4*(Epl*Mf - (E1J + E2J)*(Emi - pmi) - Mf*ppl) + 
     -    Mff*(CZ2*(Epl - ppl) + CZ3*(Epl + ppl)) + 
     -    CZ1*(-((E1J + E2J)*(Emi + pmi))
     -    + Mf*(Epl + ppl)))*SW)/Mkz**2)/deltaZ

         amp(0, 2) = amp(0, 2) + (((0,-2)*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    ((CZ2 + CZ3)*Mff*(CW*Emi*kv + (E1J + E2J + EvJ)*pmi) + 
     -    CZ1*(CW*kv*(-((E1J + E2J)*Epl) + Emi*Mf) + 
     -    (E1J + E2J + EvJ)*Mf*pmi + ((E1J + E2J)*EvJ + MB**2)*ppl) + 
     -    CZ4*(CW*kv*(-((E1J + E2J)*Epl) + Emi*Mf) + 
     -    (E1J + E2J + EvJ)*Mf*pmi + 
     -    ((E1J + E2J)*EvJ + MB**2)*ppl)))/Mkz**2)/deltaZ

         amp(0, 3) = amp(0, 3) + ((Sqrt(2.)*kv*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    ((CZ2 + CZ3)*Epl*Mff + (CZ2 - CZ3)*Mff*ppl + 
     -    CZ1*(Epl*Mf - (E1J + E2J)*(Emi - pmi) - Mf*ppl) + 
     -    CZ4*(-((E1J + E2J)*(Emi + pmi))
     -    + Mf*(Epl + ppl)))*SW)/Mkz**2)/deltaZ

c... Z FSR2
      else if (leg.eq.2) then

         CZ1 =gl(khc,f,ff)*gl(kz,f,f)*RealPart(gl(kz,kn(1),kn(1)))
         CZ2 =gl(khc,f,ff)*gr(kz,f,f)*RealPart(gl(kz,kn(1),kn(1)))
         CZ3 =gl(kz,f,f)*gr(khc,f,ff)*RealPart(gl(kz,kn(1),kn(1)))
         CZ4 =gr(khc,f,ff)*gr(kz,f,f)*RealPart(gl(kz,kn(1),kn(1)))
         
         deltaZ =dcmplx(2*E2J*EvJ + 2*CW*kJ*kv 
     -    + MB**2 - Mf**2 + Mff**2,Mf*width(f)/mx)*
     -    dcmplx(2*E1J*E2J + 2*E1J*EvJ 
     -    + 2*E2J*EvJ + 2*kJ**2 + MB**2 + Mf**2 + 
     -    Mff**2 - Mkz**2,Mkz*width(kz)/mx)

         
         amp(0, 0) = amp(0, 0) + (((0,-2)*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    ((CZ1 - CZ4)*Mf*(Emi*(E1J + E2J + EvJ) + CW*kv*pmi) - 
     -    CZ2*(Epl*((E1J + E2J)*EvJ + MB**2) 
     -    + Emi*(E1J + E2J + EvJ)*Mff + 
     -    CW*kv*(Mff*pmi + (E1J + E2J)*ppl)) + 
     -    CZ3*(Epl*((E1J + E2J)*EvJ + MB**2) 
     -    + Emi*(E1J + E2J + EvJ)*Mff + 
     -    CW*kv*(Mff*pmi + (E1J + E2J)*ppl))))/Mkz**2)/deltaZ

         amp(0, 1) = amp(0, 1) + ((Sqrt(2.)*kv*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    (-(CZ2*E1J*Emi) - CZ3*E1J*Emi
     -    + CZ2*Epl*Mff + CZ3*Epl*Mff - 
     -    CZ2*E1J*pmi + CZ3*E1J*pmi - 
     -    E2J*((CZ2 + CZ3)*Emi + (CZ2 - CZ3)*pmi) + 
     -    CZ4*Mf*(Epl - ppl) - CZ2*Mff*ppl + CZ3*Mff*ppl + 
     -    CZ1*Mf*(Epl + ppl))*SW)/Mkz**2)/deltaZ

         amp(0, 2) = amp(0, 2) + (((0,-2)*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    ((CZ1 + CZ4)*Mf*(CW*Emi*kv + (E1J + E2J + EvJ)*pmi) - 
     -    CZ2*(CW*kv*((E1J + E2J)*Epl - Emi*Mff) - 
     -    (E1J + E2J + EvJ)*Mff*pmi + ((E1J + E2J)*EvJ + MB**2)*ppl) - 
     -    CZ3*(CW*kv*((E1J + E2J)*Epl - Emi*Mff) - 
     -    (E1J + E2J + EvJ)*Mff*pmi + ((E1J + E2J)*EvJ
     -    + MB**2)*ppl)))/Mkz**2)/deltaZ
         
         amp(0, 3) = amp(0, 3) + ((Sqrt(2.)*kv*(-2*(E1J + E2J)*
     -    (E1J + EvJ) - MB**2 + Mf**2 - Mff**2 + Mkz**2)*
     -    (-(CZ2*E1J*Emi) - CZ3*E1J*Emi 
     -    + CZ2*Epl*Mff + CZ3*Epl*Mff + 
     -    CZ2*E1J*pmi - CZ3*E1J*pmi - 
     -    E2J*(CZ2*(Emi - pmi) + CZ3*
     -    (Emi + pmi)) + CZ1*Mf*(Epl - ppl) + 
     -    CZ2*Mff*ppl - CZ3*Mff*ppl + CZ4*Mf*
     -    (Epl + ppl))*SW)/Mkz**2)/deltaZ      
      
      endif

      return
      end
