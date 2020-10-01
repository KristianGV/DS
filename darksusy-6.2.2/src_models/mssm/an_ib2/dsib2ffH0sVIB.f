*******************************************************************
*** subroutine dsib2ffHisVIB add to the common block array amp(i,j), i = 0, 
*** the helicity amplitudes of the s channel VIB diagrams. 
*** 
*** Author: Francesca Calore, 2014-02-20
*******************************************************************
      subroutine dsib2ffH0sVIB(f, ih) ! ih refers to the Hi radiated off
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, ih, i
 
      real*8 Mkz, Mh3, Mkw

      complex*16 deltaH3H3, deltaH3Z, deltaZH3, deltaZZ

      complex*16 amp00H3H3, amp00H3Z, amp00ZH3, amp00ZZ,
     -     amp01H3Z, amp01ZZ, amp02H3Z, amp02ZZ, amp03H3Z, amp03ZZ
           
      integer kh(2)
      real*8 Mhk(2)
    
      kh(1) = kh2 ! hi = 1 SM Higgs
      kh(2) = kh1 ! hi = 2 heavy Higgs

c...  masses/mx mediators
      Mkz = mass(kz)/mx
      Mh3 = mass(kh3)/mx
      do i= 1,2
         Mhk(i) = mass(kh(i))/mx
      end do
      Mkw = mass(kw)/mx

c... denominators
c... H3 H3
      deltaH3H3 = (dcmplx(4*E1J**2 - Mh3**2,Mh3*width(kh3)/mx)*
     - dcmplx(4*E1J*(E1J + EvJ) - Mh3**2 + mB**2,Mh3*width(kh3)/mx))
c... H3 Z
      deltaH3Z = (dcmplx(4*E1J**2 - Mkz**2,Mkz*width(kz)/mx)*
     - dcmplx(4*E1J*(E1J + EvJ) - Mh3**2 + MB**2,Mh3*width(kh3)/mx))
c... Z H3
      deltaZH3 = (dcmplx(4*E1J**2 - Mh3**2,Mh3*width(kh3)/mx)*
     - dcmplx(4*E1J*(E1J + EvJ) - Mkz**2 + MB**2,Mkz*width(kz)/mx))
c... Z Z
      deltaZZ = (dcmplx(4*E1J**2 - Mkz**2,Mkz*width(kz)/mx)*
     - dcmplx(4*E1J*(E1J + EvJ) - Mkz**2 + MB**2,Mkz*width(kz)/mx))

c...  helicity amplitudes
      amp00H3H3  =(0,4)*E1J*Mkw*gl(kh(ih),kh3,kh3)*gr(kh3,f,f)*
     - (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))/deltaH3H3

      amp00H3Z = ((0,1)*Mf*(-4*E1J**2 + Mkz**2)*(-4 + MB**2)*gl(kz,kh3,kh(ih))*
     - (gl(kz,f,f) - gr(kz,f,f))*
     - (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1))))/(E1J*Mkz**2)/deltaH3Z

      amp00ZH3 = ((0,2)*E1J*(-4*E1J**2 + MB**2)*(4*E1J*(E1J + EvJ) - Mkz**2 +
     - MB**2)*
     - gl(kz,kh3,kh(ih))*(Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*
     - gr(kh3,f,f))/Mkz**2/deltaZH3

      amp00ZZ = ((0,1)*(2*E1J + EvJ)*Mf*Mkw*(-4*E1J**2 + Mkz**2)*
     - (-4*E1J*(E1J + EvJ) + Mkz**2 - MB**2)*
     - (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*gl(kh(ih),kz,kz)*
     - (gl(kz,f,f) - gr(kz,f,f)))/Mkz**4/deltaZZ
      
      amp01H3Z = 2*Sqrt(2.)*kv*SW*gl(kz,kh3,kh(ih))*
     - ((E1J + kJ)*gl(kz,f,f) + (E1J - kJ)*gr(kz,f,f))*
     - (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))/deltaH3Z
      
      amp02H3Z = (0,-4)*CW*kv*Mf*gl(kz,kh3,kh(ih))*(gl(kz,f,f) + gr(kz,f,f))*
     - (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))/deltaH3Z
      
      amp03H3Z = -2*Sqrt(2.)*kv*SW*gl(kz,kh3,kh(ih))*
     - ((-E1J + kJ)*gl(kz,f,f) - (E1J + kJ)*gr(kz,f,f))*
     - (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))/deltaH3Z

      amp01ZZ = -((kv*Mkw*(-4*E1J*(E1J + EvJ) + Mkz**2 - MB**2)*SW*
     - (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*gl(kh(ih),kz,kz)*
     - ((E1J + kJ)*gl(kz,f,f) + (E1J - kJ)*gr(kz,f,f)))/
     - (Sqrt(2.)*Mkz**2))/deltaZZ
      
      amp02ZZ = ((0,1)*CW*kv*Mf*Mkw*(-4*E1J*(E1J + EvJ) + Mkz**2 - MB**2)*
     - (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*gl(kh(ih),kz,kz)*
     - (gl(kz,f,f) + gr(kz,f,f)))/Mkz**2/deltaZZ
      
      amp03ZZ = -((kv*Mkw*(-4*E1J*(E1J + EvJ) + Mkz**2 - MB**2)*SW*
     - (Conjg(gl(kz,kn(1),kn(1))) + gl(kz,kn(1),kn(1)))*gl(kh(ih),kz,kz)*
     - ((E1J - kJ)*gl(kz,f,f) + (E1J + kJ)*gr(kz,f,f)))/
     - (Sqrt(2.)*Mkz**2))/deltaZZ

      amp(0, 0) = amp(0, 0) + amp00H3H3 + amp00H3Z + amp00ZH3 + amp00ZZ
      amp(0, 1) = amp(0, 1) + amp01H3Z + amp01ZZ
      amp(0, 2) = amp(0, 2) + amp02H3Z + amp02ZZ
      amp(0, 3) = amp(0, 3) + amp03H3Z + amp03ZZ

      return
      end


