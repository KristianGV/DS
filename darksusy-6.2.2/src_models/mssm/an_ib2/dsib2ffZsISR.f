*******************************************************************
*** subroutine dsib2ffZsISR add to the common block array amp(i,j) 
*** the helicity amplitudes of the s channel ISR diagrams. 
*** Z emitted from the initial legs of the  s channel diagrams.
*** H_{i}, i = 1, 2 and the sum is performed over the 2  Higges (ih index) already at the amplitude level.
*** The exchanged particle index ineu stays for neutralino.
***
*** Notice:  ISR1 = ISR2 ->  Each amplitudes is multiplied by 2
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffZsISR(f, ineu) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, ineu

      integer i, j, ih
      integer kh(2)
      real*8 Mhk(2), Mkn(4), Mh3 
      
      complex*16 deltaH3, deltaHk(2), deltaZ
      complex*16 CH31, CH32, CH33, CH34,
     &     CHk1(2), CHk2(2), CHk3(2), CHk4(2), 
     &     CZ1, CZ2, CZ3, CZ4, CZ5, CZ6
      complex*16 ampH300, ampZ00, ampHk02, ampZ02
      
      kh(1)=kh1
      kh(2)=kh2

c...  masses/mx mediators
      Mh3 = mass(kh3)/mx
      do i= 1,2
         Mhk(i) = mass(kh(i))/mx
      end do
      do j = 1, 4
         Mkn(j) = mass(kn(j))/mx
      end do  

c... denominators
c... H3
      deltaH3 = (dcmplx(4*E1J**2 - Mh3**2,Mh3*width(kh3)/mx)*
     -     dcmplx(E1J*(E1J - EvJ) + MB**2/4. - Mkn(ineu)**2,
     -     Mkn(ineu)*width(kn(ineu))/mx*0d0))
     
c... Hk (i = 1, 2)
      do i = 1, 2
         deltaHk(i) = (dcmplx(4*E1J**2 
     -        - Mhk(i)**2,Mhk(i)*width(kh(i))/mx)*
     -        dcmplx(E1J*(E1J - EvJ) + MB**2/4. - Mkn(ineu)**2,
     -        Mkn(ineu)*width(kn(ineu))/mx*0d0))
      end do 
c... Z
      deltaZ = (dcmplx(4*E1J**2 - MB**2,MB*width(kz)/mx)*
     -     dcmplx(E1J*(E1J - EvJ) + MB**2/4. - Mkn(ineu)**2,
     -     Mkn(ineu)*width(kn(ineu))/mx*0d0))

c... couplings
c... H3
      CH31 = gl(kz,kn(1),kn(ineu))*
     -     gr(kh3,f,f)*gr(kh3,kn(ineu),kn(1))
      CH32 = Conjg(gr(kh3,kn(ineu),kn(1)))*
     -     gl(kz,kn(1),kn(ineu))*gr(kh3,f,f)
      CH33 = Conjg(gl(kz,kn(1),kn(ineu)))*
     -     Conjg(gr(kh3,kn(ineu),kn(1)))*gr(kh3,f,f)
      CH34 = Conjg(gl(kz,kn(1),kn(ineu)))*
     -     gr(kh3,f,f)*gr(kh3,kn(ineu),kn(1))
c... Hk
      do i = 1,2
         CHk1(i) = gl(kz,kn(1),kn(ineu))*gr(kh(i),f,f)*
     -        gr(kh(i),kn(ineu),kn(1))
         CHk2(i) = Conjg(gr(kh(i),kn(ineu),kn(1)))*
     -        gl(kz,kn(1),kn(ineu))*gr(kh(i),f,f)
         CHk3(i) =Conjg(gl(kz,kn(1),kn(ineu)))*
     -        Conjg(gr(kh(i),kn(ineu),kn(1)))*gr(kh(i),f,f)
         CHk4(i) = Conjg(gl(kz,kn(1),kn(ineu)))*gr(kh(i),f,f)*
     -        gr(kh(i),kn(ineu),kn(1))
      end do
c... Z
      CZ1 = gl(kz,f,f)*gl(kz,kn(1),kn(ineu))**2
      CZ2  = Conjg(gl(kz,kn(1),kn(ineu)))**2*gl(kz,f,f)
      CZ3 = Conjg(gl(kz,kn(1),kn(ineu)))*
     -     gl(kz,f,f)*gl(kz,kn(1),kn(ineu))
      CZ4 = gl(kz,kn(1),kn(ineu))**2*gr(kz,f,f)
      CZ5 = Conjg(gl(kz,kn(1),kn(ineu)))**2*gr(kz,f,f)
      CZ6 = Conjg(gl(kz,kn(1),kn(ineu)))*
     -     gl(kz,kn(1),kn(ineu))*gr(kz,f,f)

c... helicity amplitudes (factor of 2* because of ISR1 = ISR2) 
c... firstly amp00, amp02 with Higgses contributions
c... H3 only amp00 contribution
      ampH300 = 0.d0
      ampH300 = + 2*(
     -     ((0,4)*E1J**2*kv*(CH31 + CH33 - (CH32 + CH34)*Mkn(ineu))
     -     )/MB)/deltaH3

c... Hk (sum over ih) only amp02 contribution
      ampHk02 = 0.d0
      do ih = 1,2
         ampHk02 = ampHk02 + 2*(
     -        ((0,-4)*E1J*kJ*kv*(CHk1(ih) + CHk3(ih)
     -        - (CHk2(ih) + CHk4(ih))*Mkn(ineu)))/MB)/deltaHk(ih)
      end do 

c... Z
      ampZ00 =0.d0
      ampZ02 = 0.d0

      ampZ00 = + 2*(
     -     ((0,2)*(CZ1 - CZ2 - CZ4 + CZ5)*
     -     kv*Mf*(-4*E1J**2 + MB**2)*Mkn(ineu))/
     -     MB**3)/deltaZ

      ampZ02 = + 2*(
     -     ((0,2)*(CZ1 - CZ2 + CZ4 - CZ5)*CW*
     -     EvJ*Mf*Mkn(ineu))/MB)/deltaZ
      
      amp(0, 0) = amp(0, 0) + ampH300 + ampZ00
      amp(0, 2) = amp(0, 2) + ampHk02 + ampZ02

      amp(0, 1) = amp(0, 1) + 2*(
     -     (Sqrt(2.)*EvJ*((CZ4 - CZ5)*(-E1J + kJ) -CZ1*(E1J + kJ)+ 
     -     CZ2*(E1J + kJ))*SW*Mkn(ineu))/MB)/deltaZ

      amp(0, 3) = amp(0, 3) + 2*(
     -     (Sqrt(2.)*EvJ*((-CZ1 + CZ2 - CZ4 + CZ5)*E1J + 
     -     (CZ1 - CZ2 - CZ4 + CZ5)*kJ)*SW*Mkn(ineu))/MB)/deltaZ

      amp(1, 1) = amp(1, 1) + 2*(
     -     (0,-1)*(1 + CW)*(-2*E1J*(CZ6*(E1J - kJ) +  
     -     CZ3*(E1J + kJ))*kv +((CZ1 - CZ2 + CZ4 - CZ5)*E1J + 
     -     (CZ1 - CZ2 - CZ4 + CZ5)*kJ)*Mkn(ineu)))/deltaZ

      amp(1, 2) = amp(1, 2) + 2*(
     -     -(Sqrt(2.)*Mf*SW*(2*(CZ3 + CZ6)*E1J*kv + 
     -     (-CZ1 + CZ2 - CZ4 + CZ5)*Mkn(ineu))))/deltaZ

      amp(1, 3) = amp(1, 3) + 2*(
     -     (0,-1)*(-1 + CW)*(-2*E1J*(CZ3*(E1J - kJ) + 
     -     CZ6*(E1J + kJ))*kv +((CZ1 - CZ2 + CZ4 - CZ5)*E1J + 
     -     (-CZ1 + CZ2 + CZ4 - CZ5)*kJ)*Mkn(ineu)))/deltaZ

      amp(-1, 1) = amp(-1, 1) + 2*(
     -     (0,1)*(-1 + CW)*(-2*E1J*(CZ6*(E1J - kJ) 
     -     + CZ3*(E1J + kJ))*kv + 
     -     ((CZ2 - CZ4 + CZ5)*E1J + (CZ2 + CZ4 - CZ5)*kJ - 
     -     CZ1*(E1J + kJ))*Mkn(ineu)))/deltaZ
       
      amp(-1, 2) = amp(-1, 2) + 2*(
     -     Sqrt(2.)*Mf*SW*(2*(CZ3 + CZ6)*E1J*kv + 
     -     (CZ1 - CZ2 + CZ4 - CZ5)*Mkn(ineu)))/deltaZ
      
      amp(-1, 3) = amp(-1, 3) + 2*(
     -     (0,1)*(1 + CW)*(-2*E1J*(CZ3*(E1J - kJ) +  
     -     CZ6*(E1J + kJ))*kv+((-CZ1 + CZ2 - CZ4 + CZ5)*E1J + 
     -     (CZ1 - CZ2 - CZ4 + CZ5)*kJ)*Mkn(ineu)))/deltaZ

      return
      end
