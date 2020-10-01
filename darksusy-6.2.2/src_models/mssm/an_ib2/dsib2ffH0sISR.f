*******************************************************************
*** subroutine dsib2ffHisISR add to the common block array amp(i,j), i = 0, 
*** the helicity amplitudes of the s channel ISR diagrams. 
*** Hi emitted from the initial legs of the  s channel diagrams.
*** The exchanged particle index ineu stays for neutralino.
***
*** Notice:  ISR1 = ISR2 ->  Each amplitudes is multiplied by 2
*** 
*** Author: Francesca Calore, 2014-02-20
*******************************************************************
      subroutine dsib2ffH0sISR(f, ineu, hi)  ! hi refers to the Hi radiated off
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, ineu, hi

      integer i, j, ih
      integer kh(2)
      real*8 Mhk(2), Mkn(4), Mh3, Mkz
      
      complex*16 deltaH3, deltaHk(2), deltaZ
      complex*16 CH31, CH32, CH33, CH34,
     &     CHk1(2), CHk2(2), CHk3(2), CHk4(2), 
     &     CZ1, CZ2, CZ3, CZ4, CZ5, CZ6, CZ7, CZ8
      complex*16 ampH300, ampZ00, ampHk02, ampZ02

      kh(1) = kh2 ! hi = 1 SM Higgs
      kh(2) = kh1 ! hi = 2 heavy Higgs

c...  masses/mx mediators
      Mh3 = mass(kh3)/mx
      do i= 1,2
         Mhk(i) = mass(kh(i))/mx
      end do
      do j = 1, 4
         Mkn(j) = mass(kn(j))/mx
      end do  
      Mkz = mass(kz)/mx

c... denominators
c... H3

c...  AG: Added missing /mx to widtgh of H3 propagator
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
      deltaZ = (dcmplx(4*E1J**2 - Mkz**2,Mkz*width(kz)/mx)*
     -     dcmplx(E1J*(E1J - EvJ) + MB**2/4. - Mkn(ineu)**2,
     -     Mkn(ineu)*width(kn(ineu))/mx*0d0))


c... couplings
c...  H3
      CH31 = Conjg(gr(kh3,kn(ineu),kn(1)))*
     - Conjg(gr(kh(hi),kn(ineu),kn(1)))*gr(kh3,f,f)      
      CH32 = Conjg(gr(kh(hi),kn(ineu),kn(1)))*gr(kh3,f,f)*
     - gr(kh3,kn(ineu),kn(1))
      CH33 = Conjg(gr(kh3,kn(ineu),kn(1)))*gr(kh3,f,f)*
     - gr(kh(hi),kn(ineu),kn(1))
      CH34 = gr(kh3,f,f)*gr(kh3,kn(ineu),kn(1))*gr(kh(hi),kn(ineu),kn(1))

c...  Hk
      do i = 1,2
         CHk1(i) =Conjg(gr(kh(hi),kn(ineu),kn(1)))*
     -    Conjg(gr(kh(i),kn(ineu),kn(1)))*gr(kh(i),f,f)
         CHk2(i) =Conjg(gr(kh(hi),kn(ineu),kn(1)))*gr(kh(i),f,f)*
     -    gr(kh(i),kn(ineu),kn(1))
         CHk3(i) =Conjg(gr(kh(i),kn(ineu),kn(1)))*gr(kh(hi),kn(ineu),kn(1))*
     -    gr(kh(i),f,f)
         CHk4(i) = gr(kh(hi),kn(ineu),kn(1))*gr(kh(i),f,f)*
     -    gr(kh(i),kn(ineu),kn(1))
      end do

c...  Z
      CZ1 = Conjg(gr(kh(hi),kn(ineu),kn(1)))*gl(kz,f,f)*
     - gl(kz,kn(1),kn(ineu))
      CZ2  = gl(kz,f,f)*gl(kz,kn(1),kn(ineu))*gr(kh(hi),kn(ineu),kn(1))
      CZ3 =Conjg(gl(kz,kn(1),kn(ineu)))*Conjg(gr(kh(hi),kn(ineu),kn(1)))*
     - gl(kz,f,f)
      CZ4 = Conjg(gl(kz,kn(1),kn(ineu)))*gl(kz,f,f)*
     - gr(kh(hi),kn(ineu),kn(1))
      CZ5 = Conjg(gr(kh(hi),kn(ineu),kn(1)))*gl(kz,kn(1),kn(ineu))*
     - gr(kz,f,f)
      CZ6 = gl(kz,kn(1),kn(ineu))*gr(kz,f,f)*gr(kh(hi),kn(ineu),kn(1))
      CZ7 = Conjg(gl(kz,kn(1),kn(ineu)))*Conjg(gr(kh(hi),kn(ineu),kn(1)))*
     - gr(kz,f,f)
      CZ8 = Conjg(gl(kz,kn(1),kn(ineu)))*gr(kz,f,f)*
     - gr(kh(hi),kn(ineu),kn(1))

c...  helicity amplitudes (factor of 2* because of ISR1 = ISR2) 
c...  firstly amp00, amp02 with Higgses contributions
c...  H3 only amp00 contribution
      ampH300 = 0.d0
      ampH300 = + 2*((0,1)*E1J*((CH32 - CH33)*(-4*E1J**2 + MB**2)
     -     + 4*(CH31 - CH34)*Mkn(ineu))/deltaH3)

c...  Hk (sum over ih) only amp02 contribution
      ampHk02 = 0.d0
      do ih = 1,2
         ampHk02 = ampHk02 + 2*(
     -    (0,-1)*kJ*((CHk2(ih) - CHk3(ih))*(-4*E1J**2 + MB**2)
     -    + 4*(CHk1(ih) - CHk4(ih))*Mkn(ineu))/deltaHk(ih))
      end do 

c... Z
      ampZ00 =0.d0
      ampZ02 = 0.d0

      ampZ00 = + 2*(((0,-0.25)*Mf*(-4*E1J**2 + Mkz**2)*
     - ((CZ2 + CZ3 - CZ6 - CZ7)*(-4 + 12*E1J**2 + MB**2) - 
     - (CZ1 + CZ4 - CZ5 - CZ8)*(-4*(1 + E1J**2) + MB**2)*Mkn(ineu)))/
     - (E1J*Mkz**2)/deltaZ)

      ampZ02 = + 2*((0,1)*CW*kv*Mf*(CZ2 + CZ3 + CZ6 + CZ7 - 
     - (CZ1 + CZ4 + CZ5 + CZ8)*Mkn(ineu))/deltaZ)
      
      amp(0, 0) = amp(0, 0) + ampH300 + ampZ00
      amp(0, 2) = amp(0, 2) + ampHk02 + ampZ02

      amp(0 , 1) = amp(0, 1) + 2*((kv*SW*(-((CZ2 + CZ3 + CZ6 + CZ7)*E1J) + 
     - (-CZ2 - CZ3 + CZ6 + CZ7)*kJ + 
     - ((CZ1 + CZ4 + CZ5 + CZ8)*E1J + 
     - (CZ1 + CZ4 - CZ5 - CZ8)*kJ)*Mkn(ineu)))/Sqrt(2.)/deltaZ)

      amp(0 , 3) = amp(0, 3) + 2*((kv*SW*(-((CZ2 + CZ3 + CZ6 + CZ7)*E1J) + 
     - (CZ2 + CZ3 - CZ6 - CZ7)*kJ + 
     - ((CZ1 + CZ4 + CZ5 + CZ8)*E1J + 
     - (-CZ1 - CZ4 + CZ5 + CZ8)*kJ)*Mkn(ineu)))/Sqrt(2.)/deltaZ)

      return
      end
