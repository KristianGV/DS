*******************************************************************
*** subroutine dsib2ffZsFSR add to the common block array amp(i,j) 
*** the helicity amplitudes of the s channel FSR diagrams. 
*** The Z is emitted from the final legs of the s channel diagrams.
*** H3 and Z diagrams are considered separately and then
*** added to the corresponding entries of amp(-1:1, 0:3). 
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffZsFSR(f, leg) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg

      real*8 q1
      complex*16 deltaH3, deltaZ
      complex*16 CH1, CH2, CZ1, CZ2, CZ3

      complex*16 ampH00, ampH01, ampH02, ampH03, 
     &     ampH10, ampH11, ampH12, ampH13, ampH20, 
     &     ampH21, ampH22, ampH23
      complex*16 ampZ00, ampZ01, ampZ02, ampZ03, 
     &     ampZ10, ampZ11, ampZ12, ampZ13, ampZ20,
     &     ampZ21, ampZ22, ampZ23

      real*8 Mh3
       
c... Defines FSR1 and FSR2 amplitudes through their simmetry relations via q1 exchange
c      if (leg.eq.1) then
         q1 = kJ
c      else 
      if (leg.eq.2) then
         q1 = -kJ
      endif

c..  masses/mx definition
      Mh3 = mass(kh3)/mx 

c... denominators
c... Higgs
      deltaH3 = (dcmplx(2*E1J*EvJ + MB**2 - 2*CW*kv*q1,Mf*width(f)/mx)*
     -     dcmplx(4*E1J*(E1J + EvJ) - Mh3**2 + MB**2,Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(2*E1J*EvJ + MB**2 - 2*CW*kv*q1,Mf*width(f)/mx)*
     -     dcmplx(4*E1J*(E1J + EvJ),MB*width(kz)/mx))


c... couplings
c... Higgs
      CH1 = gl(kz,f,f)*gr(kh3,f,f)*(Conjg(gr(kh3,kn(1),kn(1))) - 
     -     gr(kh3,kn(1),kn(1)))
      CH2 = gr(kz,f,f)*gr(kh3,f,f)*(Conjg(gr(kh3,kn(1),kn(1))) - 
     -     gr(kh3,kn(1),kn(1)))
c... Z
      CZ1 = gl(kz,f,f)**2*Dble(gl(kz,kn(1),kn(1)))
      CZ2 = gr(kz,f,f)**2*Dble(gl(kz,kn(1),kn(1)))
      CZ3 = gl(kz,f,f)*gr(kz,f,f)*Dble(gl(kz,kn(1),kn(1)))


c... helicity amplitudes (entries as for FSR1)
c... Higgs
      ampH00 = (
     -     ((0,-2)*(CH1 + CH2)*(-2*E1J**2*kv
     -     + CW*(2*E1J*EvJ + MB**2)*q1))/MB)/deltaH3

      ampH01 = (
     -     (-((Sqrt(2.)*(CH1 - CH2)*Mf*(2*E1J*EvJ + MB**2)*SW)/MB)
     -     ))/deltaH3

      ampH02 = (
     -     ((0,2)*(CH1 - CH2)*E1J*(CW*(2*E1J*EvJ + MB**2) 
     -     - 2*kv*q1))/MB)/deltaH3

      ampH03 = (
     -     -((Sqrt(2.)*(CH1 - CH2)*Mf*(2*E1J*EvJ + MB**2)*SW)/MB)
     -     )/deltaH3

      ampH10 = (
     -     Sqrt(2.)*(CH1*(-2*E1J - EvJ + kv) 
     -     - CH2*(2*E1J + EvJ + kv))*q1*SW)/deltaH3

      ampH11 = (
     -     (0,1)*(1 + CW)*(CH1*(-2*E1J - EvJ + kv) 
     -     + CH2*(2*E1J + EvJ + kv))*Mf)/deltaH3

      ampH12 = (
     -     Sqrt(2.)*E1J*(CH1*(2*E1J + EvJ - kv)
     -     - CH2*(2*E1J + EvJ + kv))*SW)/deltaH3

      ampH13 = (
     -     (0,1)*(-1 + CW)*(CH1*(-2*E1J - EvJ + kv) 
     -     + CH2*(2*E1J + EvJ + kv))*Mf)/deltaH3

      ampH20 = (
     -     Sqrt(2.)*(CH2*(-2*E1J - EvJ + kv)
     -     - CH1*(2*E1J + EvJ + kv))*q1*SW)/deltaH3

      ampH21 = (
     -     (0,-1)*(-1 + CW)*(CH2*(-2*E1J - EvJ + kv) 
     -     + CH1*(2*E1J + EvJ + kv))*Mf)/deltaH3

      ampH22 = (
     -     Sqrt(2.)*E1J*(CH2*(-2*E1J - EvJ + kv) 
     -     + CH1*(2*E1J + EvJ + kv))*SW)/deltaH3

      ampH23 = (
     -     (0,-1)*(1 + CW)*(CH2*(-2*E1J - EvJ + kv) 
     -     + CH1*(2*E1J + EvJ + kv))*Mf)/deltaH3

c... Z
      ampZ00 = (((0,8)*(CZ1 - CZ2)*E1J*(E1J + EvJ)*Mf*
     -     (kv*(2*E1J*(E1J + EvJ) + MB**2) + 
     -     CW*(-2*EvJ*(E1J + EvJ) + MB**2)*q1))/MB**3)/deltaZ
      
      ampZ01 = ((-4*Sqrt(2.)*E1J*(E1J + EvJ)*
     -     (-2*CZ3*Mf**2*(2*E1J*EvJ + MB**2) + 
     -     CZ2*(Mf**2*(2*E1J*EvJ + MB**2) - 
     -     EvJ*(-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) + 
     -     CZ1*(Mf**2*(2*E1J*EvJ + MB**2) + 
     -     EvJ*(E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)))*SW)/
     -     MB**3)/deltaZ
      
      ampZ02 = (((0,8)*E1J*(E1J + EvJ)*Mf*
     -     (2*CW*(CZ1 + CZ2)*E1J*EvJ**2 + 
     -     (-2*CZ3 + CZ1 + CZ2)*E1J*(CW*MB**2 - 2*kv*q1) + 
     -     CW*EvJ*(2*(-2*CZ3 + CZ1 + CZ2)*E1J**2 + 
     -     (CZ1 + CZ2)*(MB**2 - 2*CW*kv*q1))))/MB**3)/deltaZ
      
      ampZ03 = ((-4*Sqrt(2.)*E1J*(E1J + EvJ)*
     -     (-2*CZ3*Mf**2*(2*E1J*EvJ + MB**2) + 
     -     CZ1*(Mf**2*(2*E1J*EvJ + MB**2) - 
     -     EvJ*(-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)) + 
     -     CZ2*(Mf**2*(2*E1J*EvJ + MB**2) + 
     -     EvJ*(E1J + q1)*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1)))*SW)/
     -     MB**3)/deltaZ
      
      ampZ10 = ((4*Sqrt(2.)*E1J*(E1J + EvJ)*
     -     (-2*CZ3*kv + CZ1*(-2*E1J - EvJ + kv) + 
     -     CZ2*(2*E1J + EvJ + kv))*Mf*q1*SW)/MB**2)/deltaZ
      
      ampZ11 = (((0,4)*(1 + CW)*E1J*(E1J + EvJ)*
     -     (2*CZ3*(2*E1J + EvJ)*Mf**2 - 
     -     CZ1*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 
     -     2*kv*(E1J + (-1 + CW)*q1))) + 
     -     CZ2*((-2*E1J - EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(MB**2 + 
     -     2*(E1J*(EvJ + kv) + (kv - CW*kv)*q1)))))/MB**2)/deltaZ
      
      ampZ12 = ((4*Sqrt(2.)*E1J*(E1J + EvJ)*Mf*
     -     (-2*CZ3*E1J*(2*E1J + EvJ) + 
     -     CZ1*(E1J*(2*E1J + 3*EvJ - kv) + MB**2 - 2*CW*kv*q1) + 
     -     CZ2*(E1J*(2*E1J + 3*EvJ + kv) + MB**2 - 2*CW*kv*q1))*SW)/
     -     MB**2)/deltaZ
      
      ampZ13 = (((0,4)*(-1 + CW)*E1J*(E1J + EvJ)*
     -     (2*CZ3*(2*E1J + EvJ)*Mf**2 - 
     -     CZ2*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*(EvJ + kv) + MB**2 - 2*(1 + CW)*kv*q1))
     -     + CZ1*(-((2*E1J + EvJ + kv)*Mf**2) + 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*kv*(E1J + q1 + CW*q1))
     -     )))/MB**2)/deltaZ
      
      ampZ20 = ((-4*Sqrt(2.)*E1J*(E1J + EvJ)*
     -     (-2*CZ3*kv + CZ2*(-2*E1J - EvJ + kv) + 
     -     CZ1*(2*E1J + EvJ + kv))*Mf*q1*SW)/MB**2)/deltaZ
      
      ampZ21 = (((0,4)*(-1 + CW)*E1J*(E1J + EvJ)*
     -     (2*CZ3*(2*E1J + EvJ)*Mf**2 - 
     -     CZ1*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*(EvJ + kv) + MB**2 - 2*(1 + CW)*kv*q1))
     -     + CZ2*(-((2*E1J + EvJ + kv)*Mf**2) + 
     -     (-E1J + q1)*(2*E1J*EvJ + MB**2 - 2*kv*(E1J + q1 + CW*q1))
     -     )))/MB**2)/deltaZ
         
      ampZ22 = ((4*Sqrt(2.)*E1J*(E1J + EvJ)*Mf*
     -     (-2*CZ3*E1J*(2*E1J + EvJ) + 
     -     CZ2*(E1J*(2*E1J + 3*EvJ - kv) + MB**2 - 2*CW*kv*q1) + 
     -     CZ1*(E1J*(2*E1J + 3*EvJ + kv) + MB**2 - 2*CW*kv*q1))*SW)/
     -     MB**2)/deltaZ
      
      ampZ23 = (((0,4)*(1 + CW)*E1J*(E1J + EvJ)*
     -     (2*CZ3*(2*E1J + EvJ)*Mf**2 - 
     -     CZ2*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 
     -     2*kv*(E1J + (-1 + CW)*q1))) + 
     -     CZ1*((-2*E1J - EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(MB**2 + 
     -     2*(E1J*(EvJ + kv) + (kv - CW*kv)*q1)))))/MB**2)/deltaZ

c... exploiting symmetries
      if (leg.eq.1) then
c... FSR1
         amp(0, 0) = amp(0, 0) + ampH00 + ampZ00
         amp(0, 1) = amp(0, 1) + ampH01 + ampZ01
         amp(0, 2) = amp(0, 2) + ampH02 + ampZ02
         amp(0, 3) = amp(0, 3) + ampH03 + ampZ03
         amp(1, 0) = amp(1, 0) + ampH10 + ampZ10
         amp(1, 1) = amp(1, 1) + ampH11 + ampZ11
         amp(1, 2) = amp(1, 2) + ampH12 + ampZ12
         amp(1, 3) = amp(1, 3) + ampH13 + ampZ13
         amp(-1, 0) = amp(-1, 0) + ampH20 + ampZ20
         amp(-1, 1) = amp(-1, 1) + ampH21 + ampZ21
         amp(-1, 2) = amp(-1, 2) + ampH22 + ampZ22
         amp(-1, 3) = amp(-1, 3) + ampH23 + ampZ23
      else if (leg.eq.2) then
c... FSR2
         amp(0, 0) = amp(0, 0) - ampH00 - ampZ00
         amp(0, 1) = amp(0, 1) - ampH03 - ampZ03
         amp(0, 2) = amp(0, 2) - ampH02 - ampZ02
         amp(0, 3) = amp(0, 3) - ampH01 - ampZ01
         amp(1, 0) = amp(1, 0) - ampH20 - ampZ20
         amp(1, 1) = amp(1, 1) - ampH23 - ampZ23
         amp(1, 2) = amp(1, 2) - ampH22 - ampZ22
         amp(1, 3) = amp(1, 3) - ampH21 - ampZ21
         amp(-1, 0) = amp(-1, 0) - ampH10 - ampZ10
         amp(-1, 1) = amp(-1, 1) - ampH13 - ampZ13
         amp(-1, 2) = amp(-1, 2) - ampH12 - ampZ12
         amp(-1, 3) = amp(-1, 3) - ampH11 - ampZ11
      endif

      return
      end
      

c============== TEMPORARY =========================================
*******************************************************************
***   TEMPORARY SUBROUTINE: for IB check
***   1) w/o width propagator
***   2) isf2 -> isf1
***   3) kz -> kgamma
*******************************************************************
      subroutine dsib2ffZsFSRIB(f, leg) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg

      real*8 q1
      complex*16 deltaH3, deltaZ
      complex*16 CH1, CH2, CZ1, CZ2, CZ3

      complex*16 ampH10, ampH11, ampH12, ampH13, ampH20, 
     &     ampH21, ampH22, ampH23
      complex*16 ampZ10, ampZ11, ampZ12, ampZ13, ampZ20,
     &     ampZ21, ampZ22, ampZ23

      real*8 Mh3, Mz
       
c... Defines FSR1 and FSR2 amplitudes through their simmetry relations via q1 exchange
c      if (leg.eq.1) then
         q1 = kJ
c      else 
      if (leg.eq.2) then
         q1 = -kJ
      endif

c..  masses/mx definition
      Mh3 = mass(kh3)/mx
      Mz = mass(kz)/mx

c... denominators
c... Higgs
      deltaH3 = ((2*E1J*EvJ + MB**2 - 2*CW*kv*q1)*
     -     (4*E1J*(E1J + EvJ) - Mh3**2 + MB**2))
    
c... Z
      deltaZ = ((2*E1J*EvJ + MB**2 - 2*CW*kv*q1)*
     -     (4*E1J*(E1J + EvJ) + MB**2 - Mz**2))

c... couplings
c... Higgs
      CH1 = gl(kgamma,f,f)*gr(kh3,f,f)*(Conjg(gr(kh3,kn(1),kn(1))) - 
     -     gr(kh3,kn(1),kn(1)))
      CH2 = gr(kgamma,f,f)*gr(kh3,f,f)*(Conjg(gr(kh3,kn(1),kn(1))) - 
     -     gr(kh3,kn(1),kn(1)))
c... Z
      CZ1 = gl(kgamma,f,f)**2*Dble(gl(kz,kn(1),kn(1)))
      CZ2 = gr(kgamma,f,f)**2*Dble(gl(kz,kn(1),kn(1)))
      CZ3 = gl(kgamma,f,f)*gr(kgamma,f,f)*Dble(gl(kz,kn(1),kn(1)))

c... helicity amplitudes (entries as for FSR1)
c... Higgs

      ampH10 = (
     -     Sqrt(2.)*(CH1*(-2*E1J - EvJ + kv) 
     -     - CH2*(2*E1J + EvJ + kv))*q1*SW)/deltaH3

      ampH11 = (
     -     (0,1)*(1 + CW)*(CH1*(-2*E1J - EvJ + kv) 
     -     + CH2*(2*E1J + EvJ + kv))*Mf)/deltaH3

      ampH12 = (
     -     Sqrt(2.)*E1J*(CH1*(2*E1J + EvJ - kv)
     -     - CH2*(2*E1J + EvJ + kv))*SW)/deltaH3

      ampH13 = (
     -     (0,1)*(-1 + CW)*(CH1*(-2*E1J - EvJ + kv) 
     -     + CH2*(2*E1J + EvJ + kv))*Mf)/deltaH3

      ampH20 = (
     -     Sqrt(2.)*(CH2*(-2*E1J - EvJ + kv)
     -     - CH1*(2*E1J + EvJ + kv))*q1*SW)/deltaH3

      ampH21 = (
     -     (0,-1)*(-1 + CW)*(CH2*(-2*E1J - EvJ + kv) 
     -     + CH1*(2*E1J + EvJ + kv))*Mf)/deltaH3

      ampH22 = (
     -     Sqrt(2.)*E1J*(CH2*(-2*E1J - EvJ + kv) 
     -     + CH1*(2*E1J + EvJ + kv))*SW)/deltaH3

      ampH23 = (
     -     (0,-1)*(1 + CW)*(CH2*(-2*E1J - EvJ + kv) 
     -     + CH1*(2*E1J + EvJ + kv))*Mf)/deltaH3

c... Z
      ampZ10 = (
     -        (Sqrt(2.)*(-2*CZ3*kv + CZ1*(-4*E1J - EvJ + kv)
     -        + CZ2*(4*E1J + EvJ + kv))*
     -        Mf*(4*E1J*(E1J + EvJ) + MB**2 - Mz**2)*q1*SW)/Mz**2)/deltaZ
      ampZ11 =(
     -     ((0,1)*(1 + CW)*(-4*E1J*(E1J + EvJ) - MB**2 + Mz**2)*
     -     (-2*CZ3*(2*E1J + EvJ)*Mf**2 + 
     -     CZ1*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*kv*(E1J + (-1 + CW)*q1))) + 
     -     CZ2*((2*E1J + EvJ - kv)*Mf**2 - 
     -     (-E1J + q1)*(MB**2 + 2*(E1J*(EvJ + kv)
     -     + (kv - CW*kv)*q1)))))/Mz**2)/deltaZ
      ampZ12 = (
     -     -((Sqrt(2.)*Mf*(-4*E1J*(E1J + EvJ) - MB**2 + Mz**2)*
     -     (-2*CZ3*E1J*(2*E1J + EvJ) + 
     -     CZ1*(E1J*(2*E1J + 3*EvJ - kv) + MB**2 - 2*CW*kv*q1) + 
     -     CZ2*(E1J*(2*E1J + 3*EvJ + kv) + MB**2 - 2*CW*kv*q1))*SW)
     -     /Mz**2))/deltaZ
      ampZ13 = (
     -     ((0,1)*(-1 + CW)*(-4*E1J*(E1J + EvJ) - MB**2 + Mz**2)*
     -     (-2*CZ3*(2*E1J + EvJ)*Mf**2 + 
     -     CZ2*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*(EvJ + kv) + MB**2 - 2*(1 + CW)*kv*q1)) + 
     -     CZ1*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(-2*E1J*EvJ - MB**2 + 2*kv*(E1J + q1 + CW*q1)))))
     -     /Mz**2)/deltaZ
      ampZ20 = (
     -     -((Sqrt(2.)*(-2*CZ3*kv + CZ2*(-4*E1J - EvJ + kv)
     -     + CZ1*(4*E1J + EvJ + kv))*
     -     Mf*(4*E1J*(E1J + EvJ) + MB**2 - Mz**2)*q1*SW)/Mz**2))/deltaZ
      ampZ21 = (
     -     ((0,1)*(-1 + CW)*(-4*E1J*(E1J + EvJ) - MB**2 + Mz**2)*
     -     (-2*CZ3*(2*E1J + EvJ)*Mf**2 + 
     -     CZ1*((2*E1J + EvJ - kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*(EvJ + kv) + MB**2 - 2*(1 + CW)*kv*q1)) + 
     -     CZ2*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (-E1J + q1)*(-2*E1J*EvJ - MB**2 + 2*kv*(E1J + q1 + CW*q1)))))
     -     /Mz**2)/deltaZ
      ampZ22 = (
     -     -((Sqrt(2.)*Mf*(-4*E1J*(E1J + EvJ) - MB**2 + Mz**2)*
     -     (-2*CZ3*E1J*(2*E1J + EvJ) + 
     -     CZ2*(E1J*(2*E1J + 3*EvJ - kv) + MB**2 - 2*CW*kv*q1) + 
     -     CZ1*(E1J*(2*E1J + 3*EvJ + kv) + MB**2 - 2*CW*kv*q1))*SW)
     -     /Mz**2))/deltaZ
      ampZ23 = (
     -     ((0,1)*(1 + CW)*(-4*E1J*(E1J + EvJ) - MB**2 + Mz**2)*
     -     (-2*CZ3*(2*E1J + EvJ)*Mf**2 + 
     -     CZ2*((2*E1J + EvJ + kv)*Mf**2 + 
     -     (E1J + q1)*(2*E1J*EvJ + MB**2 - 2*kv*(E1J + (-1 + CW)*q1))) + 
     -     CZ1*((2*E1J + EvJ - kv)*Mf**2 - 
     -     (-E1J + q1)*(MB**2 + 2*(E1J*(EvJ + kv) + (kv - CW*kv)*q1)))))
     -     /Mz**2)/deltaZ
         

      if (leg.eq.1) then
c... FSR1
         amp(1, 0) = amp(1, 0) + ampH10
     -    + ampZ10
         amp(1, 1) = amp(1, 1) + ampH11 
     -    + ampZ11
         amp(1, 2) = amp(1, 2) + ampH12
     -    + ampZ12
         amp(1, 3) = amp(1, 3) + ampH13
     -    + ampZ13
         amp(-1, 0) = amp(-1, 0) + ampH20
     -    + ampZ20
         amp(-1, 1) = amp(-1, 1) + ampH21 
     -    + ampZ21
         amp(-1, 2) = amp(-1, 2) + ampH22 
     -    + ampZ22
         amp(-1, 3) = amp(-1, 3) + ampH23
     -    + ampZ23
      else if (leg.eq.2) then
c... FSR2
         amp(1, 0) = amp(1, 0) - ampH20 
     -    - ampZ20
         amp(1, 1) = amp(1, 1) - ampH23
     -    - ampZ23
         amp(1, 2) = amp(1, 2) - ampH22
     -    - ampZ22
         amp(1, 3) = amp(1, 3) - ampH21 
     -    - ampZ21
         amp(-1, 0) = amp(-1, 0) - ampH10
     -    - ampZ10
         amp(-1, 1) = amp(-1, 1) - ampH13 
     -    - ampZ13
         amp(-1, 2) = amp(-1, 2) - ampH12 
     -    - ampZ12
         amp(-1, 3) = amp(-1, 3) - ampH11
     -    - ampZ11
      endif
  

      return
      end
      
