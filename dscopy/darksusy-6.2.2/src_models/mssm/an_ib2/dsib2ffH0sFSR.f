*******************************************************************
*** subroutine dsib2ffHisFSR add to the common block array amp(i,j), i = 0, 
*** the helicity amplitudes of the s channel FSR diagrams. 
*** The Hi is emitted from the final legs of the s channel diagrams.
*** Hi and Z diagrams are considered separately and then
*** added to the corresponding entries of amp(0, 0:3). 
*** 
*** Author: Francesca Calore, 2014-02-20
*******************************************************************
      subroutine dsib2ffH0sFSR(f, leg, hi) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg, hi

      real*8 q1
      complex*16 deltaH3, deltaZ
      complex*16 CH31, CZ1, CZ2

      complex*16 ampH00, ampH01, ampH02, ampH03
      complex*16 ampZ00, ampZ01, ampZ02, ampZ03

      real*8 Mh3, Mkz
      integer kh(2)

      kh(1) = kh2 ! hi = 1 SM Higgs
      kh(2) = kh1 ! hi = 2 heavy Higgs
       
c... Defines FSR1 and FSR2 amplitudes through their simmetry relations via q1 exchange
      if (leg.eq.1) then
         q1 = kJ
      else if (leg.eq.2) then
         q1 = -kJ
      else
         write(*,*) 'dsib2ffH0sFSR: ', f, leg, hi
         stop
      endif

c..  masses/mx definition
      Mh3 = mass(kh3)/mx 
      Mkz = mass(kz)/mx

c... denominators
c... Higgs
      deltaH3 = (dcmplx(2*E1J*EvJ + MB**2 - 2*CW*kv*q1,Mf*width(f)/mx)*
     -     dcmplx(2*E1J**2 + 4*E1J*EvJ + MB**2 + 2*Mf**2 - Mh3**2 + 
     -     2*q1**2,Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(2*E1J*EvJ + MB**2 - 2*CW*kv*q1,Mf*width(f)/mx)*
     -     dcmplx(2*E1J**2 + 4*E1J*EvJ + MB**2 
     -     + 2*Mf**2 - Mkz**2 + 2*q1**2,Mkz*width(kz)/mx))

c... couplings
c... Higgs
      CH31 =gr(kh3,f,f)*(Conjg(gr(kh3,kn(1),kn(1))) - 
     - gr(kh3,kn(1),kn(1)))*gr(kh(hi),f,f)

c... Z
      CZ1 =gl(kz,f,f)*gr(kh(hi),f,f)*gl(kz,kn(1),kn(1))
      CZ2 =gr(kz,f,f)*gr(kh(hi),f,f)*gl(kz,kn(1),kn(1))
    
c... helicity amplitudes (entries as for FSR1)
c... Higgs
      ampH00 =((0,4)*CH31*(2*E1J + EvJ)*Mf)/deltaH3

      ampH01 = -2*Sqrt(2.)*CH31*kv*q1*SW/deltaH3

      ampH02 = 0d0

      ampH03 = 2*Sqrt(2.)*CH31*kv*q1*SW/deltaH3

c...  Z
      ampZ00 = (((0,2)*(CZ1 - CZ2)*(2*(2*E1J + EvJ)*Mf**2 + 
     -     E1J*(2*E1J*EvJ + MB**2 - 2*CW*kv*q1))*
     -     (MB**2 + 2*Mf**2 - Mkz**2 + 2*(E1J**2 + 2*E1J*EvJ + q1**2)))/
     -     Mkz**2)/deltaZ

      ampZ01 = ((-2*Sqrt(2.)*(CZ1 - CZ2)*kv*Mf*q1*
     -    (MB**2 + 2*Mf**2 - Mkz**2 + 2*(E1J**2 + 2*E1J*EvJ + q1**2))*SW
     -    )/Mkz**2)/deltaZ

      ampZ02 =(((0,2)*(CZ1 + CZ2)*(2*CW*kv*(-E1J + Mf)*(E1J + Mf) + 
     -     2*E1J*EvJ*q1 + MB**2*q1)*
     -     (MB**2 + 2*Mf**2 - Mkz**2 + 2*(E1J**2 + 2*E1J*EvJ + q1**2)))/
     -     Mkz**2)/deltaZ

      ampZ03 = ((2*Sqrt(2.)*(CZ1 - CZ2)*kv*Mf*q1*
     -    (MB**2 + 2*Mf**2 - Mkz**2 + 2*(E1J**2 + 2*E1J*EvJ + q1**2))*SW
     -    )/Mkz**2)/deltaZ

c... exploiting symmetries
      if (leg.eq.1) then
c... FSR1
         amp(0, 0) = amp(0, 0) + ampH00 + ampZ00
         amp(0, 1) = amp(0, 1) + ampH01 + ampZ01
         amp(0, 2) = amp(0, 2) + ampH02 + ampZ02
         amp(0, 3) = amp(0, 3) + ampH03 + ampZ03
      else if (leg.eq.2) then
c... FSR2
         amp(0, 0) = amp(0, 0) + ampH00 + ampZ00
         amp(0, 1) = amp(0, 1) + ampH03 + ampZ03
         amp(0, 2) = amp(0, 2) + ampH02 + ampZ02
         amp(0, 3) = amp(0, 3) + ampH01 + ampZ01
      endif

      return
      end
