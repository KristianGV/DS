*******************************************************************
*** subroutine dsib2ffZtuVIB add to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel VIB diagrams. 
*** Z emitted by internal legs: the sfermion - isf1 - radiates a Z by becoming a "second" sfermion - isf2 - 
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffZtuVIB(f, isf1, isf2)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf1, isf2
      
      integer sf(2), i, dsib2getsfermion 
      real*8 Msf(2)
      complex*16 delta
      complex*16 CIB1, CIB2, CIB3, CIB4
      complex*16 amp10, amp11, amp12, amp13
      

      sf(1) = dsib2getsfermion(f, isf1)
      sf(2) = dsib2getsfermion(f, isf2)
      do i=1,2
        Msf(i) = mass(sf(i))/mx ! dimensionless in order to be consistent with kinematics
      enddo
 
c TB: obsolete code after introducing dsib2getsfermion    
c      if(f.EQ.knue) sf(1) = ksnue
c      if(f.EQ.knumu) sf(1) = ksnumu
c      if(f.EQ.knutau) sf(1) = ksnutau
c      if(f.EQ.ke) sf(1) = kse(1)
c      if(f.EQ.kmu) sf(1) = ksmu(1)
c      if(f.EQ.ktau) sf(1) = kstau(1)
c      if(f.EQ.ku) sf(1) = ksu(1)
c      if(f.EQ.kd) sf(1) = ksd(1)
c      if(f.EQ.kc) sf(1) = ksc(1)
c      if(f.EQ.ks) sf(1) = kss(1)
c      if(f.EQ.kb) sf(1) = ksb(1)
c      if(f.EQ.kt) sf(1) = kst(1)
c      if(f.eq.knue.or.f.eq.knumu.or.f.eq.knutau) then
c         sf(2) = sf(1)
c      else
c         sf(2)=sf(1)+1
c      endif
c... masses of exchanged particle are divided by mx in order to be consistent with the kinematics
c      Msf(1) = mass(sf(1))/mx
c      Msf(2) = mass(sf(2))/mx

c... denominators
      delta = (
     -     dcmplx(-E1J**2 - CW*kJ*kv + Mf**2 + MB**2/4. - Msf(2)**2,
     -     Msf(2)*width(sf(2))/mx*0d0)*
     -     dcmplx(-E1J**2 + CW*kJ*kv + Mf**2 + MB**2/4. - Msf(1)**2,
     -     Msf(1)*width(sf(1))/mx*0d0)
     -     )

c... couplings
      CIB1 = Conjg(gl(sf(2),f,kn(1)))*gl(kz,sf(1),sf(2))*
     -     gl(sf(1),f,kn(1))
      CIB2 = Conjg(gl(sf(2),f,kn(1)))*gl(kz,sf(1),sf(2))*
     -     gr(sf(1),f,kn(1))
      CIB3 = Conjg(gr(sf(2),f,kn(1)))*gl(kz,sf(1),sf(2))*
     -     gl(sf(1),f,kn(1))
      CIB4 = Conjg(gr(sf(2),f,kn(1)))*gl(kz,sf(1),sf(2))*
     -     gr(sf(1),f,kn(1))


c... Helicity amplitudes (factor of 2 to get the t+u channel)
      amp(0, 0) = amp(0, 0) + 2*(
     -     ((0,-0.0625)*CW*kJ*(-16*E1J**4*(2*(CIB2 + CIB3)
     -      + (CIB1 + CIB4)*Mf) - 8*(CIB2 + CIB3)*E1J**2*(-4 + MB**2) + 
     -     (CIB1 + CIB4)*Mf*(-4 + MB**2)**2))/(E1J**2*MB))/delta
      
      amp(0, 1) = amp(0, 1) + 2*(
     -     -(CW*kJ*(CIB1*(-E1J + kJ) 
     -     + CIB4*(E1J + kJ))*(-4 + 4*E1J**2 + MB**2)*
     -     Sqrt(8*(-1 + E1J**2)**2-4*(1 + E1J**2)*MB**2 + MB**4/2.)*SW)/
     -     (16.*E1J**2*MB))/delta
      
      amp(0, 2) = amp(0, 2) + 2*(
     -     ((0,-0.0625)*CW*kJ*(-4 + 4*E1J**2 + MB**2)*
     -     (8*(CIB2 - CIB3)*E1J*kJ + 
     -     (CIB1 - CIB4)*CW*Mf*Sqrt(16*E1J**4 + (-4 + MB**2)**2 - 
     -     8*E1J**2*(4 + MB**2))))/(E1J**2*MB))/delta
      
      amp(0, 3) = amp(0, 3) + 2*(
     -     (CW*kJ*(CIB4*(-E1J + kJ) 
     -     + CIB1*(E1J + kJ))*(-4 + 4*E1J**2 + MB**2)*
     -     Sqrt(8*(-1 + E1J**2)**2-4*(1 + E1J**2)*MB**2 + MB**4/2.)*SW)/
     -     (16.*E1J**2*MB))/delta

c... exploiting simmetries
      amp10 = ((kJ*(-8*(CIB2 + CIB3)*E1J**2 + 
     -     (CIB1 + CIB4)*Mf*(-4*(1 + E1J**2) + MB**2))*SW)
     -     /(4.*Sqrt(2.)*E1J))/delta

      amp11 = (((0,0.125)*(-1 + CW**2)*kJ*(CIB1*(-E1J + kJ)
     -     + CIB4*(E1J + kJ))*Sqrt(16*E1J**4 + (-4 + MB**2)**2
     -     - 8*E1J**2*(4 + MB**2)))/E1J)/delta

      amp12 = ( (kJ*(8*(CIB2 - CIB3)*E1J*kJ + 
     -     (CIB1 - CIB4)*CW*Mf*Sqrt(16*E1J**4 + (-4 + MB**2)**2 - 
     -     8*E1J**2*(4 + MB**2)))*SW)/(4.*Sqrt(2.)*E1J))/delta

      amp13 = (((0,-0.125)*(-1 + CW**2)*kJ*(CIB4*(-E1J + kJ) 
     -     + CIB1*(E1J + kJ))*Sqrt(16*(-1 + E1J**2)**2 
     -     - 8*(1 + E1J**2)*MB**2 + MB**4))/E1J)/delta

      amp(1, 0) = amp(1, 0) + 2*amp10
      amp(1, 1) = amp(1, 1) + 2*amp11
      amp(1, 2) = amp(1, 2) + 2*amp12
      amp(1, 3) = amp(1, 3) + 2*amp13

      amp(-1, 0) = amp(-1, 0) + 2*amp10
      amp(-1, 1) = amp(-1, 1) + 2*amp11
      amp(-1, 2) = amp(-1, 2) + 2*amp12
      amp(-1, 3) = amp(-1, 3) + 2*amp13


      return
      end

c===================================================================
*******************************************************************
***   TEMPORARY: for IB check
***   1) w/o width propagator
***   2) isf2 -> isf1
***   3) kz -> kgamma
*******************************************************************
      subroutine dsib2ffZtVIBIB(f, isf1)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf1
      
      integer sf, dsib2getsfermion 
      real*8 Msf
      complex*16 delta
      complex*16 CIB1, CIB2, CIB3, CIB4
      complex*16 amp10, amp11, amp12, amp13
      
      sf = dsib2getsfermion(f, isf1)
      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics

 
c TB: obsolete code after introducing dsib2getsfermion    
c      if(f.EQ.knue) sf(1) = ksnue
c      if(f.EQ.knumu) sf(1) = ksnumu
c      if(f.EQ.knutau) sf(1) = ksnutau
c      if(f.EQ.ke) sf(1) = kse(1)
c      if(f.EQ.kmu) sf(1) = ksmu(1)
c      if(f.EQ.ktau) sf(1) = kstau(1)
c      if(f.EQ.ku) sf(1) = ksu(1)
c      if(f.EQ.kd) sf(1) = ksd(1)
c      if(f.EQ.kc) sf(1) = ksc(1)
c      if(f.EQ.ks) sf(1) = kss(1)
c      if(f.EQ.kb) sf(1) = ksb(1)
c      if(f.EQ.kt) sf(1) = kst(1)
c      if(f.eq.knue.or.f.eq.knumu.or.f.eq.knutau) then
c         sf(2) = sf(1)
c      else
c         sf(2)=sf(1)+1
c      endif
c... masses of exchanged particle are divided by mx in order to be consistent with the kinematics
c      Msf(1) = mass(sf(1))/mx
c      Msf(2) = mass(sf(2))/mx


c... denominators
      delta =  (
     -     (-E1J**2 - CW*kJ*kv + Mf**2 + MB**2/4. - Msf**2)*
     -     (-E1J**2 + CW*kJ*kv + Mf**2 + MB**2/4. - Msf**2)
     -     )

c... couplings
      CIB1 = Conjg(gl(sf,f,kn(1)))*gl(kgamma,sf,sf)*
     -     gl(sf,f,kn(1))
      CIB2 = Conjg(gl(sf,f,kn(1)))*gl(kgamma,sf,sf)*
     -     gr(sf,f,kn(1))
      CIB3 = Conjg(gr(sf,f,kn(1)))*gl(kgamma,sf,sf)*
     -     gl(sf,f,kn(1))
      CIB4 = Conjg(gr(sf,f,kn(1)))*gl(kgamma,sf,sf)*
     -     gr(sf,f,kn(1))
c$$$      CIB1 =  Conjg(gl(sf,f,kn(1)))*gl(kgamma,sf,sf)*
c$$$     -     gl(sf,f,kn(1))
c$$$      CIB2 = 0.d0
c$$$      CIB3 = 0.d0
c$$$      CIB4 = 0.d0


c... Helicity amplitudes (factor of 2 to get the t+u channel)
c... exploiting simmetries
      amp10 = ((kJ*(-8*(CIB2 + CIB3)*E1J**2 + 
     -     (CIB1 + CIB4)*Mf*(-4*(1 + E1J**2) + MB**2))*SW)
     -     /(4.*Sqrt(2.)*E1J))/delta

      amp11 = (((0,0.125)*(-1 + CW**2)*kJ*(CIB1*(-E1J + kJ)
     -     + CIB4*(E1J + kJ))*Sqrt(16*E1J**4 + (-4 + MB**2)**2
     -     - 8*E1J**2*(4 + MB**2)))/E1J)/delta

      amp12 = ( (kJ*(8*(CIB2 - CIB3)*E1J*kJ + 
     -     (CIB1 - CIB4)*CW*Mf*Sqrt(16*E1J**4 + (-4 + MB**2)**2 - 
     -     8*E1J**2*(4 + MB**2)))*SW)/(4.*Sqrt(2.)*E1J))/delta

      amp13 = (((0,-0.125)*(-1 + CW**2)*kJ*(CIB4*(-E1J + kJ) 
     -     + CIB1*(E1J + kJ))*Sqrt(16*(-1 + E1J**2)**2 
     -     - 8*(1 + E1J**2)*MB**2 + MB**4))/E1J)/delta

      amp(1, 0) = amp(1, 0) + 2*amp10
      amp(1, 1) = amp(1, 1) + 2*amp11
      amp(1, 2) = amp(1, 2) + 2*amp12
      amp(1, 3) = amp(1, 3) + 2*amp13

      amp(-1, 0) = amp(-1, 0) + 2*amp10
      amp(-1, 1) = amp(-1, 1) + 2*amp11
      amp(-1, 2) = amp(-1, 2) + 2*amp12
      amp(-1, 3) = amp(-1, 3) + 2*amp13

      return
      end

