*******************************************************************
*** subroutine dsib2ffHituVIB add to the common block array amp(i,j), i = 0,
*** the helicity amplitudes of the t+u channel VIB diagrams. 
*** Hi emitted by internal legs: the sfermion - isf1 - radiates a Z by becoming a "second" sfermion - isf2 - 
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2014-02-20
*******************************************************************
      subroutine dsib2ffH0tuVIB(f, isf1, isf2, hi)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf1, isf2, hi
      
      integer sf(2), kh(2), i
      real*8 Msf(2), Mkw
      complex*16 delta
      complex*16 CVIB1, CVIB2, CVIB3, CVIB4
      integer dsib2getsfermion

      kh(1) = kh2 ! hi = 1 SM Higgs
      kh(2) = kh1 ! hi = 2 heavy Higgs
      
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
      Mkw =1d0/mx ! H-sf-sf coupling has mass dimension (dsvertx1.f), meaning that mw is included in the definition of the coupling

c... denominators
      delta = (
     -     dcmplx(-E1J**2 - CW*kJ*kv + Mf**2 + MB**2/4. - Msf(2)**2,
     -     Msf(2)*width(sf(2))/mx*0d0)*
     -     dcmplx(-E1J**2 + CW*kJ*kv + Mf**2 + MB**2/4. - Msf(1)**2,
     -     Msf(1)*width(sf(1))/mx*0d0)
     -     )

c... couplings
      CVIB1 = Conjg(gl(sf(2),f,kn(1)))*gl(kh(hi),sf(1),sf(2))*
     - gl(sf(1),f,kn(1))
      CVIB2 =Conjg(gl(sf(2),f,kn(1)))*gl(kh(hi),sf(1),sf(2))*
     - gr(sf(1),f,kn(1))
      CVIB3 =Conjg(gr(sf(2),f,kn(1)))*gl(kh(hi),sf(1),sf(2))*
     - gl(sf(1),f,kn(1))
      CVIB4 =Conjg(gr(sf(2),f,kn(1)))*gl(kh(hi),sf(1),sf(2))*
     - gr(sf(1),f,kn(1))


c... Helicity amplitudes (factor of 2 to get the t+u channel)
      amp(0, 0) = amp(0, 0) + 2*(((0,-0.125)*Mkw*
     -     (-8*(CVIB2 + CVIB3)*E1J**2 + 
     -     (CVIB1 + CVIB4)*Mf*(-4*(1 + E1J**2) + MB**2)))/E1J/delta)

      amp(0, 1) = amp(0, 1) + 2*(-((CVIB1*(-E1J + kJ) +
     - CVIB4*(E1J + kJ))*kv*Mkw*SW)/(2.*Sqrt(2.))/delta)

      amp(0, 2) = amp(0, 2) + 2*((0,-0.5)*(2*CVIB2*kJ - 2*CVIB3*kJ +
     - (CVIB1 - CVIB4)*CW*kv*Mf)*Mkw/delta)

      amp(0, 3) = amp(0, 3) + 2*(((CVIB4*(-E1J + kJ) + 
     - CVIB1*(E1J + kJ))*kv*Mkw*SW)/(2.*Sqrt(2.))/delta)

      return
      end

