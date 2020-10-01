*******************************************************************
*** subroutine dsib2ffHituFSR add to the common block array amp(i,j), i = 0, 
*** the helicity amplitudes of the t+u channel FSR diagrams. 
*** The Hi is emitted by final legs 1 (FSR1) or 2 (FSR2) of the t+u channel diagrams.
*** isf1 is the index if the exchanged sfermion in th t-channel
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2014-02-20
*******************************************************************
      subroutine dsib2ffH0tuFSR(f, leg, isf1, hi)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf1, leg, hi
      
      integer kh(2), sf !(2)
      real*8 Msf ! (2)
      complex*16 delta
      complex*16 CFSR1L, CFSR1R, CFSR2L, CFSR2R, CFSR3L, CFSR3R
      complex*16 amp00, amp01, amp02, amp03

      real*8 q1
      integer dsib2getsfermion

      kh(1) = kh2 ! hi = 1 SM Higgs
      kh(2) = kh1 ! hi = 2 heavy Higgs

      sf  = dsib2getsfermion(f, isf1)
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

c... Defines FSR1 and FSR2 amplitudes through their symmetry relations via q1 exchange
c      if (leg.eq.1) then
         q1 = kJ
c      else 
       if (leg.eq.2) q1 = -kJ
       

c... denominators
      delta = (
     -     dcmplx(2*E1J*EvJ + MB**2 - 2*CW*kv*q1,Mf*width(f)/mx)*
     -     dcmplx(-E1J**2 + Mf**2 + MB**2/4. - CW*kv*q1 - Msf**2,
     -     Msf*width(sf)/mx*0d0))

c... couplings
      CFSR1L =Conjg(gl(sf,f,kn(1)))*gl(kh(hi),f,f)*
     - gl(sf,f,kn(1))
      CFSR2L =Conjg(gl(sf,f,kn(1)))*gl(kh(hi),f,f)*
     - gr(sf,f,kn(1))
      CFSR3L =Conjg(gr(sf,f,kn(1)))*gl(kh(hi),f,f)*
     - gr(sf,f,kn(1))
      CFSR1R =Conjg(gl(sf,f,kn(1)))*gl(sf,f,kn(1))*
     - gr(kh(hi),f,f)
      CFSR2R = Conjg(gl(sf,f,kn(1)))*gr(kh(hi),f,f)*
     - gr(sf,f,kn(1))
      CFSR3R =Conjg(gr(sf,f,kn(1)))*gr(kh(hi),f,f)*
     - gr(sf,f,kn(1))


c... helicity amplitudes (entries as for FSR1). (factor of 2 to get the t+u channel)
      amp00 =((0,-0.25)*(2*CFSR2R*Mf*(-4*(1 + E1J**2) + MB**2) + 
     - (CFSR1R + CFSR3R)*(4*E1J**4 + Mf**2*(-4 + MB**2) - 
     - E1J**2*(4 + 4*Mf**2 + MB**2 - 4*CW*kv*q1))))/E1J/delta

      amp01 = -((kv*(2*CFSR2R + (CFSR1R + CFSR3R)*Mf)*q1*SW)/Sqrt(2.))/delta

      amp02 =(0,0.25)*(CFSR1R - CFSR3R)*q1*(-MB**2 + 4*(-1 + E1J**2 + CW*kv*q1))/delta

      amp03 =(kv*(2*CFSR2R + (CFSR1R + CFSR3R)*Mf)*q1*SW)/Sqrt(2.)/delta

      if (leg.eq.1) then
c... FSR1
         amp(0, 0) = amp(0, 0) + 2*amp00
         amp(0, 1) = amp(0, 1) + 2*amp01
         amp(0, 2) = amp(0, 2) + 2*amp02
         amp(0, 3) = amp(0, 3) + 2*amp03

      else if (leg.eq.2) then
c... FSR2
         amp(0, 0) = amp(0, 0) + 2*amp00
         amp(0, 1) = amp(0, 1) + 2*amp03
         amp(0, 2) = amp(0, 2) + 2*amp02
         amp(0, 3) = amp(0, 3) + 2*amp01

      endif

      return
      end
