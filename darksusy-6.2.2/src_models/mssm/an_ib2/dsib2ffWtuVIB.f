*******************************************************************
*** subroutine dsib2ffWtuVIB add to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel VIB diagrams. 
*** W emitted by internal legs: the sfermion - isf or isff - radiates a Z by becoming a "second" sfermion - isff or isf - 
*** f stays for the fermion, ff the antifermion
*** sf: sfermion of the fermion; sff: sfermion of the antifermion
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffWtuVIB(f, isf, isff)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, isf, isff
      integer ff, dsib2getsfermion
      
      integer sf, sff  ! (2)
      real*8 Msf, Msff ! (2)
      complex*16 delta
      complex*16 CIB1, CIB2, CIB3, CIB4
      complex*16 amp10, amp11, amp12, amp13

      sf  = dsib2getsfermion(f, isf)
      Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics

      ff   = c_fbartype
      sff  = dsib2getsfermion(ff, isff)
      Msff = mass(sff)/mx ! dimensionless in order to be consistent with kinematics

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
c
c      if (f.EQ.kd.or.f.EQ.ks.or.f.EQ.kb) then
c         sf(2) = sf(1)+1
c         sff(1) = sf(1)-2
c         sff(2)= sff(1)+1
c      elseif(f.EQ.ku.or.f.EQ.kc.or.f.EQ.kt) then
c         sf(2) = sf(1)+1
c         sff(1) = sf(1)+2
c         sff(2)= sff(1)+1
c      elseif(f.EQ.ke.or.f.EQ.kmu.or.f.EQ.ktau) then
c         sf(2) = sf(1)+1
c         sff(1) = sf(1)-1
c         sff(2)= sff(1)
c      elseif(f.EQ.knue.or.f.EQ.knumu.or.f.EQ.knutau) then 
c         sf(2) = sf(1)
c         sff(1) = sf(1)+1
c         sff(2) = sff(1)+1
c      endif
c      
c      Msf(1) = mass(sf(1))/mx
c      Msf(2) = mass(sf(2))/mx
c      Msff(1) = mass(sff(1))/mx
c      Msff(2) = mass(sff(2))/mx 

c... couplings
      CIB1 = Conjg(gl(sff,ff,kn(1)))*gl(kw,sff,sf)*
     -     gl(sf,f,kn(1))
      CIB2 = Conjg(gl(sff,ff,kn(1)))*gl(kw,sff,sf)*
     -     gr(sf,f,kn(1))
      CIB3 = Conjg(gr(sff,ff,kn(1)))*gl(kw,sff,sf)*
     -     gl(sf,f,kn(1))
      CIB4 = Conjg(gr(sff,ff,kn(1)))*gl(kw,sff,sf)*
     -     gr(sf,f,kn(1))

c... denominators
      delta = (
     -     dcmplx((-2*E1J**2 - 2*E1J*E2J + 2*E1J*EvJ - 2*E2J*EvJ 
     -     - 4*CW*kJ*kv + 3*Mf**2 + Mff**2 + MB**2 - 4*Msff**2)/4.
     -     ,Msff*width(sff)/mx*0d0)*
     -     dcmplx((-2*E1J**2 - 2*E1J*E2J - 2*E1J*EvJ + 2*E2J*EvJ  
     -     + 4*CW*kJ*kv + 3*Mf**2 + Mff**2 + MB**2 - 4*Msf**2)/4.
     -     ,Msf*width(sf)/mx*0d0))

c... Helicity amplitudes
      amp(0, 0) = amp(0, 0) + 2*
     -     (((0,-0.5)*(2*CW*EvJ*kJ + (-E1J + E2J)*kv)*
     -     (2*CIB2*Epl + 2*CIB3*Epl + 
     -     (CIB1 + CIB4)*(EvJ*Emi + Epl*(Mf + Mff) 
     -     + CW*kv*pmi)))/MB)/delta
      
      amp(0, 1) = amp(0, 1) + 2*
     -     ((kv*(-2*CW*EvJ*kJ + (E1J - E2J)*kv)*(CIB1*(Epl - ppl)
     -     - CIB4*(Epl + ppl))*SW)/(2.*Sqrt(2.)*MB))/delta
      
      amp(0, 2) = amp(0, 2) + 2*
     -     (((0,0.5)*(2*CW*EvJ*kJ + (-E1J + E2J)*kv)*
     -     (2*CIB2*ppl - 2*CIB3*ppl + 
     -     (CIB1 - CIB4)*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl)))/MB)
     -     /delta
      
      amp(0, 3) = amp(0, 3) + 2*
     -     ((kv*(-2*CW*EvJ*kJ + (E1J - E2J)*kv)*
     -     (CIB4*(-Epl + ppl) + CIB1*(Epl + ppl))*SW)
     -     /(2.*Sqrt(2.)*MB))/delta

c... exploiting simmetries
      amp10 = 
     -     (-((kJ*(2*CIB2*Epl + 2*CIB3*Epl + 
     -     (CIB1 + CIB4)*(EvJ*Emi + Epl*(Mf + Mff) + CW*kv*pmi))*SW)
     -     /Sqrt(2.)))/delta

      amp11 = 
     -     ((0,-0.5)*(-1 + CW**2)*kJ*kv*(CIB1*(Epl - ppl) 
     -     - CIB4*(Epl + ppl)))/delta

      amp12 = 
     -     ((kJ*(2*CIB2*ppl - 2*CIB3*ppl + 
     -     (CIB1 - CIB4)*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl))*SW)
     -     /Sqrt(2.))/delta

       amp13 = 
     -     ((0,-0.5)*(-1 + CW**2)*kJ*kv*(CIB4*(-Epl + ppl)
     -     + CIB1*(Epl + ppl)))/delta

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
