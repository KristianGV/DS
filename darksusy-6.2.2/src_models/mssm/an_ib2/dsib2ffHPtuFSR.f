*******************************************************************
*** subroutine dsib2ffHPtuFSR add to the common block array amp(i,j), i = 0,  
*** the helicity amplitudes of the t+u channel FSR diagrams. 
*** The Hp is emitted by final legs 1 (FSR1) or 2 (FSR2) of the t+u channel diagrams.
*** isf1 is the index if the exchanged sfermion in th t-channel
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2014-02-22
***         Francesca Calore, 2015-02-15 (bug solved cf with CH)
*******************************************************************
      subroutine dsib2ffHPtuFSR(f, leg, isf1)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg, isf1
      
      integer ff, sf, sff ! (2)
      real*8 Msf, Msff

      complex*16 delta
      complex*16 CFSR1, CFSR2, CFSR3, CFSR4, CFSR5, CFSR6
      complex*16 amp00, amp01, amp02, amp03
      integer dsib2getsfermion

      ff   = c_fbartype
      
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

c... Defines FSR1 and FSR2 amplitudes 
      if (leg.eq.1) then
        sff  = dsib2getsfermion(ff, isf1)
        Msff = mass(sff)/mx ! dimensionless in order to be consistent with kinematics
      else
        goto 100
      endif  

c... couplings
         CFSR1 =Conjg(gl(sff,ff,kn(1)))*
     -    gl(khc,f,ff)*gl(sff,ff,kn(1))
         CFSR2 =Conjg(gl(sff,ff,kn(1)))*
     -    gl(khc,f,ff)*gr(sff,ff,kn(1))
         CFSR3 =Conjg(gr(sff,ff,kn(1)))*
     -    gl(khc,f,ff)*gr(sff,ff,kn(1))
         CFSR4 =Conjg(gl(sff,ff,kn(1)))*
     -    gl(sff,ff,kn(1))*gr(khc,f,ff)
         CFSR5 =Conjg(gl(sff,ff,kn(1)))*
     -    gr(khc,f,ff)*gr(sff,ff,kn(1))
         CFSR6 =Conjg(gr(sff,ff,kn(1)))*
     -    gr(khc,f,ff)*gr(sff,ff,kn(1))

c... denominators
         delta = (dcmplx(2*E1J*EvJ + Mf**2
     -    - Mff**2 + MB**2 - 2*CW*kv*kJ,
     -    Mff*width(ff)/mx)*
     -    dcmplx((-2*E1J*E2J + 2*E1J*EvJ - 2*E2J*EvJ
     -    + MB**2 + Mf**2 + Mff**2 - 
     -    4*CW*kv*kJ - 2*kJ**2 - 4*Msff**2)/4.,
     -    Msff*width(sff)/mx*0d0))

c... Helicity amplitudes (entries as for FSR1)
         
         amp00 = ((0,0.5)*(2*CFSR2*Epl*Mff + 2*CFSR5*Epl*Mff + 
     -    2*CFSR2*(Emi*EvJ + Epl*Mf + CW*kv*pmi) + 
     -    2*CFSR5*(Emi*EvJ + Epl*Mf + CW*kv*pmi) + 
     -    CFSR1*Mff*(Emi*EvJ + Epl*(Mf + Mff) + CW*kv*pmi) + 
     -    CFSR6*Mff*(Emi*EvJ + Epl*(Mf + Mff) + CW*kv*pmi) + 
     -    CFSR3*(Mff*(Emi*EvJ + CW*kv*pmi) + 
     -    Epl*(2*E1J*EvJ + MB**2 + Mf**2 + Mf*Mff - 2*CW*kv*kJ)) + 
     -    CFSR4*(Mff*(Emi*EvJ + CW*kv*pmi) + 
     -    Epl*(2*E1J*EvJ + MB**2 + Mf**2 
     -    + Mf*Mff - 2*CW*kv*kJ))))/delta

         amp01 = ((kv*(2*CFSR2*(Epl - ppl) - 2*CFSR5*(Epl + ppl) + 
     -    Mff*(CFSR1*(Epl - ppl) + CFSR3*(Epl - ppl) - 
     -    (CFSR4 + CFSR6)*(Epl + ppl)))*SW)/(2.*Sqrt(2.)))/delta
         
         amp02 = ((0,0.5)*(2*CFSR2*Mff*ppl - 2*CFSR5*Mff*ppl - 
     -     2*CFSR2*(CW*Emi*kv + EvJ*pmi + Mf*ppl) + 
     -     2*CFSR5*(CW*Emi*kv + EvJ*pmi + Mf*ppl) - 
     -     CFSR1*Mff*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl) + 
     -     CFSR6*Mff*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl) + 
     -     CFSR3*(-(CW*Emi*kv*Mff) - EvJ*Mff*pmi 
     -     + 2*E1J*EvJ*ppl + MB**2*ppl + 
     -     Mf**2*ppl - Mf*Mff*ppl - 2*CW*kv*ppl*kJ) + 
     -     CFSR4*(CW*Emi*kv*Mff + EvJ*Mff*pmi
     -        - 2*E1J*EvJ*ppl - MB**2*ppl - 
     -     Mf**2*ppl + Mf*Mff*ppl + 2*CW*kv*ppl*kJ)))/delta
         
         amp03 = ((kv*(-2*CFSR5*(Epl - ppl) + 2*CFSR2*(Epl + ppl) + 
     -    Mff*(-((CFSR4 + CFSR6)*(Epl - ppl)) + CFSR1*(Epl + ppl) + 
     -    CFSR3*(Epl + ppl)))*SW)/(2.*Sqrt(2.)))/delta


c      else 
100    if (leg.eq.2) then
         sf  = dsib2getsfermion(f, isf1)
         Msf = mass(sf)/mx ! dimensionless in order to be consistent with kinematics


c... couplings
         CFSR1 =Conjg(gl(sf,f,kn(1)))*
     -    gl(khc,f,ff)*gl(sf,f,kn(1))
         CFSR2 =Conjg(gl(sf,f,kn(1)))*
     -    gl(khc,f,ff)*gr(sf,f,kn(1))
         CFSR3 =Conjg(gr(sf,f,kn(1)))*
     -    gl(khc,f,ff)*gr(sf,f,kn(1))
         CFSR4 =Conjg(gl(sf,f,kn(1)))*
     -    gl(sf,f,kn(1))*gr(khc,f,ff)
         CFSR5 =Conjg(gl(sf,f,kn(1)))*
     -    gr(khc,f,ff)*gr(sf,f,kn(1))
         CFSR6 =Conjg(gr(sf,f,kn(1)))*
     -    gr(khc,f,ff)*gr(sf,f,kn(1))

c... denominators 
         delta = (dcmplx(2*E2J*EvJ - Mf**2 + Mff**2 + MB**2 + 2*CW*kv*kJ
     -        ,Mf*width(f)/mx)*
     -        dcmplx((-2*E1J*E2J - 2*E1J*EvJ + 2*E2J*EvJ 
     -        + MB**2 + Mf**2 + Mff**2 + 
     -        4*CW*kv*kJ - 2*kJ**2 - 4*Msf**2)/4.,
     -        Msf*width(sf)/mx*0d0))

c... Helicity amplitudes (entries as for FSR2)
         amp00 = ((0,0.5)*(2*CFSR2*Epl*Mf + 2*CFSR5*Epl*Mf + 
     -     2*CFSR2*(Emi*EvJ + Epl*Mff + CW*kv*pmi) + 
     -     2*CFSR5*(Emi*EvJ + Epl*Mff + CW*kv*pmi) + 
     -     CFSR3*Mf*(Emi*EvJ + Epl*(Mf + Mff) + CW*kv*pmi) + 
     -     CFSR4*Mf*(Emi*EvJ + Epl*(Mf + Mff) + CW*kv*pmi) + 
     -     CFSR1*(Mf*(Emi*EvJ + CW*kv*pmi) + 
     -     Epl*(2*E2J*EvJ + MB**2 + Mf*Mff + Mff**2 + 2*CW*kv*kJ)) + 
     -     CFSR6*(Mf*(Emi*EvJ + CW*kv*pmi) + 
     -     Epl*(2*E2J*EvJ + MB**2 + Mf*Mff 
     -     + Mff**2 + 2*CW*kv*kJ))))/delta
         
         amp01 = (-(kv*(-2*CFSR5*(Epl - ppl) + 2*CFSR2*(Epl + ppl) + 
     -    Mf*(-((CFSR4 + CFSR6)*(Epl - ppl)) + CFSR1*(Epl + ppl) + 
     -    CFSR3*(Epl + ppl)))*SW)/(2.*Sqrt(2.)))/delta
         
         amp02 = ((0,0.5)*(2*CFSR2*Mf*ppl - 2*CFSR5*Mf*ppl + 
     -     2*CFSR2*(CW*Emi*kv + EvJ*pmi - Mff*ppl) - 
     -     2*CFSR5*(CW*Emi*kv + EvJ*pmi - Mff*ppl) + 
     -     CFSR3*Mf*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl) - 
     -     CFSR4*Mf*(CW*Emi*kv + EvJ*pmi + Mf*ppl - Mff*ppl) + 
     -     CFSR1*(Mf*(CW*Emi*kv + EvJ*pmi - Mff*ppl) + 
     -     ppl*(2*E2J*EvJ + MB**2 + Mff**2 + 2*CW*kv*kJ)) - 
     -     CFSR6*(Mf*(CW*Emi*kv + EvJ*pmi - Mff*ppl) + 
     -     ppl*(2*E2J*EvJ + MB**2 + Mff**2 + 2*CW*kv*kJ))))/delta
         
         amp03 =((-(kv*(2*CFSR2*(Epl - ppl) - 2*CFSR5*(Epl + ppl) + 
     -    Mf*(CFSR1*(Epl - ppl) + CFSR3*(Epl - ppl) - 
     -   (CFSR4 + CFSR6)*(Epl + ppl)))*SW)/(2.*Sqrt(2.))))/delta

      endif

c... fermion = neutrinos
      if(f.eq.1.or.f.eq.3.or.f.eq.5) then
         if (leg.eq.1) then
c...  FSR1 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) + 2*amp00
            amp(0, 1) = amp(0, 1) + 2*amp01
            amp(0, 2) = amp(0, 2) + 2*amp02
            amp(0, 3) = amp(0, 3) + 2*amp03
         else if (leg.eq.2) then
c...  FSR2 (factor of 2 to get the t+u channel) 
c...  factor of 1/2. because of only one sneutrino exachanged (vs sum isf=1,2 in dsib2MsqffHP)
            amp(0, 0) = amp(0, 0) + 2*amp00/2.
            amp(0, 1) = amp(0, 1) + 2*amp01/2.
            amp(0, 2) = amp(0, 2) + 2*amp02/2.
            amp(0, 3) = amp(0, 3) + 2*amp03/2.
         endif
c... fermion = leptons
      elseif(f.eq.2.or.f.eq.4.or.f.eq.6) then
         if (leg.eq.1) then
c...  FSR1 (factor of 2 to get the t+u channel)
c...  factor of 1/2. because of only one sneutrino exachanged (vs sum isf=1,2 in dsib2MsqffHP)
            amp(0, 0) = amp(0, 0) + 2*amp00/2.
            amp(0, 1) = amp(0, 1) + 2*amp01/2.
            amp(0, 2) = amp(0, 2) + 2*amp02/2.
            amp(0, 3) = amp(0, 3) + 2*amp03/2.
         else if (leg.eq.2) then
c...  FSR2 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) + 2*amp00
            amp(0, 1) = amp(0, 1) + 2*amp01
            amp(0, 2) = amp(0, 2) + 2*amp02
            amp(0, 3) = amp(0, 3) + 2*amp03
         endif
c... quarks
      else
         if (leg.eq.1) then
c...  FSR1 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) + 2*amp00
            amp(0, 1) = amp(0, 1) + 2*amp01
            amp(0, 2) = amp(0, 2) + 2*amp02
            amp(0, 3) = amp(0, 3) + 2*amp03
         else if (leg.eq.2) then
c...  FSR2 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) + 2*amp00
            amp(0, 1) = amp(0, 1) + 2*amp01
            amp(0, 2) = amp(0, 2) + 2*amp02
            amp(0, 3) = amp(0, 3) + 2*amp03
         endif 
      endif

      return
      end

