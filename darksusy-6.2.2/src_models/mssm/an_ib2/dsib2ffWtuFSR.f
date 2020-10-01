*******************************************************************
*** subroutine dsib2ffWtuFSR add to the common block array amp(i,j) 
*** the helicity amplitudes of the t+u channel FSR diagrams. 
*** The W is emitted by final legs 1 (FSR1) or 2 (FSR2) of the t+u channel diagrams.
*** isf1 is the index if the exchanged sfermion in th t-channel
***
*** Notice: in the v->0 limit, to get the t+u amplitude a factor of 2* is required! 
***         Helicity amplitudes are multiplied by 2.
*** 
*** Author: Francesca Calore, 2013-04-10
*******************************************************************
      subroutine dsib2ffWtuFSR(f, leg, isf1)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg, isf1
      
      integer ff

      complex*16 delta
      complex*16 CFSR1, CFSR2, CFSR3

      integer tmpf, tmpff, tmpsf
      real*8 m11, m22, Mtmpsf,
     &     tmpE1J, tmpE2J, tmpq1, tmpppl

      complex*16 amp00, amp01, amp02, amp03, amp10, amp11, amp12, amp13, 
     &     amp20, amp21, amp22, amp23
      integer dsib2getsfermion


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

       ff   = c_fbartype
c... Defines FSR1 and FSR2 amplitudes through their symmetry relations
c      if (leg.eq.1) then
         tmpf = f
         tmpff = ff
         tmpsf=0
         if (leg.eq.1) tmpsf = dsib2getsfermion(ff, isf1)
c         tmpsf(1) = sff(1)
c         tmpsf(2) = sff(2)
         m11 = Mf
         m22 = Mff
c         Mtmpsf(1) = Msff(1)
c         Mtmpsf(2) = Msff(2)
         tmpE1J = E1J
         tmpE2J = E2J
         tmpq1 = kJ
         tmpppl = ppl
c      else 
       if (leg.eq.2) then
         tmpf = ff
         tmpff = f
         tmpsf = dsib2getsfermion(f, isf1)
c         tmpsf(1) = sf(1)
c         tmpsf(2) = sf(2)
         m11 = Mff
         m22 = Mf
c         Mtmpsf(1) = Msf(1)
c         Mtmpsf(2) = Msf(2)
         tmpE1J = E2J 
         tmpE2J = E1J
         tmpq1 = - kJ
         tmpppl = - ppl
      endif
      Mtmpsf = mass(tmpsf)/mx ! dimensionless in order to be consistent with kinematics



c... couplings
      CFSR1 = Conjg(gl(tmpsf,tmpff,kn(1)))*gl(kw,tmpf,tmpff)*
     -     gl(tmpsf,tmpff,kn(1))
      CFSR2 = Conjg(gl(tmpsf,tmpff,kn(1)))*gl(kw,tmpf,tmpff)*
     -     gr(tmpsf,tmpff,kn(1)) ! same as Conjg(gr(tmpsf,tmpff,kn(1)))*gl(kw,tmpf,tmpff)*gl(tmpsf,tmpff,kn(1))
      CFSR3 = Conjg(gr(tmpsf,tmpff,kn(1)))*gl(kw,tmpf,tmpff)*
     -     gr(tmpsf,tmpff,kn(1))

c... denominators
      delta = (
     -     dcmplx(m11**2 - m22**2 + MB**2 + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1,
     -     m22*width(tmpff)/mx)*
     -     dcmplx((m11**2 + m22**2 + MB**2 + 2*EvJ*tmpE1J - 
     -     2*EvJ*tmpE2J - 2*tmpE1J*tmpE2J - 4*CW*kv*tmpq1 - 2*tmpq1**2 - 
     -     4*Mtmpsf**2)/4.,Mtmpsf*width(tmpsf)/mx*0d0))

c... Helicity amplitudes (entries as for FSR1)
      amp00 = (
     -     ((0,-0.5)*(2*CFSR2*(Emi*kv*m11 - Emi*kv*m22-2*Epl*kv*tmpE1J + 
     -     CW*MB**2*tmpppl + CW*EvJ*(m11*pmi - m22*pmi + 2*Epl*tmpq1)) + 
     -     CFSR1*m22*(Emi*kv*m11 - Emi*kv*m22 - 2*Epl*kv*tmpE1J + 
     -     CW*MB**2*tmpppl + CW*EvJ*(m11*pmi - m22*pmi + 2*Epl*tmpq1)) + 
     -     CFSR3*(-(kv*(2*Epl*m22*tmpE1J + 
     -     Emi*(m11**2 - m11*m22 + MB**2 + 2*EvJ*tmpE1J))) + 
     -     2*CW**2*EvJ*kv*pmi*tmpq1 + 
     -     CW*(-(EvJ*pmi*(m11**2 + MB**2 + 2*EvJ*tmpE1J)) + 
     -     2*Emi*(EvJ**2 - MB**2)*tmpq1 + 
     -     m22*(EvJ*m11*pmi + MB**2*tmpppl + 2*EvJ*Epl*tmpq1)))))/MB)
     -     /delta

      amp01 = (
     -     (SW*(-((2*CFSR2 + (CFSR1 + CFSR3)*m22)*MB**2*(Emi - pmi)) - 
     -     2*CFSR3*EvJ**2*tmpE1J*(Epl + tmpppl) - 
     -     EvJ*(2*CFSR2*(Epl*(m11 + m22) + (-m11 + m22)*tmpppl) + 
     -     CFSR1*m22*(Epl*(m11 + m22) + (-m11 + m22)*tmpppl) + 
     -     CFSR3*(tmpppl*(m11**2 - m11*m22 + MB**2 - 2*CW*kv*tmpq1) + 
     -     Epl*(m11*(m11 + m22) + MB**2 - 2*CW*kv*tmpq1)))))
     -     /(2.*Sqrt(2.)*MB))/delta

      amp02 = (
     -     ((0,0.5)*(2*CFSR2*(kv*(m11 + m22)*pmi - 2*kv*tmpE1J*tmpppl + 
     -     CW*(Epl*MB**2 + EvJ*(Emi*(m11 + m22) + 2*tmpppl*tmpq1))) + 
     -     CFSR1*m22*(kv*(m11 + m22)*pmi - 2*kv*tmpE1J*tmpppl + 
     -     CW*(Epl*MB**2 + EvJ*(Emi*(m11 + m22) + 2*tmpppl*tmpq1))) + 
     -     CFSR3*(kv*(m11**2*pmi+m11*m22*pmi+MB**2*pmi+2*EvJ*pmi*tmpE1J- 
     -     2*m22*tmpE1J*tmpppl) - 2*CW**2*EvJ*Emi*kv*tmpq1 + 
     -     CW*(EvJ*Emi*m11*(m11 + m22) + Epl*m22*MB**2 + 
     -     MB**2*(EvJ*Emi + 2*pmi*tmpq1) + 
     -     2*EvJ*(EvJ*Emi*tmpE1J - EvJ*pmi*tmpq1 + m22*tmpppl*tmpq1)))))
     -     /MB)/delta

      amp03 = (
     -     -(SW*(2*CFSR2*(MB**2*(Emi + pmi) + 
     -     EvJ*(Epl*(m11 + m22) + (m11 - m22)*tmpppl)) + 
     -     CFSR1*m22*(MB**2*(Emi + pmi) + 
     -     EvJ*(Epl*(m11 + m22) + (m11 - m22)*tmpppl)) + 
     -     CFSR3*(m22*MB**2*(Emi + pmi) + EvJ*m11*m22*(Epl + tmpppl) + 
     -     EvJ*(Epl - tmpppl)*(m11**2 + MB**2 
     -     + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1))))
     -     /(2.*Sqrt(2.)*MB))/delta
      
      amp10 = (
     -     (SW*(CFSR1*m22*(-(m11*pmi) + m22*pmi - EvJ*tmpppl
     -     + kv*tmpppl - 2*Epl*tmpq1) - 2*CFSR2*(m11*pmi - m22*pmi
     -     + EvJ*tmpppl - kv*tmpppl + 2*Epl*tmpq1) + 
     -     CFSR3*(m11**2*pmi - m11*m22*pmi + 
     -     m22*(-(EvJ*tmpppl) + kv*tmpppl - 2*Epl*tmpq1) + 
     -     pmi*(MB**2 + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1))))
     -     /(2.*Sqrt(2.)))/delta

      amp11 = (
     -     (0,-0.25)*(1 + CW)*(2*CFSR2*(Epl*(m11 + m22) 
     -     - (-EvJ + kv)*(Emi - pmi) - m11*tmpppl + m22*tmpppl) + 
     -     CFSR1*m22*(Epl*(m11 + m22) - (-EvJ + kv)*(Emi - pmi) 
     -     - m11*tmpppl + m22*tmpppl) +
     -     CFSR3*(-((-EvJ + kv)*m22*(Emi - pmi)) + m11**2*tmpppl - 
     -     m11*m22*tmpppl + tmpppl*(MB**2 + 2*EvJ*tmpE1J  
     -     - 2*CW*kv*tmpq1) + Epl*(m11*(m11 + m22) + MB**2 
     -     + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1))))/delta

      amp12 = (
     -     (SW*(2*CFSR2*(EvJ*Epl - Epl*kv
     -     + Emi*(m11 + m22) + 2*tmpppl*tmpq1) + 
     -     CFSR1*m22*(EvJ*Epl-Epl*kv+Emi*(m11 + m22) + 2*tmpppl*tmpq1) + 
     -     CFSR3*(Emi*(m11*(m11 + m22) + MB**2 
     -     + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1) + 
     -     m22*(EvJ*Epl - Epl*kv + 2*tmpppl*tmpq1))))
     -     /(2.*Sqrt(2.)))/delta

      amp13 = (
     -     (0,-0.25)*(-1 + CW)*(2*CFSR2*(Epl*(m11 + m22)  
     -     - (-EvJ + kv)*(Emi + pmi) + m11*tmpppl - m22*tmpppl) + 
     -     CFSR1*m22*(Epl*(m11 + m22)-
     -     (-EvJ + kv)*(Emi + pmi) + m11*tmpppl - 
     -     m22*tmpppl)+CFSR3*((EvJ - kv)*m22*(Emi + pmi)-m11**2*tmpppl + 
     -     m11*m22*tmpppl-tmpppl*(MB**2 + 2*EvJ*tmpE1J-2*CW*kv*tmpq1) + 
     -     Epl*(m11*(m11 + m22) + MB**2 + 2*EvJ*tmpE1J-2*CW*kv*tmpq1)))
     -     )/delta

      amp20 = (
     -     -(SW*(2*CFSR2*(m11*pmi - m22*pmi 
     -     + (EvJ + kv)*tmpppl + 2*Epl*tmpq1) + 
     -     CFSR1*m22*(m11*pmi - m22*pmi+(EvJ + kv)*tmpppl+2*Epl*tmpq1) + 
     -     CFSR3*(-(m11**2*pmi) + m11*m22*pmi + 
     -     m22*((EvJ + kv)*tmpppl + 2*Epl*tmpq1) - 
     -     pmi*(MB**2 + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1))))
     -     /(2.*Sqrt(2.)))/delta

      amp21 = (
     -     (0,-0.25)*(-1 + CW)*(2*CFSR2*(Epl*(m11 + m22)
     -     + (EvJ + kv)*(Emi - pmi) - m11*tmpppl + m22*tmpppl) + 
     -     CFSR1*m22*(Epl*(m11 + m22)+(EvJ + kv)*
     -     (Emi - pmi)-m11*tmpppl + 
     -     m22*tmpppl)+CFSR3*((EvJ + kv)*m22*(Emi - pmi)+m11**2*tmpppl - 
     -     m11*m22*tmpppl+tmpppl*(MB**2 + 2*EvJ*tmpE1J-2*CW*kv*tmpq1) + 
     -     Epl*(m11*(m11 + m22) + MB**2 + 2*EvJ*tmpE1J-2*CW*kv*tmpq1)))
     -     )/delta
      
      amp22 = (
     -     (SW*(2*CFSR2*(Epl*(EvJ + kv)
     -     + Emi*(m11 + m22) + 2*tmpppl*tmpq1) + 
     -     CFSR1*m22*(Epl*(EvJ + kv) + Emi*(m11 + m22)+2*tmpppl*tmpq1) + 
     -     CFSR3*(Emi*(m11*(m11 + m22) + MB**2 
     -     + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1) + 
     -     m22*(Epl*(EvJ + kv) + 2*tmpppl*tmpq1))))/(2.*Sqrt(2.)))/delta
      
      amp23 = (
     -     (0,-0.25)*(1 + CW)*(2*CFSR2*(Epl*(m11 + m22) 
     -     + (EvJ + kv)*(Emi + pmi) + m11*tmpppl - m22*tmpppl) + 
     -     CFSR1*m22*(Epl*(m11 + m22)+
     -     (EvJ + kv)*(Emi + pmi) + m11*tmpppl - 
     -     m22*tmpppl)+CFSR3*((EvJ + kv)*m22*(Emi + pmi)-m11**2*tmpppl + 
     -     m11*m22*tmpppl-tmpppl*(MB**2 + 2*EvJ*tmpE1J-2*CW*kv*tmpq1) + 
     -     Epl*(m11*(m11 + m22) + MB**2 + 2*EvJ*tmpE1J-2*CW*kv*tmpq1)))
     -     )/delta
      
c... fermion = neutrinos
      if(f.eq.1.or.f.eq.3.or.f.eq.5) then
         if (leg.eq.1) then
c...  FSR1 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) + 2*amp00
            amp(0, 1) = amp(0, 1) + 2*amp01
            amp(0, 2) = amp(0, 2) + 2*amp02
            amp(0, 3) = amp(0, 3) + 2*amp03
            amp(1, 0) = amp(1, 0) + 2*amp10
            amp(1, 1) = amp(1, 1) + 2*amp11
            amp(1, 2) = amp(1, 2) + 2*amp12
            amp(1, 3) = amp(1, 3) + 2*amp13
            amp(-1, 0) = amp(-1, 0) + 2*amp20
            amp(-1, 1) = amp(-1, 1) + 2*amp21
            amp(-1, 2) = amp(-1, 2) + 2*amp22
            amp(-1, 3) = amp(-1, 3) + 2*amp23
         else if (leg.eq.2) then
c...  FSR2 (factor of 2 to get the t+u channel) 
c...  factor of 1/2. because of only one sneutrino exachanged (vs sum isf=1,2 in dsib2MsqffW)
            amp(0, 0) = amp(0, 0) - 2*amp00/2.
            amp(0, 1) = amp(0, 1) - 2*amp03/2.
            amp(0, 2) = amp(0, 2) - 2*amp02/2.
            amp(0, 3) = amp(0, 3) - 2*amp01/2.
            amp(1, 0) = amp(1, 0) - 2*amp20/2.
            amp(1, 1) = amp(1, 1) - 2*amp23/2.
            amp(1, 2) = amp(1, 2) - 2*amp22/2.
            amp(1, 3) = amp(1, 3) - 2*amp21/2.
            amp(-1, 0) = amp(-1, 0) - 2*amp10/2.
            amp(-1, 1) = amp(-1, 1) - 2*amp13/2.
            amp(-1, 2) = amp(-1, 2) - 2*amp12/2.
            amp(-1, 3) = amp(-1, 3) - 2*amp11/2.
         endif
c... fermion = leptons
      elseif(f.eq.2.or.f.eq.4.or.f.eq.6) then
         if (leg.eq.1) then
c...  FSR1 (factor of 2 to get the t+u channel)
c...  factor of 1/2. because of only one sneutrino exachanged (vs sum isf=1,2 in dsib2MsqffW)
            amp(0, 0) = amp(0, 0) + 2*amp00/2.
            amp(0, 1) = amp(0, 1) + 2*amp01/2.
            amp(0, 2) = amp(0, 2) + 2*amp02/2.
            amp(0, 3) = amp(0, 3) + 2*amp03/2.
            amp(1, 0) = amp(1, 0) + 2*amp10/2.
            amp(1, 1) = amp(1, 1) + 2*amp11/2.
            amp(1, 2) = amp(1, 2) + 2*amp12/2.
            amp(1, 3) = amp(1, 3) + 2*amp13/2.
            amp(-1, 0) = amp(-1, 0) + 2*amp20/2.
            amp(-1, 1) = amp(-1, 1) + 2*amp21/2.
            amp(-1, 2) = amp(-1, 2) + 2*amp22/2.
            amp(-1, 3) = amp(-1, 3) + 2*amp23/2.
         else if (leg.eq.2) then
c...  FSR2 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) - 2*amp00
            amp(0, 1) = amp(0, 1) - 2*amp03
            amp(0, 2) = amp(0, 2) - 2*amp02
            amp(0, 3) = amp(0, 3) - 2*amp01
            amp(1, 0) = amp(1, 0) - 2*amp20
            amp(1, 1) = amp(1, 1) - 2*amp23
            amp(1, 2) = amp(1, 2) - 2*amp22
            amp(1, 3) = amp(1, 3) - 2*amp21
            amp(-1, 0) = amp(-1, 0) - 2*amp10
            amp(-1, 1) = amp(-1, 1) - 2*amp13
            amp(-1, 2) = amp(-1, 2) - 2*amp12
            amp(-1, 3) = amp(-1, 3) - 2*amp11
         endif
c... quarks
      else
         if (leg.eq.1) then
c...  FSR1 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) + 2*amp00
            amp(0, 1) = amp(0, 1) + 2*amp01
            amp(0, 2) = amp(0, 2) + 2*amp02
            amp(0, 3) = amp(0, 3) + 2*amp03
            amp(1, 0) = amp(1, 0) + 2*amp10
            amp(1, 1) = amp(1, 1) + 2*amp11
            amp(1, 2) = amp(1, 2) + 2*amp12
            amp(1, 3) = amp(1, 3) + 2*amp13
            amp(-1, 0) = amp(-1, 0) + 2*amp20
            amp(-1, 1) = amp(-1, 1) + 2*amp21
            amp(-1, 2) = amp(-1, 2) + 2*amp22
            amp(-1, 3) = amp(-1, 3) + 2*amp23
         else if (leg.eq.2) then
c...  FSR2 (factor of 2 to get the t+u channel)
            amp(0, 0) = amp(0, 0) - 2*amp00
            amp(0, 1) = amp(0, 1) - 2*amp03
            amp(0, 2) = amp(0, 2) - 2*amp02
            amp(0, 3) = amp(0, 3) - 2*amp01
            amp(1, 0) = amp(1, 0) - 2*amp20
            amp(1, 1) = amp(1, 1) - 2*amp23
            amp(1, 2) = amp(1, 2) - 2*amp22
            amp(1, 3) = amp(1, 3) - 2*amp21
            amp(-1, 0) = amp(-1, 0) - 2*amp10
            amp(-1, 1) = amp(-1, 1) - 2*amp13
            amp(-1, 2) = amp(-1, 2) - 2*amp12
            amp(-1, 3) = amp(-1, 3) - 2*amp11
         endif 
      endif

      return
      end

