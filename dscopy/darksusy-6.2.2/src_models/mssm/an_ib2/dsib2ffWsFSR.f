*******************************************************************
*** subroutine dsib2ffWsFSR add to the common block array amp(i,j) 
*** the helicity amplitudes of the s channel FSR diagrams. 
*** The W is emitted from the final legs of the s channel diagrams.
*** H3 and Z diagrams are considered separately and then
*** added to the corresponding entries of amp(-1:1, 0:3). 
*** 
*** Author: Francesca Calore, 2013-04-10
***         Francesca Calore, 2015-24-11 (modified ff)
*******************************************************************
      subroutine dsib2ffWsFSR(f, leg) 
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, leg
      integer ff

      integer tmpf, tmpff
      real*8 m11, m22, tmpE1J, tmpE2J, tmpq1, tmpppl
      real*8 Mz, Mh3
      complex*16 CH3, CZ1, CZ2
      complex*16 deltaH3, deltaZ
      complex*16 ampH00, ampH01, ampH02, ampH03, 
     & ampH10, ampH11, ampH12, ampH13,
     &     ampH20, ampH21, ampH22, ampH23
      complex*16 ampZ00, ampZ01, ampZ02, ampZ03, 
     &     ampZ10, ampZ11, ampZ12, ampZ13,
     &     ampZ20, ampZ21, ampZ22, ampZ23

      ff = c_fbartype

c... Defines FSR1 and FSR2 amplitudes through their simmetry relations
c      if (leg.eq.1) then
         tmpf = f
         tmpff = ff
         m11 = Mf
         m22 = Mff
         tmpE1J = E1J
         tmpE2J = E2J
         tmpq1 = kJ
         tmpppl = ppl
c      else 
       if (leg.eq.2) then
         tmpf = ff
         tmpff = f
         m11 = Mff
         m22 = Mf
         tmpE1J = E2J 
         tmpE2J = E1J
         tmpq1 = - kJ
         tmpppl = - ppl
      endif

c... masses/mx exchanged particles
      Mz = mass(kz)/mx
      Mh3 = mass(kh3)/mx

c... couplings
c... H3      
      CH3 = gl(kw,tmpf,tmpff)*gr(kh3,tmpff,tmpff)*
     -     (Conjg(gr(kh3,kn(1),kn(1))) - gr(kh3,kn(1),kn(1)))

c... Z
      CZ1 = gl(kw,tmpf,tmpff)*gl(kz,tmpff,tmpff)*
     -     Dble(gl(kz,kn(1),kn(1)))
      CZ2 = gl(kw,tmpf,tmpff)*gr(kz,tmpff,tmpff)*
     -     Dble(gl(kz,kn(1),kn(1)))
c... denominators
c... H3
      deltaH3 = (
     -     dcmplx(m11**2 - m22**2 + MB**2 + 2*EvJ*tmpE1J - 2*CW*kv*tmpq1,
     -     m22*width(tmpff)/mx)*dcmplx(m11**2 + m22**2 - Mh3**2
     -     + MB**2 + 2*EvJ*tmpE1J + 2*EvJ*tmpE2J + 2*tmpE1J*tmpE2J 
     -     + 2*tmpq1**2,Mh3*width(kh3)/mx))
c... Z
      deltaZ = (dcmplx(m11**2 - m22**2 + MB**2 + 2*EvJ*tmpE1J 
     -     - 2*CW*kv*tmpq1, m22*width(tmpff)/mx)*
     -     dcmplx(m11**2 + m22**2 - Mz**2 + MB**2 + 2*EvJ*tmpE1J + 
     -     2*EvJ*tmpE2J + 2*tmpE1J*tmpE2J + 2*tmpq1**2,Mz*width(kz)/mx)) 


c... helicity amplitudes
c... H3 
      ampH00 = (
     -     ((0,-2)*CH3*(Emi*kv*m11 - Emi*kv*m22 
     -     - 2*Epl*kv*tmpE1J + CW*MB**2*tmpppl + 
     -     CW*EvJ*(m11*pmi - m22*pmi + 2*Epl*tmpq1)))/MB)/deltaH3

      ampH01 = (
     -     -((Sqrt(2.)*CH3*SW*(MB**2*(Emi - pmi) + 
     -     EvJ*(Epl*(m11 + m22) + (-m11 + m22)*tmpppl)))/MB))/deltaH3

      ampH02 = (
     -     ((0,2)*CH3*(kv*(m11 + m22)*pmi - 2*kv*tmpE1J*tmpppl + 
     -     CW*(Epl*MB**2 + EvJ*(Emi*(m11 + m22) + 2*tmpppl*tmpq1)))
     -     )/MB)/deltaH3

      ampH03  = (
     -     -((Sqrt(2.)*CH3*SW*(MB**2*(Emi + pmi) + 
     -     EvJ*(Epl*(m11 + m22) + (m11 - m22)*tmpppl)))/MB))/deltaH3

      ampH10 = (
     -     Sqrt(2.)*CH3*SW*(-(m11*pmi) + m22*pmi - EvJ*tmpppl + kv*tmpppl
     -     - 2*Epl*tmpq1))/deltaH3

      ampH11 = (
     -     (0,-1)*CH3*(1 + CW)*(Epl*(m11 + m22) - (-EvJ + kv)*(Emi - pmi) 
     -     - m11*tmpppl + m22*tmpppl))/deltaH3

      ampH12 = (
     -     Sqrt(2.)*CH3*SW*(Epl*(EvJ - kv) + Emi*(m11 + m22) 
     -     + 2*tmpppl*tmpq1))/deltaH3

      ampH13 = (
     -     (0,-1)*CH3*(-1 + CW)*(Epl*(m11 + m22) - 
     -     (-EvJ + kv)*(Emi + pmi) + m11*tmpppl - m22*tmpppl))/deltaH3

      ampH20 = (
     -     -(Sqrt(2.)*CH3*SW*(m11*pmi - m22*pmi + (EvJ + kv)*tmpppl 
     -     + 2*Epl*tmpq1)))/deltaH3

      ampH21 = (
     -     (0,-1)*CH3*(-1 + CW)*(Epl*(m11 + m22) + (EvJ + kv)*(Emi - pmi)
     -     - m11*tmpppl + m22*tmpppl))/deltaH3

      ampH22 = (
     -     Sqrt(2.)*CH3*SW*(Epl*(EvJ + kv) + Emi*(m11 + m22)
     -     + 2*tmpppl*tmpq1))/deltaH3

      ampH23 = (
     -     (0,-1)*CH3*(1 + CW)*(Epl*(m11 + m22) + (EvJ + kv)*(Emi + pmi)
     -     + m11*tmpppl - m22*tmpppl))/deltaH3

c... Z
      ampZ00 = (((0,2)*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(-(Epl*kv*(tmpE1J + tmpE2J)) + 
     -     CW*(MB**2 + EvJ*(tmpE1J + tmpE2J))*tmpppl) + 
     -     CZ1*(kv*(-(Epl*m11*(tmpE1J + tmpE2J)) + 
     -     Emi*(MB**2 + 2*tmpE1J*(EvJ + tmpE1J + tmpE2J))) - 
     -     2*CW**2*EvJ*kv*pmi*tmpq1 + 
     -     CW*(MB**2*pmi*(EvJ - tmpE1J + tmpE2J) + 
     -     m11*(MB**2 + EvJ*(tmpE1J + tmpE2J))*tmpppl - 
     -     2*EvJ*(-(EvJ*pmi*tmpE1J) + 
     -     Emi*(EvJ + tmpE1J + tmpE2J)*tmpq1)))))/(MB*Mz**2))/deltaZ
      
      ampZ01 =(-((Sqrt(2.)*SW*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (-(CZ2*m22*(Emi - pmi)*(MB**2 + EvJ*(tmpE1J + tmpE2J))) + 
     -     CZ1*(2*EvJ**2*tmpE1J*(Epl + tmpppl) + 
     -     MB**2*(m11*(Emi + pmi) - 
     -     (tmpE1J - tmpE2J)*(Epl + tmpppl)) + 
     -     EvJ*(m11*(Emi + pmi)*(tmpE1J + tmpE2J) + 
     -     (Epl + tmpppl)*(MB**2 - 2*CW*kv*tmpq1)))))/
     -     (MB*Mz**2)))/deltaZ
      
      ampZ02 = (((0,2)*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (-(CZ2*m22*(CW*Epl*(MB**2 + EvJ*(tmpE1J + tmpE2J)) - 
     -     kv*(tmpE1J + tmpE2J)*tmpppl)) + 
     -     CZ1*(kv*(MB**2*pmi + 
     -     2*pmi*tmpE1J*(EvJ + tmpE1J + tmpE2J) - 
     -     m11*(tmpE1J + tmpE2J)*tmpppl) - 
     -     2*CW**2*Emi*EvJ*kv*tmpq1 + 
     -     CW*(Emi*(2*EvJ**2*tmpE1J + 
     -     MB**2*(EvJ - tmpE1J + tmpE2J)) + 
     -     Epl*m11*(MB**2 + EvJ*(tmpE1J + tmpE2J)) - 
     -     2*EvJ*pmi*(EvJ + tmpE1J + tmpE2J)*tmpq1))))/(MB*Mz**2))/deltaZ
      
      ampZ03 = ((Sqrt(2.)*SW*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(Emi + pmi)*(MB**2 + EvJ*(tmpE1J + tmpE2J)) + 
     -     CZ1*(-(Emi*m11*(MB**2 + EvJ*(tmpE1J + tmpE2J))) + 
     -     m11*pmi*(MB**2 + EvJ*(tmpE1J + tmpE2J)) + 
     -     (Epl - tmpppl)*
     -     (MB**2*(-EvJ + tmpE1J - tmpE2J) + 
     -     2*EvJ*(-(EvJ*tmpE1J) + CW*kv*tmpq1)))))/(MB*Mz**2))/deltaZ
      
      ampZ10 = ((Sqrt(2.)*SW*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(EvJ - kv + tmpE1J + tmpE2J)*tmpppl + 
     -     CZ1*(MB**2*pmi + 
     -     m11*(EvJ + kv + tmpE1J + tmpE2J)*tmpppl - 
     -     2*CW*kv*pmi*tmpq1 + 
     -     EvJ*(pmi*(tmpE1J + tmpE2J) - 2*Emi*tmpq1) - 
     -     (tmpE1J + tmpE2J)*(kv*pmi + 2*Emi*tmpq1))))/Mz**2)/deltaZ
      
      ampZ11 = (((0,1)*(1 + CW)*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(Emi - pmi)*(EvJ - kv + tmpE1J + tmpE2J) - 
     -     CZ1*(m11*(Emi + pmi)*(EvJ + kv + tmpE1J + tmpE2J) + 
     -     (Epl + tmpppl)*
     -     (MB**2 + EvJ*(tmpE1J + tmpE2J) - 
     -     kv*(tmpE1J + tmpE2J + 2*(-1 + CW)*tmpq1)))))/Mz**2)/deltaZ
      
      ampZ12 = ((Sqrt(2.)*SW*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*Epl*m22*(-EvJ + kv - tmpE1J - tmpE2J) + 
     -     CZ1*(Epl*m11*(EvJ + kv + tmpE1J + tmpE2J) - 
     -     2*pmi*(EvJ + tmpE1J + tmpE2J)*tmpq1 + 
     -     Emi*(MB**2 + EvJ*(tmpE1J + tmpE2J) - 
     -     kv*(tmpE1J + tmpE2J + 2*CW*tmpq1)))))/Mz**2)/deltaZ
      
      ampZ13 =(((0,1)*(-1 + CW)*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(Emi + pmi)*(EvJ - kv + tmpE1J + tmpE2J) + 
     -     CZ1*(-(m11*(Emi - pmi)*(EvJ + kv + tmpE1J + tmpE2J)) - 
     -     (Epl - tmpppl)*
     -     (MB**2 + EvJ*(tmpE1J + tmpE2J) - 
     -     kv*(tmpE1J + tmpE2J + 2*(1 + CW)*tmpq1)))))/Mz**2)/deltaZ
      
      ampZ20 = ((Sqrt(2.)*SW*(CZ2*m22*(EvJ + kv + tmpE1J +tmpE2J)*
     -     tmpppl + CZ1*(MB**2*pmi + EvJ*pmi*tmpE1J + kv*pmi*tmpE1J + 
     -     EvJ*pmi*tmpE2J + kv*pmi*tmpE2J + 
     -     m11*(EvJ - kv + tmpE1J + tmpE2J)*tmpppl - 
     -     2*Emi*EvJ*tmpq1 - 2*CW*kv*pmi*tmpq1 - 
     -     2*Emi*tmpE1J*tmpq1 - 2*Emi*tmpE2J*tmpq1))*
     -     (m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J)
     -     + tmpq1**2)))/Mz**2)/deltaZ
      
      ampZ21 =(((0,1)*(-1 + CW)*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(Emi - pmi)*(EvJ + kv + tmpE1J + tmpE2J) - 
     -     CZ1*(m11*(Emi + pmi)*(EvJ - kv + tmpE1J + tmpE2J) + 
     -     (Epl + tmpppl)*
     -     (MB**2 + (EvJ + kv)*(tmpE1J + tmpE2J) - 
     -     2*(1 + CW)*kv*tmpq1))))/Mz**2)/deltaZ
      
      ampZ22 = ((Sqrt(2.)*SW*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (-(CZ2*Epl*m22*(EvJ + kv + tmpE1J + tmpE2J)) + 
     -     CZ1*(Epl*m11*(EvJ - kv + tmpE1J + tmpE2J) - 
     -     2*pmi*(EvJ + tmpE1J + tmpE2J)*tmpq1 + 
     -     Emi*(MB**2 + (EvJ + kv)*(tmpE1J + tmpE2J) - 
     -     2*CW*kv*tmpq1))))/Mz**2)/deltaZ
      
      ampZ23 = (((0,1)*(1 + CW)*(m11**2 + m22**2 + MB**2 - Mz**2 + 
     -     2*(tmpE1J*tmpE2J + EvJ*(tmpE1J + tmpE2J) + tmpq1**2))*
     -     (CZ2*m22*(Emi + pmi)*(EvJ + kv + tmpE1J + tmpE2J) + 
     -     CZ1*(m11*(Emi - pmi)*(-EvJ + kv - tmpE1J - tmpE2J) - 
     -     (Epl - tmpppl)*
     -     (MB**2 + (EvJ + kv)*(tmpE1J + tmpE2J) - 
     -     2*(-1 + CW)*kv*tmpq1))))/Mz**2)/deltaZ

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
