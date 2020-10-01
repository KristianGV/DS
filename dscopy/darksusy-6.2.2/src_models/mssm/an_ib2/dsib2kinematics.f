************************************************************************
*** subroutines dsib2kinematics writes relevant dimensionless kinematical
*** quantities to the common block (needed by helicity amplitude routines).
***
***    Input: QB                   - Boson is neutral (0) or charged (1)
***           E1 \equiv x1 = E1/mx - CMS fermion energy
***           EB \equiv xB = Eb/mx - CMS boson energy
***
***    Output: istat - error code (0 if everything OK)
***
*** Author: Francesca Calore, 2013-04-10
***         2014-02-20 (Torsten Bringmann, adapted to include scalars)
************************************************************************
      subroutine dsib2kinematics(QB, EB, E1, istat)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      real*8 EB,E1
      integer QB, istat
      real*8 E2
c-----------------------------------------------------------------------
      
      istat=0

      E2 = 2.d0 - E1 - EB

c...  BB system: dimensionless quantities
      if (QB.eq.0) then
         E1J = Sqrt(MB**2/4. + E1 + E2 - 1.d0)
         E2J = E1J
         EvJ = 1.d0/E1J - E1J - MB**2/(4.*E1J)
      elseif(QB.eq.1) then
         if (Mf.eq.0.and.Mff.eq.0) then ! massless case
            E1J = Sqrt(MB**2/4. + E1 + E2 - 1.d0)
            E2J = E1J
            EvJ = 1.d0/E1J - E1J - MB**2/(4.*E1J)
            Epl = E1J
            Emi = 0.d0
            ppl = E1J
            pmi = 0.d0
         else 
            E1J = (E1 + E2 - 1.d0 + (MB**2 - Mff**2 + Mf**2)/4.)
     -           /Sqrt(E1 + E2 - 1.d0 + MB**2/4.)
            E2J = (E1 + E2 - 1.d0 + (MB**2 - Mf**2 + Mff**2)/4.)
     -           /Sqrt(E1 + E2 - 1.d0 + MB**2/4.)
            EvJ = (2.d0 - E1 - E2 - MB**2/2.)
     -           /Sqrt(E1 + E2 - 1. + MB**2/4.)
            Epl = Sqrt((E1J*E2J + Mf*Mff + (E1J**2 - Mf**2))/2.)
            Emi = Sqrt((E1J*E2J + Mf*Mff - (E1J**2 - Mf**2))/2.)
            ppl = Sqrt((E1J*E2J - Mf*Mff + (E1J**2 - Mf**2))/2.)
            if ((Mf*E2J - Mff*E1J).gt.0.d0) then
               pmi = Sqrt((E1J*E2J - Mf*Mff - (E1J**2 - Mf**2))/2.)
            else if ((Mf*E2J - Mff*E1J).lt.0d0) then
               pmi = - Sqrt((E1J*E2J - Mf*Mff - (E1J**2 - Mf**2))/2.)
            else 
               write(*,*) 'WARNING in dsib2kinematics: 
     &              sign of pmi indeterminate', (Mf*E2J - Mff*E1J)
               istat=-1
            endif
         endif
      else
         write(*,*)'WARNING in dsib2kinematics: invalid value of QB = ',QB
         istat=-1
      endif

      kJ = Sqrt(E1J**2 - Mf**2)
      kv = Sqrt(EvJ**2 - MB**2)
    
      if (QB.eq.0) then
         CW = (E2 - E1)/(kv*kJ)               !cos (fermion,vector)
      elseif(QB.eq.1) then
         CW = (-2.*(1. - E2) + (Mf**2 + MB**2 - Mff**2)/2. + E1J*EvJ)
     -        /(kv*kJ)
      endif

      SW = Sqrt(1. - CW**2)

      if (CW.gt.1.d0.or.CW.lt.-1.d0) then
        istat=-1
        write(*,*) 'ERROR in dsib2kinematics:
     -              wrong kinematical boundaries ', E1, EB, CW 
      endif
      
      return
      end
