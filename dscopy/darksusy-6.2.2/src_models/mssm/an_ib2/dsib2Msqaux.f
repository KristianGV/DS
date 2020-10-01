************************************************************************
*** Function dsib2Msqaux returns the 3-body amplitude squared, as a 
*** function of the energy of one of the final state particles. 
*** Channel and final state particle set in common blocks
*** (auxiliary function needed for integration routines, do not call
*** directly!)
***
*** Author: Torsten Bringmann, 2014-02-20
************************************************************************
      real*8 function dsib2Msqaux(Efinal)
      implicit none
      include 'dsib2com.h'
      real*8 Efinal          ! energy/mx in the CMS

      real*8 EB, E1, sqmatr, Efsav, Mfsav, mfbarsav, E1sav,E2sav
      integer i, j, istat, cfsav, cfbarsav
      character*5 cpsav
      logical isCPsymmetric
c-----------------------------------------------------------------------

      sqmatr=0.0d0  
      dsib2Msqaux=0.0d0                                 
      isCPsymmetric=.true.

c... For charged bosons, assume CP symmetry and map everything to the case where 
c... the emitted boson is negatively charged      
      cfsav=c_ftype
      cfbarsav=c_fbartype
      cpsav=c_pfinal
      Efsav=Efinal
      E1sav=c_E1
      E2sav=c_E2      
      Mfsav=Mf
      mfbarsav=Mff

      if (isCPsymmetric.and.
     &   ((c_Btype.eq.2.or.c_Btype.eq.6).and.(mod(c_ftype,2).eq.0))) then
          c_ftype=c_ftype-1
          c_fbartype=c_fbartype+1
          Mf=mfbarsav
          Mff=Mfsav
          if (c_pfinal.eq.'B') Efinal=2.d0-c_EB-Efsav
          if (cpsav.eq.'fbar') then
            c_E1=E2sav 
            c_pfinal='f'
          endif  
         if (cpsav.eq.'f') then
            c_E2=E1sav 
            c_pfinal='fbar'
          endif  
      endif

  
c... assign energies of Boson and Fermion final states
      if (c_pfinal.eq.'B') then
        EB=c_EB
        E1=Efinal
      elseif (c_pfinal.eq.'fbar') then
        EB=Efinal
        E1=2.d0-c_E2-Efinal
      elseif (c_pfinal.eq.'f') then
        EB=Efinal
        E1=c_E1
      endif
      if (c_Btype.eq.1.or.c_Btype.eq.3.or.c_Btype.eq.4.or.c_Btype.eq.5)
     +    call dsib2kinematics(0, EB, E1, istat)      ! neutral Boson
      if (c_Btype.eq.2.or.c_Btype.eq.6)
     +    call dsib2kinematics(1, EB, E1, istat)      ! charged Boson
      if (istat.ne.0) return  


c... initialize & calculate helicity amplitudes
      do i = -1,1
         do j = 0,3
            amp(i,j) = 0d0
         end do 
      end do
C      write(*,*) 'ffZ'
      if (c_Btype.eq.1) call dsib2ffZamps(c_ftype)
C      write(*,*) 'ffW'
      if (c_Btype.eq.2) call dsib2ffWamps(c_ftype)
C      write(*,*) 'ffH01'
      if (c_Btype.eq.3) call dsib2ffH0amps(c_ftype,1)
C      write(*,*) 'ffH02'
      if (c_Btype.eq.4) call dsib2ffH0amps(c_ftype,2)
C      write(*,*) 'ffA0'
      if (c_Btype.eq.5) call dsib2ffA0amps(c_ftype)
C      write(*,*) 'ffHp'
      if (c_Btype.eq.6) call dsib2ffHPamps(c_ftype)

c... sum all helicity amplitudes
      do i = -1,1         ! NB: only i=0 is different from zero for Scalars...
         do j = 0,3
            sqmatr = sqmatr + dble(amp(i, j)*Conjg(amp(i, j))) 
         end do
      end do

c... spin average and color factors
      if(c_ftype.ge.7.and.c_ftype.le.12) sqmatr=3*sqmatr
      dsib2Msqaux = 1.d20*sqmatr/4.    ! NB: result is dimensionless; 
                                       ! divide by mx**2 to get physical result.
                                       ! Factor of 1.d20 for numerical convergence
 
 
c...restore orginal CP state
      if ((c_Btype.eq.2.or.c_Btype.eq.6).and.(mod(cfsav,2).eq.0)) then
        c_ftype=cfsav
        c_fbartype=cfbarsav
        c_pfinal=cpsav
        Efinal=Efsav
        Mf=Mfsav
        Mff=mfbarsav
      endif
 

      return
      end
