************************************************************************
*** subroutine dsIB2chinit initializes masses and other relevant 
*** parameters for a given IB2 channel and writes them to common block 
*** variables in dsib2com.h
***
***   Input: IB2ch  - 3-body channel,
***                   see dsIByieldone or dsib2sigmav for details
***
***   Output: err=1 for unknown channel IB2ch, err=0 otherwise
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2015-11-19
************************************************************************

      subroutine dsib2chinit(IB2ch,err)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer IB2ch, err
c-----------------------------------------------------------------------
c NB: ALL masses and energies are dimensionless in ib2 routines!!! 
c     Only exception: mx (everywhere), egev (input for yield routines)  
c     and Efinal (input for dsib2dnde)   
c-----------------------------------------------------------------------
      
      err = 0
      
      mx = mass(kn(1))

      c_ftype = mod(IB2ch,100) ! Fermion type (for internal use)
      c_Btype = IB2ch/100      ! Boson type (for internal use)

      if (c_ftype.lt.1.or.c_ftype.gt.12) goto 1000  
      if (c_Btype.lt.1.or.c_Btype.gt.6) goto 1000  

c... determine ANTIfermion type
      c_fbartype=c_ftype                     ! neutral boson final states
      if (c_Btype.eq.2.or.c_Btype.eq.6) then ! charged boson final states
        if (mod(c_ftype,2).EQ.0) c_fbartype = c_fbartype-1
        if (mod(c_ftype,2).EQ.1) c_fbartype = c_fbartype+1      
      endif

      Mf = mass(c_ftype)/mx      
      Mff = mass(c_fbartype)/mx      

      if (c_Btype.eq.1) mB = mass(kz)/mx     ! Z
      if (c_Btype.eq.2) mB = mass(kw)/mx     ! W
      if (c_Btype.eq.3) mB = mass(kh2)/mx    ! SM Higgs
      if (c_Btype.eq.4) mB = mass(kh1)/mx    ! heavy Higgs
      if (c_Btype.eq.5) mB = mass(kh3)/mx    ! A0
      if (c_Btype.eq.6) mB = mass(khc)/mx    ! charged Higgs
     
      return

c... channel not implemented
 1000 err=1
      write(*,*) 'ERROR in dsib2chinit: IB2ch not recognized!'
      return

      end

