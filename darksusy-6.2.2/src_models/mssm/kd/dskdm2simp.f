***********************************************************************
*** dskdm2simp returns the scattering amplitude for NEUTRALINO DM squared 
*** at zero momentum transfer, SUMMED over both initial and final spin and 
*** other internal states, and then divided by the DM internal degrees of 
*** freedom. This is only valid in the limit of relativistic scattering 
*** partners with small energies omega, where we can expand as
***
***    |M|**2 = cn*(omega/m0)**n + O( (omega/m0)**(n+1) )
***                                                                         
***  type : interface                                                       
***                                                                         
***  input: SMtype   - SM scattering partners:
***                    4,5,6 - e,m,t leptons
***                    1,2,3 - e,m,t neutrinos
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
*** updates: 2013-06-12 (made model-dependence explicit)
***          2016-12-21 (changed spin-average prescription)
***********************************************************************

      subroutine dskdm2simp(SMtype,cn,n)
      implicit none

      include 'dsmssm.h'
      include 'dskdcom.h'
      include 'dsidtag.h'

      integer n, SMtype,kf,nsf,ksf(2)
      real*8  cn
      integer i,j

      real*8 tmpres,sv,dssigmav0tot

      character*20 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/

c... This internal consistency check makes sure that the correct particle module 
c... is loaded, and should (at least) be included for all interface functions
      call dscheckmodule('MSSM','dskdm2simp')


      cn=0d0
      n=0
      tmpres=0d0



        sv=dssigmav0tot()  ! make sure all vertices 
                                        ! and masses are calculated...

c... set up sparticles
        if (SMtype.le.3) then
          nsf=1
          kf=knu(SMtype)
        endif
        if (SMtype.ge.4) then
          nsf=2
          kf=kl(SMtype-3)
        endif
        if (SMtype.ge.7) return
        do 10 i=1,nsf
          if (SMtype.le.3) ksf(i)=ksnu(SMtype)
          if (SMtype.ge.4) ksf(i)=ksl(SMtype+3*i-6)
  10    continue

        n=2

          tmpres=8*(m0/mass(kz))**4
     -            *abs(gl(kz,kn(1),kn(1)))**2
     -            *(abs(gl(kz,kf,kf))**2
     -              +abs(gr(kz,kf,kf))**2)
          do 50 i=1,nsf
           do 40 j=1,nsf
            tmpres=tmpres+2*
     -     (abs(gl(ksf(i),kf,kn(1)))**2
     -      *abs(gl(ksf(j),kf,kn(1)))**2
     -      +abs(gr(ksf(i),kf,kn(1)))**2
     -      *abs(gr(ksf(j),kf,kn(1)))**2
     -      +4*dble(conjg(gl(ksf(i),kf,kn(1)))
     -               *gr(ksf(i),kf,kn(1))
     -               *gl(ksf(j),kf,kn(1))
     -               *conjg(gr(ksf(j),kf,kn(1)))))
     -     /((mass(ksf(i))/m0)**2-1d0)
     -     /((mass(ksf(j))/m0)**2-1d0)
   40      continue
   50     continue
          cn=2*nsf*tmpres 


      if (cn.lt.0d0) then
        if (memory.ne.idtag) then
          write(*,*) 'WARNING: negative |M|^2 in dskdm2simp.f for model'
     &               ,idtag, ' !'
c          write(*,*) 'SMtype, cn = ',SMtype,cn
          memory=idtag
        endif
        cn=0d0
      endif

      return

      end





        

