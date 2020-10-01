********************************************************************************
*** subroutine dsib2getsfermion returns the sfermions corresponding to a     *** 
*** given fermion f. If the global parameter IB2sfgen (set in dsib2set)      ***
*** equals 1, only the two sfermions states are returned that have the       ***
*** largest flavour component in common with f. Otherwise (for IB2sfgen=3),  ***
*** all sfermions with identical charges are returned (without assuming      ***
*** anything about their flavour ordering).                                  ***
***                                                                          ***
***      input: f - particle code of fermion                                 ***
***             i - index of sfermion to be returned                         ***
***                                                                          ***
***      output: particle code of sfermion number i                          ***
***                                                                          ***
*** Author: Torsten Bringmann, 2019-05-31                                    ***
********************************************************************************
      integer function dsib2getsfermion(f, i)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'

      integer f, i
      
      dsib2getsfermion = 0
      
c      write(*,*) 'dsib2getsfermion : ',f,i, IB2sfgen
      
      if (IB2sfgen.eq.3) then
        if (i.lt.1.or.i.gt.6) goto 100
        if ((f.eq.knue.or.f.eq.knumu.or.f.eq.knutau).and.i.le.3) then
          dsib2getsfermion = ksnu(i)
        elseif (f.eq.ke.or.f.eq.kmu.or.f.eq.ktau) then
          dsib2getsfermion = ksl(i)
        elseif (f.eq.ku.or.f.eq.kc.or.f.eq.kt) then
          dsib2getsfermion = ksqu(i)
        elseif (f.eq.kd.or.f.eq.ks.or.f.eq.kb) then
          dsib2getsfermion = ksqd(i)
        else
          goto 100 
        endif
      elseif (IB2sfgen.eq.1) then
        if (i.lt.1.or.i.gt.2.or.f.lt.1.or.f.gt.12) goto 100
        if (f.eq.knue)   dsib2getsfermion = ksnu_flav(1,1)
        if (f.eq.knumu)  dsib2getsfermion = ksnu_flav(2,1)
        if (f.eq.knutau) dsib2getsfermion = ksnu_flav(3,1)
        if (f.eq.ke)     dsib2getsfermion = ksl_flav(1,i)
        if (f.eq.kmu)    dsib2getsfermion = ksl_flav(2,i)
        if (f.eq.ktau)   dsib2getsfermion = ksl_flav(3,i)
        if (f.eq.ku)     dsib2getsfermion = ksqu_flav(1,i)
        if (f.eq.kc)     dsib2getsfermion = ksqu_flav(2,i)
        if (f.eq.kt)     dsib2getsfermion = ksqu_flav(3,i)
        if (f.eq.kd)     dsib2getsfermion = ksqd_flav(1,i)
        if (f.eq.ks)     dsib2getsfermion = ksqd_flav(2,i)
        if (f.eq.kb)     dsib2getsfermion = ksqd_flav(3,i)
      else
        write(*,*) 'ERROR in IB2 routines: Common block variable IB2sfgen'
        write(*,*) 'should equal 1 or 3, but is set to IB2sfgen = ',IB2sfgen
        write(*,*) 'Stopping program...'
        stop
      endif 

c      write(*,*) 'dsib2getsfermion = ',dsib2getsfermion 
      
      return
 
 100  write(*,*) 'ERROR in dsib2getsfermion!'
      write(*,*) 'Unsupported set of input parameters:'
      write(*,*) 'f, i = ', f, i
      write(*,*) 'Stopping program...'
      stop

      end
