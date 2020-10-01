      real*8 function dsrmq(mscale,kpart)
      implicit none
      include 'dsmssm.h'
      include 'dsisasugra.h'

c... restructured by Torsten Bringmann, 2015-10-30


      real*8 mscale,dsrmq1loop,dsrmq4loop,cosbeta
      integer kpart


c... dsrmq should no longer stop but simply return not running masses
c... TB debug: kb should be added to list -- but this leads to presently
c...           not understood problems in the RD computations. Needs addressing !!!
c      if (kpart.ne.kc.and.kpart.ne.kb.and.kpart.ne.kt.and.kpart.ne.ktau) then
      if (kpart.ne.kc.and.kpart.ne.kt.and.kpart.ne.ktau) then
         dsrmq=mass(kpart)
c          write(*,*) 'dsrmq called for wrong particle'
c          write(*,*) 'particle =  ',pname(kpart)
c          write(*,*) 'rather than tau or c, b, t quark'
c          write(*,*) 'program stopped'
c          stop
         return
      endif       

      dsrmq=0.0d0
      if(roption.eq.'norun') then  
         dsrmq=mass(kpart)
      elseif(roption.eq.'1loop') then
        if(kpart.eq.kc.or.kpart.eq.kb.or.kpart.eq.kt) then
          dsrmq=dsrmq1loop(mscale,kpart) 
        elseif(kpart.eq.ktau) then
          dsrmq=mass(ktau)
        endif
      elseif(roption.eq.'4loop') then
        if(kpart.eq.kc.or.kpart.eq.kb.or.kpart.eq.kt) then
          dsrmq=dsrmq4loop(mscale,kpart) 
        elseif(kpart.eq.ktau) then
          dsrmq=mass(ktau)
        endif
      elseif(roption.eq.'isasu') then  
        cosbeta=1.0d0/sqrt(1.0d0+tanbe*tanbe)
        if(kpart.eq.kc) then
          dsrmq=mass(kc)
        elseif(kpart.eq.kb) then
          dsrmq=gss(5)*cosbeta*vev
        elseif(kpart.eq.kt) then
          dsrmq=gss(6)*tanbe*cosbeta*vev
        elseif(kpart.eq.ktau) then
          dsrmq=gss(4)*cosbeta*vev
        endif  
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsrmq'
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end
