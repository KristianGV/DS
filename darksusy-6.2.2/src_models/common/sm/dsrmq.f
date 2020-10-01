      real*8 function dsrmq(mscale,kpart)
      implicit none
      include 'dsparticles.h'
      include 'dssm.h'

c... restructured by Torsten Bringmann, 2015-10-30
c... TB: stripped off MSSM part, 11/2016


      real*8 mscale,dsrmq1loop,dsrmq4loop
      integer kpart


c... dsrmq should no longer stop but simply return not running masses
c... TB: added kb to list
      if (kpart.ne.kc.and.kpart.ne.kb.and.kpart.ne.kt.and.kpart.ne.ktau) then
c      if (kpart.ne.kc.and.kpart.ne.kt.and.kpart.ne.ktau) then
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
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsrmq'
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end
