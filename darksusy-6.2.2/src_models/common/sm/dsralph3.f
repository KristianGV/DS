      real*8 function dsralph3(mscale)

c... TB: stripped off MSSM part: 11/2016

      implicit none
      include 'dsparticles.h'
      include 'dssm.h'
      include 'dsmpconst.h'

      real*8 mscale,dsralph31loop,dsralph34loop

      if(roption.eq.'norun') then  
        dsralph3 = alph3mz
      elseif(roption.eq.'1loop') then
        dsralph3 = dsralph31loop(mscale) 
      elseif(roption.eq.'4loop') then
        dsralph3 = dsralph34loop(mscale) 
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsralph3'
        write(*,*) 'program stopped'
        stop
      endif  
      return
      end
