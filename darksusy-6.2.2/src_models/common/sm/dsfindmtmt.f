      subroutine dsfindmtmt
      implicit none
      include 'dsparticles.h'
      include 'dssm.h'

c...TB: stripped off MSSM part 11/2016

      if(roption.eq.'norun') then  
        mtmt=mass(kt)
      elseif(roption.eq.'1loop') then
        call dsfindmtmt1loop
      elseif(roption.eq.'4loop') then
        call dsfindmtmt4loop
      else
        write(*,*) 'invalid option roption = ',roption
        write(*,*) 'in function dsralph3'
        write(*,*) 'program stopped'
        stop
      endif

      return
      end  
