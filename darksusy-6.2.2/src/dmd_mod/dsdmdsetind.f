      subroutine dsdmdsetind(indexout,rein1)
c_______________________________________________________________________
c subroutine assigning the index 'indexout' to store in the vector
c dmdparvec the real*8 input value rein1.
c_______________________________________________________________________
      implicit none
      include 'dsdmddrvrcom.h'
      integer indexout
      real*8 rein1
ccc
ccc overwrite for temporary halo model:
ccc      
      if(ihalocur.eq.ihalotmp) then
        ndmdpartmp=ndmdpartmp+1
        dmdparvec(ndmdpartmp)=rein1
        indexout=ndmdpartmp
        return
      endif  
ccc
ccc store in the repository for all the others:
ccc      
      ndmdpar=ndmdpar+1
      if(ndmdpar.gt.ndmdparmax) then
        write(*,*) 'DS: call to dsdmdsetind exceeding # of parameters'
        stop
      endif  
      dmdparvec(ndmdpar)=rein1
      indexout=ndmdpar 
      return
      end
