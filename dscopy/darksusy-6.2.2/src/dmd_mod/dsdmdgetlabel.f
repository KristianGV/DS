      subroutine dsdmdgetlabel(labelinoutfun,labelout)
************************************************************************
*** subroutine to search in file 'labelfile' for the label corresponding
*** to the currently implemented set of npar parameters par
*** if the set of npar parameters par does not exist in labelfile, a new
*** label is added 
*** automatically created out of the suffix labelsuffix and a sequential
*** number, and label/parameters are stored in labelfile; if addlabel is
*** set to false, an error message is printed and the program stops
***
*** input:
***     labelinoutfun - external function setting internal common blocks      
*** inputs from the common block dslabelcom:
***     labelfile - in the file in which label and parameters are stored
***     npar - the number of parameters to be stored (npar.le.100)
***     par - a real*8 vector with the parameters to be stored
***     addlabel - logical variable, if set to true new label/parameters
***       can be stored in labelfile
***     labelsuffix - character*4 suffix to create automatically new 
***       labels to store new sets of parameters 
*** output:
***     labelout - the label found or autmatically generated
************************************************************************
      implicit none
      include 'dslabelcom.h'
ccc
      real*8 labelinoutfun
      external labelinoutfun
      character*10 labelout
ccc
      character*10 locallabel
      character*5 alabel
      real*8 localpar(100),diff,sum,dummy 
      integer i,ilabel,nlabsuff
      logical matchall,match
ccc
ccc load label settings
ccc        
      dummy=labelinoutfun(labset)
ccc
ccc load label parameters
ccc        
      dummy=labelinoutfun(labin)
ccc
      nlabsuff=0
      do i=1,4  ! labelsuffix is character*4
        if(labelsuffix(i:i).ne.' ') nlabsuff=nlabsuff+1
      enddo
      if(nlabsuff.eq.0) then
        write(*,*) 'DS: error in dsdmdgetlabel : '
        write(*,*) 'DS: empty labelsuffix = ',labelsuffix
        write(*,*) 'DS: program stopped'
        stop
      endif  
ccc
      open(unit=labelunit,file=labelfile,status='unknown',
     &     form='formatted')
      ilabel=0
 200  read(labelunit,1000,end=100,err=101)
     & locallabel,(localpar(i),i=1,npar)
      matchall=.true.
      call dslabcheck(10,locallabel,nlabsuff,labelsuffix,match)
      if(.not.match) matchall=.false.
      do i=1,npar
        diff=dabs(localpar(i)-par(i))
        sum=dabs(localpar(i)+par(i))
        if(diff.gt.1.d-7*sum.and.matchall) matchall=.false. 
      enddo
ccc
ccc if a match is found, assign label and return
ccc
      if(matchall) then
        labelout=locallabel
        close(labelunit)
        return
      endif
 101  continue
      ilabel=ilabel+1
      goto 200
 100  continue
      close(labelunit)
      ilabel=ilabel+1
ccc
ccc if no match is found, there are two possibilities. if addlabel is 
ccc set to true, a new label is automatically created out of the suffix 
ccc labelsuffix and ilabel:
ccc
      if(addlabel) then
        labelout='a'
        call dscharadd(labelout,labelsuffix)
        write(alabel,2000) ilabel
 2000   format(i5)
        call dscharadd(labelout,alabel)
        open(unit=labelunit,file=labelfile,status='unknown',
     &       form='formatted',access='append')
        write(labelunit,1000) labelout,(par(i),i=1,npar)
        close(labelunit)
ccc 
ccc if addlabel is set to false, write and error message and stop
ccc
      else
        write(*,*) 'DS: error in dsdmdgetlabel for file : ',labelfile
        write(*,*) 'DS: attempt to find a label for a parameter set'
        write(*,*) 'DS: which is not available'
        write(*,*) 'DS: par : ',(par(i),i=1,npar)
        write(*,*) 'DS: and possibility to write in file switched off'
        write(*,*) 'DS: program stopped'
        stop
      endif
 1000 format(1x,a10,100(1x,e14.8))
      return
      end
