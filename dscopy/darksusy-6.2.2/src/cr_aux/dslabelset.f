      subroutine dsaddlabel2(labelinoutfun,labelin)
************************************************************************
*** subroutine which adds to file 'labelfile' the new label input   
*** 'labelin' and a row of associated npar parameters par
*** if labelain already exists and it is associated to the same set of 
*** parameters, the subroutine just returns
*** if labelin already exists and it is associated to a different set
*** of parameters, an error message is printed and the program stops
***
*** inputs:
***     labelinoutfun - external function setting internal common blocks      
***     labelin - the label to be added (character*10)      
*** inputs from the common block dslabelcom:
***     labelfile - in the file in which label and parameters are stored
***     npar - the number of parameters to be stored (npar.le.100)
***     par - a real*8 vector with the parameters to be stored
************************************************************************
      implicit none
      include 'dslabelcom.h'
ccc
      real*8 labelinoutfun
      external labelinoutfun
      character*10 labelin
ccc      
      character*10 locallabel
      real*8 localpar(100),diff,sum,dummy
      integer i,j
      logical matchall,matchlab
ccc
ccc load label settings
ccc        
      dummy=labelinoutfun(labset)
ccc
ccc load label parameters
ccc        
      dummy=labelinoutfun(labin)
ccc
      open(unit=labelunit,file=labelfile,status='unknown',
     &     form='formatted')
 200  read(labelunit,1000,end=100,err=101)
     &  locallabel,(localpar(i),i=1,npar)
      matchall=.true.
      do i=1,npar
        diff=dabs(localpar(i)-par(i))
        sum=dabs(localpar(i)+par(i))
        if(diff.gt.1.d-7*sum.and.matchall) matchall=.false. 
      enddo
ccc
ccc check whether labelin and locallabel are matching or not
ccc
      matchlab=.true.
      do j=1,10
        if((labelin(j:j).eq.' ').and.(locallabel(j:j).ne.' ')) 
     &    matchlab=.false.
        if((labelin(j:j).ne.' ').and.(locallabel(j:j).ne.' ')) then
          if((labelin(j:j).ne.locallabel(j:j)).and.matchlab) 
     &      matchlab=.false.
        endif
      enddo
      if(matchall.and.matchlab) then
ccc
ccc trying to add a label/parameters which are already in the label file
ccc close file and return
ccc
        close(labelunit)
        return
      endif
ccc
      if(matchall.and.(.not.matchlab)) then
        write(*,*) 'DS: error in dsaddlabel2 for file : ',labelfile
        write(*,*) 'DS: attempt to add the new label = '
        write(*,1000) labelin 
        write(*,*) 'DS: for parameters already stored under label = '
        write(*,1000) locallabel 
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
      if((.not.matchall).and.matchlab) then
        write(*,*) 'DS: error in dsaddlabel2 for file : ',labelfile
        write(*,*) 'DS: attempt to add a new model under a label'
        write(*,*) 'DS: already tagging another model; new label = '
        write(*,1000) labelin 
        write(*,*) 'DS: label found in file = '
        write(*,1000) locallabel 
        write(*,*) 'DS: new model / stored model parameters : '
        do i=1,npar
          write(*,*) i,par(i),localpar(i)
        enddo
        write(*,*) 'DS: program stopped'
        stop
      endif
 101  continue
      goto 200
 100  close(labelunit)
ccc
ccc add to the file label and parameters
ccc
      open(unit=labelunit,file=labelfile,status='old',
     &     form='formatted',access='append')
      write(labelunit,1000) labelin,(par(i),i=1,npar)
      close(labelunit)
 1000 format(1x,a10,100(1x,e14.8))
      return
      end


      subroutine dsgetlabel2(labelinoutfun,labelout)
************************************************************************
*** subroutine to search in file 'labelfile' for the label corresponding
*** to the currently implemented set of npar parameters par
*** if the set of npar parameters par does not exist in labelfile and 
*** the logical variable addlabel is set to true, a new label is 
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
      integer i,ilabel
      logical matchall
ccc
ccc load label settings
ccc        
      dummy=labelinoutfun(labset)
ccc
ccc load label parameters
ccc        
      dummy=labelinoutfun(labin)
ccc
      open(unit=labelunit,file=labelfile,status='unknown',
     &     form='formatted')
      ilabel=0
 200  read(labelunit,1000,end=100,err=101)
     & locallabel,(localpar(i),i=1,npar)
      matchall=.true.
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
ccc if addlabel is set to true, write and error message and stop
ccc
      else
        write(*,*) 'DS: error in dsgetlabel2 for file : ',labelfile
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

      

      subroutine dsreadlabel2(labelinoutfun,labelin)
************************************************************************
*** subroutine to read in file 'labelfile' the set of npar parameters 
*** par corresponding to the label 'labelin'
***
*** inputs from the common block dslabelcom:
***     labelfile - in the file in which label and parameters are stored
***     labelin - character*10 the label to link to
************************************************************************
      implicit none
      include 'dslabelcom.h'
ccc
      real*8 labelinoutfun
      external labelinoutfun
      character*10 labelin
ccc       
      character*10 locallabel
      real*8 localpar(100),dummy
      integer i,j
      logical matchall
ccc
ccc load label settings
ccc        
      dummy=labelinoutfun(labset)
ccc
      open(unit=labelunit,file=labelfile,status='unknown',
     &     form='formatted')
 200  read(labelunit,1000,end=100,err=101)
     &  locallabel,(localpar(i),i=1,npar)
      matchall=.true.
      do j=1,10
        if((labelin(j:j).eq.' ').and.(locallabel(j:j).ne.' ')) 
     &    matchall=.false.
        if((labelin(j:j).ne.' ').and.(locallabel(j:j).ne.' ')) then
          if((labelin(j:j).ne.locallabel(j:j)).and.matchall) 
     &      matchall=.false.
        endif
      enddo
      if(matchall) then
        do i=1,npar
          par(i)=localpar(i)
        enddo
        close(labelunit)
ccc
ccc pass label parameters
ccc        
        dummy=labelinoutfun(labout)
        return
      endif
 101  continue
      goto 200
 100  continue
      write(*,*) 'DS: error in dsreadlabel2 for file : ',labelfile
      write(*,*) 'DS: there is no label matching the input'
      write(*,*) 'DS: label = ',labelin
      write(*,*) 'DS: program stopped'
      stop
 1000 format(1x,a10,100(1x,e14.8))
      end
