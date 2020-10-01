      real*8 function dsdmdsetlabel(inout)
c_______________________________________________________________________
c auxiliary function to read/write/store halo profile within the file
c dmdlabel.dat
c NOTE: since dmdlabel.dat may contain automatically generated labels
c different users may have the same halo profile stored under a
c different label, hence also linking to different tabulation files;
c to avoid this, make sure that the file dmdlabel.dat is the same and
c replace names for the corresponding linked files      
c_______________________________________________________________________
      implicit none
      include 'dslabelcom.h'
      include 'dsdmddrvrcom.h'
      integer inout
      integer ii

      dsdmdsetlabel=0.0d0 ! TB; make sure function has return value

ccc
      if(inout.ne.labin.and.inout.ne.labout.and.inout.ne.labset) then
        write(*,*) 'DS: call to dsdmdsetlabel with inout = ',inout
        write(*,*) 'DS: out of the allowed values :',labin,labout,labset
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc temporary patch: labout cannot be called now
ccc      
      if(inout.eq.labout) then
        write(*,*) 'DS: call to dsdmdsetlabel with inout = labout'
        write(*,*) 'DS: this option is temporarely inactive'
        write(*,*) 'DS: program stopped'
        stop
      endif        
ccc
ccc load settings:      
ccc
      if(inout.eq.labset) then
        npar=ndmdlab
        call dsdatafile(labelfile,'dmdlabel.dat')
        addlabel=.true.
        labelsuffix=sufflab
        return
      endif
ccc
ccc parameter list: 
ccc      
      if(inout.eq.labin) then
        do ii=1,ndmdlab
           par(ii)=dmdparvec(ihalotag(idmdlab(ii),ihalocur))
        enddo
c      elseif(inout.eq.labout) then  
c        do ii=1,ndmdlab
c          dmdparvec(ihalotag(idmdlab(ii),ihalocur))=par(ii)
c        enddo
      endif
      return
      end
