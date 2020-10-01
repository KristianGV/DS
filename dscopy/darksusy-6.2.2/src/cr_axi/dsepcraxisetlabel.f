************************************************************************
***  function setting the correspondence between a label and a set of
***  parameters involved in the ep axisymmetric green function
***  tabulation
***  
************************************************************************
      real*8 function dsepcraxisetlabel(inout)
      implicit none
      include 'dslabelcom.h'
      include 'dscraxicom.h' ! for pbcraxiaddlab & pbcraxisuf + diff. par.
      integer inout
ccc
      if(inout.ne.labin.and.inout.ne.labout.and.inout.ne.labset) then
        write(*,*) 'DS: call to dsepcraxisetlabel with inout = ',inout
        write(*,*) 'DS: out of the allowed values :',labin,labout,labset
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc load settings:      
ccc
      if(inout.eq.labset) then
        npar=2
        call dsdatafile(labelfile,'axieplabel.dat')
        addlabel=epcraxiaddlab
        labelsuffix=epcraxisuf
        goto 100
      endif
ccc
ccc parameter list:
ccc
      if(inout.eq.labin) then
        par(1)=diffhh
        par(2)=diffrcep
      elseif(inout.eq.labout) then  
        diffhh=par(1)
        diffrcep=par(2)
      endif
ccc some fake initialization for the compiler not to complain:
 100  dsepcraxisetlabel=0.d0
      return
      end

