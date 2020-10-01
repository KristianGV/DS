************************************************************************
***  function setting the correspondence between a label and a set of
***  paramters defining to the position in the galaxy where craxi
***  tabulations are computed
***  
************************************************************************
      real*8 function dscraxirzsetlabel(inout)
      implicit none
      include 'dslabelcom.h'
      include 'dscraxicom.h'
      integer inout
ccc
      real*8 Rcraxi,zcraxi
      common/dscraxirzckcom/Rcraxi,zcraxi      
ccc
      if(inout.ne.labin.and.inout.ne.labout.and.inout.ne.labset) then
        write(*,*) 'DS: call to dscraxirzsetlabel with inout = ',inout
        write(*,*) 'DS: out of the allowed values :',labin,labout,labset
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc load settings:      
ccc
      if(inout.eq.labset) then
        npar=2
        call dsdatafile(labelfile,'axirzlabel.dat')
        addlabel=rzaxiaddlab
        labelsuffix=rzaxisuf
        goto 100
      endif
ccc
ccc parameter list:
ccc      
      if(inout.eq.labin) then
        par(1)=Rcraxi
        par(2)=zcraxi
      elseif(inout.eq.labout) then  
        Rcraxi=par(1)
        zcraxi=par(2)
      endif
ccc some fake initialization for the compiler not to complain:
 100  dscraxirzsetlabel=0.d0
      return
      end
