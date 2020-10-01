************************************************************************
***  function setting the correspondence between a label and a set of
***  parameters involved in the pb and db diffusion time tabulations
***
***  NOTE: this version is connected to the specific choice of the
***  diffusion coefficient function dskdiff
***      
***  type : REPLACEABLE                                                   
***  
************************************************************************
      real*8 function dspbcraxisetlabel(inout)
      implicit none
      include 'dslabelcom.h'
      include 'dscraxicom.h'  ! for pbcraxiaddlab & pbcraxisuf
      integer inout
ccc
      if(inout.ne.labin.and.inout.ne.labout.and.inout.ne.labset) then
        write(*,*) 'DS: call to dspbcraxisetlabel with inout = ',inout
        write(*,*) 'DS: out of the allowed values :',labin,labout,labset
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
ccc load settings:      
ccc
      if(inout.eq.labset) then
        npar=14
        call dsdatafile(labelfile,'axipblabel.dat')
        addlabel=pbcraxiaddlab
        labelsuffix=pbcraxisuf
        goto 100
      endif
ccc
ccc parameter list:
ccc
      if(inout.eq.labin) then
        par(1)=diffhh
        par(2)=diffRh
        par(3)=diffhg
        par(4)=diffng
        par(5)=diffnh
        par(6)=diffcvel
        par(7)=diffrcpb
ccc
ccc up to here are needed for any choice of the diffusion coefficient
ccc function dskdiff; the following are specific for dskdiff and would
ccc a different setup if such function is replaced
ccc        
        par(8)=k0halo
        par(9)=k0gasdisc
        par(10)=kdiffrig0
        par(11)=kdiffdelta
        par(12)=kdiffdeltalow
        par(13)=kdiffeta
        par(14)=dble(nkdiff)
      elseif(inout.eq.labout) then  
        diffhh=par(1)
        diffRh=par(2)
        diffhg=par(3)
        diffng=par(4)
        diffnh=par(5)
        diffcvel=par(6)
        diffrcpb=par(7)
ccc
ccc up to here are needed for any choice of the diffusion coefficient
ccc function dskdiff; the following are specific for dskdiff and would
ccc a different setup if such function is replaced
ccc        
        k0halo=par(8)
        k0gasdisc=par(9)
        kdiffrig0=par(10)
        kdiffdelta=par(11)
        kdiffdeltalow=par(12)
        kdiffeta=par(13)
        nkdiff=int(par(14))
      endif
ccc some fake initialization for the compiler not to complain:
 100  dspbcraxisetlabel=0.d0
      return
      end

