************************************************************************
***  subroutine checking whether you need to reload the label associated 
***  to diffusion parameters involved in the ep axisymmetric green
***  function tabulation. If they are the integer dsepcraxino is 
***  increased by one unit
***  
************************************************************************
      subroutine dsepcraxick
      implicit none
      include 'dslabelcom.h'
      include 'dscraxicom.h' ! for dsepcraxino
      integer i
      real*8 diff,sum,dummy,dsepcraxisetlabel
      logical nomatch
ccc
      real*8 parstore(100)
      integer nparstore
      common/dsepcraxickcom/parstore,nparstore      
ccc      
      logical first
      data first/.true./
      save first
ccc
      dummy=dsepcraxisetlabel(labset) ! for npar
      dummy=dsepcraxisetlabel(labin) ! for par(1) -- par(npar)
ccc
      nomatch=.false.
      if(nparstore.ne.npar) nomatch=.true. 
      do i=1,npar
        diff=dabs(parstore(i)-par(i))
        sum=dabs(parstore(i)+par(i))
        if(diff.gt.1.d-7*sum) nomatch=.true.
      enddo
      if(first) then
        dsepcraxino=0
        first=.false.
        nomatch=.true.
      endif
      if(nomatch) then
        dsepcraxino=dsepcraxino+1
        nparstore=npar
        do i=1,npar
          parstore(i)=par(i)
        enddo
      endif
      return
      end


      integer function dsepcraxinumber()
************************************************************************
***  function returning the integer dsepcraxino associated to diffusion
***  parameters involved in the ep green functions
***  
************************************************************************
      implicit none
      include 'dscraxicom.h'
      dsepcraxinumber=dsepcraxino
      return
      end
      
