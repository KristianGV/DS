************************************************************************
***  subroutine checking whether you need to reload the label associated 
***  to diffusion parameters involved in the pb and db diffusion time
***  tabulations. If they are the integer dspbcraxino is increased by
***  one unit
***  
************************************************************************
      subroutine dspbcraxick
      implicit none
      include 'dslabelcom.h'
      include 'dscraxicom.h' ! for dspbcraxino
      integer i
      real*8 diff,sum,dummy,dspbcraxisetlabel
      logical nomatch
ccc
      real*8 parstore(100)
      integer nparstore
      common/dspbcraxickcom/parstore,nparstore      
ccc      
      logical first
      data first/.true./
      save first
ccc
      dummy=dspbcraxisetlabel(labset) ! for npar
      dummy=dspbcraxisetlabel(labin) ! for par(1) -- par(npar)
ccc
      nomatch=.false.
      if(nparstore.ne.npar) nomatch=.true. 
      do i=1,npar
        diff=dabs(parstore(i)-par(i))
        sum=dabs(parstore(i)+par(i))
        if(diff.gt.1.d-7*sum) nomatch=.true.
      enddo
      if(first) then
        dspbcraxino=0
        first=.false.
        nomatch=.true.
      endif
      if(nomatch) then
        dspbcraxino=dspbcraxino+1
        nparstore=npar
        do i=1,npar
          parstore(i)=par(i)
        enddo
      endif
      return
      end


      integer function dspbcraxinumber()
************************************************************************
***  function returning the integer dspbcraxino associated to diffusion
***  parameters involved in the pb and db diffusion time tabulations.
***  
************************************************************************
      implicit none
      include 'dscraxicom.h'
      dspbcraxinumber=dspbcraxino
      return
      end
