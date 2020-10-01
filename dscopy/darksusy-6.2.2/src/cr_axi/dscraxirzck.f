************************************************************************
***  subroutine checking whether you need to reload the label associated 
***  to the position in the galaxy where craxi tabulations are computed
***  
***  output: noout -> an integer which is increased for any new set of
***          parameters      
***  
************************************************************************
      subroutine dscraxirzck(R,z,noout)
      implicit none
      include 'dscraxicom.h'
      real*8 R,z
      integer noout
      real*8 diff,sum
      logical nomatch
ccc
      real*8 Rcraxi,zcraxi
      common/dscraxirzckcom/Rcraxi,zcraxi      
ccc
      integer dscraxirzno
      common/dscraxirznocom/dscraxirzno
ccc
      logical first
      data first/.true./
      save first
ccc
      nomatch=.false.
      diff=dabs(Rcraxi-R)
      sum=dabs(Rcraxi+R)
      if(diff.gt.1.d-7*sum) nomatch=.true. 
      diff=dabs(zcraxi-z)
      sum=dabs(zcraxi+z)
      if(diff.gt.1.d-7*sum) nomatch=.true. 
      if(first) then
        dscraxirzno=0
        first=.false.
        nomatch=.true.
      endif
      if(nomatch) then
        dscraxirzno=dscraxirzno+1
        Rcraxi=R
        zcraxi=z
      endif
      noout=dscraxirzno
      return
      end
