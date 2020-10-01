      real*8 function dstaufunsph(xx,ilosi)
************************************************************************
*** function needed for absorption along the l.o.s.i. as implemented in
*** the function dslosisph and dsomlosisph
***
*** this function is linked after a call to dstausetsph; if in that call
*** locilosi has been set equal to ilosiabscomp then this function calls
*** dslosidriver and get the value of dstaufunsph, otherwise use the       
*** tabulated version loaded in dstausetsph
***
*** xx is a line of sight variable, as measured from the closest point
*** to the center of the spherical system
*** ilosi = index for (eventually) linking the driver dslosidriver
***
*** dstaufunsph is dimensionless  
************************************************************************
      implicit none
      include 'dsdmdintcom.h' ! for taureschoweq1,locilosi
      include 'dslosidrvrcom.h'
      include 'dsdvcom.h'
      real*8 xx
      integer ilosi
      real*8 out,result
      integer iwin
c
      real*8 xxvec(1000),yyvec(1000),yyvec2(1000)
      integer nintxx,howtau
      common/dstauxxtabcom/xxvec,yyvec,yyvec2,nintxx,howtau
c
      if(locilosi.eq.ilosiabscomp) then
ccc
ccc dstaufunsph computed in dslosidriver
ccc        
        iwin=ilositaufun+ilosi
        dvsph(idvrsph)=xx ! this is conventional for 1 dependent variable
        call dslosidriver(iwin,ntysph,dvsph,out)
        dstaufunsph=out
        return
      endif
ccc
ccc otherwise use the tabulation:
ccc      
      if(howtau.eq.1) then
        if(xx.ge.xxvec(1).and.xx.le.xxvec(nintxx)) then
          call dssplint(xxvec,yyvec,yyvec2,nintxx,xx,result)
          dstaufunsph=result
        else
          dstaufunsph=0.d0
        endif
      elseif(howtau.eq.2) then
        if(xx.ge.dexp(xxvec(1)).and.xx.le.dexp(xxvec(nintxx))) then
          call dssplint(xxvec,yyvec,yyvec2,nintxx,dlog(xx),result)
          dstaufunsph=dexp(result)
        else
          dstaufunsph=0.d0
        endif
      else
        write(*,*) 'DS: in dstaufunsph, howtau not properly',
     &    ' initialized = ',howtau
        stop
      endif
c      write(*,*) 'xx, dstaufunsph ',xx,dstaufunsph
      return
      end
