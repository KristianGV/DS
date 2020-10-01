      subroutine dslinint(xvec,yvec,n,x,y)
c simple linear interpolation
      implicit none
      integer i,n
      real*8 xvec(n),x,yvec(n),y
      real*8 x1,x2,y1,y2,coeffa,coeffb
      call dslocate(xvec,n,x,i)
      if(i.eq.n) then
        y=yvec(n)
        return
      endif
      x1=xvec(i)
      x2=xvec(i+1)
      y1=yvec(i)
      y2=yvec(i+1)
      coeffa=(y1-y2)/(x1-x2)
      coeffb=y1-coeffa*x1
      y=coeffa*x+coeffb
      return
      end     



      subroutine dslocate(xx,n,x,j)
c fast interval finding
      integer j,n
      real*8 x,xx(n)
      integer jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end
