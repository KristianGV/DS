!Midpoint integration for power-law power spectra
!Adapted from Numerical Recipes.
!Author: Pat Scott (p.scott@imperial.ac.uk)
!Date: Some time in 2010/11.
!Added to DarkSUSY: Oct 13 2014

      subroutine midpnt(func,n_spec,alpha,k,a,b,s,n)
      implicit none
      integer n
      double precision n_spec,alpha,k,a,b,s,func
      external func
      integer it,j
      double precision ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5d0*(a+b),n_spec,alpha,k)
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.d0*tnm)
        ddel=del+del
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x,n_spec,alpha,k)
          x=x+ddel
          sum=sum+func(x,n_spec,alpha,k)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.d0
      endif
      return
      end

