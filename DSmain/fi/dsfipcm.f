
c_______________________________________________________________________
c  Function dsfipcm.
c       Finds momentum to particle 1 and 2 in particle X's 
c       reference frame
c
c  input:
c    mcm - mass for stationary particle X
c    m1  - mass for particle 1
c    m2  - mass for particle 2
c
c  author: Kristian Gjestad Vangsnes (kristgva@uio.no)   2020-08-14
c=======================================================================
      real*8 function dsfipcm(mcm,m1,m2)
      implicit none
      real*8 mcm,m1,m2,a,b,c,lambda
      a=mcm*mcm
      b=m1*m1
      c=m2*m2
c     Källen function:      
      lambda=a*a-2*a*(b+c)+(b-c)**2
      dsfipcm=sqrt(lambda)/(2*mcm)
      return
      end 