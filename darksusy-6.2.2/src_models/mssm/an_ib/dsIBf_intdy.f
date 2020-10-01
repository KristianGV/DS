*****************************************************************************
*** integrates a function f(x,y) over the kinematic variable 
***
*** y = (p+k)^2/(4 mx^2)
***
*** that appears for 3-body final states. (k is the photon momentum and p
*** that of one of the other final sates). The result is the
*** (dimensionless) differential photon yield per annihilation
***
*** called by dsIByield
***
*** author: Torsten Bringmann (bringman@sissa.it)
*** date:   2007-04-18
*** update: 2008-03-10 - better numerical convergence
*****************************************************************************

      real*8 function dsIBf_intdy(IBch,x,mwimp,mp1,mp2)

      implicit none
      include 'dsibcom.h'

      real*8 x,mwimp,mp1,mp2,eta,diffeta
      integer IBch

c... these parameters are needed by dqagse
      real*8 aint,bint,result,abserr
      integer neval,ier, limit
      parameter (limit=20)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last

      real*8 dsIBintsel
      external dsIBintsel

      dsIBf_intdy=0.d0

c... set parameters needed by dsIBintsel
      if ((IBch.ge.1.and.IBch.le.12).or.
     &    (IBch.ge.101.and.IBch.le.112)) then
          intch = IBch
      else
          return
      endif
      ibcom_x=x
      ibcom_mx=mwimp
      ibcom_mp1=mp1
      ibcom_mp2=mp2


c...  if kinematical limits satisfied, determine integration limits 
      if (((1.-(mp1+mp2)**2/4./mwimp**2)).gt.x) then   
         eta=(mp1/mwimp)**2
         diffeta=(mp2/mwimp)**2
         diffeta=(eta-diffeta)/4.
         aint=eta/4. + diffeta*x/2./(1.-x)
     &    + (1. - sqrt( (1.+diffeta/(1.-x))**2 - eta/(1.-x) ) )*x/2d0
         bint=eta/4. + diffeta*x/2./(1.-x)
     &    + (1. + sqrt( (1.+diffeta/(1.-x))**2 - eta/(1.-x) ) )*x/2d0
c         aint=aint*0.9999d0       ! this can be used to ensure numerical stability
c         bint=bint*0.9999d0       ! should not be necessary anymore
      else
         return
      endif

      call dqagse(dsIBintsel,aint,bint,1d-10,IBacc,20,result,abserr,
     &            neval,ier,alist,blist,rlist,elist,iord,last)

      if (ier.ne.0) iberr=ibset(iberr,1)

      dsIBf_intdy=result
     

      end
