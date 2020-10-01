**********************************************************************
*** function dsanthreshold returns the correction factor to a 2-body
*** rate close to a kinematical threshold. We consider a process
*** DM DM -> A B, and assume that just above threshold the cross 
*** section scales as 
*** 
***   (sigma v)_AB \propto [\lambda(1, m_A^2/s, m_B^2/s)]^(n+1/2)
***
*** For s-wave annihilation, the cross section including the 
*** possibility to produce a virtual B (assumed to decay exclusively 
*** to much lighter particles) is then given by
***
***   (sigma v)_AB^* = (sigma v)_red * dsanthreshold
***
*** where the reduced cross section is given by
***
***   (sigma v)_red = S (sigma v)_AB / [\lambda(1, m_A^2/s, m_B^2/s)]^(n+1/2)
***
*** which by construction is non-zero at threshold. The symmetry factor 
*** is S=2 if A=B, otherwise S=1. [more complicated symmetry factors arise 
*** if A is the same as one of the decay products of B, see 1705.03466]
*** 
***    Input:
***
***          s      - CMS energy squared
***          mA, mB - mass of particle A, B
***          width  - total decay width of B
***          n      - integer (definition see above)
***
***  The output is dimension-less.
***
*** author: Torsten.Bringmann.fys.uio.no
*** date: 2017-05-11
**********************************************************************

      real*8 function dsanthreshold(s,mA,mB,width,n)
      implicit none
      include 'dsmpconst.h'

      real*8 s, mA,mB,width
      integer n

      real*8 mumax, tmpres, thrint
      external thrint 

c... common block needed for integration     
      real*8 muA, muB, gammaB, nn
      common /thrintcom/ muA, muB, gammaB, nn
      save /thrintcom/

     
      dsanthreshold=0.0d0
      
c... this only implements 3-body final states, so we need at least
c... enough energy to produce A
      if (mA.gt.sqrt(s)) return      

c... far below the threshold, we also do not expect a sizeable effect
c... so we simply set the cross section to zero. This effective threshold must
c... be consistently implemented also in dsrdparticles
      if (sqrt(s).lt.(mA+mB-16.0d0*width)) return      
      
c... far above threshold, we expect to return a factor of 1
c... (and the numerical integartion becomes unstable)
      if (sqrt(s).gt.(mA+100.d0*mB)) then
        dsanthreshold=1.0d0
        return
      endif

c... preparing integration of thrint (defined below)
      gammaB = width/mB
      muA = mA**2/s
      muB = mB**2/s
      nn = 0.5d0 + dble(n)

      mumax = ((sqrt(s) - mA)/mB)**2

      call dgadap(0.0d0,mumax,thrint,1.d-7,tmpres)


      dsanthreshold= tmpres/pi

      dsanthreshold=min(1.0d0,dsanthreshold) ! JE TMP
      
      return
      end

******************************************      
      real*8 function lambda(x,y,z)
      implicit none
      real*8 x,y,z
      lambda =  x**2 + y**2 + z**2 - 2.d0*(x*y + x*z + y*z)
      return
      end
******************************************      
      
******************************************      
      real*8 function thrint(mu)
      implicit none
      real*8 mu, lambda
      real*8 muA, muB, gammaB, nn
      common /thrintcom/ muA, muB, gammaB, nn
      save /thrintcom/
            
      thrint = gammaB*mu/((mu-1.d0)**2+gammaB**2) 
     &         * lambda(1.d0,muA,mu*muB)**nn

      return
      end
******************************************      
   
      
