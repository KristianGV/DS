************************************************************************
*** function dsIB2intres returns the energy integral over the supplied 
*** function, assuming that this function features the same relevant 
*** resonances as the amplitude. 
***
***   Input:  intfunction - function to be integrated over
***           ptype       - type of final state particle for which the 
***                         energy integral is performed
***                         ('B','f' or 'fbar')
***           xmin,xmax   - energy range considered [in units of m_chi]
***           acc         - required accuracy
***   Output: ier         - error returned by dqagse
***                         (corresponding contributions to value of
***                          integral are always set to 0)
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-03-06
************************************************************************

      real*8 function dsib2intres(intfunction,ptype,xmin,xmax,acc,ier)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'


c------------------------ variables ------------------------------------
      character*(*) ptype
      real*8 intfunction,xmin,xmax,acc
      external intfunction

      integer nresmax
      parameter (nresmax=20)
      real*8 x_res(nresmax),dx_res(nresmax), dx_resp(nresmax),xint(nresmax) 
      real*8 dx_restmp(nresmax),tmpres, result
      integer nres,i,j, jmin,kres(nresmax),nint

      integer ier
c-----------------------------------------------------------------------


      result=0d0
      tmpres=0d0
      ier=0

      nint=1
      do 10 i=1,nresmax 
        xint(i)=0d0
        kres(i)=0
 10   continue
      xint(1)=xmin
      xint(2)=xmax


      if (ptype.eq.'B') then      
c there must be a smarter way to do this (i.e. select only relevant resonances...)
         nres=4
         kres(1)=kh1
         kres(2)=kh2
         kres(3)=kh3
         kres(4)=kz
         if (c_Btype.eq.2.or.c_Btype.eq.6) then
           nres=nres+2
           kres(5)=khc
           kres(6)=kw           
         endif
         
c... something like this *should* work better, but doesn't...!?
c         nres=0
c         minres=0.01d0
c         sv=dssigmav(13+c_ftype)
c         if ((c_Btype.eq.1.and.dssigmav(8)/sv.gt.minres).or.  
c     &       (c_Btype.eq.3.and.dssigmav(2)/sv.gt.minres).or.
c     &       (c_Btype.eq.4.and.dssigmav(1)/sv.gt.minres).or.
c     &       (c_Btype.eq.5.and.dssigmav(5)/sv.gt.minres)) then
c           nres=nres+1
c           kres(nres)=kh1
c         endif
c         if ((c_Btype.eq.1.and.dssigmav(9)/sv.gt.minres).or.  
c     &       (c_Btype.eq.3.and.dssigmav(3)/sv.gt.minres).or.
c     &       (c_Btype.eq.4.and.dssigmav(2)/sv.gt.minres).or.
c     &       (c_Btype.eq.5.and.dssigmav(6)/sv.gt.minres)) then
c           nres=nres+1
c           kres(nres)=kh2
c         endif
c         if ((c_Btype.eq.1.and.dssigmav(10)/sv.gt.minres).or.  
c     &       (c_Btype.eq.3.and.dssigmav(6)/sv.gt.minres).or.
c     &       (c_Btype.eq.4.and.dssigmav(5)/sv.gt.minres).or.
c     &       (c_Btype.eq.5.and.dssigmav(4)/sv.gt.minres)) then
c           nres=nres+1
c           kres(nres)=kh3
c         endif
c         if ((c_Btype.eq.1.and.dssigmav(12)/sv.gt.minres).or.  
c     &       (c_Btype.eq.3.and.dssigmav(9)/sv.gt.minres).or.
c     &       (c_Btype.eq.4.and.dssigmav(8)/sv.gt.minres).or.
c     &       (c_Btype.eq.5.and.dssigmav(10)/sv.gt.minres)) then
c           nres=nres+1
c           kres(nres)=kz
c         endif
c         if ((c_Btype.eq.2.and.dssigmav(11)/sv.gt.minres).or.  
c     &       (c_Btype.eq.6.and.dssigmav(7)/sv.gt.minres)) then
c           nres=nres+1
c           kres(nres)=khc
c         endif
c         if ((c_Btype.eq.2.and.dssigmav(13)/sv.gt.minres).or.  
c     &       (c_Btype.eq.6.and.dssigmav(11)/sv.gt.minres)) then
c           nres=nres+1
c           kres(nres)=kw
c         endif


c... resonance at z = 1 from FSR decay of top quarks
         if (c_Btype.eq.2.and.(c_ftype.eq.11.or.c_fbartype.eq.11)) then
            nres = nres+1
            if (ptype.eq.'fbar') x_res(nres) = 1d0 - c_E1
            if (ptype.eq.'f') x_res(nres) = 1d0 - c_E2
c            if (Mff.gt.MB) x_res(nres) = 1d0 - c_E1
c            if (Mf.gt.MB) x_res(nres) = 1d0 - c_E2
            kres(nres) = kt
         endif

         do 50 i=1,nres
            if (kres(i).ne.kt)
     &        x_res(i) = 1.d0+mB**2/4.- mass(kres(i))**2/(4.*mx**2)
            dx_res(i) = intres*mass(kres(i))*width(kres(i))/mx**2 !  int range
            if (kres(i).eq.kh2) 
     &           dx_res(i)=10.*dx_res(i) ! need 'larger' range for very narrow SM H!
 50      continue
         
         do 70 i=1,nres-1 
            jmin=i
            do 60 j=i+1,nres
               if (x_res(j).lt.x_res(jmin)) jmin=j
 60         continue
            if (jmin.gt.i) then               ! NB: use (nres+1) as dummy variable
               x_res(nres+1)=x_res(i)
               x_res(i)=x_res(jmin)
               x_res(jmin)=x_res(nres+1)
               dx_res(nres+1)=dx_res(i)
               dx_res(i)=dx_res(jmin)  
               dx_res(jmin)=dx_res(nres+1)
               kres(nres+1)=kres(i)
               kres(i)=kres(jmin)
               kres(jmin)=kres(nres+1)
            endif
 70      continue
         
c.. check for overlapping resonances (note asymmetry: change x+dx, but not x-dx)
         do i = 1, nres-1
            dx_resp(i)=dx_res(i)
            if(x_res(i + 1) - x_res(i).lt.dx_res(i)) then
               dx_resp(i) = (x_res(i + 1) - x_res(i))/2.
               dx_restmp(i) = x_res(i+1) - dx_res(i+1) - x_res(i) 
               if(dx_restmp(i).gt.0d0.and.dx_restmp(i).lt.dx_res(i)) then
                  dx_resp(i)= dx_restmp(i)      ! (~if i+1 is more narrow than i)
               endif
            endif
         enddo
         dx_resp(nres)=dx_res(nres)
         
c... split integral around potentially important resonances
         do 100 i=1,nres
            if (x_res(i)-dx_res(i).gt.xint(nint).and.
     &           x_res(i)-dx_res(i).lt.xmax) then
               nint=nint+1
               xint(nint)=x_res(i)-dx_res(i)
            endif
            if (x_res(i)+dx_resp(i).gt.xint(nint).and.
     &           x_res(i)+dx_resp(i).lt.xmax) then
               nint=nint+1
               xint(nint)=x_res(i)+dx_resp(i)
            endif
 100      continue

      elseif ((ptype.eq.'fbar'.and.Mf.gt.MB).or.
     &        (ptype.eq.'f'.and.Mff.gt.MB)) then      
c... only possible resonance: at z = 1 from FSR of top quarks

         nres=1
         dx_res(nres) = 10.*intres*mass(kt)*width(kt)/mx**2  
         if (1.d0-dx_res(nres).gt.xint(nint).and.
     &       1.d0-dx_res(nres).lt.xmax) then
                nint=nint+1
                xint(nint)=1.d0-dx_res(nres)
         endif
         if (1.d0+dx_res(nres).gt.xint(nint).and.
     &       1.d0+dx_res(nres).lt.xmax) then
                nint=nint+1
                xint(nint)=1.d0+dx_res(nres)
         endif

      endif
      xint(nint+1)=xmax

      do  i=1,nint
        call dgadap(xint(i),xint(i+1),intfunction,acc,tmpres)
        result=result + tmpres
      enddo


      dsib2intres=result

      return
      end

