************************************************************************
*** Function dsIB2dsde_aux returns the full SUSY differential cross
*** section, as a function of one of the 3-body final state particles P.
*** WARNING: If you want to call this function, you need to call 
***          dsib2chinit() first to specify the annihilation channel
***          and set common block c_pfinal to choose the particle P
***          ('B', 'f' and 'fbar' currently implemented)
***
***   Input: zint  - energy of particle c_pfinal, divided by mx
***
***   Output: (d(sv)_{3-body}/dz) as function of z = E_P/mx.
***           [units: 1/GeV**2]
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-03-05
************************************************************************
      real*8 function dsIB2dsde_aux(zint) 
      implicit none
      include 'dsib2com.h'
      include 'dsmpconst.h'
      include 'dsmssm.h'

c------------------------ functions ------------------------------------
      real*8   dsib2Msqaux,dsib2intres
      external dsib2Msqaux
c------------------------ variables ------------------------------------
      real*8   zint, result, zmax, Emin, Emax, Etmp1, Etmp2
      real*8   mfinal1, mfinal2, mfinal3
      integer  ier
c-----------------------------------------------------------------------


      dsIB2dsde_aux=0.d0
      result=0.d0

      if (c_pfinal.eq.'B') then
         c_EB=zint
         mfinal1=MB
         mfinal2=Mf
         mfinal3=Mff
      elseif (c_pfinal.eq.'f') then
         c_E1=zint
         mfinal1=Mf
         mfinal2=MB
         mfinal3=Mff
      elseif (c_pfinal.eq.'fbar') then
         c_E2 = zint
         mfinal1=Mff
         mfinal2=MB
         mfinal3=Mf
      else
         write(*,*) 'ERROR in dsIB2dsde_aux: c_pfinal = ',c_pfinal
         return   
      endif

      zmax= 1. + (mfinal1**2 - (mfinal2 + mfinal3)**2)/4.
      if (zint.lt.mfinal1.or.zint.gt.zmax) return 
      if (2.d0.le.(mfinal1+mfinal2+mfinal3)) return
 
      Etmp1=(2-zint)*(-4+4*zint-mfinal2**2+mfinal3**2-mfinal1**2)/
     -      (-8 + 8*zint - 2.*mfinal1**2)
      Etmp2=Sqrt((zint - mfinal1)*(zint + mfinal1)*
     -      (-4 + 4*zint + (mfinal2 - mfinal3)**2 - mfinal1**2)*
     -      (-4 + 4*zint + (mfinal2 + mfinal3)**2 - mfinal1**2))/
     -      (-8 + 8*zint - 2.*mfinal1**2)
      Emin=Etmp1+Etmp2
      Emax=Etmp1-Etmp2
c... WARNING: result is very sensitive to integration limits. 
c      Emin=1.00001*Emin
c      Emax=0.99999*Emax     
      if (Emax.lt.Emin) return

      if (c_pfinal.eq.'B') 
     &      result=dsib2intres(dsib2Msqaux,'f',Emin,Emax,IB2acc/3.,ier)
          
      if (c_pfinal.eq.'f'.or.c_pfinal.eq.'fbar') 
     &      result=dsib2intres(dsib2Msqaux,'B',Emin,Emax,IB2acc/3.,ier)

      if (ier.ne.0) write(*,*) 'dsIB2dsde_aux: WARNING - integration error:',
     &                          ier, c_pfinal, c_Btype, c_ftype

      
c...  overal normalization Msq -> sv
      dsIB2dsde_aux=result/((2.*pi)**3*16.*mx**2)
      
      return 
      end


