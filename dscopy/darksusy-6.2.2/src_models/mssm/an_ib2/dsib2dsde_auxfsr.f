************************************************************************
*** Function dsIB2dndz_fsr returns the energy distribution of one of 
*** the 3-body final state particles, implementing the model-
*** independent result for final-state radiation from 1104.2996 (1009.0224)
*** WARNING: If you want to call this function, you need to call 
***          dsib2chinit() first to specify the channel.
***
***   Input: pconv - final state particle type
***                  ('B', 'f' and 'fbar' currently implemented)
***          zint  - energy of particle pconv, divided by mx
***          Vtype - 1 for Z, 2 for W
***          ch    - fermion pair as defined in dsIB2yieldone
***
***   Output: (d(sv)_{3-body}/dz)
***
*** Author: Francesca Calore, 2012-09-12
************************************************************************
      real*8 function dsIB2dsde_auxfsr(zint,pconv, Vtype, ch) 
      implicit none
      include 'dsib2com.h'
      include 'dsmssm.h'
      include 'dsmpconst.h'
      

***************
*
* FIXME: NOT TESTED SUFFICIENTLY
*
***************





c------------------------ functions ------------------------------------
      real*8 PartonDistr
c------------------------ variables ------------------------------------
      real*8 zint
      character*(*) pconv 
      integer Vtype,ch
      character*10 pfin 
      real*8 tmpres,tmpres_trans,tmpres_long
      integer ch_up, ch_down
      real*8 s0
c-----------------------------------------------------------------------

      dsib2dsde_auxfsr=0.d0
      tmpres = 0.d0
      tmpres_trans =0.d0
      tmpres_long = 0.d0

      if (pconv.eq.'B') then
         pfin = 'B'
      
         if (c_Btype.eq.1) then ! factor of two since Z is a real particle and average over f fbar polarizations          
            tmpres_trans =2*(PartonDistr(ch, 0, pfin, 0, zint)
     -           +  PartonDistr(ch, 1, pfin, 0, zint))/2.
            
            tmpres_long = 2*PartonDistr(ch, 0, pfin, 1, zint) ! longitudinal only from unpolarized top
           
         elseif(c_Btype.eq.2) then ! average over f fbar polarizations
            if(mod(ch,2).EQ.1)then  ! up-type fermion
               ch_up = ch
               ch_down = ch+1  
            elseif(mod(ch,2).EQ.0)then  ! down-type fermion. With current dsib2chinit this is never realized
               ch_up = ch+1
               ch_down = ch
            endif
            tmpres_trans = (PartonDistr(ch_up, 0, pfin, 0, zint)
     -           +  PartonDistr(ch_up, 1, pfin, 0, zint) 
     -           +  PartonDistr(ch_down, 0, pfin, 0, zint)
     -           +  PartonDistr(ch_down, 1, pfin, 0, zint))/2.
            
            tmpres_long = (PartonDistr(ch_up, 0, pfin, 1, zint) ! longitudinal only from unpolarized t or b
     -           +  PartonDistr(ch_down, 0, pfin, 1, zint)) 
         endif

         tmpres=(2*tmpres_trans + tmpres_long)/3. ! unpolarized spectrum from unpolarized f fbar initial 2-body final state
c$$$      elseif (pconv.eq.'f'.or.pconv.eq.'fbar') then ! factor of two for average also for neutrinos.
c$$$         if (c_Btype.eq.1)then
c$$$            pfin = 'f'
c$$$            tmpres=(PartonDistr(ch, 0, pfin, 0, zint)
c$$$     -           +  PartonDistr(ch, 1, pfin, 0, zint))/2.
c$$$
c$$$         elseif(c_Btype.eq.2)then
c$$$            pfin = 'f_partner' ! from splitting function D(f -> f')
c$$$            tmpres= PartonDistr(ch, 0, pfin, 0, zint)/2. 
c$$$          
c$$$         endif
      endif

      s0 = 1.73272168d-31       !tmpres normalized accordingly to 1104.2996
      dsib2dsde_auxfsr=tmpres*s0*1.d-3 ! normalized accordingly to dsib2dndz 

      return 
      end

******************************************************************************
***   Auxiliary interface functions for dsib2dndz_fsr (all with _fsr)
******************************************************************************
***   Universal kinematical function L(z), z = E_{P}/mx = energy fraction
******************************************************************************
      real*8 function LFunction(E_frac)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'
 
      real*8 E_frac
      real*8 E2_cm
      real*8 Lfunc_tmp, arg
      
      LFunction =  0.d0

      Lfunc_tmp = 0.d0
      arg= 0.d0

      E2_cm = 4*mx**2
      arg = (E2_cm*E_frac**2)/(4.*MB**2)
   
      if (E_frac.le.(MB/mx))then
         Lfunc_tmp = 0.d0       
      elseif (E_frac.gt.(MB/mx)) then
         Lfunc_tmp = Log(arg) + 2*Log(1. + Sqrt(1. - 1/arg)) 
      else
         write(*, *)'ERROR: universal kinematical function
     &  has wrong argument ! ' ! due to kinematical boundaries
      endif

      LFunction =  Lfunc_tmp

      return 
      end
******************************************************************************
*** Generalized splitting functions for massive partons
******************************************************************************
      real*8 function Splitting(type_pOUT, x, Vir, Yuk)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'
      integer type_pOUT ! type of particle: 0 = F; 1 = V
      real*8 x
      integer Vir, Yuk

      real*8 Splitting_tmp
      real*8 Lfunc_std ! for massive parton w/o corrections
      real*8 LFunction
      real*8 s

      Splitting = 0.d0
      Splitting_tmp = 0.d0

      s = 4*mx**2
      Lfunc_std = Log(s/MB**2) 
c...  P_out = F = f or fbar
      if(type_pOUT.eq.0)then
c...  Vir = 0 : real splitting function
         if (Vir.eq.0) then 
            if (Yuk.eq.0)then
               Splitting_tmp = ((1. + x**2)*Lfunction(1. - x))/(1. - x)
            elseif(Yuk.eq.1)then
               Splitting_tmp = (1. - x)*lfunc_std
            endif
c...  Vir = 1 : virtual splitting function
         elseif (Vir.eq.1) then  
            if (Yuk.eq.0)then
               Splitting_tmp = (3*Lfunc_std)/2. - Lfunc_std**2/2.
            elseif(Yuk.eq.1)then
               Splitting_tmp = - Lfunc_std/2.
            endif
         endif

c...  P_out = V
      elseif(type_pOUT.eq.1) then
         if (Vir.eq.0.and.Yuk.eq.0) then 
            Splitting_tmp = ((1.+(1.-x)**2)*Lfunction(x))/x
         elseif(Vir.eq.1.and.Yuk.eq.0)then
            Splitting_tmp = (3*Lfunc_std)/2. - Lfunc_std**2/2.
         endif

c...  P_OUT = S (only top, bottom quarks)
      elseif(type_pOUT.eq.2) then   
         if (Vir.eq.0.and.Yuk.eq.1) then 
            Splitting_tmp = x*Lfunc_std
         endif
      else
         write(*,*) 'WARNING: missing case in splitting functions! '
      endif

      Splitting = Splitting_tmp
 
      return
      end
******************************************************************************
***   Parton distribution functions D(I -> J)(x). 
***   Tree-level + O(alpha_ew) corrections
******************************************************************************
      real*8 function PartonDistr(F, chir, pfin, polar, x)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'
      include 'dsmpconst.h'
      integer F, chir    ! F = initial 2-body fermion with chirality chir
      character*10 pfin
      integer polar ! polarization of vector bosons
      real*8 x

      real*8 alph2ew,y_top, alphatop, const_ew, wiso
      real*8 prefactor,PDistr_tmp

      real*8 Zcoupling, Splitting

      PartonDistr = 0.d0
      PDistr_tmp = 0.
      prefactor = 0.

c...  define D(I->J) constants
      alph2ew = alphem/s2thw
      const_ew = alph2ew/(2.*pi)
      y_top = mass(kt)/174.     ! as defined in the reference paper
      alphatop = y_top**2/(4.*pi)
      
      if (chir.eq.0) wiso =  wiso3(F) !left handed 
      if (chir.eq.1) wiso =  0  !right handed, wiso3(F) = 0

      if(pfin.eq.'B') then
         
         if(polar.eq.0) then ! transversal polarization
            
            if (c_Btype.eq.1) then
               prefactor = const_ew*(Zcoupling(F, wiso)**2/costhw**2) ! prefactor of the Splitting function
               PDistr_tmp = prefactor*Splitting(1, x, 0, 0)

            elseif(c_Btype.eq.2)then
               prefactor = const_ew*2*wiso**2
               PDistr_tmp = prefactor*Splitting(1, x, 0, 0)

            endif
            
         elseif(polar.eq.1) then ! longitudinal polarization only due to Yukawa, unpolarized
            prefactor = alphatop/(8.*pi)
                  
            if (c_Btype.eq.1) then
               if(F.eq.kt)
     -              PDistr_tmp = prefactor*Splitting(2, x, 0, 1)

            elseif(c_Btype.eq.2)then
               if (F.eq.kt.or.F.eq.kb) 
     -              PDistr_tmp = prefactor*Splitting(2, x, 0, 1)               
            endif
         else
            write(*,*)'ERROR in PartonDistr: 
     -                    Undefined vector polarization!'
         endif

      elseif(pfin.eq.'f')then
         polar = 0
         prefactor = const_ew*(Zcoupling(F, wiso)**2/costhw**2 + 
     -           echarg(F)**2*s2thw) ! prefactor of the Splitting function
         PDistr_tmp = prefactor*Splitting(0, x, 0, 0)
            
         if(abs(x - 1.d0).lt.1.d-10)then
            prefactor = const_ew*(2*wiso**2 + Zcoupling(F, wiso)**2/
     -           costhw**2 + echarg(F)**2*s2thw) ! prefactor of Pvir(F->F)
            PDistr_tmp = PDistr_tmp + (1 + prefactor*
     -           Splitting(0, x, 1, 0))
         endif

c...  for t, b quarks extra splittings due to the top Yukawa interaction are present
c...  for unpolarized quark

         if(F.eq.kt)then
            prefactor = alphatop/(4.*pi) 
            PDistr_tmp = PDistr_tmp + 
     -           prefactor*Splitting(0, x, 0, 1)
            if(abs(x - 1d0).lt.1d-10)then
               prefactor = (3*alphatop)/(8.*pi) 
               PDistr_tmp = PDistr_tmp + 
     -              prefactor*Splitting(0, x, 1, 1)
            endif
         elseif(F.eq.kb)then
            if(abs(x - 1d0).lt.1d-10)then
               prefactor = alphatop/(8.*pi) 
               PDistr_tmp = PDistr_tmp + 
     -              prefactor*Splitting(0, x, 1, 1)
            endif
         endif

      elseif(pfin.eq.'f_partner')then
         polar = 0
         prefactor = const_ew*2*wiso**2
         PDistr_tmp = prefactor*Splitting(0, x, 0, 0)
c...  for t, b quarks extra splittings due to the top Yukawa interaction are present
c...  for unpolarized quark
         if (F.eq.kt.or.F.eq.kb) then ! up-type fermion splitting in bottom t -> b or down-type b -> t
            prefactor = alphatop/(8.*pi)    
            PDistr_tmp = prefactor*Splitting(0, x, 0, 1)
         endif

      endif

      PartonDistr = PDistr_tmp

      return
      end
******************************************************************************
*** It returns g_{f} \equiv T3 - s2thw*Q 
******************************************************************************
      real*8 function Zcoupling(f, wiso)
      implicit none
      include 'dsmssm.h'
      include 'dsib2com.h'
      include 'dsmpconst.h'
      integer f
      real*8 wiso
      
      Zcoupling =  wiso - s2thw*echarg(f)
   
      return 
      end      
