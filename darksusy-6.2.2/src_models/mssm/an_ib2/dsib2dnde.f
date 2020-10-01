************************************************************************
*** function dsIB2dnde returns the normalized energy spectrum of a given 
*** final state particle:
c BeginTex
c   \begin{displaymath}
c    \frac{dN}{dE} \equiv 
c    \frac{1}{(\sigma v)_f} \frac{d(\sigma v)_f}{dE_P}
c   \end{displaymath}
c EndTex
***
***   Input: IB2chinput  - selects channel f for 3-body final state:
***                        1XX for ffZ
***                        2XX for fFW
***                        3XX for ffh (XX like for ffZ)
***                        4XX for ffH (XX like for ffZ)
***                        5XX for ffA (XX like for ffZ)
***                        6XX for ffH+- (XX like for fFW)
***
***                    XX   | \bar f f Z      |   \bar f f W
***                   ------+-----------------+-------------------
***                    01   | nu_e   nu_e   Z |   e+     nu_e   W-
***                    02   | e      e      Z |   nu_e   e-     W+
***                    03   | nu_mu  nu_mu  Z |   mu+    nu_mu  W-
***                    04   | mu     mu     Z |   nu_mu  mu-    W+
***                    05   | nu_tau nu_tau Z |   tau+   nu_tau W-
***                    06   | tau    tau    Z |   nu_tau tau-   W+
***                    07   | u      u      Z |   dbar   u      W-
***                    08   | d      d      Z |   ubar   d      W+
***                    09   | c      c      Z |   sbar   c      W-
***                    10   | s      s      Z |   cbar   s      W+
***                    11   | t      t      Z |   bbar   t      W-
***                    12   | b      b      Z |   tbar   b      W+
***                    (note that this numbering must be consistent 
***                     with dmssm.h !!! )
***         
***          pfinal - type of considered final state particle:
***                   'f' for fermion
***                   'fbar' for antifermion
***                   'B' for vectorboson or scalar
***          Efinal - energy of considered final state particle [in GeV]
***
***   Output: differential energy spectrum [1/GeV]
***
*** Author: Torsten Bringmann (torsten.bringmann@fys.uio.no) 
*** Date:   2014-12-07
************************************************************************

      real*8 function dsib2dnde(IB2chinput,pfinal,Efinal)
      implicit none
      include 'dsib2com.h'
      include 'dsmpconst.h'

c------------------------ variables ------------------------------------
      integer IB2chinput
      real*8  Efinal, svib2, dnde, xfinal
      character*5 pfinal
c------------------------ functions ------------------------------------
      real*8   dsIB2dsde_aux,dsib2sigmav      
c-----------------------------------------------------------------------


      dsib2dnde=0d0
      svib2=dsib2sigmav(IB2chinput,1) ! this also sets all necessary 
                                      ! common blocks

      if (pfinal.ne.'f'.and.pfinal.ne.'fbar'.and.pfinal.ne.'B') then
        write(*,*) 'ERROR in dsib2dnde: unknown particle identifier: ',pfinal
        return
      endif

      c_pfinal=pfinal
      xfinal=Efinal/mx
      
      dnde=dsIB2dsde_aux(xfinal)*gev2cm3s/mx/svib2

      dsib2dnde=dnde/1.d20  !correct for factor introduced in dsib2Msqaux

      return
      end

