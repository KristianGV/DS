**********************************************************************
*** function dsanyieldline returns the branching fraction for some
*** line features.
*** The calculation is done at v=0 and does NOT include gamma rays,
*** these are instead handled by dsanyieldgammaline.
*** Input:
***   pdg = PDG code for final state particle (annihilation to pdg and -pdg
***         is assumed).
***         11 e+ e-
***         12 nu_e nu_e-bar
***         14 nu_mu nu_nu-bar
***         16 nu_tau nu_tau-bar
*** Output: branching fraction to channel
*** Units: unitless
**********************************************************************

      real*8 function dsanyieldline(pdg)
      implicit none
      integer pdg
      real*8 dssigmav0tot,dssigmav0
      real*8 sigv0

      sigv0=dssigmav0tot()
      if (pdg.eq.11) then ! e+ e-
         dsanyieldline=dssigmav0(11,-11)/sigv0
      elseif (pdg.eq.12) then ! nu_e nu_e-bar
         dsanyieldline=dssigmav0(12,-12)/sigv0
      elseif (pdg.eq.14) then ! nu_mu nu_mu-bar
         dsanyieldline=dssigmav0(14,-14)/sigv0
      elseif (pdg.eq.16) then ! nu_tau nu_tau-bar
         dsanyieldline=dssigmav0(16,-16)/sigv0
      else
         write(*,*) 'DS Warning in dsline: incorrect pdg code ',
     &   'pdg = ',pdg
         dsanyieldline=0.d0
      endif

      return
      end
