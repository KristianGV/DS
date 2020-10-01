
***********************************************************************
*** Function dsib3svqqgratio returns the ratio between (FSR-subtracted)
*** qqg and tree-level cross section of a neutralino pair annihilating 
*** via an s-wave into a quark pair. See Bringmann, Galea & Walia (2015).
*** 
***  input:   kpn -- particle code for neutralino
***           kpq -- particle code for quark 
***
*** author: torsten bringmann
*** date: 2015-10-30
***********************************************************************
      real*8 function dsib3svqqgratio(kpn,kpq)
        implicit none
        include 'dsmssm.h'
        include 'dsmpconst.h'

        integer kpn,kpq, istat, iq

        real*8 mn, mq, alphas, VIBem, Emin

        integer dsidnumber
        integer idold
        data idold/-123456789/
        save idold

        real*8 rsave(42:45,7:12)
        data rsave/24*-1.2345d50/
        save rsave


c... functions
        real*8 dsralph3, dsIByieldone

c-----------------------------------------------------------------------
        dsib3svqqgratio=0.0d0
        
        if (kpn.lt.42.or.kpn.gt.45.or.kpq.lt.7.or.kpq.gt.12) return

        if (kpn.ne.42) return ! FIXME: IB routines currently cannot handle
                              ! other neutralinos

        if (idold.ne.dsidnumber()) then

          mn=mass(kpn)
          alphas=dsralph3(2.d0*mn)
          Emin=1.d-3*mn ! result does not change significantly for smaller Emin

          do iq=7,12
            mq=mass(iq)
            if (mn.gt.mq) then
              VIBem=dsIByieldone(Emin,iq,52,istat) 
              if (iq.eq.7.or.iq.eq.9.or.iq.eq.11) then        ! u, c, t quarks
                rsave(kpn,iq) = VIBem*(4./3.d0 * 9./4.d0)*alphas/alphem
              elseif (iq.eq.8.or.iq.eq.10.or.iq.eq.12) then  ! d, s, b quarks
                rsave(kpn,iq) = VIBem*(4./3.d0 * 9.)*alphas/alphem
              endif  
            endif
          enddo

          idold=dsidnumber()

        endif
 
        dsib3svqqgratio = rsave(kpn,kpq)

        return

      end
