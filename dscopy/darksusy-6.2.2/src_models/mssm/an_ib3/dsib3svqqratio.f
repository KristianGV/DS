
***********************************************************************
*** Function dsib3svqqratio returns the ratio between QCD-corrected and
*** tree-level cross section of a neutralino pair annihilating via an
*** s-wave into a quark pair. See Eqs. (A17, A21) of 
*** Bringmann, Galea & Walia (2015).
*** 
***  input:   Ecm   -- total CMS energy 
***           kpn -- particle code for neutralino 
***           kpq -- particle code for quark 
***
*** author: torsten bringmann & parampreet walia 
*** date: 2015-10-30
***********************************************************************
      real*8 function dsib3svqqratio(Ecm,kpn,kpq)
        implicit none
        include 'dsmssm.h'
        include 'dsmpconst.h'

        real*8 Ecm
        integer kpq, kpn

        real*8 fact1,fact2,b0,brt, mq, alphas
        real*8 tf1,tf2,tf3,tf4,tf5,tf6,tf7,tf8,tf9,abeta

c... functions
        real*8 dsralph3, dsrmq, dsdilog

c-----------------------------------------------------------------------
        dsib3svqqratio=1.0d0
        if (kpn.lt.42.or.kpn.gt.45.or.kpq.lt.7.or.kpq.gt.12) return

        mq=mass(kpq)
        alphas=dsralph3(Ecm)
        if (Ecm.le.2.*mq) return


c...use re-summed cross section, (A21)
        if (0.925d0*mq.lt.Ecm) then  
        
            fact1 = 3.d0  !uncomment the below for first-order corrections
c           fact1=3.d0+(mq/Ecm)**2*(14.-8.*log(mq/Ecm)) 
            fact2= (dsrmq(Ecm,kpq)/dsrmq(2.d0*mq,kpq))**2
            
            dsib3svqqratio = fact2*(1. + (alphas/pi)*fact1)

c... O(alpha) cross section, see also Drees & Hikasa (1990), doi:10.1016/0370-2693(90)91130-4
        else   
            b0= sqrt(1.d0 - (2.*mq/Ecm)**2) ! this should be Ecm-> 2mchi, but we 
                                            ! are NR here and this prescription ensures 
                                            ! a smooth transition to the re-summed case above
            brt= (1.d0-b0)/(1.d0+b0)
            tf1=1.d0+b0**2
            tf2=dsdilog(brt)
            tf3=dsdilog(-brt)
            tf4=log(2/(1.d0+b0))*log(1/brt)
            tf5=log(b0)*log(1/brt)
            tf6=b0*log(4/(1.d0-b0**2))
            tf7=b0*log(b0)
            abeta=tf1*(4*tf2+2*tf3-3*tf4-2*tf5)-3*tf6-4*tf7
            tf8=1/(16.d0*b0)*(19.d0+2.d0*b0**2+3.d0*b0**4)*log(1/brt)
            tf9=3/8.d0*(7.d0-b0**2)
            dsib3svqqratio = 1+4/3.d0*alphas/pi*(abeta/b0+tf8+tf9)
        endif

        return

      end
