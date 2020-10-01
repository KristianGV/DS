*******************************************************************************
*** Subroutine dsddDMCRquenching returns the nuclear recoil rate Tr that    ***
*** corresponds to a given quenched energy Tquench, as wee as the           ***
*** derivative (dTr/dTquench) at the same quenched energy Tquench.          ***
*** A call to dsddDMCRquenching_set determines which quenching prescription ***
*** is used. See there for available options.                               ***
***                                                                         ***
***  Input:                                                                 ***
***    Tquench - quenched energy as observable at the detector              ***
***                                                                         ***
***  Output:                                                                ***
***    Tr - nuclear recoil energy that corresponds to Tquench               ***
***    dTrdTq - (dTr/dTquench) evalauted at Tquench                         ***
***                                                                         ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2018-06-27                                                         ***
*******************************************************************************
      subroutine dsddDMCRquenching(Tquench, Tr, dTrdTq)
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      real*8 Tquench, Tr, dTrdTq
      real*8 lnTq
      integer k,dk

c... this is the no quenching case
      Tr = Tquench   
      dTrdTq = 1.0d0      

      if (quenchhow.eq.2.and.(quenching_set)) then ! interpolate from table

        lnTq = log(Tquench)
        if ((lnTq.gt.lntqdat(nqdat)).or.(lnTq.lt.lntqdat(1))) then
           khi = nqdat
           khi = 1
           if (prtlevel.ge.1) then
             write(*,*) 'WARNING: dsddDMCRquenching has been called'
             write(*,*) 'outside the tabulated ranch of quenched energies.'
             write(*,*) 'NO QUENCHING will be used.'
           endif
           return         
        endif
        dk=1
 100    if (lntqdat(khi).lt.lnTq) then 
           khi=khi+dk
           dk=2*dk
           if (khi.lt.nqdat) goto 100
           khi=nqdat
        endif
 110    if (lntqdat(klo).gt.lnTq) then
           klo=klo-dk
           dk=2*dk
           if (klo.gt.1) goto 110
           klo=1
        endif
 120    if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if (lntqdat(k).gt.lnTq) then
              khi=k
           else
              klo=k
           endif
           goto 120
        endif
         Tr = exp(lnTrecdat(klo)+(lnTrecdat(khi)-lnTrecdat(klo))
     &          *(lnTq-lntqdat(klo))/(lntqdat(khi)-lntqdat(klo)))
         dTrdTq = exp(lndTdTdat(klo)+(lndTdTdat(khi)-lndTdTdat(klo))
     &          *(lnTq-lntqdat(klo))/(lntqdat(khi)-lntqdat(klo)))

      else ! no quenching
         if (quenchhow.ne.1.and.prtlevel.ge.1) then
            write(*,*) 'WARNING: dsddDMCRquenching has been called'
            write(*,*) 'without first calling dsddDMCRquenching_set with'
            write(*,*) 'a supported option. NO QUENCHING will be used.'
          endif 
      endif

      return
      end

