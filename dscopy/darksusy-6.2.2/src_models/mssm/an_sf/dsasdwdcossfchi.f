************************************************************************
*** SUBROUTINE dsasdwdcossfchi                                       ***
*** computes dW_{ij}/dcostheta                                       ***
*** for sfermion(1) + neutralino(2) (or chargino(2))                 ***
*** plus sfermion(1) + neutralino(2) (or chargino(2))                ***
***                                                                  *** 
*** AUTHOR: Piero Ullio, ullio@sissa.it                              ***
*** Date: 01-11-04                                                   ***
*** Modified: Joakim Edsjo, Mia Schelke                              ***
***   to include gluon final states, 2002-03-21                      ***
*** Modified: Piero Ullio                                            ***
***   to switch to ampl2 with physical polarizations, 02-07-01       ***
*** Modified: Mia Schelke                                            ***
***   to allow for initial state squarks of different family than    ***
***   final state quarks in selected cases, June 2006, Jan 2007      ***
*** modified: Piero Ullio to include a new labelling of states,      ***
***   08-05-30                                                       ***
*** modified: Joakim Edsjo to properly include inter-family mixing   ***
***   2016-10-29                                                     ***
************************************************************************

      real*8 function dsasdwdcossfchi(p,costhe,kp1,kp2)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      real*8 p,costhe
      integer kp1,kp2,kp3,kp4,kp2int,i,f,icase
      real*8 w,par
      logical wok
************************************************************************
      if (p.lt.0.0d0) then
         dsasdwdcossfchi=0.0d0
      return
      endif

***** tollerance on the position of a threshold
      thstep=1.d-15

      if(kp2.ge.kn1.and.kp2.le.kn4) kp2int=kn1
      if(kp2.ge.kcha1.and.kp2.le.kcha2) kp2int=kcha1
************************************************************************
*****
***** set the askin variables
      p12=p
      costheta=costhe
***** masses in initial state:
      mass1=mass(kp1)
      mass2=mass(kp2)
***** kinematic variables:
      call dsaskinset1
***** family indices:
      iifam(1)=ivfam(kp1)
      iifam(2)=0
***** internal degrees of freedom of the initial particles:
***** note: here we do not include particle anti-particle degrees of
*****       freedom. Only spin and colour are included here.
*****       The particle anti-particle degrees of freedom are
*****       taken into account when w is summed below.
      gg1c=1.d0
      if(iifam(1).ge.ivfam(ksu1).and.iifam(1).le.ivfam(ksqd_flav(3,2))) then
        gg1c=3.d0
      endif
      gg2c=2.d0  ! both for neutralinos and charginos since particle
                 ! antiparticle states are not included here.
***** 
*****    2 cases:
***** 1) sfermion + neutralino
      if(kp2int.eq.kn1) then
        do i=1,54
          prtial(i)=0.0d0
        enddo
        if(kp1.ne.kp1c.or.kp2int.ne.kp2c) then ! to save CPU, don't recalc
        kp1c=kp1
        kp2c=kp2int
        cgammain=.false.
        cgluonin=.false.
        ciaux=iifam(1)/10
        ciaux=iifam(1)-ciaux*10
c        do f=knue,kb ! not needed anymore (JE)
c          if(iifam(1).eq.ivfam(f)) then 
c            fsave=f
c          endif
c        enddo
c        kcfers=fsave
c...Add virtual s-channel fermions (made array by 3)
        ncfers=3
        if (iifam(1).eq.12.or.iifam(1).eq.22.or.iifam(1).eq.32) then
           kcfersv(1)=ke
           kcfersv(2)=kmu
           kcfersv(3)=ktau
        elseif (iifam(1).eq.11.or.iifam(1).eq.21.or.iifam(1).eq.31) then
           kcfersv(1)=knue
           kcfersv(2)=knumu
           kcfersv(3)=knutau
        elseif (iifam(1).eq.42.or.iifam(1).eq.52.or.iifam(1).eq.62) then
           kcfersv(1)=kd
           kcfersv(2)=ks
           kcfersv(3)=kb
        elseif (iifam(1).eq.41.or.iifam(1).eq.51.or.iifam(1).eq.61) then
           kcfersv(1)=ku
           kcfersv(2)=kc
           kcfersv(3)=kt
        else
           write(*,*) 'ERROR in DS: dsasdwdcossfchi:'
           write(*,*) 'Illegal iifam= ',iifam(1)
           write(*,*) 'Something is wrong.'
           stop
        endif
        
        if(kcfersv(1).le.ktau.and.kcfersv(1).ge.knue) then ! slepton case
          if(ciaux.eq.1) then ! sneutrino
c            ncsfert=1 ! add all three?
c            kcsfertn(1)=ksnu((fsave+1)/2)
             ncsfert=3 ! JE CORRECTION these 4 lines
             kcsfertn(1)=ksnu(1)
             kcsfertn(2)=ksnu(2)
             kcsfertn(3)=ksnu(3)
c            ncferd=1 ! add all three
c            kcferd(1)=fsave+1
             ncferd=3 ! JE CORRECTION these 4 lines
             kcferd(1)=ke
             kcferd(2)=kmu
             kcferd(3)=ktau
          endif
          if(ciaux.eq.2) then ! charged slepton
c            ncsfert=2 ! add all six?
c            kcsfertn(1)=ksl(fsave/2)
c            kcsfertn(2)=ksl(fsave/2+3)
            ncsfert=6 ! JE CORRECTION these 7 lines
            kcsfertn(1)=ksl(1)
            kcsfertn(2)=ksl(2)
            kcsfertn(3)=ksl(3)
            kcsfertn(4)=ksl(4)
            kcsfertn(5)=ksl(5)
            kcsfertn(6)=ksl(6)
c            ncferd=1 ! add all three?
c            kcferd(1)=fsave-1
            ncferd=3 ! JE CORRECTION these 4 lines
            kcferd(1)=knue
            kcferd(2)=knumu
            kcferd(3)=knutau
            cgammain=.true.
          endif
        endif
        if(kcfersv(1).le.kb.and.kcfersv(1).ge.ku) then ! squarks
          if(ciaux.eq.1) then ! up type 
c            ncsfert=2 ! add all 6?
c            kcsfertn(1)=ksqu((fsave-6+1)/2)
c            kcsfertn(2)=ksqu((fsave-6+1)/2+3)
            ncsfert=6 ! JE CORRECTION these 7 lines
            kcsfertn(1)=ksqu(1)
            kcsfertn(2)=ksqu(2)
            kcsfertn(3)=ksqu(3)
            kcsfertn(4)=ksqu(4)
            kcsfertn(5)=ksqu(5)
            kcsfertn(6)=ksqu(6)
            ncferd=3
            kcferd(1)=kd
            kcferd(2)=ks
            kcferd(3)=kb
            cgammain=.true.
            cgluonin=.true.
          endif
          if(ciaux.eq.2) then ! down type
c            ncsfert=2 ! add all 6?
c            kcsfertn(1)=ksqd((fsave-6)/2)
c            kcsfertn(2)=ksqd((fsave-6)/2+3)
            ncsfert=6 ! JE CORRECTION these 7 lines
            kcsfertn(1)=ksqd(1)
            kcsfertn(2)=ksqd(2)
            kcsfertn(3)=ksqd(3)
            kcsfertn(4)=ksqd(4)
            kcsfertn(5)=ksqd(5)
            kcsfertn(6)=ksqd(6)
            ncferd=3
            kcferd(1)=ku
            kcferd(2)=kc
            kcferd(3)=kt
            cgammain=.true.
            cgluonin=.true.
          endif
        endif
********
        endif
******************************************************** sf chi^0 -> Z f
        kp3=kz
        icase=1
        do f=1,3
c           kcfers=kcfersv(f)
           kp4=kcfersv(f)
           call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
           prtial(1)=prtial(1)+par
        enddo

**************************************************** sf chi^0 -> gamma f
        if(cgammain) then
          kp3=kgamma
          icase=7
          do f=1,3
c             kcfers=kcfersv(f)
             kp4=kcfersv(f)
             call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
             prtial(2)=prtial(2)+par
          enddo
        endif
******************************************************** sf chi^0 -> h f
        kp3=kh2
        icase=1
        do f=1,3
c           kcfers=kcfersv(f)
           kp4=kcfersv(f)
           call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
           prtial(3)=prtial(3)+par
        enddo
******************************************************** sf chi^0 -> H f
        kp3=kh1
        icase=1
        do f=1,3
c           kcfers=kcfersv(f)
           kp4=kcfersv(f)
           call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
           prtial(4)=prtial(4)+par
        enddo
******************************************************** sf chi^0 -> a f
        kp3=kh3
        icase=2
        do f=1,3
c           kcfers=kcfersv(f)
           kp4=kcfersv(f)
           call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
           prtial(5)=prtial(5)+par
        enddo
**************************************************** sf chi^0 -> gluon f
        if(cgluonin) then
          kp3=kgluon
          icase=9
          do f=1,3
c             kcfers=kcfersv(f)
             kp4=kcfersv(f)
             call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
             prtial(8)=prtial(8)+par
          enddo
        endif
********
c        kcfers=fsave ! FIX ME is this needed?
        if(ciaux.eq.1) then ! up type
          prtial(6)=0.d0
          prtial(7)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            if(kp4.le.ktau.and.kp4.ge.knue) then
c              ncsfertc=2
c              kcsfertc(1)=ksl(kp4/2)
c              kcsfertc(2)=ksl(kp4/2+3)
               ncsfertc=6 ! JE CORRECTION these 7 lines
               kcsfertc(1)=ksl(1)
               kcsfertc(2)=ksl(2)
               kcsfertc(3)=ksl(3)
               kcsfertc(4)=ksl(4)
               kcsfertc(5)=ksl(5)
               kcsfertc(6)=ksl(6)
            else
c              ncsfertc=2
c              kcsfertc(1)=ksqd((kp4-6)/2)
c              kcsfertc(2)=ksqd((kp4-6)/2+3)
               ncsfertc=6 ! JE CORRECTION these 7 lines
               kcsfertc(1)=ksqd(1)
               kcsfertc(2)=ksqd(2)
               kcsfertc(3)=ksqd(3)
               kcsfertc(4)=ksqd(4)
               kcsfertc(5)=ksqd(5)
               kcsfertc(6)=ksqd(6)
            endif
*********************************************** sfu^* chi^0 -> W^- fdbar
            kp3=kw
            icase=2
c            kcfers=kp4 ! JE CORRECTION
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par   
*********************************************** sfu^* chi^0 -> h^- fdbar
            kp3=khc
            icase=3
c            kcfers=kp4 ! JE CORRECTION
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
        elseif(ciaux.eq.2) then  ! down type
********
          prtial(6)=0.d0
          prtial(7)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            if(kp4.le.ktau.and.kp4.ge.knue) then
c              ncsfertc=1
c              kcsfertc(1)=ksnu((kp4+1)/2)
               ncsfertc=3 ! JE CORRECTION these 4 lines
               kcsfertc(1)=ksnu(1)
               kcsfertc(2)=ksnu(2)
               kcsfertc(3)=ksnu(3)
            else
c              ncsfertc=2
c              kcsfertc(1)=ksqu((kp4-6+1)/2)
c              kcsfertc(2)=ksqu((kp4-6+1)/2+3)
              ncsfertc=6 ! JE CORRECTION these 7 lines
              kcsfertc(1)=ksqu(1)
              kcsfertc(2)=ksqu(2)
              kcsfertc(3)=ksqu(3)
              kcsfertc(4)=ksqu(4)
              kcsfertc(5)=ksqu(5)
              kcsfertc(6)=ksqu(6)
            endif
**************************************************** sfd chi^0 -> W^- fu
            kp3=kw
            icase=2
c            kcfers=kp4+1 ! JE CORRECTION
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par
**************************************************** sfd chi^0 -> h^- fu
            kp3=khc
            icase=4
c            kcfers=kp4+1 ! JE CORRECTION
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
        else
          write(*,*) 'DS: wrong ciaux value =',ciaux
          write(*,*) 'DS: in dsasdwdcossfchi. program stopped.'
          stop
        endif
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sf chi^0 -> Z f
        w=w+prtial(2)            ! sf chi^0 -> gamma f
        w=w+prtial(3)            ! sf chi^0 -> h f
        w=w+prtial(4)            ! sf chi^0 -> H f
        w=w+prtial(5)            ! sf chi^0 -> a f
        w=w+prtial(6)            ! sfu^* chi^0 -> W^- fdbar
                                 ! or  sfd chi^0 -> W^- + fu
        w=w+prtial(7)            ! sfu^* chi^0 -> h^- fdbar
                                 ! or  sfd chi^0 -> h^- fu
        w=w+prtial(8)            ! sf chi^0 -> gluon f
        dsasdwdcossfchi = w      ! The weighting factor is 1
c        write(*,*) 'sfchi, prtial(1)=',prtial(1)
c        write(*,*) 'sfchi, prtial(2)=',prtial(2)
c        write(*,*) 'sfchi, prtial(3)=',prtial(3)
c        write(*,*) 'sfchi, prtial(4)=',prtial(4)
c        write(*,*) 'sfchi, prtial(5)=',prtial(5)
c        write(*,*) 'sfchi, prtial(6)=',prtial(6)
c        write(*,*) 'sfchi, prtial(7)=',prtial(7)
c        write(*,*) 'sfchi, prtial(8)=',prtial(8)
c        write(*,*) 'sfchi, w=',w

************************************************************************
*****
***** 2) sfermion + chargino
      elseif(kp2int.eq.kcha1) then
        do i=1,54
          prtial(i)=0.0d0
        enddo
        if(kp1.ne.kp1c.or.kp2int.ne.kp2c) then ! only recalcualte when needed
        kp1c=kp1
        kp2c=kp2int
        cgluonin=.false.
        ciaux=iifam(1)/10
        ciaux=iifam(1)-ciaux*10
c        do f=knue,kb! not needed anymore
c          if(iifam(1).eq.ivfam(f)) then 
c            fsave=f
c          endif
c        enddo
c        kcfers=fsave
c JE FIX. Is this enough with fsave and kcfers? Do we need to make kcfers an
c array of dimension 3
c...JE CORRECTION
c...Add final state fermions to be used in this routine
        ncfers=3
        if (iifam(1).eq.12.or.iifam(1).eq.22.or.iifam(1).eq.32) then
           kcfersv(1)=ke
           kcfersv(2)=kmu
           kcfersv(3)=ktau
        elseif (iifam(1).eq.11.or.iifam(1).eq.21.or.iifam(1).eq.31) then
           kcfersv(1)=knue
           kcfersv(2)=knumu
           kcfersv(3)=knutau
        elseif (iifam(1).eq.42.or.iifam(1).eq.52.or.iifam(1).eq.62) then
           kcfersv(1)=kd
           kcfersv(2)=ks
           kcfersv(3)=kb
        elseif (iifam(1).eq.41.or.iifam(1).eq.51.or.iifam(1).eq.61) then
           kcfersv(1)=ku
           kcfersv(2)=kc
           kcfersv(3)=kt
        else
           write(*,*) 'ERROR in DS: dsasdwdcossfchi:'
           write(*,*) 'Illegal iifam= ',iifam(1)
           write(*,*) 'Something is wrong.'
           stop
        endif

        if(kcfersv(1).le.ktau.and.kcfersv(1).ge.knue) then
          if(ciaux.eq.1) then ! sneutrino
c            ncsfert=1
c            kcsfertn(1)=ksnu((fsave+1)/2)
            ncsfert=3 ! JE CORRECTION these 4 lines
            kcsfertn(1)=ksnu(1)
            kcsfertn(2)=ksnu(2)
            kcsfertn(3)=ksnu(3)
c            ncfers=1
c            kcfersv(1)=kcfers
c            ncfers=3 ! JE CORRECTION these 4 lines ! done already
c            kcfersv(1)=kcfinal(1)
c            kcfersv(2)=kcfinal(2)
c            kcfersv(3)=kcfinal(3)
c            ncferd=1
c            kcferd(1)=fsave+1
            ncferd=3 ! JE CORRECTION these 4 lines
            kcferd(1)=ke
            kcferd(2)=kmu
            kcferd(3)=ktau
c            ncsfertc=2
c            kcsfertc(1)=ksl((fsave+1)/2)
c            kcsfertc(2)=ksl((fsave+1)/2+3)
           ncsfertc=6 ! JE CORRECTION these 7 lines
           kcsfertc(1)=ksl(1)
           kcsfertc(2)=ksl(2)
           kcsfertc(3)=ksl(3)
           kcsfertc(4)=ksl(4)
           kcsfertc(5)=ksl(5)
           kcsfertc(6)=ksl(6)
          endif
          if(ciaux.eq.2) then ! charged slepton
c            ncsfert=2
c            kcsfertn(1)=ksl(fsave/2)
c            kcsfertn(2)=ksl(fsave/2+3)
            ncsfert=6 ! JE CORRECTION these 7 lines
            kcsfertn(1)=ksl(1)
            kcsfertn(2)=ksl(2)
            kcsfertn(3)=ksl(3)
            kcsfertn(4)=ksl(4)
            kcsfertn(5)=ksl(5)
            kcsfertn(6)=ksl(6)
c            ncfers=1
c            kcfersv(1)=kcfers
c            ncfers=3 ! JE CORRECTION these 4 lines ! done already
c            kcfersv(1)=kcfinal(1)
c            kcfersv(2)=kcfinal(2)
c            kcfersv(3)=kcfinal(3)
c            ncferd=1
c            kcferd(1)=fsave-1
            ncferd=3 ! JE CORRECTION these 4 lines
            kcferd(1)=knue
            kcferd(2)=knumu
            kcferd(3)=knutau
c            ncsfertc=1
c            kcsfertc(1)=ksnu(fsave/2)
            ncsfertc=3 ! JE CORRECTION these 4 lines
            kcsfertc(1)=ksnu(1)
            kcsfertc(2)=ksnu(2)
            kcsfertc(3)=ksnu(3)
          endif
        endif
        if(kcfersv(1).le.kb.and.kcfersv(1).ge.ku) then ! squarks
          if(ciaux.eq.1) then ! up type
c            ncsfert=2
c            kcsfertn(1)=ksqu((fsave-6+1)/2)
c            kcsfertn(2)=ksqu((fsave-6+1)/2+3)
            ncsfert=6 ! JE CORRECTION these 7 lines
            kcsfertn(1)=ksqu(1)
            kcsfertn(2)=ksqu(2)
            kcsfertn(3)=ksqu(3)
            kcsfertn(4)=ksqu(4)
            kcsfertn(5)=ksqu(5)
            kcsfertn(6)=ksqu(6)
c            ncfers=3 ! done already
c            kcfersv(1)=ku
c            kcfersv(2)=kc
c            kcfersv(3)=kt
            ncferd=3
            kcferd(1)=kd
            kcferd(2)=ks
            kcferd(3)=kb
            ncsfertc=6
            do i=1,6
              kcsfertc(i)=ksqd(i)
            enddo
            cgluonin=.true.
          endif
          if(ciaux.eq.2) then ! down type
c            ncsfert=2
c            kcsfertn(1)=ksqd((fsave-6)/2)
c            kcsfertn(2)=ksqd((fsave-6)/2+3)
            ncsfert=6 ! JE CORRECTION these 7 lines
            kcsfertn(1)=ksqd(1)
            kcsfertn(2)=ksqd(2)
            kcsfertn(3)=ksqd(3)
            kcsfertn(4)=ksqd(4)
            kcsfertn(5)=ksqd(5)
            kcsfertn(6)=ksqd(6)
c            ncfers=3 ! done already
c            kcfersv(1)=kd
c            kcfersv(2)=ks
c            kcfersv(3)=kb
            ncferd=3
            kcferd(1)=ku
            kcferd(2)=kc
            kcferd(3)=kt
            ncsfertc=6
            do i=1,6
              kcsfertc(i)=ksqu(i)
            enddo
            cgluonin=.true.
          endif
        endif
********
        endif
********
        if(ciaux.eq.1) then ! up type
************************************************* sfu^* chi^+ -> Z fdbar
          kp3=kz
          prtial(1)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=1
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(1)=prtial(1)+par
          enddo
********************************************* sfu^* chi^+ -> gamma fdbar
          kp3=kgamma
          prtial(2)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=5
            call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
            prtial(2)=prtial(2)+par
          enddo
************************************************* sfu^* chi^+ -> h fdbar
          kp3=kh2
          prtial(3)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=1
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(3)=prtial(3)+par
          enddo
************************************************* sfu^* chi^+ -> H fdbar
          kp3=kh1
          prtial(4)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=1
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(4)=prtial(4)+par
          enddo
************************************************* sfu^* chi^+ -> a fdbar
          kp3=kh3
          prtial(5)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=2
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(5)=prtial(5)+par
          enddo
**************************************************** sfu chi^+ -> W^+ fu
          kp3=kw
          prtial(6)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=4
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(6)=prtial(6)+par
          enddo
*********************************************** sfu^* chi^+ -> W^+ fubar
          kp3=kw
          prtial(7)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=3
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
**************************************************** sfu chi^+ -> H^+ fu
          kp3=khc
          prtial(8)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=5
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
            prtial(8)=prtial(8)+par
          enddo
*********************************************** sfu^* chi^+ -> H^+ fubar
          kp3=khc
          prtial(9)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=4
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(9)=prtial(9)+par
          enddo
*********************************************** sfu^* chi^+ -> g fdbar
          if(cgluonin) then
            kp3=kgluon
            prtial(10)=0.d0
            do i=1,ncferd
              kp4=kcferd(i)
              icase=6
              call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
              prtial(10)=prtial(10)+par
            enddo
          else
            prtial(10)=0.0d0
          endif
**************************************************** sum partial results
          w=0.d0
          w=w+prtial(1)            ! sfu^* chi^+ -> Z fdbar
          w=w+prtial(2)            ! sfu^* chi^+ -> gamma fdbar
          w=w+prtial(3)            ! sfu^* chi^+ -> h fdbar
          w=w+prtial(4)            ! sfu^* chi^+ -> H fdbar
          w=w+prtial(5)            ! sfu^* chi^+ -> a fdbar
          w=w+prtial(6)            ! sfu chi^+ -> W^+ fu
          w=w+prtial(7)            ! sfu^* chi^+ -> W^+ fubar
          w=w+prtial(8)            ! sfu chi^+ -> H^+ fu
          w=w+prtial(9)            ! sfu^* chi^+ -> H^+ fubar
          w=w+prtial(10)           ! sfu^* chi^+ -> g fdbar
          dsasdwdcossfchi = w*0.5d0 ! 0.5d0 <- we should return weff for sf chi+ ann
                                    ! with sf and chi+ combined part. and anti-particle states
          
********
        elseif(ciaux.eq.2) then ! down type
****************************************************** sfd chi^+ -> Z fu
          kp3=kz
          prtial(1)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
            prtial(1)=prtial(1)+par
          enddo
************************************************** sfd chi^+ -> gamma fu
          kp3=kgamma
          prtial(2)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=8
            call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
            prtial(2)=prtial(2)+par
          enddo
****************************************************** sfd chi^+ -> h fu
          kp3=kh2
          prtial(3)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)          
            prtial(3)=prtial(3)+par
          enddo
****************************************************** sfd chi^+ -> H fu
          kp3=kh1
          prtial(4)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)          
            prtial(4)=prtial(4)+par
          enddo  
****************************************************** sfd chi^+ -> a fu
          kp3=kh3
          prtial(5)=0.d0
          do i=1,ncferd
            kp4=kcferd(i)
            icase=3
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)          
            prtial(5)=prtial(5)+par
          enddo
**************************************************** sfd chi^+ -> W^+ fd
          kp3=kw
          prtial(6)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=5  
            call dsaschicasea(kp1,kp2,kp3,kp4,icase,par)
c            write(*,*) 'sfd chi+ -> W+ fd: ',kp1,kp2,kp3,kp4,par
            prtial(6)=prtial(6)+par
          enddo
**************************************************** sfd chi^+ -> h^+ fd
          kp3=khc
          prtial(7)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=6
            call dsaschicasec(kp1,kp2,kp3,kp4,icase,par)
            prtial(7)=prtial(7)+par
          enddo
*********************************************** sfd^* chi^+ -> W^+ fdbar
          kp3=kw
          prtial(8)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=4
            call dsaschicaseb(kp1,kp2,kp3,kp4,icase,par)
            prtial(8)=prtial(8)+par
          enddo
*********************************************** sfd^* chi^+ -> h^+ fdbar
          kp3=khc
          prtial(9)=0.d0
          do i=1,ncfers
            kp4=kcfersv(i)
            icase=5
            call dsaschicased(kp1,kp2,kp3,kp4,icase,par)
            prtial(9)=prtial(9)+par
          enddo
*********************************************** sfd chi^+ -> g fu
          if(cgluonin) then
            kp3=kgluon
            prtial(10)=0.d0
            do i=1,ncferd
              kp4=kcferd(i)
              icase=10
              call dsaschizero(kp1,kp2,kp3,kp4,icase,par)
              prtial(10)=prtial(10)+par
            enddo 
          else
            prtial(10)=0.0d0
          endif
**************************************************** sum partial results
          w=0.d0
          w=w+prtial(1)            ! sfd chi^+ -> Z fu
          w=w+prtial(2)            ! sfd chi^+ -> gamma fu
          w=w+prtial(3)            ! sfd chi^+ -> h fu
          w=w+prtial(4)            ! sfd chi^+ -> H fu
          w=w+prtial(5)            ! sfd chi^+ -> a fu
          w=w+prtial(6)            ! sfd chi^+ -> W^+ fd
          w=w+prtial(7)            ! sfd chi^+ -> h^+ fd
          w=w+prtial(8)            ! sfd^* chi^+ -> W^+ fdbar
          w=w+prtial(9)            ! sfd^* chi^+ -> h^+ fdbar
          w=w+prtial(10)           ! sfd chi^+ -> g fu
          dsasdwdcossfchi = w*0.5d0 
               ! 0.5d0 <- we should return weff for sf chi+ ann
               ! with sf and chi+ combined part. and anti-particle states
c          write(*,*) 'sf-cha: prtial(1) = ',prtial(1)
c          write(*,*) 'sf-cha: prtial(2) = ',prtial(2)
c          write(*,*) 'sf-cha: prtial(3) = ',prtial(3)
c          write(*,*) 'sf-cha: prtial(4) = ',prtial(4)
c          write(*,*) 'sf-cha: prtial(5) = ',prtial(5)
c          write(*,*) 'sf-cha: prtial(6) = ',prtial(6)
c          write(*,*) 'sf-cha: prtial(7) = ',prtial(7)
c          write(*,*) 'sf-cha: prtial(8) = ',prtial(8)
c          write(*,*) 'sf-cha: prtial(9) = ',prtial(9)
c          write(*,*) 'sf-cha: prtial(10)= ',prtial(10)
********
        else
          write(*,*) 'DS: wrong ciaux value =',ciaux
          write(*,*) 'DS: in dsasdwdcossfchi. program stopped.'
          stop
        endif
      else 
        write(*,*) 'DS: in dsasdwdcossfchi wrong part. in init. state:'
        write(*,*) pname(kp2),' instead of neutralino or chargino.'
        write(*,*) 'DS: program stopped.'
        stop
      endif 


      if(aszeroprint) then
      wok=.true.
      do i=1,10
        if (prtial(i).lt.0.0d0) wok=.false.
      enddo
      if (.not.wok) then
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: dsasdwdcossfchi called with kp1=',kp1,'kp2=',kp2
        write(*,*) 'DS:  p=',p,' costh=',costhe,' w=',w
        do i=1,10
          write(*,*) 'DS: i= ',i,' prtial= ',prtial(i)
        enddo
      endif
      endif
  
      if (w.lt.0.0d0) w=abs(w)
      return
      end








