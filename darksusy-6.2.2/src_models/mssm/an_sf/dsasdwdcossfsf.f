************************************************************************
*** SUBROUTINE dsasdwdcossfsf                                        ***
*** computes dW_{ij}/dcostheta                                       ***
*** for sfermion(1) + antisfermion(2) plus sfermion(1) + sfermion(2) ***
***                                                                  *** 
*** AUTHOR: Piero Ullio, ullio@sissa.it                              ***
*** Date: 01-08-10                                                   ***
*** modified: Joakim Edsjo, Mia Schelke to include squarks with      ***
***   gauge and Higgs boson final states and gluons                  ***
***   02-05-22                                                       ***
***   corrected switching of initial states fixed 020613 (edsjo)     ***
*** modified: Piero Ullio                                            ***
***   02-03-22                                                       ***
*** modified: Piero Ullio                                            ***
***   02-07-01                                                       ***
*** modified: Mia Schelke to include squark mixing:                  ***
***   diff family i.s. sq's with bosons in final state               ***
***   Jan 2007                                                       ***
*** modified: Piero Ullio to include mixing in fermion final states  ***
***   and a new labelling of states                                  ***
***   08-05-30                                                       ***
***   Joakim Edsjo, correction in slepton-squark part (case 3)       ***
***   2015-04-08                                                     ***
***   Joakim Edsjo, made more general (full 6x6 mixing),             ***
***   Cleaned up a bit (not as many if's needed with full 6x6        ***
***   mixing). 2016-11-16                                            ***
************************************************************************

      real*8 function dsasdwdcossfsf(p,costhe,kp1,kp2)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      include 'dsandwcom.h'
      include 'dsidtag.h'
      real*8 p,costhe
      integer kp1,kp2,kp3,kp4,itmp,icase
      integer kp1i,kp2i
      integer f,f1,f2,i,fsave,f2save
      real*8 w,result,par
      logical wok

      data fsave,f2save/0,0/ ! to avoid compiler warnings
************************************************************************
      if (p.lt.0.0d0) then
         dsasdwdcossfsf=0.0d0
      return
      endif

***** tollerance on the position of a threshold and on a negative
***** fermion channel
      thstep=1.d-15
      fertoll=1.d-10

c...Since we want to be able to switch kp1 and kp2 without having
c...them changed on return, let's use some temporary variables
      kp1i=kp1
      kp2i=kp2

************************************************************************
*****
***** set the askin variables
      p12=p
      costheta=costhe

 100  continue
***** family indices:
      iifam(1)=ivfam(kp1i)
      iifam(2)=ivfam(kp2i)

***** decide on type for kp1i and kp2i, i.e.:
*****   itype(iii)=ivfam(ku) for up-type (s)quark
*****   itype(iii)=ivfam(kd) for down-type (s)quark
*****   itype(iii)=ivfam(iii) for (s)leptons, JE CORR ivfam(ke) for (s)leptons
      itype(1)=ivtype(kp1i)
      itype(2)=ivtype(kp2i)

***** NOTE: boson final states are computed just for kp1i up-type and 
***** kp2i down-type, not viceversa, so in the latter case switch them
***** reordering this squark first if it appears in case 3
      if(itype(1).eq.(itype(2)+1)) then
        itmp=kp1i
        kp1i=kp2i
        kp2i=itmp
        goto 100
      endif 
***** reordering this squark first if it appears in case 3
      if(itype(2).ge.ivtype(ku).and.itype(1).lt.ivtype(ku)) then
        itmp=kp1i
        kp1i=kp2i
        kp2i=itmp
        goto 100
      endif 

***** masses in initial state:
      mass1=mass(kp1i)
      mass2=mass(kp2i)

***** kinematic variables:
      call dsaskinset1
***** 
*****    3 cases:
c======================================================================
c======================================================================
***** 1) sfermions with same type index, i.e.:
*****    slepton + slepton* (both charged or both neutral) 
*****    up-type-squark1 + up-type-squark2  or 
*****    down-type-squark1 + down-type-squark2 
c======================================================================
c======================================================================
      if(itype(1).eq.itype(2)) then
        do i=1,66
          prtial(i)=0.0d0
        enddo
ccc
ccc check whether to reload final/intermediate state arrays/variables
ccc
        if(kp1i.ne.kp1s.or.kp2i.ne.kp2s) then ! recalculate when needed
        kp1s=kp1i
        kp2s=kp2i
ccc
ccc do the resetting in the various cases:
ccc
        do i=1,66
          chon(i)=.false.
        enddo
        gluonin=.false.
        gammain=.false.
        neutcurr=.false.
        nosneutrinov=.true.
ccc
        if(itype(1).ge.ivtype(knue).and.itype(1).le.ivtype(ktau)) then
c============================================================
ccc input:   \tilde{l}_1,2(i) \tilde{l}^*_1,2(i)
c============================================================
ccc
ccc fermion channels:
ccc
        i=0
***** f fbar 
        do f=knue,kb
          i=i+1
          kp3in(i)=f
          kp4in(i)=f
          chon(i)=.true.
       enddo                    ! i=1-12 now filled

c..JE CORR following states added
***** up-type quark(i), up-type antiquark(j) i,j matching initial state
        i=12
        do f1=1,3 ! i=13-18 filled here
        do f2=1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
c          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4)) then ! JE CORR
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
c          endif ! JE CORR
          endif         
        enddo
        enddo
c...JE CORR following states added        
***** down-type quark(i), down-type antiquark(j) i.ne.j 
        do f1=1,3 ! i=19-24 filled here
        do f2=1,3 
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo
       
***** lepton1,  lepton2
c...JE CORR Added all families as final states, ch 25-30
        i=24
        if (itype(1).eq.ivtype(knue)) then ! sneutrino
           do f1=1,3
              do f2=f1,3
                 i=i+1
                 chon(i)=.true.
                 kp3in(i)=knu(f1)
                 kp4in(i)=knu(f2)
              enddo
           enddo
           fsave=knue           ! just for type testing later
        else                    ! charged slepton
           do f1=1,3
              do f2=f1,3
                 i=i+1
                 chon(i)=.true.
                 kp3in(i)=kl(f1)
                 kp4in(i)=kl(f2)
              enddo
           enddo
           fsave=ke ! just for type testing
        endif
        
c..Old one family code. Replaced by code above        
c        i=25 ! lepton_1 lepton_1 (same type)
c        chon(i)=.true.
c        do f=knue,ktau
c          if(itype(1).eq.ivtype(f)) then 
c            kp3in(i)=f
c            kp4in(i)=f
c            fsave=f
c          endif
c        enddo      

c...Mixed final state lepton states (OK even w/o flavour violation
c...for initial states with different flavours)        
        i=54
        do f1=1,3 ! i=55-60 filled here, nu_i nu_j 
           do f2=1,3
              kp3=knu(f1)
              kp4=knu(f2)
              if (kp3.ne.kp4) then
                 i=i+1
                 if (itype(1).eq.itype(2)) then ! same type (up/down)
                    kp3in(i)=kp3 ! sl sl* ->
                    kp4in(i)=kp4
                    chon(i)=.true.
                 endif
              endif
           enddo
        enddo

        do f1=1,3 ! i=61-66 filled here, l_i l_j
           do f2=1,3
              kp3=kl(f1)
              kp4=kl(f2)
              if (kp3.ne.kp4) then
                 i=i+1
                 if (itype(1).eq.itype(2)) then ! same type (up/down)
                    kp3in(i)=kp3 ! sl sl* ->
                    kp4in(i)=kp4
                    chon(i)=.true.
                 endif
              endif
           enddo
        enddo
        
        
ccc
ccc gauge boson channels:
ccc
        if(dabs(echarg(fsave)).gt.1.d-15) gammain=.true.
        neutcurr=.true.
ccc
ccc t- and u-channels:
ccc
        if(dabs(echarg(fsave)).gt.1.d-15) then ! for charged sleptons:
c          ksfertc(1)=ksnu(fsave/2)
c          nsfertc=1
          nsfertc=3 ! JE CORRECTION 4 lines
          ksfertc(1)=ksnu(1)
          ksfertc(2)=ksnu(2)
          ksfertc(3)=ksnu(3)
          nsferuc=0
c          ksfertn(1)=ksl(fsave/2)
c          ksfertn(2)=ksl(fsave/2+3)
c          nsfertn=2
          nsfertn=6 ! JE CORRECTION 7 lines
          ksfertn(1)=ksl(1)
          ksfertn(2)=ksl(2)
          ksfertn(3)=ksl(3)
          ksfertn(4)=ksl(4)
          ksfertn(5)=ksl(5)
          ksfertn(6)=ksl(6)
c          ksferun(1)=ksl(fsave/2)
c          ksferun(2)=ksl(fsave/2+3)
c          nsferun=2
          nsferun=6 ! JE CORRECTION 7 lines
          ksferun(1)=ksl(1)
          ksferun(2)=ksl(2)
          ksferun(3)=ksl(3)
          ksferun(4)=ksl(4)
          ksferun(5)=ksl(5)
          ksferun(6)=ksl(6)
        else ! for sneutrinos:
          nsfertc=0
c          ksferuc(1)=ksl((fsave+1)/2)
c          ksferuc(2)=ksl((fsave+1)/2+3)
c          nsferuc=2
          nsferuc=6 ! JE CORRECTION 7 lines
          ksferuc(1)=ksl(1)
          ksferuc(2)=ksl(2)
          ksferuc(3)=ksl(3)
          ksferuc(4)=ksl(4)
          ksferuc(5)=ksl(5)
          ksferuc(6)=ksl(6)
c          ksfertn(1)=ksnu((fsave+1)/2)
c          nsfertn=1
          nsfertn=3 ! JE CORRECTION 4 lines
          ksfertn(1)=ksnu(1)
          ksfertn(2)=ksnu(2)
          ksfertn(3)=ksnu(3)
c          ksferun(1)=ksnu((fsave+1)/2)
c          nsferun=1
          nsferun=3 ! JE CORRECTION 4 lines
          ksferun(1)=ksnu(1)
          ksferun(2)=ksnu(2)
          ksferun(3)=ksnu(3)
          nosneutrinov=.false.
        endif
ccc
        gg1=1.d0 
        gg2=1.d0
        chcol=1
        cfactini=1.d0
        cfactfin=3.d0
        endif
ccc
c        if(iifam(1).eq.iifam(2).and.itype(1).eq.ivtype(ku)) then
        if(itype(1).eq.ivtype(ku)) then ! JE CORR
c============================================================
ccc input:   \tilde{u}_1,2(i) \tilde{u}^*_1,2(j) (all i,j)
c============================================================
ccc
ccc fermion channels:
ccc
        i=0
***** f fbar
        do f=knue,kb
          i=i+1
          kp3in(i)=f
          kp4in(i)=f
          chon(i)=.true.
        enddo

c..JE CORR following states added
***** up-type quark(i), up-type antiquark(j) i,j matching initial state
        i=12
        do f1=1,3 ! i=13-18 filled here
        do f2=1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
c          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4)) then ! JE CORR
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
c          endif ! JE CORR
          endif         
        enddo
        enddo
        
***** down-type quark(i), down-type antiquark(j) i.ne.j 
        do f1=1,3 ! i=19-24 filled here
        do f2=1,3 
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo

***** up-type quark1, up-type quark2 
        do f1=1,3 ! i=25-30 filled here
        do f2=f1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          i=i+1  
c          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4).or. ! JE CORR
c     &       iifam(1).eq.ivfam(kp4).and.iifam(2).eq.ivfam(kp3)) then
            chon(i)=.true.
            kp3in(i)=kp3
            kp4in(i)=kp4
c            fsave=f1 ! JE CORR
c          endif
        enddo
        enddo
ccc
ccc boson channels:
ccc
        gluonin=.true.
        gammain=.true.
        neutcurr=.true.
ccc
ccc t- and u-channels:
ccc
        nsfertc=0 
        do i=1,6
          ksferuc(i)=ksqd(i)
        enddo
        nsferuc=6
c        ksfertn(1)=ksqu(fsave)
c        ksfertn(2)=ksqu(fsave+3)
c        nsfertn=2
        nsfertn=6 ! JE CORRECTION 7 lines
        ksfertn(1)=ksqu(1)
        ksfertn(2)=ksqu(2)
        ksfertn(3)=ksqu(3)
        ksfertn(4)=ksqu(4)
        ksfertn(5)=ksqu(5)
        ksfertn(6)=ksqu(6)
c        ksferun(1)=ksqu(fsave)
c        ksferun(2)=ksqu(fsave+3)
c        nsferun=2
        nsferun=6 ! JE CORRECTION 7 lines
        ksferun(1)=ksqu(1)
        ksferun(2)=ksqu(2)
        ksferun(3)=ksqu(3)
        ksferun(4)=ksqu(4)
        ksferun(5)=ksqu(5)
        ksferun(6)=ksqu(6)
ccc
ccc degrees of freedom and color factors
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
ccc
        endif
ccc
c        if(iifam(1).eq.iifam(2).and.itype(1).eq.ivtype(kd)) then
        if(itype(1).eq.ivtype(kd)) then ! JE CORR
c============================================================
ccc input:   \tilde{d}_1,2(i) \tilde{d}^*_1,2(j) (all i,j)
c============================================================
        i=0
***** f fbar
        do f=knue,kb
          i=i+1
          kp3in(i)=f
          kp4in(i)=f
          chon(i)=.true.
        enddo
***** up-type quark(i), up-type antiquark(j) i.ne.j 
        do f1=1,3
        do f2=1,3
          kp3=ku+2*(f1-1)
          kp4=ku+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
          endif         
        enddo
        enddo

c...JE CORR added d_i d_j final states to be more general
***** down-type quark(i), down-type antiquark(j) i,j matching initial state 
        do f1=1,3
        do f2=1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          if(kp3.ne.kp4) then
          i=i+1
c          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4)) then ! JE CORR
          chon(i)=.true.
          kp3in(i)=kp3
          kp4in(i)=kp4
c          endif ! JE CORR
          endif         
        enddo
        enddo
      
***** down-type quark1, down-type quark2 
        do f1=1,3
        do f2=f1,3
          kp3=kd+2*(f1-1)
          kp4=kd+2*(f2-1)
          i=i+1
c          if(iifam(1).eq.ivfam(kp3).and.iifam(2).eq.ivfam(kp4).or. ! JE CORR
c     &       iifam(1).eq.ivfam(kp4).and.iifam(2).eq.ivfam(kp3)) then
            chon(i)=.true.
            kp3in(i)=kp3
            kp4in(i)=kp4
c            fsave=f1 ! JE COEE
c          endif ! JE CORR
        enddo
        enddo  
ccc
ccc boson channels:
ccc
        gluonin=.true.
        gammain=.true.
        neutcurr=.true.
ccc
ccc t- and u-channels:
ccc
        do i=1,6
          ksfertc(i)=ksqu(i)
        enddo
        nsfertc=6
        nsferuc=0
c        ksfertn(1)=ksqd(fsave)
c        ksfertn(2)=ksqd(fsave+3)
c        nsfertn=2
        nsfertn=6 ! JE CORRECTION 7 lines
        ksfertn(1)=ksqd(1)
        ksfertn(2)=ksqd(2)
        ksfertn(3)=ksqd(3)
        ksfertn(4)=ksqd(4)
        ksfertn(5)=ksqd(5)
        ksfertn(6)=ksqd(6)
c        ksferun(1)=ksqd(fsave)
c        ksferun(2)=ksqd(fsave+3)
c        nsferun=2
        nsferun=6 ! JE CORRECTION 7 lines
        ksferun(1)=ksqd(1)
        ksferun(2)=ksqd(2)
        ksferun(3)=ksqd(3)
        ksferun(4)=ksqd(4)
        ksferun(5)=ksqd(5)
        ksferun(6)=ksqd(6)
ccc
ccc degrees of freedom and color factors
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
ccc
        endif
ccc
ccc
        endif
*****
******************************************************* sf sf* -> f fbar
        do i=1,12 
          if(chon(i)) then
          kp3=kp3in(i)
          kp4=kp4in(i)
          if(chcol.eq.1) then
            icase=1
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            if(ncolor(i).gt.2.d0) then
              result=result*cfactfin 
            endif
          elseif(chcol.eq.2) then
             icase=1
             call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
          else
            write(*,*) 'chcol not set correctly, = ',chcol
            stop
          endif  
          prtial(i)=result
          endif
        enddo  

c...Add mixed slepton final states sl sl* -> ...
        do i=55,66
          if(chon(i)) then
          kp3=kp3in(i)
          kp4=kp4in(i)
          if(chcol.eq.1) then
            icase=3
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
          elseif(chcol.eq.2) then ! only for sleptons
c             icase=1
c             call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
             result=0.d0
          else
            write(*,*) 'chcol not set correctly, = ',chcol
            stop
          endif  
          prtial(i)=result
          endif
        enddo  
        
***************************************************** sf sf* -> q1 q2bar
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state, not in the final state: (JE CORR: this is more general now)
*****
        if(chcol.eq.2) then  ! JE FIX ME. Do we want this for sleptons too?
        do i=13,24
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=1
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result
          endif
        enddo
        endif
************************************************* sl sl -> lepton lepton
        if(chcol.eq.1) then  
          do i=25,30
             if(chon(i)) then
                kp3=kp3in(i)
                kp4=kp4in(i)
                icase=1
                call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
                prtial(i)=result
             endif
          enddo
        endif 
********************************************************* sf sf -> q1 q2
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state, not in the final state: (JE CORR this is more general now)
*****
        if(chcol.eq.2) then  
        do i=25,30
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=1
            call dsasferecol(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result
          endif
        enddo
        endif
******************************************************* sf sf* ->  w+ w-
        kp3=kw
        kp4=kw
        if(nsferuc.ne.0) then
          icase=2
        else
          icase=3
        endif
        call dsasgbgb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(31)=par
******************************************************** sf sf* -> h- w+
        kp3=kw
        kp4=khc
        if(nsferuc.ne.0) then
          icase=5
        else
          icase=6
        endif
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(47)=2.d0*par    !factor of 2 to include sf sf* -> h+ w-
******************************************************** sf sf* -> h+ h-
        kp3=khc
        kp4=khc
        if(nsferuc.ne.0) then
          icase=6
        else
          icase=7
        endif
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(48)=par
****************************************************** sf sf* -> z gamma
        if(gammain) then 
          kp3=kz
          kp4=kgamma
          icase=2
          call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(33)=par
************************************************** sf sf* -> gamma gamma
          kp3=kgamma
          kp4=kgamma
          icase=1
          call dsasgbgb2exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(34)=par
****************************************************** sf sf* -> gamma h
          kp3=kgamma
          kp4=kh2
          icase=7
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(38)=par
****************************************************** sf sf* -> gamma H
          kp3=kgamma
          kp4=kh1
          icase=7
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(39)=par
****************************************************** sf sf* -> gamma a
          kp3=kgamma
          kp4=kh3
          icase=7
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(40)=par
        endif  
********************************************************** sf sf* -> z z
        if(neutcurr) then
          kp3=kz
          kp4=kz
          icase=4
          call dsasgbgb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(32)=par
********************************************************** sf sf* -> z h
          kp3=kz
          kp4=kh2
          icase=8
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(35)=par
********************************************************** sf sf* -> z H
          kp3=kz
          kp4=kh1
          icase=8
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(36)=par
********************************************************** sf sf* -> z a
          kp3=kz
          kp4=kh3
          icase=9
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(37)=par
********************************************************** sf sf* -> h h
          kp3=kh2
          kp4=kh2
          icase=3
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(41)=par
********************************************************** sf sf* -> h a
          kp3=kh2
          kp4=kh3
          icase=4
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(42)=par
********************************************************** sf sf* -> H a
          kp3=kh1
          kp4=kh3
          icase=4
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(43)=par
********************************************************** sf sf* -> a a
          kp3=kh3
          kp4=kh3
          icase=5
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(44)=par
********************************************************** sf sf* -> H h
          kp3=kh1
          kp4=kh2
          icase=3
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(45)=par
********************************************************** sf sf* -> H H
          kp3=kh1
          kp4=kh1
          icase=3
          call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(46)=par
        endif
************************************************** sf sf* -> gluon gluon
        if(gluonin) then
          kp3=kgluon
          kp4=kgluon
          icase=2
          call dsasgbgb2exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(49)=par  
************************************************** sf sf* -> Z gluon
          kp3=kz
          kp4=kgluon
          icase=3
          call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(50)=par
************************************************** sf sf* -> gamma gluon
          kp3=kgamma
          kp4=kgluon
          icase=3
          call dsasgbgb2exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(51)=par
************************************************** sf sf* -> gluon H1
          kp3=kgluon
          kp4=kh1
          icase=10
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(52)=par
************************************************** sf sf* -> gluon H2
          kp3=kgluon
          kp4=kh2
          icase=10
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(53)=par
************************************************** sf sf* -> gluon H3
          kp3=kgluon
          kp4=kh3
          icase=10
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(54)=par
        endif  
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sf sf* -> nue nuebar
        w=w+prtial(2)            ! sf sf* -> e+ e-
        w=w+prtial(3)            ! sf sf* -> numu numubar
        w=w+prtial(4)            ! sf sf* -> mu+ mu-
        w=w+prtial(5)            ! sf sf* -> nutau nutaubar
        w=w+prtial(6)            ! sf sf* -> tau+ tau-
        w=w+prtial(7)            ! sf sf* -> u ubar
        w=w+prtial(8)            ! sf sf* -> d dbar
        w=w+prtial(9)            ! sf sf* -> c cbar
        w=w+prtial(10)           ! sf sf* -> s sbar
        w=w+prtial(11)           ! sf sf* -> t tbar
        w=w+prtial(12)           ! sf sf* -> b bbar
        w=w+prtial(13)           ! sf sf* -> u cbar
        w=w+prtial(14)           ! sf sf* -> u tbar
        w=w+prtial(15)           ! sf sf* -> c ubar
        w=w+prtial(16)           ! sf sf* -> c tbar
        w=w+prtial(17)           ! sf sf* -> t ubar
        w=w+prtial(18)           ! sf sf* -> t cbar
        w=w+prtial(19)           ! sf sf* -> d sbar
        w=w+prtial(20)           ! sf sf* -> d bbar
        w=w+prtial(21)           ! sf sf* -> s dbar
        w=w+prtial(22)           ! sf sf* -> s bbar
        w=w+prtial(23)           ! sf sf* -> b dbar
        w=w+prtial(24)           ! sf sf* -> b sbar
        w=w+prtial(25)           ! sf sf -> lepton_1 lepton_1 or u u or d d
        w=w+prtial(26)           ! sf sf -> lepton_1 lepton 2 or u c  or  d s 
        w=w+prtial(27)           ! sf sf -> lepton_1 lepton 3 or u t  or  d b
        w=w+prtial(28)           ! sf sf -> lepton_2 lepton_2 or c c  or  s s
        w=w+prtial(29)           ! sf sf -> lepton_2 lepton_3 or  c t  or  s b
        w=w+prtial(30)           ! sf sf -> lepton_3 lepton_3 or t t  or  b b
        w=w+prtial(31)           ! sf sf* ->  w+ w-
        w=w+prtial(32)           ! sf sf* -> z z
        w=w+prtial(33)           ! sf sf* -> z gamma
        w=w+prtial(34)           ! sf sf* -> gamma gamma
        w=w+prtial(35)           ! sf sf* -> z h
        w=w+prtial(36)           ! sf sf* -> z H
        w=w+prtial(37)           ! sf sf* -> z a
        w=w+prtial(38)           ! sf sf* -> gamma h
        w=w+prtial(39)           ! sf sf* -> gamma H
        w=w+prtial(40)           ! sf sf* -> gamma a
        w=w+prtial(41)           ! sf sf* -> h h
        w=w+prtial(42)           ! sf sf* -> h a
        w=w+prtial(43)           ! sf sf* -> H a
        w=w+prtial(44)           ! sf sf* -> a a
        w=w+prtial(45)           ! sf sf* -> H h
        w=w+prtial(46)           ! sf sf* -> H H
        w=w+prtial(47)           ! sf sf* -> h- w+ and h+ w-
        w=w+prtial(48)           ! sf sf* -> h+ h-
        w=w+prtial(49)           ! sf sf* -> gluon gluon
        w=w+prtial(50)           ! sf sf* -> Z gluon
        w=w+prtial(51)           ! sf sf* -> gamma gluon
        w=w+prtial(52)           ! sf sf* -> gluon H
        w=w+prtial(53)           ! sf sf* -> gluon h
        w=w+prtial(54)           ! sf sf* -> gluon A
        w=w+prtial(55)           ! sf sf* -> nu_e nu_mu-bar
        w=w+prtial(56)           ! sf sf* -> nu_e nu_tau-bar
        w=w+prtial(57)           ! sf sf* -> nu_mu nu_e-bar
        w=w+prtial(58)           ! sf sf* -> nu_mu nu_tau-bar
        w=w+prtial(59)           ! sf sf* -> nu_tau nu_e-bar
        w=w+prtial(60)           ! sf sf* -> nu_tau nu_mu-bar
        w=w+prtial(61)           ! sf sf* -> e- mu+
        w=w+prtial(62)           ! sf sf* -> e- tau+
        w=w+prtial(63)           ! sf sf* -> mu- e+
        w=w+prtial(64)           ! sf sf* -> mu- tau+
        w=w+prtial(65)           ! sf sf* -> tau- e+
        w=w+prtial(66)           ! sf sf* -> tau- mu+

c        do i=1,66 ! JE TMP
c           write(*,'(A,I2,1x,I2,1x,I3,1x,E12.6,1x,E12.6,1x,L1)')
c     &       'sfsf case 1: ',kp1,kp2,i,mass1,prtial(i),chon(i)
c        enddo
        
c
c check for large negative terms in the fermion final states
c
        wok=.true.
        do i=1,30
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo

        do i=55,66
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo
c
c write error message:
c
        if(aszeroprint) then
        if (.not.wok) then
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: large negative term in dsasdwdcossfsf'
          write(*,*) 'DS: called with kp1i=',kp1i,' kp2i=',kp2i
          write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
          do i=1,66
            write(*,*) 'DS: i=',i,' prtial=',prtial(i)
          enddo
        endif
        endif
c
        dsasdwdcossfsf = w*0.5d0 ! 0.5d0 <- we should return weff for 
              ! sf sf ann with sf combined part. and anti-particle state
c======================================================================
c======================================================================
*****
***** 2) sneutrino + slepton or
*****    up-type squark + down-type squark
c======================================================================
c======================================================================
      elseif(abs(itype(1)-itype(2)).eq.1) then
        do i=1,39
          prtial(i)=0.0d0
        enddo
ccc
ccc check whether to reload final/intermediate state arrays/variables
ccc
        if(kp1i.ne.kp1s.or.kp2i.ne.kp2s) then
        kp1s=kp1i
        kp2s=kp2i
*****
***** check whether to include gluon final states
***** 
        gluonin=.false.
ccc

        if(itype(1).ge.ivtype(knue).and.itype(1).le.ivtype(ktau)) then
c============================================================
ccc input:   \tilde{nu}(i) \tilde{l}^*_1,2(j)
c============================================================
ccc
ccc fermion channels: 
ccc    neutrino(j) lepton^+(j) -- all allowed 
ccc    up-type quark(j) down-type antiquark(k) -- all allowed
ccc
ccc    neutrino(i) lepton^-(j) -- all three families
        i=12
        do f1=1,3 ! JE CORR all families
           kp3=knu(f1)
           do f2=1,3
              i=i+1
              kp4=kl(f2)
              kp3in(i)=kp3
              kp4in(i)=kp4
           enddo
        enddo

c...Need to add mixed final states snu sl* -> nu' l-bar         
        i=33
        do f1=1,3 ! JE CORR all families
           kp3=knu(f1)
           do f2=1,3
              kp4=kl(f2)
              if ((ivfam(kp3)+1).ne.ivfam(kp4)) then ! new mixed state
                 i=i+1
                 kp3in(i)=kp3
                 kp4in(i)=kp4
              endif
           enddo
        enddo

              
c...Old code below
c        do f=1,3
c          kp3=knue+2*(f-1)
c          kp4=ke+2*(f-1)
c          if(ivtype(kp3).eq.itype(1)) then
c            kp3in(i)=kp3
c            kp4in(i)=kp4
c            fsave=f 
c          endif
c        enddo

ccc
ccc gauge boson channels:
ccc
ccc
ccc t- and u-channels:
ccc
c        ksfertc(1)=ksl(fsave)
c        ksfertc(2)=ksl(fsave+3)
c        nsfertc=2
        ksfertc(1)=ksl(1) ! JE CORR all 6 below
        ksfertc(2)=ksl(2)
        ksfertc(3)=ksl(3)
        ksfertc(4)=ksl(4)
        ksfertc(5)=ksl(5)
        ksfertc(6)=ksl(6)
        nsfertc=6

c        ksferuc(1)=ksnu(fsave)
c        nsferuc=1
        ksferuc(1)=ksnu(1) ! JE CORR all 3 below
        ksferuc(2)=ksnu(2)
        ksferuc(3)=ksnu(3)
        nsferuc=3

c        ksfertn(1)=ksnu(fsave)
c        nsfertn=1
        ksfertn(1)=ksnu(1) ! JE CORR all 3 below
        ksfertn(2)=ksnu(2)
        ksfertn(3)=ksnu(3)
        nsfertn=3

c        ksferun(1)=ksl(fsave)
c        ksferun(2)=ksl(fsave+3)
c        nsferun=2
        ksferun(1)=ksl(1) ! JE CORR all 6 below
        ksferun(2)=ksl(2)
        ksferun(3)=ksl(3)
        ksferun(4)=ksl(4)
        ksferun(5)=ksl(5)
        ksferun(6)=ksl(6)
        nsferun=6
ccc
        gg1=1.d0 
        gg2=1.d0
        chcol=1
        cfactini=1.d0
        cfactfin=3.d0
        endif
ccc
        if(itype(1).ge.ivtype(ku).and.(iifam(2)-iifam(1)).eq.1) then
c============================================================
ccc input:   \tilde{u}_1,2(i) \tilde{d}^*_1,2(i)
c============================================================
ccc
ccc fermion channels: 
ccc    neutrino(j) lepton^+(j) -- same family j \in (1,2,3)
ccc    up-type quark(j) down-type antiquark(k) -- all allowed
ccc
ccc    up-type quark(i) down-type antiquark(i) -- same family 
ccc       as in the initial state

c...JE CORR. No need to set fsave here           
c        do f=1,3
c          kp3=ku+2*(f-1) ! up-type quark
c          if(ivfam(kp3).eq.iifam(1)) then
c            fsave=f
c            f2save=f
c          endif
c        enddo

ccc
ccc gauge boson channels:
ccc
        gluonin=.true.
ccc
ccc t- and u-channels:
ccc
c        ksfertc(1)=ksqd(fsave)
c        ksfertc(2)=ksqd(fsave+3)
c        nsfertc=2
        ksfertc(1)=ksqd(1) ! JE CORR all 6
        ksfertc(2)=ksqd(2)
        ksfertc(3)=ksqd(3)
        ksfertc(4)=ksqd(4)
        ksfertc(5)=ksqd(5)
        ksfertc(6)=ksqd(6)
        nsfertc=6

c        ksferuc(1)=ksqu(fsave)
c        ksferuc(2)=ksqu(fsave+3)
c        nsferuc=2
        ksferuc(1)=ksqu(1) ! JE CORR all 6
        ksferuc(2)=ksqu(2)
        ksferuc(3)=ksqu(3)
        ksferuc(4)=ksqu(4)
        ksferuc(5)=ksqu(5)
        ksferuc(6)=ksqu(6)
        nsferuc=6

c        ksfertn(1)=ksqu(fsave)
c        ksfertn(2)=ksqu(fsave+3)
c        nsfertn=2
        ksfertn(1)=ksqu(1) ! JE CORR all 6
        ksfertn(2)=ksqu(2)
        ksfertn(3)=ksqu(3)
        ksfertn(4)=ksqu(4)
        ksfertn(5)=ksqu(5)
        ksfertn(6)=ksqu(6)
        nsfertn=6

c        ksferun(1)=ksqd(fsave)
c        ksferun(2)=ksqd(fsave+3)
c        nsferun=2
        ksferun(1)=ksqd(1) ! JE CORR all 6
        ksferun(2)=ksqd(2)
        ksferun(3)=ksqd(3)
        ksferun(4)=ksqd(4)
        ksferun(5)=ksqd(5)
        ksferun(6)=ksqd(6)
        nsferun=6
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
        endif
ccc
        if(itype(1).ge.ivtype(ku).and.(iifam(2)-iifam(1)).ne.1) then
c============================================================
ccc input:   \tilde{u}(i) \tilde{d}^*_1,2(j)  i.ne.j
c============================================================
ccc
ccc fermion channels: 
ccc    neutrino(j) lepton^+(j) -- same family j \in (1,2,3)
ccc    up-type quark(j) down-type antiquark(k) -- all allowed
ccc
ccc    up-type quark(i) down-type antiquark(i) -- same family 
ccc       as in the initial state
c...JE CORR. This is not needed           
c        do f1=1,3
c          kp3=ku+2*(f1-1) ! up-type quark
c          if(ivfam(kp3).eq.iifam(1)) then
c            fsave=f1
c          endif
c        enddo
c        do f2=1,3
c          kp4=kd+2*(f2-1) ! down-type quark
c          if(ivfam(kp4).eq.iifam(2)) then
c            f2save=f2
c          endif
c        enddo

ccc
ccc gauge boson channels:
ccc
        gluonin=.true.
ccc
ccc t- and u-channels:
ccc
c        ksfertc(1)=ksqd(f2save)
c        ksfertc(2)=ksqd(f2save+3)
c        nsfertc=2
        ksfertc(1)=ksqd(1) ! JE CORR all 6
        ksfertc(2)=ksqd(2)
        ksfertc(3)=ksqd(3)
        ksfertc(4)=ksqd(4)
        ksfertc(5)=ksqd(5)
        ksfertc(6)=ksqd(6)
        nsfertc=6
        
c        ksferuc(1)=ksqu(fsave)
c        ksferuc(2)=ksqu(fsave+3)
c        nsferuc=2
        ksferuc(1)=ksqu(1) ! JE CORR all 6
        ksferuc(2)=ksqu(2)
        ksferuc(3)=ksqu(3)
        ksferuc(4)=ksqu(4)
        ksferuc(5)=ksqu(5)
        ksferuc(6)=ksqu(6)
        nsferuc=6
        
c        ksfertn(1)=ksqu(fsave)
c        ksfertn(2)=ksqu(fsave+3)
c        nsfertn=2
        ksfertn(1)=ksqu(1) ! JE CORR all 6
        ksfertn(2)=ksqu(2)
        ksfertn(3)=ksqu(3)
        ksfertn(4)=ksqu(4)
        ksfertn(5)=ksqu(5)
        ksfertn(6)=ksqu(6)
        nsfertn=6

c        ksferun(1)=ksqd(f2save)
c        ksferun(2)=ksqd(f2save+3)
c        nsferun=2
        ksferun(1)=ksqd(1) ! JE CORR all 6
        ksferun(2)=ksqd(2)
        ksferun(3)=ksqd(3)
        ksferun(4)=ksqd(4)
        ksferun(5)=ksqd(5)
        ksferun(6)=ksqd(6)
        nsferun=6
ccc
        gg1=3.d0
        gg2=3.d0
        chcol=2
        call dsascolset(chcol)
        cfactini=1.d0  ! set inside the routine dsasfercol
        cfactfin=1.d0
        endif
ccc
        endif
ccc

*************************************************** sfu sfd* -> fu fdbar
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state and in the final state; for sleptons flavour mixing for 
***** quarks in the final state:
*****
        if(chcol.eq.1) then ! sleptons in initial state
***** first the lepton final states:
          do f=1,3
            kp3=knue+2*(f-1)
            kp4=ke+2*(f-1)
            icase=2
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(f)=result
          enddo

c...Now the mixed ones
          do i=34,39
             kp3=kp3in(i)
             kp4=kp4in(i)
             icase=2
             call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
             prtial(i)=result
          enddo

***** then the quark final states:
          i=0
          do f1=1,3
          do f2=1,3
            kp3=ku+2*(f1-1) ! up-type quark
            kp4=kd+2*(f2-1) ! down-type quark
            icase=2
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            result=result*cfactfin
            i=i+1  
            prtial(3+i)=result
          enddo
          enddo
        else ! squarks in initial state
***** first the lepton final states:
          do f=1,3
            kp3=knue+2*(f-1)
            kp4=ke+2*(f-1)
            icase=2
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(f)=result
          enddo
***** then the quark final states:
          i=0
          do f1=1,3
          do f2=1,3
            kp3=ku+2*(f1-1) ! up-type quark
            kp4=kd+2*(f2-1) ! down-type quark
            icase=2
            call dsasfercol(kp1i,kp2i,kp3,kp4,icase,result)
            i=i+1  
            prtial(3+i)=result
          enddo
          enddo
        endif
******************************************************* sfu sfd -> fu fd
*****
***** for squark initial states, possible flavour mixing in the initial 
***** state and in the final state
*****
***** first the lepton final states:
        if(chcol.eq.1) then     ! slepton initial state
          do i=13,21
             kp3=kp3in(i)
             kp4=kp4in(i)
             icase=2
             call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)

             prtial(i)=result
          enddo
*****
***** then the quark final states:
        else ! squark initial state
        i=0
        do f1=1,3
        do f2=1,3
          kp3=ku+2*(f1-1) ! up-type quark
          kp4=kd+2*(f2-1) ! down-type quark
          icase=2
          call dsasferecol(kp1i,kp2i,kp3,kp4,icase,result)
          i=i+1  
          prtial(12+i)=result
        enddo
        enddo
        endif
******************************************************* sfu sfd* -> w+ z
        kp3=kw
        kp4=kz
        icase=1
        call dsasgbgb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(22)=par
*************************************************** sfu sfd* -> w+ gamma
        kp3=kw
        kp4=kgamma
        icase=1
        call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(23)=par
******************************************************* sfu sfd* -> w+ h
        kp3=kw
        kp4=kh2
        icase=3
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(24)=par
******************************************************* sfu sfd* -> w+ H
        kp3=kw
        kp4=kh1
        icase=3
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(25)=par
******************************************************* sfu sfd* -> w+ a
        kp3=kw
        kp4=kh3
        icase=4
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(26)=par
*************************************************** sfu sfd* -> h+ gamma
        kp3=kgamma
        kp4=khc
        icase=2
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(27)=par
******************************************************* sfu sfd* -> h+ z
        kp3=kz
        kp4=khc
        icase=1
        call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(28)=par
******************************************************* sfu sfd* -> h+ h
        kp3=khc
        kp4=kh2
        icase=1
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(29)=par
******************************************************* sfu sfd* -> h+ H
        kp3=khc
        kp4=kh1
        icase=1
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(30)=par
******************************************************* sfu sfd* -> h+ a
        kp3=khc
        kp4=kh3
        icase=2
        call dsashbhb(kp1i,kp2i,kp3,kp4,icase,par)
        prtial(31)=par
        if(gluonin) then
******************************************************* sfu sfd* -> W+ g
          kp3=kw
          kp4=kgluon
          icase=4
          call dsasgbgb1exp(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(32)=par
******************************************************* sfu sfd* -> g H+
          kp3=kgluon
          kp4=khc
          icase=11
          call dsasgbhb(kp1i,kp2i,kp3,kp4,icase,par)
          prtial(33)=par
        endif
**************************************************** sum partial results
        w=0.d0
        w=w+prtial(1)            ! sfu sfd* -> nue e+
        w=w+prtial(2)            ! sfu sfd* -> numu mu+
        w=w+prtial(3)            ! sfu sfd* -> nutau tau+
        w=w+prtial(4)            ! sfu sfd* -> u dbar
        w=w+prtial(5)            ! sfu sfd* -> u sbar
        w=w+prtial(6)            ! sfu sfd* -> u bbar
        w=w+prtial(7)            ! sfu sfd* -> c dbar
        w=w+prtial(8)            ! sfu sfd* -> c sbar
        w=w+prtial(9)            ! sfu sfd* -> c bbar
        w=w+prtial(10)           ! sfu sfd* -> t dbar
        w=w+prtial(11)           ! sfu sfd* -> t sbar
        w=w+prtial(12)           ! sfu sfd* -> t bbar
        w=w+prtial(13)           ! sfu sfd -> u d or nu_e e-
        w=w+prtial(14)           ! sfu sfd -> u s or nu_e mu-
        w=w+prtial(15)           ! sfu sfd -> u b or nu_e tau-
        w=w+prtial(16)           ! sfu sfd -> c d or nu_mu e-
        w=w+prtial(17)           ! sfu sfd -> c s or nu_mu mu-
        w=w+prtial(18)           ! sfu sfd -> c b or nu_mu tau-
        w=w+prtial(19)           ! sfu sfd -> t d or nu_tau e-
        w=w+prtial(20)           ! sfu sfd -> t s or nu_tau mu-
        w=w+prtial(21)           ! sfu sfd -> t b or nu_tau tau-
        w=w+prtial(22)           ! sfu sfd* -> w+ z
        w=w+prtial(23)           ! sfu sfd* -> w+ gamma
        w=w+prtial(24)           ! sfu sfd* -> w+ h
        w=w+prtial(25)           ! sfu sfd* -> w+ H
        w=w+prtial(26)           ! sfu sfd* -> w+ a
        w=w+prtial(27)           ! sfu sfd* -> h+ gamma
        w=w+prtial(28)           ! sfu sfd* -> h+ z
        w=w+prtial(29)           ! sfu sfd* -> h+ h
        w=w+prtial(30)           ! sfu sfd* -> h+ H
        w=w+prtial(31)           ! sfu sfd* -> h+ a
        w=w+prtial(32)           ! sfu sfd* -> W+ g
        w=w+prtial(33)           ! sfu sfd* -> g H+
        w=w+prtial(34)           ! sfu sfd* -> nu_e mu+
        w=w+prtial(35)           ! sfu sfd* -> nu_e tau+
        w=w+prtial(36)           ! sfu sfd* -> nu_mu e+
        w=w+prtial(37)           ! sfu sfd* -> nu_mu tau+
        w=w+prtial(38)           ! sfu sfd* -> nu_tau e+
        w=w+prtial(39)           ! sfu sfd* -> nu_tau mu+

c        do i=1,39               ! JE TMP
c           write(*,*) 'sfsf case 2: i=',i,prtial(i)
c        enddo
        
        
C
C Check for large negative terms in the fermion final states
c
        wok=.true.
        do i=1,21
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo

        do i=34,39
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo

c
c write error message:
c
        if(aszeroprint) then
        if (.not.wok) then
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: large negative term in dsasdwdcossfsf'
          write(*,*) 'DS: called with kp1i=',kp1i,' kp2i=',kp2i
          write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
          do i=1,39
            write(*,*) 'DS: i=',i,' prtial=',prtial(i)
          enddo
        endif
        endif
c
        dsasdwdcossfsf = w*0.5d0 ! 0.5d0 <- we should return weff for 
              ! sf sf ann with sf combined part. and anti-particle state
c======================================================================
c======================================================================
************************************************************************
*****
***** 3) one squark + one slepton
***** JE CORR. Note. Before this case was used also for sleptons in
***** different families, but this is handled above instead.
***** This part is almost entirely rewritten by J. Edsjo
c======================================================================
c======================================================================
      else
c...The final state particles can be different for both the first
c...and second outgoing particles.
c...JE CORR. Changed 9 final states to 27.         
c        write(*,*) 'Now in case 3...'
        do i=1,27
          prtial(i)=0.0d0
        enddo
ccc
ccc check whether to reload final/intermediate state arrays/variables
ccc
        if(kp1i.ne.kp1s.or.kp2i.ne.kp2s) then
        kp1s=kp1i
        kp2s=kp2i
        do i=1,27
          chon(i)=.false.
        enddo
***** identify particles in initial and final state
c============================================================
c=== up squark - slepton
c============================================================        
***** the first only can be a up-type squark:
        if(itype(1).eq.ivtype(ku)) then
          i=0
          do f1=1,3
             kp3=ku+2*(f1-1) ! the different kinds of up quarks
             do f2=1,3
                i=i+1
                if (itype(2).eq.ivtype(ksnu1)) then ! sneutrino
                   kp4=knu(f2)
                else            ! charged slepton
                   kp4=kl(f2)
                endif
                kp3in(i)=kp3 ! squ sl* -> u l-bar
                kp4in(i)=kp4
                chon(i)=.true.
                kp3in(i+9)=kp3 ! squ sl -> u l
                kp4in(i+9)=kp4
                chon(i+9)=.true.
             enddo
          enddo

          i=18
          do f1=1,3
             kp3=kd+2*(f1-1)
             do f2=1,3
                i=i+1
                if (itype(2).eq.ivtype(ksnu1)) then ! sneutrino
                   kp4=kl(f2)
                else            ! charged slepton
                   kp4=knu(f2)
                endif
                kp3in(i)=kp3
                kp4in(i)=kp4
                chon(i)=.true.
             enddo
          enddo
          ick1=1
***** 1 squark and 1 slepton
          gg1=3.d0
          gg2=1.d0
          chcol=1
          cfactini=3.d0
          cfactfin=1.d0
c============================================================
c=== down squark - slepton
c============================================================        
***** ... or a down-type squark:
        elseif(itype(1).eq.ivtype(kd)) then
          i=0
          do f1=1,3
             kp3=kd+2*(f1-1)
             do f2=1,3
                i=i+1
                if (itype(2).eq.ivtype(ksnu1)) then ! sneutrino
                   kp4=knu(f2)
                else            ! charged slepton
                   kp4=kl(f2)
                endif
                kp3in(i)=kp3 ! sqd sl* -> d l-bar
                kp4in(i)=kp4
                chon(i)=.true.
                kp3in(i+9)=kp3 ! sqd sl -> d l
                kp4in(i+9)=kp4
                chon(i+9)=.true.
             enddo
          enddo

          i=18
          do f1=1,3             
             kp3=ku+2*(f1-1)
             do f2=1,3
                i=i+1
                if (itype(2).eq.ivtype(ksnu1)) then ! sneutrino
                   kp4=kl(f2)
                else            ! charged slepton
                   kp4=knu(f2)
                endif
                kp3in(i)=kp3
                kp4in(i)=kp4
                chon(i)=.true.
             enddo
          enddo
          ick1=2
***** 1 squark and 1 slepton
          gg1=3.d0
          gg2=1.d0
          chcol=1
          cfactini=3.d0
          cfactfin=1.d0
        endif
       
***** identify the second slepton.
        if (itype(2).eq.ivtype(knue)) then
           ick2=1
        else
           ick2=2
        endif
ccc
        endif
*****
*************************************************** sf1 sf2* -> f1 f2bar
        do i=1,9
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=3
            call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result*cfactini
          endif
        enddo  
******************************************************* sf1 sf2 -> f1 f2
        do i=10,18
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            icase=3
            call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
            prtial(i)=result*cfactini
          endif
        enddo  
*************************************************** sf1 sf2* -> f3 f4bar
**************************************************** or sf1 sf2 -> f3 f4
*****  f1 and f3 switched up/down, f2 and f4 switched up/down
        do i=19,27
          if(chon(i)) then
            kp3=kp3in(i)
            kp4=kp4in(i)
            if(ick1.eq.ick2) then
              icase=3
              call dsasfer(kp1i,kp2i,kp3,kp4,icase,result)
            else
              icase=3
              call dsasfere(kp1i,kp2i,kp3,kp4,icase,result)
            endif
            prtial(i)=result*cfactini
          endif  
        enddo
**************************************************** sum partial results
        w=0.d0 ! JE CORR, more states below (9->27)
        w=w+prtial(1)           ! sq1 sl2* -> q1_1 lbar_1
        w=w+prtial(2)           ! sq1 sl2* -> q1_1 lbar_2
        w=w+prtial(3)           ! sq1 sl2* -> q1_1 lbar_3
        w=w+prtial(4)           ! sq1 sl2* -> q1_2 lbar_1
        w=w+prtial(5)           ! sq1 sl2* -> q1_2 lbar_2
        w=w+prtial(6)           ! sq1 sl2* -> q1_2 lbar_3
        w=w+prtial(7)           ! sq1 sl2* -> q1_3 lbar_1
        w=w+prtial(8)           ! sq1 sl2* -> q1_3 lbar_2
        w=w+prtial(9)           ! sq1 sl2* -> q1_3 lbar_3
        w=w+prtial(10)          ! sq1 sl2 -> q1_1 l_1
        w=w+prtial(11)          ! sq1 sl2 -> q1_1 l_2
        w=w+prtial(12)          ! sq1 sl2 -> q1_1 l_3
        w=w+prtial(13)          ! sq1 sl2 -> q1_2 l_1
        w=w+prtial(14)          ! sq1 sl2 -> q1_2 l_2
        w=w+prtial(15)          ! sq1 sl2 -> q1_2 l_3
        w=w+prtial(16)          ! sq1 sl2 -> q1_3 l_1
        w=w+prtial(17)          ! sq1 sl2 -> q1_3 l_2
        w=w+prtial(18)          ! sq1 sl2 -> q1_3 l_3
c...States with final states switched w.r.t initial state regarding up/down
        w=w+prtial(19)          ! sq1 sl2* -> q1'_1 l'_1
        w=w+prtial(20)          ! sq1 sl2* -> q1'_1 l'_2
        w=w+prtial(21)          ! sq1 sl2* -> q1'_1 l'_3
        w=w+prtial(22)          ! sq1 sl2* -> q1'_2 l'_1
        w=w+prtial(23)          ! sq1 sl2* -> q1'_2 l'_2
        w=w+prtial(24)          ! sq1 sl2* -> q1'_2 l'_3
        w=w+prtial(25)          ! sq1 sl2* -> q1'_3 l'_1
        w=w+prtial(26)          ! sq1 sl2* -> q1'_3 l'_2
        w=w+prtial(27)          ! sq1 sl2* -> q1'_3 l'_3

c
c check for large negative terms in the fermion final states
c
        wok=.true.
        do i=1,27
          if (prtial(i).lt.0.0d0) then
            if (dabs(prtial(i)/w).gt.fertoll) then
              wok=.false.
            else
              w=w-prtial(i)
              prtial(i)=0.0d0
            endif  
          endif  
        enddo
c
c write error message:
c
        if(aszeroprint) then
        if (.not.wok) then
          write(*,*) 'DS: Model: ',idtag
          write(*,*) 'DS: large negative term in dsasdwdcossfsf'
          write(*,*) 'DS: called with kp1i=',kp1i,' kp2i=',kp2i
          write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
          do i=1,12
            write(*,*) 'DS: i=',i,' prtial=',prtial(i)
          enddo
        endif
        endif
c
        dsasdwdcossfsf = w*0.5d0 ! 0.5d0 <- we should return weff for 
              ! sf sf ann with sf combined part. and anti-particle state
************************************************************************
***** if statement corresponding to the 3 cases closed
      endif

      if (w.lt.0.0d0) then
        write(*,*) 'DS: Model: ',idtag
        write(*,*) 'DS: dsasdwdcossfsf called with kp1i=',kp1i,
     &    ' kp2i=',kp2i
        write(*,*) 'DS: p=',p,' costh=',costhe,' w=',w
        write(*,*) 'DS: negative w, program stopped!'
        stop
      endif  

      return
      end








