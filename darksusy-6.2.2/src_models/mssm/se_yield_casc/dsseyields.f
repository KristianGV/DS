*****************************************************************************
*** function dsseyields calculates the yield above threshold (kind=1) or the
*** differential yield (kind=2) from a given scalar
*** boson decaying in flight, the energy of the scalar boson should be given
*** in eh.
*** scalars hno = 1-4 are supported (S10, S20, S30 and S+/S-)
*** units: 1.0e-30 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dsseyields(eh,emuth,thmax,hno,wh,kind,
     &  type,istat,seerror)
      implicit none
      include 'dsseyieldcom.h'
      include 'dsanyieldmodelcom.h'

c------------------------ variables ------------------------------------

      real*8 eh,emuth,thmax
      integer hno,istat,kind,type,chi,seerror
      character*2 wh

      real*8 e1,e2,yield

c------------------------ functions ------------------------------------

      real*8 dsseyields2,dsseemean,dsseyieldfth

c-----------------------------------------------------------------------

c...loop through different decay channels
      yield=0.0d0

c---------- neutral scalar bosons ----------
      if (hno.le.3) then

c..."fundamental" channels
        do 100 chi=1,14
          if (ans0br(chi2ch(chi),hno).gt.0.0d0) then
            yield=yield+ans0br(chi2ch(chi),hno)*2.0d0*
     &        dsseyieldfth(eh,ans0m(hno),msim(chi),msim(chi),emuth,thmax,
     &          chi,wh,kind,type,istat,seerror)
          endif
  100   continue

c..."complex" channels
        if (ans0br(1,hno).gt.0.0d0) then    ! S10 S10 channel
          e1=dsseemean(eh,ans0m(hno),ans0m(1),ans0m(1))
          yield=yield+2.0*ans0br(1,hno)*dsseyields2(e1,emuth,thmax,
     &      1,wh,kind,type,istat,seerror)
        endif

        if (ans0br(2,hno).gt.0.0d0) then    ! S10 S20 channel
          e1=dsseemean(eh,ans0m(hno),ans0m(1),ans0m(2))
          yield=yield+ans0br(2,hno)*dsseyields2(e1,emuth,thmax,
     &      1,wh,kind,type,istat,seerror)
          e2=dsseemean(eh,ans0m(hno),ans0m(2),ans0m(1))
          yield=yield+ans0br(2,hno)*dsseyields2(e2,emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
        endif

        if (ans0br(3,hno).gt.0.0d0) then    ! S20 S20 channel
          e1=dsseemean(eh,ans0m(hno),ans0m(2),ans0m(2))
          yield=yield+2.0*ans0br(3,hno)*dsseyields2(e1,emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
        endif

        if (ans0br(4,hno).gt.0.0d0) then    ! S30 S30 channel
          e1=dsseemean(eh,ans0m(hno),ans0m(3),ans0m(3))
          yield=yield+2.0*ans0br(4,hno)*dsseyields2(e1,emuth,thmax,
     &      3,wh,kind,type,istat,seerror)
        endif

        if (ans0br(5,hno).gt.0.0d0) then    ! S10 S30 channel
          e1=dsseemean(eh,ans0m(hno),ans0m(1),ans0m(3))
          yield=yield+ans0br(5,hno)*dsseyields2(e1,emuth,thmax,
     &      1,wh,kind,type,istat,seerror)
          e2=dsseemean(eh,ans0m(hno),ans0m(3),ans0m(1))
          yield=yield+ans0br(5,hno)*dsseyields2(e2,emuth,thmax,
     &      3,wh,kind,type,istat,seerror)
        endif

        if (ans0br(6,hno).gt.0.0d0) then    ! S20 S30 channel
          e1=dsseemean(eh,ans0m(hno),ans0m(2),ans0m(3))
          yield=yield+ans0br(6,hno)*dsseyields2(e1,emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
          e2=dsseemean(eh,ans0m(hno),ans0m(3),ans0m(2))
          yield=yield+ans0br(6,hno)*dsseyields2(e2,emuth,thmax,
     &      3,wh,kind,type,istat,seerror)
        endif

        if (ans0br(7,hno).gt.0.0d0) then    ! S+ S- channel
          e1=dsseemean(eh,ans0m(hno),anscm,anscm)
          yield=yield+2.0*
     &      ans0br(7,hno)*dsseyields2(e1,emuth,thmax,4,wh,
     &        kind,type,istat,seerror)
        endif

        if (ans0br(8,hno).gt.0.0d0) then    ! z0 S10 channel
          e2=dsseemean(eh,ans0m(hno),ans0m(1),msim(9))
          yield=yield+ans0br(8,hno)*
     +    dsseyieldfth(eh,ans0m(hno),msim(9),ans0m(1),emuth,thmax,
     &      9,wh,kind,type,istat,seerror)
          yield=yield+ans0br(8,hno)*dsseyields2(e2,emuth,thmax,
     &      1,wh,kind,type,istat,seerror)
        endif

        if (ans0br(9,hno).gt.0.0d0) then    ! z0 S20 channel
          e2=dsseemean(eh,ans0m(hno),ans0m(2),msim(9))
          yield=yield+ans0br(9,hno)*
     +    dsseyieldfth(eh,ans0m(hno),msim(9),ans0m(2),emuth,thmax,
     &      9,wh,kind,type,istat,seerror)
          yield=yield+ans0br(9,hno)*dsseyields2(e2,emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
        endif

        if (ans0br(10,hno).gt.0.0d0) then    ! z0 S30 channel
          e2=dsseemean(eh,ans0m(hno),ans0m(3),msim(9))
          yield=yield+ans0br(10,hno)*
     +    dsseyieldfth(eh,ans0m(hno),msim(9),ans0m(3),emuth,thmax,
     &      9,wh,kind,type,istat,seerror)
          yield=yield+ans0br(10,hno)*dsseyields2(e2,emuth,thmax,
     &      3,wh,kind,type,istat,seerror)
        endif

        if (ans0br(11,hno).gt.0.0d0) then    ! w-S+ w+S- channel
          e2=dsseemean(eh,ans0m(hno),anscm,msim(8))
          yield=yield+ans0br(11,hno)*
     +    dsseyieldfth(eh,ans0m(hno),msim(8),anscm,emuth,thmax,
     &      8,wh,kind,type,istat,seerror)
          yield=yield+ans0br(11,hno)*dsseyields2(e2,emuth,thmax,
     &      4,wh,kind,type,istat,seerror)
        endif

        if (ans0br(29,hno).gt.0.0d0) then    ! Z0 gamma
          yield=yield+ans0br(29,hno)*
     +    dsseyieldfth(eh,ans0m(hno),msim(9),0.0d0,emuth,thmax,
     &      9,wh,kind,type,istat,seerror)
        endif
c        write(*,*) '   ...after complex: ',yield

c---------- charged scalar bosons ----------
      else

        if (anscbr(1).gt.0.0d0) then ! u d-bar
           yield=yield+anscbr(1)*
     &    dsseyieldfth(eh,anscm,msim(1),msim(2),emuth,thmax,
     &      1,wh,kind,type,istat,seerror)
          yield=yield+anscbr(1)*
     &    dsseyieldfth(eh,anscm,msim(2),msim(1),emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
        endif

        if (anscbr(2).gt.0.0d0) then ! u s-bar
           yield=yield+anscbr(2)*
     &    dsseyieldfth(eh,anscm,msim(3),msim(2),emuth,thmax,
     &      3,wh,kind,type,istat,seerror)
          yield=yield+anscbr(2)*
     &    dsseyieldfth(eh,anscm,msim(2),msim(3),emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
        endif

        if (anscbr(3).gt.0.0d0) then ! u b-bar
          yield=yield+anscbr(3)*
     &    dsseyieldfth(eh,anscm,msim(5),msim(2),emuth,thmax,
     &      5,wh,kind,type,istat,seerror)
          yield=yield+anscbr(3)*
     &    dsseyieldfth(eh,anscm,msim(2),msim(5),emuth,thmax,
     &      2,wh,kind,type,istat,seerror)
        endif

        if (anscbr(4).gt.0.0d0) then ! c d-bar
          yield=yield+anscbr(4)*
     &    dsseyieldfth(eh,anscm,msim(4),msim(1),emuth,thmax,
     &      4,wh,kind,type,istat,seerror)
          yield=yield+anscbr(4)*
     &    dsseyieldfth(eh,anscm,msim(1),msim(4),emuth,thmax,
     &      1,wh,kind,type,istat,seerror)
        endif

        if (anscbr(5).gt.0.0d0) then ! c s-bar
          yield=yield+anscbr(5)*
     &    dsseyieldfth(eh,anscm,msim(4),msim(3),emuth,thmax,
     &    4,wh,kind,type,istat,seerror)
          yield=yield+anscbr(5)*
     &    dsseyieldfth(eh,anscm,msim(3),msim(4),emuth,thmax,
     &    3,wh,kind,type,istat,seerror)
        endif

        if (anscbr(6).gt.0.0d0) then ! c b-bar
          yield=yield+anscbr(6)*
     &    dsseyieldfth(eh,anscm,msim(4),msim(5),emuth,thmax,
     &    4,wh,kind,type,istat,seerror)
          yield=yield+anscbr(6)*
     &    dsseyieldfth(eh,anscm,msim(5),msim(4),emuth,thmax,
     &    5,wh,kind,type,istat,seerror)
        endif

        if (anscbr(7).gt.0.0d0) then ! t d-bar
          yield=yield+anscbr(7)*
     &    dsseyieldfth(eh,anscm,msim(6),msim(1),emuth,thmax,
     &    6,wh,kind,type,istat,seerror)
          yield=yield+anscbr(7)*
     &    dsseyieldfth(eh,anscm,msim(1),msim(6),emuth,thmax,
     &    1,wh,kind,type,istat,seerror)
        endif

        if (anscbr(8).gt.0.0d0) then ! t s-bar
          yield=yield+anscbr(8)*
     &    dsseyieldfth(eh,anscm,msim(6),msim(3),emuth,thmax,
     &    6,wh,kind,type,istat,seerror)
          yield=yield+anscbr(8)*
     &    dsseyieldfth(eh,anscm,msim(3),msim(6),emuth,thmax,
     &    3,wh,kind,type,istat,seerror)
        endif

        if (anscbr(9).gt.0.0d0) then ! t b-bar
          yield=yield+anscbr(9)*
     &    dsseyieldfth(eh,anscm,msim(6),msim(5),emuth,thmax,
     &    6,wh,kind,type,istat,seerror)
          yield=yield+anscbr(9)*
     &    dsseyieldfth(eh,anscm,msim(5),msim(6),emuth,thmax,
     &    5,wh,kind,type,istat,seerror)
        endif

        if (anscbr(12).gt.0.0d0) then ! tau nu_tau
          yield=yield+anscbr(12)*
     &    dsseyieldfth(eh,anscm,msim(11),msim(14),emuth,thmax,
     &    11,wh,kind,type,istat,seerror)
          yield=yield+anscbr(12)*
     &    dsseyieldfth(eh,anscm,msim(14),msim(11),emuth,thmax,
     &    14,wh,kind,type,istat,seerror)
        endif

        if (anscbr(13).gt.0.0d0) then ! w+ h1
          yield=yield+anscbr(13)*
     &    dsseyieldfth(eh,anscm,msim(8),ans0m(1),emuth,thmax,
     &    8,wh,kind,type,istat,seerror)
          e2=dsseemean(eh,anscm,ans0m(1),msim(8))
          yield=yield+anscbr(13)*
     &      dsseyields2(e2,emuth,thmax,1,wh,kind,type,istat,seerror)
        endif

        if (anscbr(14).gt.0.0d0) then ! w+ h2
          yield=yield+anscbr(14)*
     &    dsseyieldfth(eh,anscm,msim(8),ans0m(2),emuth,thmax,
     &    8,wh,kind,type,istat,seerror)
          e2=dsseemean(eh,anscm,ans0m(2),msim(8))
          yield=yield+anscbr(14)*
     &      dsseyields2(e2,emuth,thmax,2,wh,kind,type,istat,seerror)
        endif

        if (anscbr(15).gt.0.0d0) then ! w+ h3
          yield=yield+anscbr(15)*
     &    dsseyieldfth(eh,anscm,msim(8),ans0m(3),emuth,thmax,
     &    8,wh,kind,type,istat,seerror)
          e2=dsseemean(eh,anscm,ans0m(3),msim(8))
          yield=yield+anscbr(15)*
     &      dsseyields2(e2,emuth,thmax,3,wh,kind,type,istat,seerror)
        endif

      endif

      dsseyields = yield

      end




