*****************************************************************************
*** function dsanyields3 calculates the yield above threshold (yieldk=1) or the
*** differntial yield (yieldk=2) from a given scalar
*** boson decaying in flight, the energy of the scalar boson should be given
*** in eh.
*** scalars hno = 1-4 are supported (S10, S20, S30 and S+/S-)
*** units: 1.0e-30 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dsanyields3(eh,egev,hno,yieldk,istat)
      implicit none
      include 'dsanyieldcom.h'
      include 'dsanyieldmodelcom.h'

c------------------------ variables ------------------------------------

      real*8 eh,egev
      integer hno,istat,yieldk

      real*8 e1,e2,yield

c------------------------ functions ------------------------------------

      real*8 dsanyields4,dsanemean,dsanyieldfth,dsanyielddec

c-----------------------------------------------------------------------


c...loop through different decay channels
      yield=0.0d0

c---------- neutral scalar bosons ----------
      if (hno.le.3) then

c..."fundamental" channels
        yield=yield+dsanyielddec(eh,hno,egev,yieldk,istat)

c..."complex" channels
        if (ans0br(1,hno).gt.0.0d0) then    ! S10 S10 channel
          e1=dsanemean(eh,ans0m(hno),ans0m(1),ans0m(1))
          yield=yield+2.0*ans0br(1,hno)*
     &      dsanyields4(e1,egev,1,yieldk,istat)
        endif

        if (ans0br(3,hno).gt.0.0d0) then    ! S20 S20 channel
          e1=dsanemean(eh,ans0m(hno),ans0m(2),ans0m(2))
          yield=yield+2.0*ans0br(3,hno)*
     &      dsanyields4(e1,egev,2,yieldk,istat)
        endif

        if (ans0br(4,hno).gt.0.0d0) then    ! S30 S30 channel
          e1=dsanemean(eh,ans0m(hno),ans0m(3),ans0m(3))
          yield=yield+2.0*ans0br(4,hno)*
     &      dsanyields4(e1,egev,3,yieldk,istat)
        endif

        if (ans0br(7,hno).gt.0.0d0) then    ! S+ S- channel
          e1=dsanemean(eh,ans0m(hno),anscm,anscm)
          yield=yield+2.0*
     &      ans0br(7,hno)*dsanyields4(e1,egev,4,yieldk,istat)
        endif

        if (ans0br(8,hno).gt.0.0d0) then    ! z0 S10 channel
          e2=dsanemean(eh,ans0m(hno),ans0m(1),msim(9))
          yield=yield+ans0br(8,hno)*
     +    dsanyieldfth(eh,ans0m(hno),msim(9),
     &      ans0m(1),egev,9,yieldk,istat)
          yield=yield+ans0br(8,hno)*dsanyields4(e2,egev,1,yieldk,istat)
        endif

        if (ans0br(9,hno).gt.0.0d0) then    ! z0 S20 channel
          e2=dsanemean(eh,ans0m(hno),ans0m(2),msim(9))
          yield=yield+ans0br(9,hno)*
     +    dsanyieldfth(eh,ans0m(hno),msim(9),
     &       ans0m(2),egev,9,yieldk,istat)
          yield=yield+ans0br(9,hno)*dsanyields4(e2,egev,2,yieldk,istat)
        endif

        if (ans0br(10,hno).gt.0.0d0) then    ! z0 S30 channel
          e2=dsanemean(eh,ans0m(hno),ans0m(3),msim(9))
          yield=yield+ans0br(10,hno)*
     +    dsanyieldfth(eh,ans0m(hno),msim(9),
     &       ans0m(3),egev,9,yieldk,istat)
          yield=yield
     &      +ans0br(10,hno)*dsanyields4(e2,egev,3,yieldk,istat)
        endif

        if (ans0br(11,hno).gt.0.0d0) then    ! w-S+ w+S- channel
          e2=dsanemean(eh,ans0m(hno),anscm,msim(8))
          yield=yield+ans0br(11,hno)*
     +    dsanyieldfth(eh,ans0m(hno),msim(8),anscm,egev,8,yieldk,istat)
          yield=yield
     &      +ans0br(11,hno)*dsanyields4(e2,egev,4,yieldk,istat)
        endif

c---------- charged scalar bosons ----------
      else

        if (anscbr(3).gt.0.0d0) then ! u b-bar
          yield=yield+anscbr(3)*
     &    dsanyieldfth(eh,anscm,msim(5),msim(2),egev,5,yieldk,istat)
          yield=yield+anscbr(3)*
     &    dsanyieldfth(eh,anscm,msim(2),msim(5),egev,2,yieldk,istat)
        endif

        if (anscbr(4).gt.0.0d0) then ! c d-bar
          yield=yield+anscbr(4)*
     &    dsanyieldfth(eh,anscm,msim(4),msim(1),egev,4,yieldk,istat)
          yield=yield+anscbr(4)*
     &    dsanyieldfth(eh,anscm,msim(1),msim(4),egev,1,yieldk,istat)
        endif

        if (anscbr(5).gt.0.0d0) then ! c s-bar
          yield=yield+anscbr(5)*
     &    dsanyieldfth(eh,anscm,msim(4),msim(3),egev,4,yieldk,istat)
          yield=yield+anscbr(5)*
     &    dsanyieldfth(eh,anscm,msim(3),msim(4),egev,3,yieldk,istat)
        endif

        if (anscbr(6).gt.0.0d0) then ! c b-bar
          yield=yield+anscbr(6)*
     &    dsanyieldfth(eh,anscm,msim(4),msim(5),egev,
     &    4,yieldk,istat)
          yield=yield+anscbr(6)*
     &    dsanyieldfth(eh,anscm,msim(5),msim(4),egev,
     &    5,yieldk,istat)
        endif

        if (anscbr(7).gt.0.0d0) then ! t d-bar
          yield=yield+anscbr(7)*
     &    dsanyieldfth(eh,anscm,msim(6),msim(1),egev,
     &    6,yieldk,istat)
          yield=yield+anscbr(7)*
     &    dsanyieldfth(eh,anscm,msim(1),msim(6),egev,
     &    1,yieldk,istat)
        endif

        if (anscbr(8).gt.0.0d0) then ! t s-bar
          yield=yield+anscbr(8)*
     &    dsanyieldfth(eh,anscm,msim(6),msim(3),egev,
     &    6,yieldk,istat)
          yield=yield+anscbr(8)*
     &    dsanyieldfth(eh,anscm,msim(3),msim(6),egev,
     &    3,yieldk,istat)
        endif

        if (anscbr(9).gt.0.0d0) then ! t b-bar
          yield=yield+anscbr(9)*
     &    dsanyieldfth(eh,anscm,msim(6),msim(5),egev,
     &    6,yieldk,istat)
          yield=yield+anscbr(9)*
     &    dsanyieldfth(eh,anscm,msim(5),msim(6),egev,
     &    5,yieldk,istat)
        endif

        if (anscbr(12).gt.0.0d0) then ! tau nu_tau
          yield=yield+anscbr(12)*
     &    dsanyieldfth(eh,anscm,msim(11),0.0d0,egev,
     &    11,yieldk,istat)
        endif

        if (anscbr(13).gt.0.0d0) then ! w+ h1
          yield=yield+anscbr(13)*
     &    dsanyieldfth(eh,anscm,msim(8),ans0m(1),egev,
     &    8,yieldk,istat)
          e2=dsanemean(eh,anscm,ans0m(1),msim(8))
          yield=yield+anscbr(13)*
     &      dsanyields4(e2,egev,1,yieldk,istat)
        endif

        if (anscbr(14).gt.0.0d0) then ! w+ h2
          yield=yield+anscbr(14)*
     &    dsanyieldfth(eh,anscm,msim(8),ans0m(2),egev,
     &    8,yieldk,istat)
          e2=dsanemean(eh,anscm,ans0m(2),msim(8))
          yield=yield+anscbr(14)*
     &      dsanyields4(e2,egev,2,yieldk,istat)
        endif

        if (anscbr(15).gt.0.0d0) then ! w+ h3
          yield=yield+anscbr(15)*
     &    dsanyieldfth(eh,anscm,msim(8),ans0m(3),egev,
     &    8,yieldk,istat)
          e2=dsanemean(eh,anscm,ans0m(3),msim(8))
          yield=yield+anscbr(15)*
     &      dsanyields4(e2,egev,3,yieldk,istat)
        endif

      endif

      dsanyields3 = yield

      end


