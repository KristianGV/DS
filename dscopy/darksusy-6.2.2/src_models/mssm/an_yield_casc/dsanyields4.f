*****************************************************************************
*** function dsanyields4 calculates the yield above threshold (yieldk=1) or the
*** differential yield (yieldk=2) from a given scalar
*** boson decaying in flight, the energy of the scalar boson should be given
*** in eh.
*** scalars hno = 1-4 are supported (S10, S20, S30 and S+/S-)
*** units: 1.0e-30 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dsanyields4(eh,egev,hno,yieldk,istat)
      implicit none
      include 'dsanyieldcom.h'
      include 'dsanyieldmodelcom.h'
      include 'dsidtag.h'

c------------------------ variables ------------------------------------

      real*8 eh,egev
      integer hno,istat,yieldk
      integer dch

      real*8 yield

c------------------------ functions ------------------------------------

      real*8 dsanyieldfth,dsanyielddec

c-----------------------------------------------------------------------


c...loop through different decay channels
      yield=0.0d0

c---------- neutral scalar bosons ----------
      if (hno.le.3) then

c..."fundamental" channels
         yield=yield+dsanyielddec(eh,hno,egev,yieldk,istat)

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

      endif

      dsanyields4 = yield

c...check for inconsistent declarations in ans0br and anscbr
      if (hno.le.3) then  ! neutral scalars
        do dch=1,11
          if (ans0br(dch,hno).gt.0.0d0) then
          write(6,*) 'error in dsanyields4: inconsistent decay widths',
     +      ' declared in ans0br'
          write(6,*) 'model: ',idtag
          endif
        enddo
      else                ! charged scalar
        do dch=13,15
          if (anscbr(dch).gt.0.0d0) then
          write(6,*) 'error in dsanyields4: inconsistent decay widths',
     +      ' declared in h0scbr'
          write(6,*) 'model: ',idtag
          endif
        enddo
      endif

      end




