*****************************************************************************
*** function dsseyields4 calculates the yield above threshold (kind=1) or the
*** differential yield (kind=2) from a given scalar
*** boson decaying in flight, the energy of the scalar boson should be given
*** in eh.
*** scalars hno = 1-4 are supported (S10, S20, S30 and S+/S-)
*** units: 1.0e-30 m**-2 (annihilation)**-1
*****************************************************************************

      real*8 function dsseyields4(eh,emuth,thmax,hno,wh,
     &  kind,type,istat,seerror)
      implicit none
      include 'dsseyieldcom.h'
      include 'dsanyieldmodelcom.h'
      include 'dsidtag.h'

c------------------------ variables ------------------------------------

      real*8 eh,emuth,thmax
      integer hno,istat,kind,type,seerror
      character*2 wh

      real*8 yield
      integer chi,ch

c------------------------ functions ------------------------------------

      real*8 dsseyieldfth

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

      endif

      dsseyields4 = yield

c...check for inconsistent declarations in hdw
      if (hno.le.3) then ! neutral scalars
        do ch=1,11
          if (ans0br(ch,hno).gt.0.0d0) then
          write(6,*) 
     &      'DS ERROR in dsseyields4: inconsistent decay widths',
     +      ' declared in ans0br.'
          write(6,*) 'model: ',idtag
          endif
        enddo
      else               ! charged scalar
        do ch=13,15
          if (anscbr(ch).gt.0.0d0) then
          write(6,*) 'error in dsseyields4: inconsistent decay widths',
     +      ' declared in anscbr.'
          write(6,*) 'model: ',idtag
          endif
        enddo
      endif

      end




