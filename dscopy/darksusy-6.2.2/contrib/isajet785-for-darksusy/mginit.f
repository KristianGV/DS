CDECK  ID>, MGINIT.
      SUBROUTINE MGINIT
C
C          Initialize common blocks for MadGraph code in ISAJET
C          Note the QCD coupling constant is g=1.
C
      IMPLICIT NONE
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
C=====     Begin common blocks used by MadGraph
      REAL*8            GW, GWWA, GWWZ
      COMMON /COUP1/    GW, GWWA, GWWZ
      SAVE /COUP1/
      REAL*8            GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/   GAL,   GAU,   GAD,   GWF
      SAVE /COUP2A/
      REAL*8            GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/   GZN,   GZL,   GZU,   GZD,   G1
      SAVE /COUP2B/
      REAL*8            GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/    GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      SAVE /COUP3/
      COMPLEX*16        GCHF(2,12)
      COMMON /COUP4/    GCHF
      SAVE /COUP4/
      REAL*8            WMASS,WWIDTH,ZMASS,ZWIDTH
      COMMON /VMASS1/   WMASS,WWIDTH,ZMASS,ZWIDTH
      SAVE /VMASS1/
      REAL*8            AMASS,AWIDTH,HMASS,HWIDTH
      COMMON /VMASS2/   AMASS,AWIDTH,HMASS,HWIDTH
      SAVE /VMASS2/
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
      SAVE /FERMIONS/
      REAL*8            GG(2), G
      COMMON /COUPQCD/  GG,    G
      SAVE /COUPQCD/
C=====     End common blocks used by MadGraph
C
      INTEGER I
      REAL AMGMW
      REAL*8 SW2
C
C          Fermion masses and widths
C
      FMASS(1) = AMGMW(IDE,1)
      FMASS(2) = 0D0
      FMASS(3) = AMGMW(IDUP,1)
      FMASS(4) = AMGMW(IDDN,1)
      FMASS(5) = AMGMW(IDMU,1)
      FMASS(6) = 0D0
      FMASS(7) = AMGMW(IDCH,1)
      FMASS(8) = AMGMW(IDST,1)
      FMASS(9) = AMGMW(IDTAU,1)
      FMASS(10)= 0D0
      FMASS(11)= AMGMW(IDTP,1)
      FMASS(12)= AMGMW(IDBT,1)
      DO 100 I=1,12
        FWIDTH(I)=0D0
100   CONTINUE
C
C          Boson masses and widths
C
      AMASS=0D0
      AWIDTH=0D0
      WMASS=AMGMW(IDW,1)
      WWIDTH=AMGMW(IDW,2)
      ZMASS=AMGMW(IDZ,1)
      ZWIDTH=AMGMW(IDZ,2)
      HMASS=AMGMW(IDH,1)
      HWIDTH=AMGMW(IDH,2)
      SW2=AMGMW(1,3)
C
C          Calls to Helas routines to set couplings
C
      CALL COUP1X(SW2,GW,GWWA,GWWZ)
      CALL COUP2X(SW2,GAL,GAU,GAD,GWF,GZN,GZL,GZU,GZD,G1)
      CALL COUP3X(SW2,ZMASS,HMASS,GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH)
      DO 110 I=1,12
         CALL COUP4X(SW2,ZMASS,FMASS(I),GCHF(1,I))
110   CONTINUE
C
C          QCD couplings
C
      G = 1D0
      GG(1)=-G
      GG(2)=-G
      RETURN
      END
