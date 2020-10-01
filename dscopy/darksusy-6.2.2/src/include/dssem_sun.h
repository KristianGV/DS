*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                          dssem_sun.h                             ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2003-11-24, 2005-11-24

c...JE FIX: this solution doesn't work, leads to multiply defined vars
c      include 'dssem_sun-const.h' ! def of zmax and isomax
      
      integer zmax,isomax
      parameter(zmax=92,isomax=10)
      integer sdmax,sdn,sdnne
      parameter(sdmax=3000)
      integer sdelmax
      parameter (sdelmax=500)

      real*8 sdr(sdmax),sdm(sdmax),sdrho(sdmax),sdphi(sdmax),
     &  sdmfr(zmax,0:isomax,sdmax),sdcdens(sdmax,0:2)
      real*8 sdrne(sdmax),sdne(sdmax)
      real*8 sdabund(zmax),sdaa(zmax,0:isomax),sdma(zmax,0:isomax),
     &  r_sun,m_sun,cd_sun(0:2),sdaberr(zmax),sdsp(zmax,0:isomax),
     &  sdfr(zmax,0:isomax)

      character*2 sdname(1:zmax)
      integer sdzsi(sdelmax),sdisosi(sdelmax),sdelsi
      integer sdzsd(sdelmax),sdisosd(sdelmax),sdelsd

      common /dssun/sdr,sdm,sdrho,sdphi,sdmfr,sdcdens,
     &  sdabund,sdaberr,sdaa,sdma,sdsp,sdfr,
     &  r_sun,m_sun,cd_sun,sdrne,sdne,sdn,sdnne,
     &  sdzsi,sdisosi,sdelsi,
     &  sdzsd,sdisosd,sdelsd,sdname
      save /dssun/

      character cdt
      common /dssuncd/cdt
      save /dssuncd/


      character*12 sunid
      character*200 sunfile,sunnefile,sunabundfile,sunisofile
      character*2 absrc
      integer sunfiletype
      logical sdread
      common /dssungen/sunfiletype,sdread,sunfile,sunnefile,
     &  sunabundfile,sunisofile,sunid,absrc
      save /dssungen/


c...stored capture information
c...JE FIX. Want to move out of here, but need zmax and isomax
      real*8 capsunsi(zmax,0:isomax),capsunsd(zmax,0:isomax)
      common /dscapsun/capsunsi,capsunsd
      save /dscapsun/

************************* end of dssem_sun.h ***************************






