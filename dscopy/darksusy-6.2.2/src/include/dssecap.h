*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dssecap.h                             ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2003-11-27


c...common block with available tables for capture rates (w/o num FF int)
      integer nc,ntea,ntsu
      parameter(nc=2500,  ! number of entries in the tables
     &          ntea=6,   ! number of earth and sun tables that can be
     &          ntsu=6)   ! loaded simultaneously
      real*8 ctabea(0:nc,ntea),ctabsusi(0:nc,ntsu),ctabsusd(0:nc,ntsu)
      common /dssecap/ctabea,ctabsusi,ctabsusd
      save /dssecap/

c...common block with information about what is loaded and with what
c...(w/o num FF integration)      
      integer nealoaded,nsuloaded
      character*200 fileea(ntea),filesu(ntsu)
      common /dssecap2/fileea,filesu,nealoaded,nsuloaded
      save /dssecap2/

c...common block with available tables for capture rates
c...with numerical FF integration
      integer ncff,ntsuff
      parameter(ncff=1000,  ! number of entries in the tables
     &          ntsuff=6)   ! number of tables that can be loaded simultaneously
      real*8 ctabffsu(0:ncff,ntsuff,6) ! six combinations of couplings
      common /dssecapff/ctabffsu
      save /dssecapff/

c...common block with information about what is loaded and with what
c...(witn numerical FF integration)      
      integer nsuffloaded
      character*200 filesuff(ntsuff)
      common /dssecapff2/filesuff,nsuffloaded
      save /dssecapff2/

************************* end of dssecap.h *****************************






