      function dsrdfuncs(wrate,u)
c_______________________________________________________________________
c  10^15 * dsrdfunc.
c  input:
c    u - integration variable
c  uses dsrdfunc
c  used for gaussian integration with gadap.f
c  author: joakim edsjo (edsjo@fysik.su.se)
c  date: 97-01-17
c  modified 2018-04-26 torsten bringmann: external wrate as input 
c                                         (instead of hardcoded dsrdwintp)
c=======================================================================
      implicit none
      real*8 dsrdfuncs,wrate,u
      real*8 dsrdfunc,dsrdwintp
      external dsrdfunc,dsrdwintp,wrate

c-----------------------------------------------------------------------
      real*8 rdx
      common /gadint2/ rdx
c-----------------------------------------------------------------------

c      dsrdfuncs=dsrdfunc(u,rdx,dsrdwintp)*1.0d15
      dsrdfuncs=dsrdfunc(u,rdx,wrate)*1.0d15

      end








































