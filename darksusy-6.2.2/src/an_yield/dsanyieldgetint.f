*****************************************************************************
*** function dsanyieldgetint returns the integrated yields for the given
*** choice of arguments. Compared to reading the arrays directly, this
*** routine interpolates in the coalescence momenta for dbars.
*** The routine is called by dsanyieldget.f that also smoothes the arrays.
*** The coalescence momentum used is determined by dbp0bar in common block.
*** If zero, the best-fit p0 dbp0fit will be used. If non-zero, it indicates
*** the deviation in units of the 1sigma uncertainty on p0.
*** For the actual p0 values, they are in MeV and tabulations exist
*** up to 300 MeV.
*** kind is used (currently) for dbar yields. 1=yield, 2=error on yield
*** Author: Joakim Edsjo, edsjo@fysik.su.se
*** Date: November, 2014
*****************************************************************************

      real*8 function dsanyieldgetint(zi,mxi,ch,fi,kind)
      implicit none
      include 'dsanyieldcom.h'

c------------------------ variables ------------------------------------

      integer zi,mxi,ch,fi,dpi,kind
      real*8 p0,dpl

c-----------------------------------------------------------------------

      if (fi.lt.9.or.fi.gt.13) then ! all yields but dbar

         dsanyieldgetint=phiint(zi,mxi,ch,fi)

      else ! dbar yields

         if (dbp0bar.gt.0.d0) then
            p0=dbp0fit(fi)+dbp0bar*(dbp0high(fi)-dbp0fit(fi))
         else
            p0=dbp0fit(fi)+dbp0bar*(dbp0fit(fi)-dbp0low(fi))
         endif

        call dsanifind(p0,dbdpindex(-1),dpl,dpi,-1,dpn-2)

        if (dpi.ge.-1) then
           dsanyieldgetint=(1.0-dpl)*phiintdb(zi,dpi,mxi,ch,fi,kind)
     &       +dpl*phiintdb(zi,dpi+1,mxi,ch,fi,kind)
        else
           dsanyieldgetint=0.0d0
        endif
      endif

      return
      end
