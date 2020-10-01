c      block data  dsrdcom
      subroutine  dsrdcom
      implicit none
      include 'dsrdcom.h'

      data cosmin,waccd,dpminr,dpthr,wdiffr,wdifft,
     &  hstep,hmin,compeps,xinit,xfinal,umax,cfr
     &  /0.996195d0,0.005d0,1d-4,5d-4,0.05d0,0.02d0,
     &  0.001d0,1.0d-9,0.001d0,2.0d0,200.0d0,10.0d0,0.5d0/
c...0.996195d0 - 5 degrees, 0.999048d0 - 2.5 degrees
      data thavint/1/  ! dgadap
      data rdluerr,rdlulog/6,6/
      data rdprt/0/
      data pdivr,dpres/2.0d0,0.5d0/
      data nlow,nhigh,npres,nthup,cthtest,spltest
     &     /20,10,4,4,0,1/
c...rdt_max is the maximum time (in seconds) to spend on the relic density
c...calculation. It is checked against before a new point is added to the
c...W_eff tabulation, so in reality it can spend a bit longer than the limit.
c...The time is evaluated with the function CPU_TIME meaning that it is the
c...total CPU_TIME (not wall time) that is used. E.g. if four cores are used
c...the CPU time is four times the wall time.      
      data rdt_max/1.0d30/      ! max time for RD calculation (seconds)
      data rdinit/0/ ! flag to keep track of if initialization is done

      return
      end

