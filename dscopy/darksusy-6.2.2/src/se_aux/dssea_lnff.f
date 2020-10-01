      real*8 function dssea_lnff(x)
      implicit none
      real*8 dssea_ff,x
      dssea_lnff=dssea_ff(exp(x))*exp(x)
      return
      end
