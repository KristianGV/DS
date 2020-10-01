      real*8 function dssea_ff3(x)
      implicit none
      real*8 dssea_ff2,x
      dssea_ff3=dssea_ff2(exp(x))*exp(x)*1.0d15
      return
      end
