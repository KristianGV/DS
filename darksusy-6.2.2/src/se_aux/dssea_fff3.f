      real*8 function dssea_fff3(x)
      implicit none
      real*8 dssea_fff2,x
      dssea_fff3=dssea_fff2(exp(x))*exp(x)*1.0d15
      return
      end
