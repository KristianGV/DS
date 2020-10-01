      real*8 function dsmqpole1loop(mqmq)
      implicit none
      include 'dsmpconst.h'
      real*8 mqmq,dsralph31loop
c
      dsmqpole1loop=mqmq*(1.d0+4.d0/3.d0/pi*dsralph31loop(mqmq))
      return
      end
