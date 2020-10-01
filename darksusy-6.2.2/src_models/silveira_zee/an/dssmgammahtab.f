      function dssmgammahtab(sqrts)
c... Tabulation of standard model total Higgs decay width in GeV
c... From Dittmaier et al 1101.0593
c... Author: Paolo Gondolo 2016
      implicit none
      real*8 dssmgammahtab,sqrts
      real*8 mh(52),wh(52)
      integer i
      data mh/
     &      90.d0,  95.d0, 100.d0, 105.d0, 110.d0, 115.d0, 120.d0, 125.d0, 130.d0, 135.d0,
     &     140.d0, 145.d0, 150.d0, 155.d0, 160.d0, 165.d0, 170.d0, 175.d0, 180.d0, 185.d0,
     &     190.d0, 195.d0, 200.d0, 210.d0, 220.d0, 230.d0, 240.d0, 250.d0, 260.d0, 270.d0,
     &     280.d0, 290.d0, 300.d0, 310.d0, 320.d0, 330.d0, 340.d0, 350.d0, 360.d0, 370.d0,
     &     380.d0, 390.d0, 400.d0, 410.d0, 420.d0, 430.d0, 440.d0, 450.d0, 460.d0, 470.d0,
     &     480.d0, 490.d0/
      data wh/
     &     2.20d-3, 2.32d-3, 2.46d-3, 2.62d-3, 2.82d-3, 3.09d-3, 3.47d-3, 4.03d-3, 4.87d-3, 6.14d-3,
     &     8.12d-3, 1.14d-2, 1.73d-2, 3.02d-2, 8.29d-2, 2.46d-1, 3.80d-1, 5.00d-1, 6.31d-1, 8.32d-1,
     &     1.04d0, 1.24d0, 1.43d0, 1.85d0, 2.31d0, 2.82d0, 3.40d0, 4.04d0, 4.76d0, 5.55d0,
     &     6.43d0, 7.39d0, 8.43d0, 9.57d0, 10.8d0, 12.1d0, 13.5d0, 15.2d0, 17.6d0, 20.2d0,
     &     23.1d0, 26.1d0, 29.2d0, 32.5d0, 35.9d0, 39.4d0, 43.1d0, 46.9d0, 50.8d0, 54.9d0,
     &     59.1d0, 63.5d0/
      if (sqrts.lt.mh(1).or.sqrts.gt.mh(52)) then
        write (*,*) 'WARNING in dssmgammahtab: sqrts out of range. Return dssmgammahtab=0'
        dssmgammahtab=0.d0
        return
      endif
      i=1
 1000 continue
      if (sqrts.gt.mh(i)) then
        i=i+1
        goto 1000
      endif
      dssmgammahtab=wh(i-1)+(sqrts-mh(i-1))/(mh(i)-mh(i-1))*(wh(i)-wh(i-1))
      return
      end

