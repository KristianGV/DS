      function dssmgammah_hdecay_tab(ichannel,sqrts)
***
*** Returns the standard model partial Higgs decay width for channel ichannel
*** This tabulation is made with HDECAY.      
***
*** List of channels:
***   ichallen = 0  :  total Width      
***   ichannel = 1  :  nue + anti-nue (outside 90GeV-300GeV only)
***   ichannel = 2  :  e+ + e- (outside 90GeV-300GeV only)
***   ichannel = 3  :  numu + anti-numu (outside 90GeV-300GeV only)
***   ichannel = 4  :  mu+ + mu-
***   ichannel = 5  :  nutau + anti-nutau (outside 90GeV-300GeV only)
***   ichannel = 6  :  tau+ + tau-
***   ichannel = 7  :  u + ubar (outside 90GeV-300GeV only)
***   ichannel = 8  :  d + dbar (outside 90GeV-300GeV only)
***   ichannel = 9  :  c + cbar
***   ichannel = 10  :  s + sbar
***   ichannel = 11  :  t + tbar
***   ichannel = 12  :  b + bbar
***   ichannel = 13  :  gamma + gamma (90GeV-300GeV only)
***   ichannel = 14  :  W+ + W-
***   ichannel = 15  :  Z + Z
***   ichannel = 16  :  g + g (90GeV-300GeV only)
***   ichannel = 17  :  H + H
***   ichannel = 18  :  Z + gamma (90GeV-300GeV only)
***
*** Author: Paolo Gondolo 2016
*** Modified: Joakim Edsjö, to use hdecay tables instead of results from
*** the literature      
      implicit none
      include 'dsio.h'
      
      real*8 dssmgammah_hdecay_tab,sqrts
      integer ichannel,i
      character*128 file1,file2,scr

      integer ni
      parameter(ni=1020)
      real*8 mh(ni+1),wh(ni+1,0:18)
      save mh,wh
      real*8 w1,w2
      
      logical first
      data first/.true./
      save first

      include 'dsdir.h'

      if (first) then
c...Note: the files contain the mass and branching fractions and
c...the total width. As we want this routine to return the total width
c...we multiply the branching fractions with the total width
c...when interpolating the table         
         file1=dsdatapath
         call dscharadd(file1,'br-hdecay.sm1')
         file2=dsdatapath
         call dscharadd(file2,'br-hdecay.sm2')
         if (prtlevel.ge.2) then
           write(*,*) 'dssmgammah_hdecay_tab: reading in files...'
           write(*,*) file1
           write(*,*) file2
         endif
         open(unit=81,file=file1,status='old',form='formatted')
         open(unit=82,file=file2,status='old',form='formatted')
         do i=1,3
            read(81,*) scr ! header line
            read(82,*) scr ! header line
         enddo
         
         do i=1,ni
            wh(i,1)=0.d0 ! nu_e nu_e-bar
            wh(i,2)=0.d0 ! e+ e-
            wh(i,3)=0.d0 ! nu_mu nu_mu-bar
            wh(i,5)=0.d0 ! nu_tau nu_tau-bar
            wh(i,7)=0.d0 ! u u-bar
            wh(i,8)=0.d0 ! d d-bar
            wh(i,17)=0.d0 ! H H
            read(81,*) mh(i),wh(i,12),wh(i,6),wh(i,4),wh(i,10),
     &           wh(i,9),wh(i,11)
            read(82,*) mh(i),wh(i,16),wh(i,13),wh(i,18),wh(i,14),
     &           wh(i,15),wh(i,0)
         enddo

c...Make sure we have one extra element for interpolation at end         
         mh(ni+1)=mh(ni)+mh(ni)-mh(ni-1)
         do i=0,18
            wh(ni+1,i)=wh(ni,i)
         enddo
         
         close(81)
         close(82)
         if (prtlevel.ge.3) write(*,*) 'done.'
         first=.false.
      endif

c...Start interpolating      
      dssmgammah_hdecay_tab=0.d0
      if (sqrts.lt.mh(1).or.sqrts.gt.mh(ni)) then
         write (*,*) 'WARNING in dssmgammah_hdecay_tab:',
     &        'sqrts out of range.',
     &        ' Returning dssmgammah_hdecay_tab=0.'
        write(*,*) '  sqrts = ',sqrts
        dssmgammah_hdecay_tab=0.d0
        return
      endif
      
      if (ichannel.lt.1.or.ichannel.gt.18) then
         write (*,*) 'WARNING in dssmgammah_hdecay_tab:',
     &        ' untabulated channel=',ichannel,
     &        '. Returning dssmgammah_hdecay_tab=0.'
        dssmgammah_hdecay_tab=0.d0
        return
      endif

      call dshunt(mh,ni,sqrts,i)

      w1=wh(i,ichannel)*wh(i,0)
      w2=wh(i+1,ichannel)*wh(i+1,0)
      dssmgammah_hdecay_tab=w1
     &  +(sqrts-mh(i))/(mh(i+1)-mh(i))*(w2-w1)

      return
      end

