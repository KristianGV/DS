      subroutine dsddffsd(q,a,z,l2jjpp,l2jjnn,l2jjpn,ierr)
c_______________________________________________________________________
c  Spin-dependent structure functions for direct detection.
c  input:
c    q : real*8  : momentum transfer in GeV ( q=sqrt(M*E/2/mu^2) )
c    a : integer : mass number
c    z : integer : atomic number
c  output:
c    l2jjnn, l2jjpp, l2jjpn : the factors \lambda^2 J (J+1) made functions
c       of the momentum transfer q, for pp, nn, and pn; they are
c       related to Spp, Snn, and Spn through
c           S = \lambda^2 J (J+1) (2J+1)/\pi
c       In turn, Spp, Snn, and Spn are related to the spin-dependent 
c       structure functions S_{00}(q), S_{01}(q), S_{11}(q) defined in
c       Engel (PLB264,114,1991) via
c           Spp = S00+S11+S01
c           Snn = S00+S11-S01
c           Spn = 2*(S00-S11)
c  author: paolo gondolo (paolo@physics.utah.edu) 2004
c  modified: pg 040605 switched s01 and s11 in Na-23
c  modified: pg 080217 changed to \lambda^2 J (J+1)
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsio.h'
      include 'dsmpconst.h'

      real*8 q
      integer a,z,ierr,jerr
      real*8 l2jjpp,l2jjnn,l2jjpn
      real*8 myq

      ierr=0
      jerr=0
      
c set q=0 if requested
      myq=q
      if (ddsf(2)(1:1).eq.'0') myq=0.d0

      if (a.eq.1.and.z.eq.1) then ! proton
        l2jjpp = 0.75d0
        l2jjnn = 0.d0
        l2jjpn = 0.d0
        
      elseif (a.eq.1.and.z.eq.0) then ! neutron
        l2jjpp = 0.d0
        l2jjnn = 0.75d0
        l2jjpn = 0.d0
        
c best available form factor
        
      elseif (ddsf(2)(2:).eq.'best') then
        call dsddffism(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)
        jerr=1                  ! JE TMP
        if (jerr.ne.0) then
          if (prtlevel.gt.1) write (*,*) 
     &         'dsddffsd: switching to OddG....'
          call dsddffoddg(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)
          if (jerr.ne.0) then
            if (prtlevel.gt.1) write (*,*) 
     &           'dsddffsd: switching to SPSM....'
            call dsddffspsm(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)
            if (jerr.ne.0) then
              if (prtlevel.gt.1) write (*,*) 
     &             'dsddffsd: switching to Simple....'
              call dsddffsimsd(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)
              if (jerr.ne.0) then
                if (prtlevel.gt.1) write (*,*) 
     &               'dsddffsd: giving up.... sigma set to zero'
                l2jjpp=0.d0
                l2jjnn=0.d0
                l2jjpn=0.d0
                if (prtlevel.ge.0) write(*,'(a,i3,a,i3)') 'WARNING dsddffsd: no spin-dependent form factor found for a=',a,' z=',z
              endif
            endif
          endif
        endif
        
c interacting shell model

      elseif (ddsf(2)(2:).eq.'ISM'.or.ddsf(2)(2:).eq.'ISMR') then
        call dsddffism(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)

c single particle shell model

      elseif (ddsf(2)(2:).eq.'SPSM') then
        call dsddffspsm(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)

c odd-group model

      elseif (ddsf(2)(2:).eq.'OddG') then
        call dsddffoddg(myq,a,z,l2jjpp,l2jjnn,l2jjpn,jerr)

c unrecognized form factor label

      else
        write (*,*) 
     &       'dsddffsd : Unrecognized spin-dependent form factor',ddsf(2)
        ierr=1
        return
      endif
      
      if (abs(l2jjpp).gt.2.d2.or.(abs(l2jjnn).gt.2.d2).or.(abs(l2jjpn).gt.2.d2)) then
        write (*,*) 
     &       'dsddffsd : SD form factor too large. Suspect numerical extrapolation problem.'
        ierr=2
        return
      endif
      
      if (jerr.ne.0) ierr=jerr
      
      return
      end
