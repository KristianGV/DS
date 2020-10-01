      subroutine dsdmddriver_choice(iwin,nrein,rein,nchin,chin,labin,
     &  reout)
c_______________________________________________________________________
c subroutine looping over available halo model drivers; this is the
c default version; if you need to add other profiles you need to replace
c this subroutine      
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      include 'dsdvcom.h'
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
      logical match
c
c if you are calling this function to create a profile you need first
c to identify which profile 
c      
      if(iwin.eq.idmdcreate) then
c
c is it a spherical nfw? if 'labin' contains the string 'nfwsuff' it
c is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchnfwsuff,nfwsuff,match)
        if(match) then
c store that it is a spherical nfw:         
          ihalotag(iwhichprof,ihalocur)=infwsuff
c store the coordinate type that the spherical nfw assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
c
c is it a spherical burkert? if 'labin' contains the string 'bursuff' it
c is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchbursuff,bursuff,match)
        if(match) then
c store that it is a spherical burkert:         
          ihalotag(iwhichprof,ihalocur)=ibursuff
c store the coordinate type that the spherical burkert assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
c
c is it a spherical einasto? if 'labin' contains the string 'einsuff' it
c is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,ncheinsuff,einsuff,match)
        if(match) then
c store that it is a spherical einasto:         
          ihalotag(iwhichprof,ihalocur)=ieinsuff
c store the coordinate type that the spherical einasto assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
c
c is it a spherical model to be interpolated from a table? if 'labin'
c contains the string 'numsuff' it is assumed to be so:
c        
        call dslabcheck(nchhalotag,labin,nchnumsuff,numsuff,match)
        if(match) then
c store that it is a radially tabulated profile:         
          ihalotag(iwhichprof,ihalocur)=inumsuff
c store the coordinate type that this radial profile assumes:         
          ihalotag(itypeprof,ihalocur)=itydvsph
          goto 100
        endif
ccc
ccc you called this function with a label for which we cannot find a match
ccc a message is written and the program is stopped      
ccc      
        write(*,*) 'DS: linking to dsdmddriver_choice with option'
        write(*,*) 'DS: idmdcreate and labin = ',labin
        write(*,*) 'DS: not containing any of the defined suffices'
        write(*,*) 'DS: namely:',nfwsuff,' & ',bursuff,' & ',einsuff
        write(*,*) 'DS: program stopped'
        stop
 100    continue
      endif
ccc
      if(ihalotag(iwhichprof,ihalocur).eq.infwsuff) then
        call dsdmddriver_nfw(iwin,nrein,rein,nchin,chin,labin,reout)
      elseif(ihalotag(iwhichprof,ihalocur).eq.ibursuff) then
        call dsdmddriver_bur(iwin,nrein,rein,nchin,chin,labin,reout)
      elseif(ihalotag(iwhichprof,ihalocur).eq.ieinsuff) then
        call dsdmddriver_ein(iwin,nrein,rein,nchin,chin,labin,reout)
      elseif(ihalotag(iwhichprof,ihalocur).eq.inumsuff) then
        call dsdmddriver_num(iwin,nrein,rein,nchin,chin,labin,reout)
      endif
ccc
      return
      end
