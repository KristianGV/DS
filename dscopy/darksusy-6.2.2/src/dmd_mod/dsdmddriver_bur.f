      subroutine dsdmddriver_bur(iwin,nrein,rein,nchin,chin,labin,reout)
c_______________________________________________________________________
c this is the sample driver to initialize, print or link corresponding
c functions in case of a dark matter density profile described by the
c spherical Burkert profile, i.e.:
c
c     rho(r) = rhosbur/(1+(r/rsbur))/(1+(r/rsbur)^2)
c      
c inputs:
c   - iwin: integer setting the action of the driver; it can be equal
c     to:
c       idmdcreate -> profile creation is completed here (first settings
c         are in the dsdmddriver subroutine)      
c       idmdparprint -> the current halo profile is printed
c         depending on the input labin, see below      
c       idmddensph -> the driver links to the spherical density profile
c   - labin: external label of the current halo model: you got here in
c       case of labin containing the string 'bur'; if it contains also
c       the string 'mw' the model is intepreted as referring to the
c       Milky Way, you can define it assigning the local halo density
c       rather than the reference density 'rhosbur' and an internal
c       label for tabulation files is automatically created      
c   - rein(nrein) real*8 vector for different inputs according to
c      different iwin values       
c   - chin(nchin) character*10 vector for input parameter labels
c output:         
c   - reout: real*8 output for the different function linking as
c     specified by the input value iwin 
c_______________________________________________________________________
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmddrvrcom.h'
      include 'dsdvcom.h'
      include 'dsmpconst.h'
      include 'dsio.h'
      integer iwin
      character(*) labin
      integer nrein,nchin
      real*8 rein(nrein),reout
      character(*) chin(nchin)
      character*10 labelout
      real*8 dsdmdsetlabel
      external dsdmdsetlabel
ccc      
      integer ii,jj,indexout,ndmdstart
      real*8 r,rs,rhos,x,dsdmdprof_bur,radouttr,radintr,normx,rho0,
     &  dsI1_bur,units
      logical mwlab,match,setrho
ccc
      if(iwin.eq.idmdcreate.or.iwin.eq.idmdparprint) then
c
c set first the list of parameters which needs to be associated to
c an internal tag (via dsdmdsetlabel.f & dsdmdgetlabel.f):        
c the reference density needs to be the first entry since in case of
c Milky Way profiles you can set instead the value of the density at
c the local position in the Galaxy, see below.
c        
        idmdlab(1)=irhosbur
        chdmdlab(1)=chrhosbur
        idmdlab(2)=irsbur
        chdmdlab(2)=chrsbur
        idmdlab(3)=iradintr
        chdmdlab(3)=chradintr
        idmdlab(4)=iradouttr
        chdmdlab(4)=chradouttr
        ndmdlab=4
        sufflab=bursuff
c
c add then extra parameters to be also stored for this halo profile:
c
        idmdlab(5)=iobjdist
        chdmdlab(5)=chobjdist
        ndmdlabtot=5
      endif
c
c then go thorugh the list of iwin options:
c      
      if(iwin.eq.idmdcreate) then
c
c load whether it is Milky Way or not
c        
        mwlab=.false.
        if(ihalotag(iwhichobject,ihalocur).eq.imwsuff) mwlab=.true. 
c        
c complete the creation of the dark matter density profile: if it is
c not the Milky Way load all label parameters, if it is the Milky Way
c skip rhosbur at first        
c
        ndmdstart=1
        if(mwlab) ndmdstart=2 
        do ii=ndmdstart,ndmdlabtot
          match=.false. 
          do jj=1,nchin
            if(chin(jj).eq.chdmdlab(ii)) then
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(idmdlab(ii),ihalocur)=indexout
              match=.true.
            endif
          enddo
          if(.not.match) then
            write(*,*) 'DS: in dsdmddriver_bur missing initialization'
            write(*,*) 'DS: of the halo parameter : ',chdmdlab(ii)
            write(*,*) 'DS: when setting the halo model : ',labin
            write(*,*) 'DS: program stopped '
            stop
          endif  
        enddo
        if(mwlab) then
c
c adding rhosbur, possibly in terms of the local halo density
c
          setrho=.false.
          do jj=1,nchin
            if(chin(jj).eq.chrho0) then      
              r=dmdparvec(ihalotag(iobjdist,ihalocur))
              rs=dmdparvec(ihalotag(irsbur,ihalocur))
              x=r/rs
              normx=dsdmdprof_bur(x)
              normx=rein(jj)/normx
              call dsdmdsetind(indexout,normx)
              ihalotag(irhosbur,ihalocur)=indexout
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(irho0,ihalocur)=indexout
              setrho=.true.
            endif
          enddo
          do jj=1,nchin
            if(chin(jj).eq.chrhosbur) then
              if(setrho) then
                write(*,*) 'DS: inconsistency !!! In dsdmddriver_bur'
                write(*,*) 'DS: you are setting both rho0 and rhos'
                stop
              endif
              call dsdmdsetind(indexout,rein(jj))
              ihalotag(irhosbur,ihalocur)=indexout
              r=dmdparvec(ihalotag(iobjdist,ihalocur))
              rs=dmdparvec(ihalotag(irsbur,ihalocur))
              rhos=dmdparvec(ihalotag(irhosbur,ihalocur))
              x=r/rs
              rho0=rhos*dsdmdprof_bur(x)
              call dsdmdsetind(indexout,rho0)
              ihalotag(irho0,ihalocur)=indexout
              setrho=.true.
            endif  
          enddo
          if(.not.setrho) then
            write(*,*) 'DS: Error: I n dsdmddriver_bur you are'
            write(*,*) 'DS: not setting neither rho0 nor rhos'
            write(*,*) 'DS: program stopped '
            stop
          endif
        else
c
c if this is not the Milky Way, rho0 is set to zero: 
c
          rho0=0.d0
          call dsdmdsetind(indexout,rho0)
          ihalotag(irho0,ihalocur)=indexout
        endif  
c 
c if it is Milky Way and a non-temporarily-defined halo model load
c a label tag for CR green function tabulations:
c
        if(mwlab.and.(ihalocur.ne.ihalotmp)) then
          call dsdmdgetlabel(dsdmdsetlabel,labelout)
          halotabtag(ihalocur)=labelout
c
c NOTE: in the future you may need tabulations also for models
c different from the Milky Way, add them below  
c          
        else
          halotabtag(ihalocur)='none'
        endif  
        return
      endif
ccc
      if (iwin.eq.idmdparprint) then
c         if (prtlevel.gt.2) then
            write(*,*)
            write(*,*) 'DS: the halo model with access label: '
     &        ,halotag(ihalocur)
            write(*,*) 'DS: has a corresponding automatically generated'
            write(*,*) 'DS: label for tabulation files equal to: '
     &        ,halotabtag(ihalocur)
            write(*,*) 'DS: the list of parameters for this model is: '
            do ii=1,ndmdlabtot
               write(*,*) 'DS: par & value : ',chdmdlab(ii)
     &                    ,dmdparvec(ihalotag(idmdlab(ii),ihalocur)) 
            enddo
c        endif
        return
      endif
ccc
      if(iwin.eq.idmddensph) then
c        
c link to spherical dark matter density profile with inner density
c cutoff and outer truncation radius:
c
        r=rein(idvrsph)
        radouttr=dmdparvec(ihalotag(iradouttr,ihalocur))
        if(r.gt.radouttr) then
          reout=0.d0
          return
        endif  
        radintr=dmdparvec(ihalotag(iradintr,ihalocur))
        r=max(r,radintr)
        rs=dmdparvec(ihalotag(irsbur,ihalocur))
        rhos=dmdparvec(ihalotag(irhosbur,ihalocur))
        x=r/rs
        reout=rhos*dsdmdprof_bur(x)  ! GeV cm^-3
        return
      endif
ccc
      if(iwin.eq.idmdmassph) then
c        
c link to spherical dark matter mass profile with outer truncation
c radius:
c
        r=rein(idvrsph)
        radouttr=dmdparvec(ihalotag(iradouttr,ihalocur))
        rs=dmdparvec(ihalotag(irsbur,ihalocur))
        rhos=dmdparvec(ihalotag(irhosbur,ihalocur))
        x=(min(r,radouttr))/rs
        reout=rhos*4.d0*pi*rs**3*dsI1_bur(x) ! GeV cm^-3 kpc^3
        units=(cmperpc*1.d3)**3/GeVperSolarMass
        reout=reout*units ! M_sun
        return
      endif
ccc
ccc you called this function with an iwin which is not set, an error
ccc message is written and the program is stopped      
ccc      
      write(*,*) 'DS: linking to dsdmddriver_bur with iwin = ',iwin
      write(*,*) 'DS: out of the defined range of values'
      write(*,*) 'DS: program stopped'
      stop
      end
