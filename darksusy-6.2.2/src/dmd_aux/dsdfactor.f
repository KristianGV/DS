      real*8 function dsdfactor(labhalo,psi,theta)
**********************************************************************
***   the function dsdfactor gives the so-called Dfactor, namely the 
***   line of sight integral of the square of a spherical dark matter
***   density profile, to be implemented in gamma-ray and neutrino
***   fluxes from particle decays in DM halos; it links to dsomlosisph
***   and saves values under subsequent calls.
***   
***   inputs:
***     labhalo = halo model access label, to load the halo model via 
***         dsdmdselect_halomodel
***     psi = the angular offset between the pointing direction
***         and the direction of the center of the distributionin
***         (in rad)
***     theta = aperture of the acceptance cone (defines the
***         solid angle within which line-of-sight-integrals are 
***         performed assuming a step-function response from the 
***         instrument; more options are available in dsomlosisph, use
***         that to compute line-of-sight-integrals in your main file)
***
***   type : commonly used
***   desc : D-factor (l.o.s. integral for decaying DM)
***
***      
***   unit of return value: kpc sr * GeV cm^-3
**********************************************************************
      implicit none
      include 'dsdmdcom.h'  ! psdecay,dmdradouttr
      include 'dslosidrvrcom.h' ! ilosigasph
      character(*) labhalo
      real*8 psi,theta
      integer power
      integer nthetain,npsfin
      real*8 thetain(2),dsomlosisph
      logical loadlab
      real*8 parvec(3),dfactor
      logical check,match
ccc
ccc loading the halo model corresponding to the label 'halolab' and
ccc transfer corresponding parameters    
ccc
      call dsdmdselect_halomodel(labhalo)
      loadlab=.false.  ! do not load again the halo in dsomlosisph
ccc
ccc check whether this dfactor has already been computed:
ccc      
      check=.true.
      parvec(1)=psi
      parvec(2)=theta
      call dsdfactorhdb(dmdihalotag,parvec,dmdtmp,check,match)
      if(match) then
        dsdfactor=parvec(3)
        return
      endif
ccc        
ccc compute the dfactor - value of line-of sight integration over
ccc a solid angle:      
ccc    \int d\Omega\int ds rho(s) [kpc sr GeV cm^-3]
ccc        
      npsfin=1 ! integral over solid angle assuming step-function 
               ! response from the instrument within a cone of  
               ! aperture thetain(2) (in rad) in the direction
               ! thetain(1) (in rad)
      nthetain=2 ! 2 thetain entries
      thetain(1)=psi
      thetain(2)=theta
      power=psdecay
c      write(*,*) 'computing dfactor ...'
      dfactor=dsomlosisph(thetain,nthetain,npsfin,dmdradouttr,
     &         ilosigasph,power,labhalo,loadlab)
ccc
ccc store the value:
ccc      
      check=.false.
      parvec(3)=dfactor
      call dsdfactorhdb(dmdihalotag,parvec,dmdtmp,check,match)
      dsdfactor=dfactor
      return
      end



      subroutine dsdfactorhdb(dmdindex,parvec,tmp,check,match)
ccc
ccc tmp.eq.true/false -> temporary halo model/halo model in database
ccc check.eq.true/false -> check match/add new result
ccc      
      implicit none
      integer dmdindex
      real*8 parvec(3)
      logical tmp,check,match
ccc      
      integer ndatdmax,nparvec
      parameter(ndatdmax=100,nparvec=3)
      integer ntmpdatd,ndatd,dmdihatmpdatd,dmdihadatd(ndatdmax)
      real*8 tmpdatd(nparvec,ndatdmax),datd(nparvec,ndatdmax)
      common/dfactordatcom/tmpdatd,datd,dmdihatmpdatd,dmdihadatd,
     &  ntmpdatd,ndatd
ccc
      integer ii,jj
      real*8 dif,sum
      logical first
      data first/.true./
      save first
ccc
      if(first) then
        ndatd=0
        first=.false.
      endif
ccc
      if(check) then
        if(tmp) then
          if(dmdindex.ne.dmdihatmpdatd) then
            ntmpdatd=0
          else
            do ii=1,ntmpdatd
              match=.true.
              do jj=1,nparvec-1 
                dif=dabs(tmpdatd(jj,ii)-parvec(jj))
                sum=dabs(tmpdatd(jj,ii)+parvec(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo  
              if(match) then
                parvec(nparvec)=tmpdatd(nparvec,ii)
                return
              endif
            enddo
          endif
          match=.false.
          return
        else
          do ii=1,ndatd
            if(dmdindex.eq.dmdihadatd(ii)) then 
              match=.true.
              do jj=1,nparvec-1 
                dif=dabs(datd(jj,ii)-parvec(jj))
                sum=dabs(datd(jj,ii)+parvec(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo  
              if(match) then
                parvec(nparvec)=datd(nparvec,ii)
                return
              endif
            endif
          enddo
          match=.false.
          return
        endif
      else
        if(tmp) then
          ntmpdatd=ntmpdatd+1
          if(ntmpdatd.gt.ndatdmax) ntmpdatd=1
          dmdihatmpdatd=dmdindex
          do jj=1,nparvec 
            tmpdatd(jj,ntmpdatd)=parvec(jj)
          enddo
        else
          ndatd=ndatd+1
          if(ndatd.gt.ndatdmax) ndatd=1
          dmdihadatd(ndatd)=dmdindex
          do jj=1,nparvec 
            datd(jj,ndatd)=parvec(jj)
          enddo
        endif  
      endif      
      return
      end
        
