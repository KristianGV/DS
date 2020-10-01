      real*8 function dsjfactor(labhalo,psi,theta)
**********************************************************************
***   the function dsjfactor gives the so-called Jfactor, namely the 
***   line of sight integral of the square of a spherical dark matter
***   density profile, to be implemented in gamma-ray and neutrino
***   fluxes from pair annihilations in DM halos (in the limit of zero
***   relative velocity); it links to dsomlosisph and saves values
***   under subsequent calls.
***   
***   inputs:
***     labhalo = halo model access label, to load the halo model via 
***         dsdmdselect_halomodel; NOTE: in case it is set to 'none'
***         no model is loaded and jfactor and dfactor are NOT
***         computed here, but given as an input
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
***   desc : J-factor (l.o.s. integral for annihilating DM)
***
***      
***   unit of return value: kpc sr * GeV^2 cm^-6
**********************************************************************
      implicit none
      include 'dsdmdcom.h'  ! psannihi,dmdradouttr
      include 'dslosidrvrcom.h' ! ilosigasph
      character(*) labhalo
      real*8 psi,theta
      integer power
      integer nthetain,npsfin
      real*8 thetain(2),dsomlosisph
      logical loadlab
      real*8 parvec(3),jfactor
      logical check,match
ccc
ccc loading the halo model corresponding to the label 'halolab' and
ccc transfer corresponding parameters    
ccc
      call dsdmdselect_halomodel(labhalo)
      loadlab=.false.  ! do not load again the halo in dsomlosisph
ccc
ccc check whether this jfactor has already been computed:
ccc      
      check=.true.
      parvec(1)=psi
      parvec(2)=theta
      call dsjfactorhdb(dmdihalotag,parvec,dmdtmp,check,match)
      if(match) then
        dsjfactor=parvec(3)
        return
      endif
ccc        
ccc compute the jfactor - value of line-of sight integration over
ccc a solid angle:      
ccc    \int d\Omega\int ds (rho(s))^2 [kpc sr GeV^2 cm^-6]
ccc        
      npsfin=1 ! integral over solid angle assuming step-function 
               ! response from the instrument within a cone of  
               ! aperture thetain(2) (in rad) in the direction
               ! thetain(1) (in rad)
      nthetain=2 ! 2 thetain entries
      thetain(1)=psi
      thetain(2)=theta
      power=psannihi
c      write(*,*) 'computing jfactor ...'
      jfactor=dsomlosisph(thetain,nthetain,npsfin,dmdradouttr,
     &         ilosigasph,power,labhalo,loadlab)
ccc
ccc store the value:
ccc      
      check=.false.
      parvec(3)=jfactor
      call dsjfactorhdb(dmdihalotag,parvec,dmdtmp,check,match)
      dsjfactor=jfactor
      return
      end



      subroutine dsjfactorhdb(dmdindex,parvec,tmp,check,match)
ccc
ccc tmp.eq.true/false -> temporary halo model/halo model in database
ccc check.eq.true/false -> check match/add new result
ccc      
      implicit none
      integer dmdindex
      real*8 parvec(3)
      logical tmp,check,match
ccc      
      integer ndatjmax,nparvec
      parameter(ndatjmax=100,nparvec=3)
      integer ntmpdatj,ndatj,dmdihatmpdatj,dmdihadatj(ndatjmax)
      real*8 tmpdatj(nparvec,ndatjmax),datj(nparvec,ndatjmax)
      common/jfactordatcom/tmpdatj,datj,dmdihatmpdatj,dmdihadatj,
     &  ntmpdatj,ndatj
ccc
      integer ii,jj
      real*8 dif,sum
      logical first
      data first/.true./
      save first
ccc
      if(first) then
        ndatj=0
        first=.false.
      endif
ccc
      if(check) then
        if(tmp) then
          if(dmdindex.ne.dmdihatmpdatj) then
            ntmpdatj=0
          else
            do ii=1,ntmpdatj
              match=.true.
              do jj=1,nparvec-1 
                dif=dabs(tmpdatj(jj,ii)-parvec(jj))
                sum=dabs(tmpdatj(jj,ii)+parvec(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo  
              if(match) then
                parvec(nparvec)=tmpdatj(nparvec,ii)
                return
              endif
            enddo
          endif
          match=.false.
          return
        else
          do ii=1,ndatj
            if(dmdindex.eq.dmdihadatj(ii)) then 
              match=.true.
              do jj=1,nparvec-1 
                dif=dabs(datj(jj,ii)-parvec(jj))
                sum=dabs(datj(jj,ii)+parvec(jj))
                if(dif.gt.1.d-7*sum) match=.false.
              enddo  
              if(match) then
                parvec(nparvec)=datj(nparvec,ii)
                return
              endif
            endif
          enddo
          match=.false.
          return
        endif
      else
        if(tmp) then
          ntmpdatj=ntmpdatj+1
          if(ntmpdatj.gt.ndatjmax) ntmpdatj=1
          dmdihatmpdatj=dmdindex
          do jj=1,nparvec 
            tmpdatj(jj,ntmpdatj)=parvec(jj)
          enddo
        else
          ndatj=ndatj+1
          if(ndatj.gt.ndatjmax) ndatj=1
          dmdihadatj(ndatj)=dmdindex
          do jj=1,nparvec 
            datj(jj,ndatj)=parvec(jj)
          enddo
        endif  
      endif      
      return
      end
        
