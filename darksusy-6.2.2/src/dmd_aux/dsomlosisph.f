      real*8 function dsomlosisph(thetain,nthetain,npsfin,outradius,
     &  ilosi,power,labhalo,loadlab)
************************************************************************
*** integral over the solid angle of the line-of-sight-integral of the 
*** product:              
***   dslosifunsph [which includes eventually exp(-absorption factor)]
***        * (eventually) gaussian (elliptical) instrument response 
***          function
*** for a given spherically symmetric configuration, i.e.:
***                                                          
***  = int_{PSF} dOmega 
***    int_{l.o.s.} dl dslosifunsph * dsomlosiweight         
***                                                          
*** inputs:                                                  
***   thetain(1) (in rad) the angular offset from the line of sight in 
***     the direction of the center of the distribution (except in case 
***     npsfin = 6, see below)
***
***   nthetain = number of entries and dimension of the vector thetain  
***
***   npsfin = 1 -> integral over solid angle assuming step-function 
***          response from the instrument within a cone of aperture 
***          thetain(2) (in rad)
***   npsfin = 2 -> integral over solid angle assuming gaussion response  
***          from the instrument with angular resolution thetain(2) 
***          (in rad)
***   npsfin = 3 -> integral over solid angle assuming gaussion  
***          elliptical response from the instrument with angular   
***          resolution thetain(2) (in rad) in the direction of the 
***          center of the system, angular resolution thetain(3) (in  
***          rad) perpendicular to it 
***   npsfin = 4 -> the same as npsfin = 2 except for thetain(2) (in 
***          rad) being the full width at half maximum 
***   npsfin = 5 -> the same as npsfin = 3 except for thetain(2,3) 
***          (in rad) being the full widths at half maximum 
***
***   npsfin = 6 -> integral over solid angle, assuming step-function 
***          response from the instrument within an angular window in
***          galactic coordinates: 
***                  b_min = thetain(1) in rad 
***                  b_max = thetain(2) in rad
***                  l_min = thetain(3) in rad
***                  l_max = thetain(4) in rad
***          for the moment you must restrict to: 
***             0 < b < pi/2                        
***             0 < l < pi                        
***
***   outradius (in kpc) is the maximal radius within which the l.o.s. 
***     integral has to be performed            
***
***   ilosi = option selecting which dslosifunsph to integrate; in
***     case it exceeds ilosiabssft (set in dslosidrvrcom.h) steps
***     to include absorption are performed       
***
***   power = 1 or 2 selecting annihilation or decay
***
***   labhalo = halo model access label, within the set of models 
***     defined in the halo model repository  
***   loadlab = logical variable: in case it is set to true the halo
***     labhalo is selected within this function, otherwise it is
***     assumed that it has been loaded before linking to this function
***     such as in the function dsgadphidEsph     
***
*** output:
***   result in kpc sr * [dslosifunsph]
************************************************************************
      implicit none
      include 'dsdmdcom.h' ! to interface dmd halo parameters
      include 'dsdmdintcom.h' ! setting maxradius,thetav(n),sphilosi
ccc
      integer nthetain,npsfin,ilosi,power
      real*8 thetain(nthetain),outradius
      character(*) labhalo
      logical loadlab
ccc
      real*8 thetasup,thetainf
      real*8 onemcosinf,onemcossup,tollf2,tollab,prec,resint,pi,partial,
     &  dsomlosisph1i,eps,onemcosthetacut,thetacut,sinthetacut,tinter,
     &  tsuptemp,check1,check2
      real*8 loctheta,pio2,costhetasup,costhetainf
      integer nmax,ihow
      external dsomlosisph1i,dsomlosisph1
ccc
      real*8 tinf,tsup
      common/losintextrcom/tinf,tsup
ccc      
      integer powerloc
      common/dslosifunsphcom/powerloc
ccc
      powerloc=power
ccc
ccc some consistency checks in the inputs
ccc
      if(nthetain.lt.2.or.nthetain.gt.10) then
        write(*,*) 'DS: in dsomlosisph, nthetain = ',nthetain
        write(*,*) 'DS: not properly set, program stopped'
        stop
      endif
      if(npsfin.lt.1.or.npsfin.gt.6) then
        write(*,*) 'DS: in dsomlosisph, npsf = ',npsfin
        write(*,*) 'DS: not properly set, program stopped'
        stop
      endif
      if(loadlab) then
ccc
ccc loading the halo model corresponding to the label 'halolab' and
ccc transfer corresponding parameters    
ccc
        call dsdmdselect_halomodel(labhalo)
      endif  
      objdistance=dmdobjdist
      radinnertr=dmdradintr
      radoutertr=dmdradouttr
ccc      
      maxradius=outradius
ccc
ccc if maxradius is larger than the outer truncation radius replaced it
ccc with the latter:
ccc
      if(maxradius.gt.radoutertr) maxradius=radoutertr
      pi=4.d0*datan(1.d0)
      sphilosi=ilosi
      npsf=npsfin
ccc
      if(npsf.eq.1) then    ! step function angular acceptance
        thetav(1)=thetain(1)  ! angle between direction of observation
                              ! and direction of the center of the 
                              ! distribution 
        thetav(2)=thetain(2)  ! angle defining the aperture of the cone
                            ! within which the angular integral is 
                            ! performed with weight function equal to 1 
      elseif(npsf.eq.2.or.npsf.eq.4) then ! gaussian angular acceptance
        thetav(1)=thetain(1)  ! angle between direction of observation
                              ! and direction of the center of the 
                              ! distribution 
        thetav(3)=thetain(2)  ! angular resolution or fwhm
        if(npsf.eq.4) then
          thetav(3)=thetav(3)/(2.d0*dsqrt(dlog(4.d0))) 
                                            ! fwhm -> angular resolution
          npsf=npsf-2
        endif
        thetav(2)=min(10.d0*thetav(3),pi/2.d0) ! 3 sigma: angle defining  
               ! the aperture of the cone within which the angular 
               ! integral is performed with gaussian weight function
      elseif(npsf.eq.3.or.npsf.eq.5) then ! gaussian elliptical angular 
                                          ! acceptance
        thetav(1)=thetain(1)  ! angle between direction of observation
                              ! and direction of the center of the 
                              ! distribution 
        thetav(3)=thetain(2)  ! angular resolution or fwhm, center dir.
        thetav(4)=thetain(3)  ! angular resolution or fwhm, perpend. dir.
        if(npsf.eq.5) then
          thetav(3)=thetav(3)/(2.d0*dsqrt(dlog(4.d0)))
                                            ! fwhm -> angular resolution
          thetav(4)=thetav(4)/(2.d0*dsqrt(dlog(4.d0)))
                                            ! fwhm -> angular resolution
          npsf=npsf-2
        endif
        thetav(2)=max(thetav(3),thetav(4))
        thetav(2)=min(10.d0*thetav(2),pi/2.d0) ! 3 sigma: angle defining  
               ! the aperture of the cone within which the angular 
               ! integral is performed with gaussian weight function
      endif
ccc
      if(npsf.ge.2.and.npsf.le.5.and.thetav(3).lt.4.8d-7) then
        ! warning that the result may not be accurate for such small 
        ! thetav(3)
        write(*,*) 'DS: you have calle dsomlosisph with an angular'
        write(*,*) 'DS: resolution = ',thetav(3)*2.06265d5,' arcmin,'
        write(*,*) 'DS: below the 0.1 arcmin resolution for which'
        write(*,*) 'DS: dsomlosisph is fully tested.'
        write(*,*) 'DS: be aware that the result might not be accurate'
        write(*,*) 'DS: and consider using dslosisph instead.'
      endif
ccc
      sinthetacut=radinnertr/objdistance
      thetacut=dasin(sinthetacut)
      if(thetacut.lt.1.d-4) then
        onemcosthetacut=thetacut**2/2.d0-thetacut**4/24.d0
      else
        onemcosthetacut=1.d0-dcos(thetacut)
      endif
      if(onemcosthetacut.le.0.d0) then
        write(*,*) 'DS: in dsomlosisph, onemcosthetacut = ',
     &                onemcosthetacut
        write(*,*) 'DS: not properly set, program stopped'
        stop
      endif
ccc
      if(npsf.ge.1.and.npsf.le.5) then
 100    continue
        thetasup=min(pi,thetav(1)+thetav(2))
        thetainf=max(0.d0,thetav(1)-thetav(2))
ccc
        if(thetav(2).lt.pi/2.d0*0.9d0.and.(thetav(1)-thetav(2)).ge.0.d0
     &   .and.(npsf.eq.2.or.npsf.eq.3)) then 
        ! extra check: in case a very extreme dsnpsph the contribution
        ! from the central region could still be relevant even if it is
        ! outside the 3 sigma cone; reset thetav(2) in such a case:  
          tinf=0.d0
          tsup=2.d0
          check1=dsomlosisph1i(onemcosthetacut)
          check2=dsomlosisph1i(1.d0-dcos(thetav(1)))
          if(check1.gt.1.d-8*check2) then
            thetav(2)=min(thetav(1)*1.1d0,pi/2.d0)
            goto 100
          endif
        endif
ccc
ccc compute the range of 1-cos(theta)
ccc
        if(thetainf.lt.1.d-4) then
          onemcosinf=thetainf**2/2.d0-thetainf**4/24.d0
        else
          onemcosinf=1.d0-dcos(thetainf)
        endif
        if(thetasup.eq.pi) then
          onemcossup=2.d0
        elseif(thetasup.lt.1.d-4) then
          onemcossup=thetasup**2/2.d0-thetasup**4/24.d0
        else
          onemcossup=1.d0-dcos(thetasup)
        endif
ccc
ccc
ccc
      elseif(npsf.eq.6) then
        thetav(1)=thetain(1)  ! minumum latitude b in rad
        thetav(2)=thetain(2)  ! maximum latitude b in rad
        if(thetav(1).gt.thetav(2)) then
          loctheta=thetav(1)
          thetav(1)=thetav(2)
          thetav(2)=loctheta
        endif
        thetav(3)=thetain(3)  ! minumum longitude l in rad
        thetav(4)=thetain(4)  ! maximum longitude l in rad
        if(thetav(3).gt.thetav(4)) then
          loctheta=thetav(3)
          thetav(3)=thetav(4)
          thetav(4)=loctheta
        endif
ccc
ccc check that: 0 < b < + pi/2   &   0 < l < + pi
ccc
        pio2=2.d0*datan(1.d0)
        if(thetav(1).gt.pio2.or.thetav(1).lt.0.d0) then
          write(*,*) 'DS: calling dsomlosisph with npsf = 6'
          write(*,*) 'DS: and minimum b out of range = ',thetav(1)
          write(*,*) 'DS: rather than [0, ',pio2,']'
          write(*,*) 'DS: program stopped'
          stop
        endif
        if(thetav(2).gt.pio2.or.thetav(2).lt.0.d0) then
          write(*,*) 'DS: calling dsomlosisph with npsf = 6'
          write(*,*) 'DS: and maximum b out of range = ',thetav(2)
          write(*,*) 'DS: rather than [0, ',pio2,']'
          write(*,*) 'DS: program stopped'
          stop
        endif
        if(thetav(3).gt.2.d0*pio2.or.thetav(3).lt.0.d0) then
          write(*,*) 'DS: calling dsomlosisph with npsf = 6'
          write(*,*) 'DS: and minimum b out of range = ',thetav(3)
          write(*,*) 'DS: rather than [0, ',2.d0*pio2,']'
          write(*,*) 'DS: program stopped'
          stop
        endif
        if(thetav(4).gt.2.d0*pio2.or.thetav(4).lt.0.d0) then
          write(*,*) 'DS: calling dsomlosisph with npsf = 6'
          write(*,*) 'DS: and minimum b out of range = ',thetav(4)
          write(*,*) 'DS: rather than [0, ',2.d0*pio2,']'
          write(*,*) 'DS: program stopped'
          stop
        endif
ccc
ccc compute the range of 1-cos(theta)
ccc
        if(dcos(thetain(4)).ge.0.d0) then
          costhetasup=dcos(thetain(1))*dcos(thetain(3))
          costhetainf=dcos(thetain(2))*dcos(thetain(4))
        elseif(dcos(thetain(3)).le.0.d0) then
          costhetainf=dcos(thetain(1))*dcos(thetain(4))
          costhetasup=dcos(thetain(2))*dcos(thetain(3))
        else
          costhetainf=dcos(thetain(1))*dcos(thetain(4))
          costhetasup=dcos(thetain(1))*dcos(thetain(3))
        endif
        onemcosinf=1.d0-costhetasup
        onemcossup=1.d0-costhetainf
      else
        write(*,*) 'DS: in dsomlosisph npsf = ',npsf,' not properly'
        write(*,*) 'DS: initialized, program stopped'
        stop
      endif
ccc
ccc     
ccc
      partial=0.d0
      if(onemcosthetacut.gt.onemcosinf) then
        tinf=onemcosinf
        tsup=min(onemcosthetacut,onemcossup)
        prec=1.d-4
        eps=(tsup-tinf)*prec*dsomlosisph1i(tsup)/10.d0
        call dsfun_intb(dsomlosisph1i,tinf,tsup,eps,prec,resint)
        partial=partial+resint
c        write(*,*) 'partial 1 : ',tinf,tsup,eps,resint,partial
      endif
      if(onemcossup.gt.onemcosthetacut) then
        tinf=max(onemcosthetacut,onemcosinf)
        tsup=onemcossup
c        write(*,*) 'tinf,tsup : ',tinf,tsup
        if(npsf.eq.2.or.npsf.eq.3) then
          tinter=1.d0-dcos(max(0.d0,thetav(1)-thetav(3)))
          if(tinter.gt.tinf) then            
          ihow=2
          tsuptemp=tsup
          tsup=1.d0-dcos(thetav(1))
          tollf2=2.d0
          tollab=(dlog(tsup)-dlog(tinf))/100.d0
          nmax=1000
          prec=1.d-4
          call dsfun_intparb(dsomlosisph1i,dsomlosisph1,tinf,tsup,
     &      tollf2,tollab,nmax,prec,ihow,resint)
          partial=partial+resint
c          write(*,*) 'partial 2 : ',tinf,tsup,resint,partial
          tinf=tsup
          tsup=tsuptemp
          endif
        endif
        ihow=2
        tollf2=2.d0
        tollab=(dlog(tsup)-dlog(tinf))/100.d0
        ihow=1
        tollf2=2.d0
        tollab=((tsup)-(tinf))/100.d0
        nmax=1000
        prec=1.d-4
        call dsfun_intparb(dsomlosisph1i,dsomlosisph1,tinf,tsup,
     &     tollf2,tollab,nmax,prec,ihow,resint)
        partial=partial+resint
c        write(*,*) 'partial 3 : ',tinf,tsup,resint,partial
      endif
      dsomlosisph=partial
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosisph1i(onemcostheta)
      implicit none
      include 'dsdmdintcom.h' ! npsf
      real*8 onemcostheta
      real*8 dsomlosisph1,dsomlosiweight,phifact
ccc
      real*8 tinf,tsup
      common/losintextrcom/tinf,tsup
ccc
ccc avoid troubles outside the interval of integration
ccc
      if(onemcostheta.gt.tsup.or.onemcostheta.lt.tinf) then
        dsomlosisph1i=0.d0
        return
      endif
ccc
      phifact=dsomlosiweight(onemcostheta)
      dsomlosisph1i=phifact*dsomlosisph1(onemcostheta)
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosiweight(onemcostheta)
ccc
ccc weight function for angular acceptance
ccc
      implicit none
      include 'dsdmdintcom.h' ! npsf
      real*8 onemcostheta,dsomlosicosweight,dsomlosigaweight,
     &  dsomlosiblweight
ccc
      if(npsf.eq.1) then
        dsomlosiweight=dsomlosicosweight(onemcostheta)
      elseif(npsf.eq.2.or.npsf.eq.3) then
        dsomlosiweight=dsomlosigaweight(onemcostheta)
      elseif(npsf.eq.6) then
        dsomlosiweight=dsomlosiblweight(onemcostheta)
      else
        write(*,*) 'DS: in dsomlosiweight npsf not properly set = ',
     &             npsf
        stop
      endif
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosicosweight(onemcostheta)
ccc
ccc step function angular acceptance and constant = 1 weight
ccc
      implicit none
      include 'dsdmdintcom.h' ! for thetav(1),thetav(2)
      real*8 onemcostheta
      real*8 sintheta,cosphi,phifact,sintheta1,theta,onemcosphi
ccc
      if(thetav(1).lt.1.d-15) then
        phifact=8.d0*datan(1.d0)
      else
        sintheta=dsqrt(2.d0-onemcostheta)
        sintheta=sintheta*dsqrt(onemcostheta)
        theta=dacos(1.d0-onemcostheta)
        if(dabs(thetav(1)).lt.1.d-4.and.dabs(thetav(2)).lt.1.d-4) then
          cosphi=(thetav(1)+thetav(2))/2.d0
          cosphi=cosphi*(1.d0-(thetav(1)**2+thetav(2)**2)/12.d0)
          cosphi=cosphi*(-thetav(2)+thetav(1))
          cosphi=cosphi+onemcostheta*dcos(thetav(1))
          if(sintheta.gt.1.d-16) then
            cosphi=cosphi/sintheta
          else
            cosphi=1.d16*cosphi
          endif
          sintheta1=dsin(thetav(1))
          if(sintheta1.gt.1.d-16) then
            cosphi=cosphi/sintheta1
          else
            cosphi=1.d16*cosphi
          endif
          goto 100
        endif
        if(dabs(theta-thetav(1)).lt.1.d-4
     &     .and.dabs(thetav(2)).lt.1.d-4) then
          onemcosphi=(thetav(2)-theta+thetav(1))/2.d0
          onemcosphi=onemcosphi
     &               *(1.d0+(thetav(2)**2+(theta-thetav(1))**2)/12.d0)
          onemcosphi=onemcosphi*(thetav(2)+theta-thetav(1))
          if(sintheta.gt.1.d-16) then
            onemcosphi=onemcosphi/sintheta
          else
            onemcosphi=1.d16*onemcosphi
          endif
          sintheta1=dsin(thetav(1))
          if(sintheta1.gt.1.d-16) then
            onemcosphi=onemcosphi/sintheta1
          else
            onemcosphi=1.d16*onemcosphi
          endif
          cosphi=1.d0-onemcosphi
          goto 100
        endif
        cosphi=dcos(thetav(2))-(1.d0-onemcostheta)*dcos(thetav(1))
        if(sintheta.gt.1.d-16) then
          cosphi=cosphi/sintheta
        else
          cosphi=1.d16*cosphi
        endif
        sintheta1=dsin(thetav(1))
        if(sintheta1.gt.1.d-16) then
          cosphi=cosphi/sintheta1
        else
          cosphi=1.d16*cosphi
        endif
 100    continue
        if(cosphi.ge.1.d0.or.cosphi.lt.-1.d0) then
          phifact=8.d0*datan(1.d0)
        else
          phifact=dabs(2.d0*dacos(cosphi))
        endif
      endif
      dsomlosicosweight=phifact
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosiblweight(onemcostheta)
ccc
ccc step function angular acceptance and constant = 1 weight
ccc within b,l window
ccc
      implicit none
      include 'dsdmdintcom.h' ! for thetav(1) -> thetav(4)
      real*8 onemcostheta
      real*8 cosbmin,cosbmax,cosbmin1,cosbmax1
      real*8 alpha,onemcosb,sin2alphamin,sin2alphamax,inf,sup
ccc
      if(dcos(thetav(3))*dcos(thetav(4)).ge.0.d0) then
        cosbmin1=min((1.d0-onemcostheta)/dcos(thetav(3)),
     &    (1.d0-onemcostheta)/dcos(thetav(4)))
        cosbmax1=max((1.d0-onemcostheta)/dcos(thetav(3)),
     &    (1.d0-onemcostheta)/dcos(thetav(4)))
      else
        if((1.d0-onemcostheta).ge.0.d0) then
          cosbmin1=(1.d0-onemcostheta)/dcos(thetav(3))
          cosbmax1=(1.d0-onemcostheta)/1.d-16
        else
          cosbmin1=(1.d0-onemcostheta)/dcos(thetav(4))
          cosbmax1=-(1.d0-onemcostheta)/1.d-16
        endif
      endif
      cosbmin=max(cosbmin1,dcos(thetav(2)))
      cosbmax=min(cosbmax1,dcos(thetav(1)))
ccc
      alpha=0.d0
      if(cosbmax.gt.cosbmin) then
        onemcosb=1.d0-cosbmax
        sin2alphamin=onemcosb*(2.d0-onemcosb)
     &    /onemcostheta/(2.d0-onemcostheta)
c        write(*,*) 'sin2 min = ',sin2alphamin
        if(sin2alphamin.gt.1.d0) sin2alphamin=1.d0
        onemcosb=1.d0-cosbmin
        sin2alphamax=onemcosb*(2.d0-onemcosb)
     &   /onemcostheta/(2.d0-onemcostheta)
c        write(*,*) 'sin2 max = ',sin2alphamax
        if(sin2alphamax.gt.1.d0) sin2alphamax=1.d0
        inf=dsqrt(sin2alphamin)
        sup=dsqrt(sin2alphamax)
        if(inf.lt.sup) alpha=dasin(sup)-dasin(inf) 
      endif
      dsomlosiblweight=alpha
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosigaweight(onemcostheta)
ccc
ccc gaussian (elliptical) angular acceptance
ccc
      implicit none
      include 'dsdmdintcom.h' ! for thetav(1),thetav(2)
      real*8 onemcostheta
      real*8 sintheta,cosphi,dsomlosigaweint,onemcosphi
      real*8 prec,eps,resint,sintheta1,theta
      external dsomlosigaweint
c
      real*8 onemctheta,phiinf,phisup
      common/omlosigawecom/onemctheta,phiinf,phisup
c
      onemctheta=onemcostheta
      if(onemcostheta.eq.0.d0.or.thetav(1).lt.1.d-15) then
        phisup=4.d0*datan(1.d0)        
      else
        sintheta=dsqrt(2.d0-onemcostheta)
        sintheta=sintheta*dsqrt(onemcostheta)
        theta=dacos(1.d0-onemcostheta)
        if(dabs(thetav(1)).lt.1.d-4.and.dabs(thetav(2)).lt.1.d-4) then
          cosphi=(thetav(1)+thetav(2))/2.d0
          cosphi=cosphi*(1.d0-(thetav(1)**2+thetav(2)**2)/12.d0)
          cosphi=cosphi*(-thetav(2)+thetav(1))
          cosphi=cosphi+onemcostheta*dcos(thetav(1))
          if(sintheta.gt.1.d-16) then
            cosphi=cosphi/sintheta
          else
            cosphi=1.d16*cosphi
          endif
          sintheta1=dsin(thetav(1))
          if(sintheta1.gt.1.d-16) then
            cosphi=cosphi/sintheta1
          else
            cosphi=1.d16*cosphi
          endif
          goto 100
        endif
        if(dabs(theta-thetav(1)).lt.1.d-4
     &     .and.dabs(thetav(2)).lt.1.d-4) then
          onemcosphi=(thetav(2)-theta+thetav(1))/2.d0
          onemcosphi=onemcosphi
     &               *(1.d0-(thetav(2)**2+(theta-thetav(1))**2)/12.d0)
          onemcosphi=onemcosphi*(thetav(2)+theta-thetav(1))
          if(sintheta.gt.1.d-16) then
            onemcosphi=onemcosphi/sintheta
          else
            onemcosphi=1.d16*onemcosphi
          endif
          sintheta1=dsin(thetav(1))
          if(sintheta1.gt.1.d-16) then
            onemcosphi=onemcosphi/sintheta1
          else
            onemcosphi=1.d16*onemcosphi
          endif
          cosphi=1.d0-onemcosphi
          goto 100
        endif
        cosphi=dcos(thetav(2))-(1.d0-onemcostheta)*dcos(thetav(1))
        if(sintheta.gt.1.d-16) then
          cosphi=cosphi/sintheta
        else
          cosphi=1.d16*cosphi
        endif
        sintheta1=dsin(thetav(1))
        if(sintheta1.gt.1.d-16) then
          cosphi=cosphi/sintheta1
        else
          cosphi=1.d16*cosphi
        endif
 100    continue
        if(cosphi.ge.1.d0.or.cosphi.le.-1.d0) then
          phisup=4.d0*datan(1.d0)
        else
          phisup=dabs(dacos(cosphi))
        endif
      endif
      phiinf=-phisup
      prec=1.d-3
      eps=(phisup-phiinf)*prec*dsomlosigaweint(0.d0)/100.d0
      call dsfun_int(dsomlosigaweint,phiinf,phisup,eps,prec,resint)
      dsomlosigaweight=resint
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosigaweint(phi)
      implicit none
      real*8 phi
      include 'dsdmdintcom.h' ! for thetav(1)
ccc  
      real*8 sintheta,onemcosthetanew,sinphinew2,tanphinew
      real*8 dsomlosigaweint2
ccc
      real*8 onemctheta,phiinf,phisup
      common/omlosigawecom/onemctheta,phiinf,phisup
ccc
      if(phi.gt.phisup.or.phi.lt.phiinf) then
        dsomlosigaweint=0.d0
        return
      endif
      sintheta=dsqrt(2.d0-onemctheta)
      sintheta=sintheta*dsqrt(onemctheta)
      if(thetav(1).lt.1.d-4) then
        onemcosthetanew=onemctheta
     &    +(thetav(1)**2/2.d0-thetav(1)**4/24.d0)*(1.d0-onemctheta)
     &    -dsin(thetav(1))*sintheta*dcos(phi)
      else
        onemcosthetanew=(1.d0-dcos(thetav(1)))
     &    +dcos(thetav(1))*onemctheta
     &    -dsin(thetav(1))*sintheta*dcos(phi)
        if(onemcosthetanew.lt.0.d0.and.onemcosthetanew.gt.-1.d-15)
     &  onemcosthetanew=0.d0
      endif 
      if(onemcosthetanew.gt.2.d0.or.onemcosthetanew.lt.0.d0) then
        write(*,*) 'DS: in dsomlosigaweint wrong onemcosthetanew',
     &    onemcosthetanew
        stop
      endif
      tanphinew=dsin(thetav(1))*dsin(phi)/
     &     (dcos(phi)*sintheta*dcos(thetav(1))
     &      +dsin(thetav(1))*(1.d0-onemctheta))  
                                           !!!! or is it a minus sign???
      if(dabs(tanphinew).lt.1.d0) then
        sinphinew2=tanphinew**2/(1.d0+tanphinew**2)
      else
        sinphinew2=1.d0/(1.d0+1.d0/tanphinew**2)
      endif
      dsomlosigaweint=dsomlosigaweint2(onemcosthetanew,sinphinew2)
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosigaweint2(onemcostheta,sinphi2)
      implicit none
      real*8 onemcostheta,sinphi2
      include 'dsdmdintcom.h' ! for npsf, thetav(3), thetav(4)
ccc  
      real*8 tantheta2,cosphi2,phifact
ccc
      if(onemcostheta.gt.2.d0.or.onemcostheta.lt.0.d0) then
        dsomlosigaweint2=1.d0
        return
      endif
      tantheta2=onemcostheta*(2.d0-onemcostheta)/
     &          (1.d0-onemcostheta*(2.d0-onemcostheta))
      cosphi2=1.d0-sinphi2
      if(npsf.eq.2) then ! gaussion response 
        phifact=dexp(-0.5d0*tantheta2/dtan(thetav(3))**2)
      elseif(npsf.eq.3) then ! gaussian elliptical response
        phifact=dexp(-0.5d0*(tantheta2/dtan(thetav(3))**2*cosphi2+
     &                       tantheta2/dtan(thetav(4))**2*sinphi2))
      else
        write(*,*) 'DS: in dsomlosigaweint2 npsf not properly set = ',
     &             npsf
        stop
      endif
      dsomlosigaweint2=phifact
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosisph1(onemcostheta)
ccc
ccc kpc * [dslosifunsph]
ccc
      implicit none
      include 'dsdmdintcom.h' ! for maxradius,sphilosi
                              ! setting: sinobjdist2,npairs_norm,locilosi 
      include 'dslosidrvrcom.h' ! for ilosiabssft,ilosinsph
ccc
      real*8 onemcostheta
      real*8 xxglobmax,sinobjdist,cosobjdist,eps,dsomlosisph2
ccc
      real*8 xxmax,xxmin,tollf2,tollab,prec,partial,result,theta
      real*8 xxminvec(3),xxmaxvec(3),weight(3),checktau,dstausph,
     &  dstaufunsph
      integer nstep,nstepnopar
      integer nmax,ihow,i,nstepi
      external dsomlosisph2
ccc
ccc check on whether absoption needs to be included or not
ccc
      taucut=1.d-4
      taumax=100.d0
      if(sphilosi.gt.ilosiabssft
     &   .and.sphilosi.le.(ilosiabssft+ilosinsph)) then
        theta=dasin(dsqrt(onemcostheta*(2.d0-onemcostheta)))
        checktau=dstausph(theta,maxradius,sphilosi)
        if(checktau.lt.taucut*10.d0) then
          sphilosii=sphilosi-ilosiabssft
        else
          sphilosii=sphilosi
        endif
      else
        sphilosii=sphilosi
      endif
ccc
      cosobjdist=(1.d0-onemcostheta)*objdistance
      sinobjdist2=onemcostheta*(2.d0-onemcostheta)
      sinobjdist2=sinobjdist2*objdistance
      sinobjdist2=sinobjdist2*objdistance
      sinobjdist=dsqrt(onemcostheta*(2.d0-onemcostheta))
      sinobjdist=sinobjdist*objdistance
      if(maxradius.lt.sinobjdist) then
        dsomlosisph1=0.d0
        return
      endif 
ccc
ccc in case of costheta greater than zero make the integration into 2 
ccc   steps: first from the miniminum (alias the closest point to the 
ccc   observer) to the maximum of dsnpsph, then from the maximum to the 
ccc   furthest point with respect to the observer
ccc
      if(onemcostheta.lt.1.d0-radinnertr/objdistance) then
        xxglobmax=dsqrt(maxradius**2-sinobjdist**2)
        xxmax=min(cosobjdist,xxglobmax)
        xxmin=radinnertr
        if(xxmax.lt.xxmin) then
          nstep=1
          nstepnopar=1
          xxminvec(1)=-xxmax
          xxmaxvec(1)=xxmax
          weight(1)=1.d0
        else
          if(sphilosii.gt.ilosiabssft
     &       .and.sphilosii.le.(ilosiabssft+ilosinsph)) then
            nstep=3
            nstepnopar=2
            xxminvec(1)=xxmin
            xxmaxvec(1)=xxmax
            weight(1)=1.d0
            xxminvec(2)=-xxmin
            xxmaxvec(2)=xxmin
            weight(2)=1.d0
            xxminvec(3)=xxmin
            xxmaxvec(3)=xxglobmax
            weight(3)=1.d0
          else
            xxminvec(1)=xxmin
            xxmaxvec(1)=xxmax
            weight(1)=2.d0
            nstepnopar=2
            xxminvec(2)=-xxmin
            xxmaxvec(2)=xxmin
            weight(2)=1.d0
            if(cosobjdist.gt.xxglobmax) then
              nstep=2
            else
              nstep=3
              xxminvec(3)=xxmax
              xxmaxvec(3)=xxglobmax
              weight(3)=1.d0
            endif
          endif
        endif
        partial=0.d0
        do i=1,nstep
          xxmin=xxminvec(i)
          xxmax=xxmaxvec(i)
          if(sphilosii.gt.ilosiabssft
     &       .and.sphilosii.le.(ilosiabssft+ilosinsph)) then
            nstepi=i
            if(i.eq.nstepnopar.and.i.ne.2) nstepi=4
            call dstausetsph(xxmin,xxmax,nstepi,sphilosii)
            if(dstaufunsph(xxmin,sphilosii).gt.taumax
     &         .and.dstaufunsph(xxmax,sphilosii).gt.taumax) then
              result=0.d0
              npairs_norm=1.d0
              goto 100
            endif
          endif
          npairs_norm=1.d0
          locilosi=sphilosii
          npairs_norm=max(dsomlosisph2(xxmin),dsomlosisph2(xxmax))
          if(npairs_norm.le.0.d0) then
            result=0.d0
          else
            if(i.eq.nstepnopar) then
              prec=1.d-4
              eps=2.d0*xxmax*prec*dsomlosisph2(xxmax)/10.d0
              call dsfun_int(dsomlosisph2,xxmin,xxmax,eps,prec,result)
            else
              tollf2=100.d0
              tollab=(dlog(xxmax)-dlog(xxmin))/1000.d0
              nmax=1000
              prec=1.d-5
              ihow=2
              call dsfun_intpar(dsomlosisph2,dsomlosisph2,xxmin,xxmax,
     &          tollf2,tollab,nmax,prec,ihow,result)
            endif
          endif
 100      partial=partial+result*npairs_norm*weight(i)
        enddo
        dsomlosisph1=partial
ccc
ccc in case of costheta costheta lower than or equal to zero make the 
ccc   integration in 1 step
ccc
      else ! if(onemcostheta.ge.1.d0-radinnertr/objdistance) then 
        partial=0.d0
        xxmin=dabs(cosobjdist)
        if(xxmin.lt.radinnertr) then
          if(onemcostheta.le.1.d0) then
            xxmin=-xxmin
            xxmax=radinnertr
          else
            xxmax=radinnertr
          endif     
          npairs_norm=1.d0
          locilosi=sphilosii
          if(locilosi.gt.ilosiabssft
     &       .and.locilosi.le.(ilosiabssft+ilosinsph))
     &        locilosi=locilosi-ilosiabssft ! neglect absoption in this
                                           ! small interval
          npairs_norm=max(dsomlosisph2(xxmin),dsomlosisph2(xxmax))
          if(npairs_norm.le.0.d0) then
            dsomlosisph1=0.d0
            return
          endif
          prec=1.d-4
          eps=(xxmax-xxmin)*prec*dsomlosisph2(xxmin)/10.d0
          call dsfun_int(dsomlosisph2,xxmin,xxmax,eps,prec,result)
          partial=partial+result*npairs_norm 
          xxmin=radinnertr
        endif
        xxmax=dsqrt(maxradius**2-sinobjdist2)
        if(sphilosii.gt.ilosiabssft
     &     .and.sphilosii.le.(ilosiabssft+ilosinsph)) then
          nstepi=5
          call dstausetsph(xxmin,xxmax,nstepi,sphilosii)
          if(dstaufunsph(xxmin,sphilosii).gt.taumax
     &       .and.dstaufunsph(xxmax,sphilosii).gt.taumax) then
            dsomlosisph1=0.d0
            return
          endif
        endif
        npairs_norm=1.d0
        locilosi=sphilosii
        npairs_norm=max(dsomlosisph2(xxmin),dsomlosisph2(xxmax))
        if(npairs_norm.le.0.d0) then
          dsomlosisph1=0.d0
          return
        endif
        partial=0.d0
        tollf2=10.d0
        tollab=(dlog(xxmax)-dlog(xxmin))/1000.d0
        nmax=1000
        prec=1.d-5
        ihow=2
        call dsfun_intpar(dsomlosisph2,dsomlosisph2,xxmin,xxmax,
     &    tollf2,tollab,nmax,prec,ihow,result)
        partial=partial+result*npairs_norm
        dsomlosisph1=partial
      endif
      return
      end
ccc
ccc
ccc
      real*8 function dsomlosisph2(xx)
ccc
ccc dimensionless
ccc
      implicit none
      include 'dsdmdintcom.h' ! for maxradius,sinobjdist2,npairs_norm,locilosi
      include 'dslosidrvrcom.h' ! for ilosiabssft,ilosinsph
      real*8 xx,rr,rr2,dslosifunsph,dstaufunsph
ccc
      rr2=sinobjdist2
      rr2=rr2+xx**2
      rr=dsqrt(rr2)
      if(rr.gt.maxradius*1.00000001d0) then
        dsomlosisph2=0.d0
        return
      endif
      dsomlosisph2=dslosifunsph(rr,locilosi)/npairs_norm
      if(dsomlosisph2.lt.0.d0) then
        write(*,*) 'DS: dsomlosisph2 = ',dsomlosisph2
        stop
      endif  
      if(locilosi.gt.ilosiabssft
     &   .and.locilosi.le.(ilosiabssft+ilosinsph)) then
        dsomlosisph2=dsomlosisph2*dexp(-dstaufunsph(xx,locilosi))
      endif
      return
      end

