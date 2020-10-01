      real*8 function dslosisph(theta0in,outradius,ilosi,power,labhalo,
     &  loadlab)
************************************************************************
*** integral over the line of sight of the product:              
***   dslosifunsph * (eventually) exp(-absorption factor)
*** for a given spherically symmetric configuration, i.e.:
***                                                          
***  = int_{l.o.s.} dl dslosifunsph * (eventually) dexp(-dstaufunsph)         
***
*** inputs:                                                  
***   thetain(1) (in rad) the angular offset from the line of sight in 
***     the direction of the center of the distribution 
***   outradius (in kpc) is the maximal radius within which the l.o.s. 
***     integral has to be performed            
***   ilosi = option selecting which dslosifunsph to integrate; in
***     case it exceeds ilosiabssft (set in dslosidrvrcom.h) steps
***     to include absorption are performed       
***   power = 1 or 2 selecting annihilation or decay
***   labhalo = halo model access label, within the set of models defined
***     in the halo model repository  
***   loadlab = logical variable: in case it is set to true the halo
***     labhalo is selected within this function, otherwise it is
***     assumed that it has been loaded before linking to this function
***     such as in the function dsgadphidEsph     
***
*** output:
***   result in kpc * [dslosifunsph]
************************************************************************
      implicit none
      include 'dsdmdcom.h' ! to interface dmd halo parameters
      include 'dsdmdintcom.h' ! setting maxradius,sphilosi
      include 'dslosidrvrcom.h'
      real*8 theta0in,outradius
      integer ilosi,power
      character(*) labhalo
      logical loadlab
      real*8 dslosisph1,onemcostheta0,pi,checktau,dstausph
ccc
ccc temporary FIXME!
ccc      
      integer powerloc
      common/dslosifunsphcom/powerloc
      powerloc=power
ccc
ccc loading the halo model corresponding to the label 'halolab' and
ccc transfer corresponding parameters    
ccc
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
      taucut=1.d-4
      taumax=100.d0
cccc
cccc check on whether absoption needs to be included or not
cccc
      if(sphilosii.gt.ilosiabssft
     &       .and.sphilosii.le.(ilosiabssft+ilosinsph)) then
c        write(*,*) 'DS: computing optical depth'
        checktau=dstausph(theta0in,maxradius,sphilosi)
        if(checktau.lt.taucut*10.d0) then
c          write(*,*) 'DS: theta, optical depth : ',theta0in,checktau
c          write(*,*) 'DS: switching off absorption'
          sphilosii=sphilosi-ilosiabssft
        else
          sphilosii=sphilosi
        endif
c        write(*,*) 'DS: optical depth = ',checktau
      else
        sphilosii=sphilosi
      endif
ccc
      if(theta0in.eq.0.d0) then
        onemcostheta0=0.d0
      elseif(theta0in.eq.pi) then
        onemcostheta0=2.d0
      else
        onemcostheta0=1.d0-dcos(theta0in)
      endif
      dslosisph=dslosisph1(onemcostheta0)
      return
      end
ccc
ccc
ccc
      real*8 function dslosisph1(onemcostheta)
      implicit none
      include 'dsdmdintcom.h' ! for maxradius,sphilosii
                              ! setting: sinobjdist2,npairs_norm,locilosi 
      include 'dslosidrvrcom.h' ! for ilosiabssft,ilosinsph
      real*8 onemcostheta
ccc
      real*8 xxglobmax,sinobjdist,cosobjdist,
     &  eps,dsomlosisph2
      real*8 xxmax,xxmin,tollf2,tollab,prec,partial,result,dstaufunsph
      real*8 xxminvec(3),xxmaxvec(3),weight(3)
      integer nstep,nstepnopar
      integer nmax,ihow,i,nstepi
      external dsomlosisph2
ccc
      cosobjdist=(1.d0-onemcostheta)*objdistance
      sinobjdist2=onemcostheta*(2.d0-onemcostheta)
      sinobjdist2=sinobjdist2*objdistance
      sinobjdist2=sinobjdist2*objdistance
      sinobjdist=dsqrt(onemcostheta*(2.d0-onemcostheta))
      sinobjdist=sinobjdist*objdistance
      if(maxradius.lt.sinobjdist) then
        dslosisph1=0.d0
        return
      endif
ccc
ccc in case of costheta greater than zero make the integration into 3 
ccc   steps: first from the miniminum (alias the closest point to the 
ccc   observer) to the maximum of dsnpaxi, then from the maximum to the 
ccc   furthest point with respect to the observer; integrals are in log 
ccc   scales, hence a small integration innterval around $x=0$ needs to be 
ccc   added.
ccc
      if(onemcostheta.lt.1.d0-radinnertr/objdistance) then
      xxmin=radinnertr
      xxglobmax=dsqrt(maxradius**2-sinobjdist**2)
      xxmax=min(cosobjdist,xxglobmax)
      if(xxmax.lt.xxmin) then
        nstep=1
        nstepnopar=1
        xxminvec(1)=-xxmax
        xxmaxvec(1)=xxmax
        weight(1)=1.d0
ccc  absorption neglected in this small interval
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
     &       .and.dstaufunsph(xxmax,sphilosii).gt.taumax) then
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
     &      tollf2,tollab,nmax,prec,ihow,result)
          endif
        endif
 100    partial=partial+result*npairs_norm*weight(i)
      enddo
      dslosisph1=partial
ccc
ccc in case of costheta lower than or equal to zero make the 
ccc   integration in 1 or 2 steps
ccc
      else ! if(onemcostheta.ge.1.d0-radinnertr/objdistance) then
      partial=0.d0
c
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
     &     .and.locilosi.le.(ilosiabssft+ilosinsph))
     &      locilosi=locilosi-ilosiabssft ! neglect absoption in this
                                         ! small interval
        npairs_norm=max(dsomlosisph2(xxmin),dsomlosisph2(xxmax))
        if(npairs_norm.le.0.d0) then
          dslosisph1=0.d0
          return
        endif
ccc
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
     &     .and.dstaufunsph(xxmax,sphilosii).gt.taumax) then
          dslosisph1=0.d0
          return
        endif
      endif
      npairs_norm=1.d0
      locilosi=sphilosii
      npairs_norm=max(dsomlosisph2(xxmin),dsomlosisph2(xxmax))
      if(npairs_norm.le.0.d0) then
        dslosisph1=0.d0
        return
      endif
ccc
      tollf2=100.d0
      tollab=(dlog(xxmax)-dlog(xxmin))/1000.d0
      nmax=1000
      prec=1.d-5
      ihow=2
      call dsfun_intpar(dsomlosisph2,dsomlosisph2,xxmin,xxmax,
     &  tollf2,tollab,nmax,prec,ihow,result)
      partial=partial+result*npairs_norm
      dslosisph1=partial
      endif
      return
      end


