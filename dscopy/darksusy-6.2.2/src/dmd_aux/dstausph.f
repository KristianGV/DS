      real*8 function dstausph(theta0in,outerradius,ilosi)
************************************************************************
*** function giving the total optical depth as needed for spherical
*** signals including absorption
***
*** a call to dslosidriver selects whether an explicit function is set
*** therein or dstausph is computed here via an integral over the line
*** of sight of an appropriate entry in the function dslosifunsph for a 
*** given spherically symmetric configuration, i.e.:   
***                                                          
***  = int_l.o.s. dl dslosifunsph
***
*** as it applies in case of self-absorption      
***                                                          
*** inputs:                                                  
***   theta0 (in rad) the angular offset from the line of sight in the 
***     direction of the center of the distribution 
***   outerradius (in kpc) is the maximal radius within which the l.o.s. 
***     integral has to be performed            
***   ilosi = option selecting which dslosifunsph to integrate; you
***     should link with ilosi in the interval:
***       (ilosiabssft+1,ilosiabssft+ilosinsph)
***     as set in dslosidrvrcom.h
***
*** inputs from common blocks:                               
***   objdistance (in kpc) is the distance between the observer and the 
***     center of the spherical system          
***   radinnertr (in kpc) is the innner truncation radius below which 
***     the integrated function is assumed to be constant
***                                                          
*** output: dimensionless (setting a proper link in dslosifunsph)
************************************************************************
      implicit none
      include 'dsdmdcom.h' ! for objdistance,radinnertr,radoutertr
      include 'dsdmdintcom.h' ! for taureschoweq1,maxradius,locilosi
      include 'dslosidrvrcom.h'
      include 'dsdvcom.h'
      real*8 theta0in,outerradius
      integer ilosi
ccc
      real*8 out,dstausph1,onemcostheta0,pi
      integer iwin
      iwin=ilositautot+ilosi
      dvsph(idvrsph)=theta0in  ! this is conventional for 1 dependent variable
      call dslosidriver(iwin,ntysph,dvsph,out)
ccc      
      if(locilosi.eq.ilosiabscomp) then
ccc
ccc dstausph computed in dslosidriver
ccc        
        dstausph=out
        return
      endif
ccc
ccc otherwise compute it as a losi:
ccc
      maxradius=outerradius
ccc
ccc if maxradius is larger than the outer truncation radius replaced it
ccc with the latter:
ccc
      if(maxradius.gt.radoutertr) maxradius=radoutertr
      pi=4.d0*datan(1.d0)
      if(theta0in.eq.0.d0) then
        onemcostheta0=0.d0
      elseif(theta0in.eq.pi) then
        onemcostheta0=2.d0
      else
        onemcostheta0=1.d0-dcos(theta0in)
      endif
      dstausph=dstausph1(onemcostheta0) ! kpc * [dslosifunsph] 
      dstausph=dstausph*3.0856d21       ! cm * [dslosifunsph] 
      dstausph=dstausph*taureschoweq1
      return
      end


c_______________________________________________________________________

      real*8 function dstausph1(onemcostheta)
      implicit none
      include 'dsdmdcom.h'
      include 'dsdmdintcom.h'
      real*8 onemcostheta
ccc
      real*8 xxglobmax,sinobjdist,cosobjdist,eps,dsomlosisph2
ccc
      real*8 xxmax,xxmin,tollf2,tollab,prec,partial,result
      real*8 xxminvec(3),xxmaxvec(3),weight(3)
      integer nstep,nstepnopar
      integer nmax,ihow,i
      external dsomlosisph2
ccc
ccc in case of costheta greater than zero make the integration into 2 
ccc   steps: first from the miniminum (alias the closest point to the 
ccc   observer) to the maximum of dsnpsph, then from the maximum to the 
ccc   furthest point with respect to the observer
ccc
      if(onemcostheta.lt.1.d0*(1.d0+1.d-15)) then
        cosobjdist=(1.d0-onemcostheta)*objdistance
        sinobjdist2=onemcostheta*(2.d0-onemcostheta)
        sinobjdist2=sinobjdist2*objdistance
        sinobjdist2=sinobjdist2*objdistance
        sinobjdist=dsqrt(onemcostheta*(2.d0-onemcostheta))
        sinobjdist=sinobjdist*objdistance
        if(maxradius.lt.sinobjdist) then
          dstausph1=0.d0
          return
        endif 
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
        partial=0.d0
        do i=1,nstep
          xxmin=xxminvec(i)
          xxmax=xxmaxvec(i)
          npairs_norm=1.d0
          npairs_norm=dsomlosisph2(0.d0)
          if(npairs_norm.le.0.d0) then
            result=0.d0
          else
            if(i.eq.nstepnopar) then
              prec=1.d-4
              eps=(xxmax-xxmin)*prec*dsomlosisph2(xxmax)/10.d0
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
          partial=partial+result*npairs_norm*weight(i)
c          write(*,*) i,xxmin,xxmax,result,weight(i),partial
        enddo
        dstausph1=partial
c
c in case of costheta costheta lower than or equal to zero make the 
c   integration in 1 step
c
      else ! if(onemcostheta.ge.1.d0) then
        if(maxradius.le.objdistance) then
          dstausph1=0.d0
          return
        endif 
        cosobjdist=(1.d0-onemcostheta)*objdistance
        sinobjdist2=onemcostheta*(2.d0-onemcostheta)
        sinobjdist2=sinobjdist2*objdistance
        sinobjdist2=sinobjdist2*objdistance
        xxmax=dsqrt(maxradius**2-sinobjdist2)
        xxmin=max(1.d-15,dabs(cosobjdist))
        npairs_norm=1.d0
        npairs_norm=dsomlosisph2(0.d0)
        if(npairs_norm.le.0.d0) then
          dstausph1=0.d0
          return
        endif
c
        partial=0.d0
c
        tollf2=100.d0
        tollab=(dlog(xxmax)-dlog(xxmin))/1000.d0
        nmax=1000
        prec=1.d-5
        ihow=2
        call dsfun_intpar(dsomlosisph2,dsomlosisph2,xxmin,xxmax,
     &    tollf2,tollab,nmax,prec,ihow,result)
        partial=partial+result*npairs_norm
        dstausph1=partial
      endif
      return
      end
