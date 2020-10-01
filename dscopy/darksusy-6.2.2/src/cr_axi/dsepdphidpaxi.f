      real*8 function dsepdphidpaxi(ppin,phiin,how,labhalo)
************************************************************************
*** function which computes the local galactic differential flux of 
*** positrons at the momentum pp as a result of dark matter annihilation  
*** or decay in the halo.      
*** inputs:
***     ppin - positron momentum in GeV
***     phiin - solar modulation parameter in GV, assuming that solar 
***             modulation can be treated with the force-field method
***             - if phiin.le.1.d-3, the interstellar flux is returned
***     how = 2 - a table is tabulated on first  call, and then 
***               interpolated
***           3 - as 2, but also write the table to disk at the 
***               first call
***           4 - read table from disk on first call, and use the 
***               subsequent calls. If the file does not exist, it 
***               will be created (as in 3). (use as default)
***     labhalo = halo model access label, within the set of models 
***               defined in the halo model repository  
*** output: cm^-2 s^-1 GeV^-1 sr^-1
***
***   type : commonly used
***   desc : Local galactic differential positron flux from dark matter
***
***
***  author: Piero Ullio
************************************************************************
      implicit none
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      include 'dscraxicom.h' ! for ivopt
      include 'dsmpconst.h'
      real*8 ppin,phiin
      integer iv,how
      character(*) labhalo
      real*8 pp,tpin,tp,smrf,R,z,par,dsepdndpaxi,vep
      real*8 source,eecheck,dsmwimp,dscrsource,dscrsource_line
      real*8 eline,widthline
      integer power,pdg2line,istat
      logical loadlab
ccc
ccc loading the halo model corresponding to the label 'halolab' and
ccc transfer corresponding parameters    
ccc
      call dsdmdselect_halomodel(labhalo)      
ccc
ccc check whether you have called this function for a halo model which
ccc can be used for Milky Way rates or not
ccc      
      if(.not.dmdmw) then
        write(*,*) 'DS: call to dsdbdphidtaxi with the halo label: ',
     &     labhalo    
        write(*,*) 'DS: which has not been initialized as suitable for'
        write(*,*) 'DS: Milky Way rates. Program stopped'
        stop
      endif  
      loadlab=.false.
ccc
      R=dmdobjdist
      z=0.d0
      iv=ivopt
ccc
ccc eventually add in solar modulation, pp is the interstellar momentum 
ccc
      tpin=ppin*dsqrt(1.d0+(m_e/ppin)**2)-m_e/ppin
      call dsffsolmodaxi(tpin,phiin,m_e,1.d0,1.d0,tp,smrf)
      pp=tp*dsqrt(1.d0+2.d0*m_e/tp)
      par=0.d0      
ccc
ccc include (eventually) a decay term decay and annihilation:        
ccc
      power=psdecay
      eecheck=pp
      source=dscrsource(eecheck,1,-11,power,0.0d0,istat)
      eecheck=dsmwimp()*0.95d0
      source=source+dscrsource(eecheck,1,-11,power,0.0d0,istat)
      source=source
     &   +dscrsource_line(-11,1,power,0.0d0,eline,widthline,pdg2line,
     &        istat)
      if(source.gt.0.d0) then
        par=par+dsepdndpaxi(R,z,pp,power,iv,how,labhalo,loadlab) ! cm^-3 GeV^-1
      endif
ccc
      power=psannihi
      eecheck=pp
      source=dscrsource(eecheck,1,-11,power,0.0d0,istat)
      eecheck=dsmwimp()*0.95d0
      source=source+dscrsource(eecheck,1,-11,power,0.0d0,istat)
      source=source
     &   +dscrsource_line(-11,1,power,0.0d0,eline,widthline,pdg2line,
     &        istat)
      if(source.gt.0.d0) then
        par=par+dsepdndpaxi(R,z,pp,power,iv,how,labhalo,loadlab) ! cm^-3 GeV^-1
      endif
c...convert density to flux
      vep=c_light*1d5 ! cm s^-1
      dsepdphidpaxi=smrf*par*vep/(4.0d0*pi) ! cm^-2 s^-1 GeV^-1 sr^-1
      return
      end

