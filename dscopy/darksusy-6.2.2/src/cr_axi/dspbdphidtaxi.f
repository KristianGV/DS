      real*8 function dspbdphidtaxi(tpin,phiin,how,labhalo)
************************************************************************
*** function which computes the local galactic differential flux of 
*** antiprotons at kinetic energy tpin as a result of dark matter 
*** annihilation or decay in the halo.      
***
*** type : commonly used
*** desc : Local galactic differential antiproton flux from dark matter
***
*** inputs:
***     tpin = antiproton kinetic energy tp in GeV
***     phiin - solar modulation parameter in GV, assuming that solar 
***             modulation can be treated with the force-field method
***             - if phiin.le.1.d-3, the interstellar flux is returned
***     how = 1 - calculate the flux only for requested tpin
***           2 - a table is tabulated on first  call, and then 
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
***  author: Piero Ullio
************************************************************************
      implicit none
      include 'dsdmdcom.h' ! to interface dmd halo parameters
      include 'dscraxicom.h' ! for isopt
      include 'dsmpconst.h' 
      real*8 tpin,phiin
      integer is,how,istat,power
      character(*) labhalo
      real*8 tp,smrf,R,z,par,dspbtdaxitab,dspbdmasaxi,vpb,dscrsource,
     &     source,td
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
        write(*,*) 'DS: call to dspbdphidtaxi with the halo label: ',
     &     labhalo    
        write(*,*) 'DS: which has not been initialized as suitable for'
        write(*,*) 'DS: Milky Way rates. Program stopped'
        stop
      endif
      loadlab=.false.
ccc
      R=dmdobjdist
      z=0.d0
      is=isopt
ccc
ccc eventually add in solar modulation, tp is the interstellar kinetic 
ccc energy: 
ccc
      call dsffsolmodaxi(tpin,phiin,m_p,1.d0,1.d0,tp,smrf)
      par=0.d0
ccc
ccc include (eventually) a decay term decay and annihilation:        
ccc
      power=psdecay
      source=dscrsource(tp,1,-2212,power,0.d0,istat)
      if(source.gt.0.d0) then
        td=dspbtdaxitab(R,z,tp,power,is,how,labhalo,loadlab) ! 10^15 s
        par=par+source*td*dspbdmasaxi(R,z,power)
      endif
ccc
ccc include (eventually) an annihilation term:        
ccc
      power=psannihi
      source=dscrsource(tp,1,-2212,power,0.d0,istat)
      if(source.gt.0.d0) then
        td=dspbtdaxitab(R,z,tp,power,is,how,labhalo,loadlab) ! 10^15 s
        par=par+source*td*dspbdmasaxi(R,z,power)
      endif
ccc
      vpb=dsqrt(1.d0+2.d0*m_p/tp)/(1.d0+m_p/tp)*c_light*1.d5 ! cm s^-1
      dspbdphidtaxi=smrf*1.d15*par*vpb/(4.0d0*pi) 
                      ! cm^-2 s^-1 GeV^-1 sr^-1
      return
      end
