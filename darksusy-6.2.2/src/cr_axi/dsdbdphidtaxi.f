      real*8 function dsdbdphidtaxi(tdin,phiin,how,labhalo)
************************************************************************
*** function which computes the local galactic differential flux of 
*** antideuterons at the kinetic energy per nucleon tdin as a result of 
*** dark matter pair annihilations or decay in the halo.      
*** inputs:
***     tdin = antideuteron kinetic energy per nucleon td in GeV
***     phiin - solar modulation parameter in GV, assuming that solar 
***             modulation can be treated with the force-field method
***             - if phiin.le.1.d-3, the interstellar flux is returned
***     how = 1 - calculate the flux only for requested tdin
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
***   type : commonly used
***   desc : Local galactic differential antideuteron flux from dark matter
***
***
***  author: Piero Ullio
***  modified: Joakim Edsjo, to use new yields and crsource directly    
************************************************************************
      implicit none
      include 'dsdmdcom.h'  ! to interface dmd halo parameters
      include 'dscraxicom.h' ! for isopt
      include 'dsmpconst.h'
      real*8 tdin,phiin
      integer is,how,istat,power
      character(*) labhalo
      real*8 td,smrf,R,z,par,dsdbtdaxitab,dspbdmasaxi,vdb,tdd,
     &  dscrsource,source
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
      is=isopt
ccc
ccc eventually add in solar modulation, td is the interstellar kinetic 
ccc energy: 
ccc
      call dsffsolmodaxi(tdin,phiin,m_d,2.d0,1.d0,td,smrf)
ccc
      par=0.d0
ccc
ccc include (eventually) a decay term decay and annihilation:        
ccc
      power=psdecay
      ! dscrsource uses per particle energies, we here want per nucleon
      ! this gives a factor of two in argument and result
      source=dscrsource(td*2.d0,1,-1000010020,power,0.d0,istat)*2.d0
      if(source.gt.0.d0) then
        tdd=dsdbtdaxitab(R,z,td,power,is,how,labhalo,loadlab) ! 10^15 s
        par=par+source*tdd*dspbdmasaxi(R,z,power)
      endif
ccc
ccc include (eventually) an annihilation term:        
ccc
      power=psannihi
      ! dscrsource uses per particle energies, we here want per nucleon
      ! this gives a factor of two in argument and result
      source=dscrsource(td*2.d0,1,-1000010020,power,0.d0,istat)*2.d0
      if(source.gt.0.d0) then
        tdd=dsdbtdaxitab(R,z,td,power,is,how,labhalo,loadlab) ! 10^15 s
        par=par+source*tdd*dspbdmasaxi(R,z,power)
      endif
      
      vdb=dsqrt(1.d0+m_d/td)/(1.d0+m_d/td/2.d0)*c_light*1d5 ! cm s^-1
      dsdbdphidtaxi=smrf*1.d15*par*vdb/(4.0d0*pi) 
                      ! cm^-2 s^-1 GeV^-1 sr^-1      

      return
      end
