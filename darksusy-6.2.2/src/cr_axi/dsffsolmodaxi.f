      subroutine dsffsolmodaxi(tpin,phiin,mass,A,modZ,tp,smrf)
************************************************************************
***  subroutine which computes the shift in kinetic energy per nucleon
***  to account for solar modulation in the local flux of cosmic rays.
***  the 1-parameter force-field method is assumed.  
***
***  inputs:
***     tpin = kinetic energy per nucleon at earth position [GeV]
***     phiin = solar modulation parameter in GV, assuming that solar 
***             modulation can be treated with the force-field method
***             - if phiin.le.1.d-3, tp = tpin is returned
***     mass = mass of the cr particle
***     A = number of nucleons of the cr particle
***     modZ = absolute value of the charge of the cr particle      
***           
***  output:
***     tp  = kinetic energy per nucleon in the local interstellar (LIS)
***           medium [GeV]
***     smrf = flux rescaling factor: flux_earth = smrf * flux_interst. 
***
***  NOTE: this routine assumes the rigid shift:
***          E_LIS = E_earth + modZ * phiin
***  and not, as sometimes done, a rigid shift in rigidity (p/modZ):
***          R_LIS = R_earth + phiin
***  this second option is commented out within the file.      
***  author: Piero Ullio
************************************************************************
      implicit none
      real*8 tpin,phiin,mass,A,modZ,tp,smrf
      real*8 eein,ee
c      real*8 ppin,pp
ccc
      if(phiin.le.1.d-3) then
        tp=tpin
        smrf=1.d0
        return
      endif  
      eein=A*tpin+mass
      ee=eein+modZ*phiin
      tp=(ee-mass)/A
      smrf=(eein**2-mass**2)/(ee**2-mass**2)
c      eein=A*tpin+mass
c      ppin=dsqrt(eein**2-mass**2)
c      pp=ppin+modZ*phiin
c      ee=dsqrt(pp**2+mass**2)
c      tp=(ee-mass)/A
c      smrf=ppin**2/pp**2
      return
      end
