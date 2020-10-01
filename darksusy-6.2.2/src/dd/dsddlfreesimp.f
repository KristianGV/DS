*******************************************************************************
*** Function dsddlfreesimp returns a simple estimate for the mean free      ***
*** path of strongly interacting DM through the Earth crust. We include     ***
*** contributions from the 11 most abundant elements, with densities as     ***
*** returned from dssem_earthdenscomp.                                      ***
***                                                                         ***
***  Input:                                                                 ***
***    sigsi - spin-independent cross section per nucleon [cm**2]           ***
***    depth - detector location below surface [cm]                         ***
***    how   - assume equal couplings of DM to protons and neutrons (how=1) ***
***            or assume DM only couples to protons (how=2)                 ***
***                                                                         ***
***  Output: *average* mean free path in cm, between surface and depth      ***
***          specified as input.                                            ***
***                                                                         ***
*** author: Torsten.Bringmann.fys.uio.no                                    ***
*** date 2019-02-06                                                         ***
*******************************************************************************
      real*8 function dsddlfreesimp(sigsi, depth, how)
      implicit none
      include 'dsio.h'
      include 'dsmpconst.h'
      include 'dsddcom.h'

      real*8 sigsi, depth
      integer how
      real*8 linv, dsmwimp, mdm, r
      integer i
      real*8 nav(Nelements), mN(Nelements)
      integer acom
      common /earthdensaux/ acom 
      
      real*8 dsddearthdens_aux, dsf_int, dssem_earthdenscomp
      external dsddearthdens_aux

c... This is Table 2 in 1611.05453 
c      data nav/345.0,177.0,117.0,611.0,7.94,1.47,23.3,10.9/   ! in units of 10^20 cm^-3
 
      linv=1.0d-50
      if (sigsi.lt.0.0d0) goto 100
      if (sigsi.lt.1.d-35) goto 50 ! no attenuation

      mdm = dsmwimp()
      do i=1,Nelements
         mn(i)  = mNaU(i)*atomicmassunit ! convert masses to GeV
         acom   = an(i)
c... calculate average density -- this would be needed only if the density 
c... changes on the way to the detector location
c         nav(i) = dsf_int(dsddearthdens_aux,r_earth-depth*1.0d-2,r_earth,1.0d-3)
c         nav(i) = 1.0d20*nav(i) / (1.0d-2*depth) ! average density in cm**-3
         r = r_earth-5.d2 ! any depth just below surface...
         nav(i) = dssem_earthdenscomp(r,acom)
c         write(*,*) i, an(i), mn(i), nav(i)
      enddo
      
      if (how.eq.1) then

        linv=0.0
        do i=1,Nelements
           linv = linv + nav(i)*AN(i)**2*mN(i)**3/(mn(i)+mdm)**4
        enddo
        linv = linv*2*sigsi*(mdm+m_p)*(mdm+m_n)*mdm/m_p/m_n

      elseif (how.eq.2) then

        linv=0.0
        do i=i,Nelements
           linv = linv + nav(i)*ZN(i)**2*mN(i)**3/(mn(i)+mdm)**4
        enddo
        linv = linv*2*sigsi*(mdm+m_p)**2*mdm/m_p**2

      else
        if (prtlevel.gt.1) then
          write(*,*) 'ERROR in dsddlfreesimp:'
          write(*,*) 'unknown option how = ',how
          write(*,*) 'Will return 1.0d40 cm for mean free path...'
        endif  
      endif      

 50   dsddlfreesimp=1.0/linv            
      return

 100  write(*,*) 'FATAL ERROR in dsddlfreesimp:'
      write(*,*) 'You have supplied a negative cross section!'
      stop

      end



*******************************************************************************
*** auxiliary routine just for integration
*******************************************************************************
      real*8 function dsddearthdens_aux(r) ! earth radius in m
      implicit none
      real*8 r, dssem_earthdenscomp
      integer acom
      common /earthdensaux/ acom 
      
c... rescale to improve convergence
      dsddearthdens_aux = 1.0d-20*dssem_earthdenscomp(r,acom)

      return
      end
     
