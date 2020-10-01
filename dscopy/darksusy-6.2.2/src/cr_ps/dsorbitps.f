************************************************************************
*** subroutine dsorbitps sets the orbit of a point source.  
***                                                                       
***  type : REPLACEABLE                                                   
***                                                                       
*** inputs:
***     t0 - time at which the position along the orbit is computed
***          (10^15 s)
***     it assumes that the observer is located at x,y,z=0
*** inputs through common blocks:
***     horbit - which type of orbit you are considering; currently
***     the available option is:
***     horbit=1 - motion on a straight line in the reference frame
***       in which vy=0 and y0 is constant = y0ini (kpc); 
***       the trajectory is set by 3 extra parameters, which we define 
***       to be the velocity in the x and z direction, vxin and 
***       vzin (km s^-1) and the distance of the trajectory to the 
***       observer, dminin (kpc) (of course,  dminin >= y0in);
***       time is normalized assuming t0=0 when the source is at 
***       the distance dminin
***     horbit=2 - motion on a circular orbit in the galactic plane 
***       with circular velocity vcirc and radius of the orbit radorb
***       assuming that the observer is, in galactic coordinates, at
***       x=rosb, y=z=0 (and that t0=0 corresponds to the source at the
***       position x0=radorb, y0=z0=0). the initial time ti is set 
***       equal the time when the source is located at x0 = - radorb, 
***       eventually a number of full orbits "norbit" before the 
***       current orbit. for the output a coordinate change is made and
***       the observer is assumed to be at x,y,z=0 
***
*** outputs:
***     x0,y0,z0 - position of the source (kpc) at t0
***     ti - time at which the source enters the diffusion region
***     tf - time at which the source gets out of the diffusion region
***          (10^15 s)
***
***  author: Piero Ullio
***  modified: Torsten Bringmann, 11/06/2015 (adapted to replaceable   
***                                           function concept)
************************************************************************
      subroutine dsorbitps(t0,x0,y0,z0,ti,tf)
      implicit none
      include 'dscraxicom.h'
      include 'dsmpconst.h' 
      real*8 t0,x0,y0,z0,ti,tf
      real*8 vx,vy,vz,v2,xi,yi,zi,dmin,deltaxi,deltazi
      real*8 vcirc,period,omega
ccc
      if(horbit.eq.1) then
ccc
ccc parametrization: \vec{r}_0 = \vec{r}_i + \vec{v} * (t0-ti)
ccc in case t0=0 corresponds to the source being the closest point to
ccc the origin, i.e. \vec{r}=0, along the trajectory. you find
ccc   ti= \vec{v}\cdot\vec{r}_i / |\vec{v}|^2
ccc   dmin^2 = |\vec{r}_i|^2 - (\vec{v}\cdot\vec{r}_i)^2 / |\vec{v}|^2
ccc we choose the reference frame fixing v_y=0 and y_i=y0in; we also
ccc assume that the source comes in from the vertical boundary, namely
ccc at: 
ccc    z_i = - sign(v_z) diffhh
ccc and find:
ccc    x_i = (v_x/v_z)*z_i +/- v2/v_z*dsqrt(dmin^2-y0in^2) 
ccc (+ or - is equivalent, we choose it so that |ti|<|tf|)
ccc
ccc note that radial boundary is not implemented in our treatment,
ccc orbits that have too low v_z would formally start at very large 
ccc x_i. to compensate for this and to avoid having to refer to 
ccc specify the reference frame in which v_y is 0, a sharp cut at
ccc |x0| > diffRh (for large diffRh and the observer not too close to
ccc the radial boundary this should be ok). note also that when too
ccc far away points starts to be included, it is unlikely that the
ccc straight line is a fair approximation.
ccc
ccc velocities in km/s
        vx=vxin
        vy=0.d0
        vz=vzin
ccc velocities changed to units of kpc/(10^15 s)  
        vx=vx/30.8567758d0
        vy=vy/30.8567758d0
        vz=vz/30.8567758d0
        v2=(vx**2+vy**2+vz**2)
ccc lengths in kcp
        yi=y0in
        dmin=dminin
ccc
ccc vertical or horizonthal boundary?
ccc
        if(diffhh*dabs(vx).le.diffRh*dabs(vz)) then
            ! the sharp cut at diffRh is introduced here to avoid 
            ! accounting for the specific orbit
          if(vz.lt.0.d0) then
            zi=diffhh
          else
            zi=-diffhh
          endif
          xi=zi*(vx/vz)
          deltaxi=v2/vz*dsqrt((dmin+y0in)*dabs(dmin-y0in))
          if(xi.gt.0.d0) then
            xi=xi-deltaxi
          else
            xi=xi+deltaxi
          endif
          ti=(vx*xi+vy*yi+vz*zi)/v2   ! times in 10^15 s
          tf=2.d0*diffhh/dabs(vz)+ti
        else
          if(vx.lt.0.d0) then
            xi=diffRh
          else
            xi=-diffRh
          endif
          zi=xi*(vz/vx)
          deltazi=v2/vx*dsqrt((dmin+y0in)*dabs(dmin-y0in))
          if(zi.gt.0.d0) then
            zi=zi-deltazi
          else
            zi=zi+deltazi
          endif
          ti=(vx*xi+vy*yi+vz*zi)/v2   ! times in 10^15 s
          tf=2.d0*diffRh/dabs(vx)+ti
        endif
ccc
ccc position at time t0
ccc
        x0=xi+vx*(t0-ti)
        y0=yi+vy*(t0-ti)
        z0=zi+vz*(t0-ti)
ccc
      elseif(horbit.eq.2) then
ccc
ccc circular orbit in the galactic plane with circular velocity vcirc: 
ccc     xp = radorb * cos(vcirc/radorb*t0)
ccc     yp = radorb * sin(vcirc/radorb*t0)
ccc as seen from the observer, which in galactic coordinates is at:
ccc     xp = rosb
ccc     yp = 0.d0
ccc while in the coordinates implemented in the positron routines is 
ccc at the origin, i.e. \vec{r}=0. It is also assumed t0=0, 
ccc corresponding to the source being the closest distance 
ccc     |radorb-rosb|
ccc the initial time ti is set equal the time when the source is 
ccc located at xp = - radorb, eventually a number of full orbits 
ccc     norbit
ccc before the current orbit. there is no time tf at which the source
ccc exits the diffusion region.
ccc
ccc circular velocity in km/s
        vcirc=vcircin
ccc circular velocity changed to units of kpc/(10^15 s)  
        vcirc=vcirc/30.8567758d0
ccc angular velocity in 10^-15 s^-1
        omega=vcirc/radorb
ccc period in 10^15 s:
        period=2.d0*pi/omega
        ti=-period/2.d0-norbit*period
        tf=1.d40
ccc
ccc position at time t0
ccc
        x0=-rosb+radorb*dcos(omega*t0)
        y0=radorb*dsin(omega*t0)
        z0=0.d0
      else
        write(*,*) 'DS: dsorbitps called with wrong horbit = ',
     &             horbit
        write(*,*) 'DS: program stopped'
        stop
      endif
ccc
      return
      end

